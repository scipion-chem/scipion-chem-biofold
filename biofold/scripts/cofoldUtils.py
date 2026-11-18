#!/usr/bin/env python3
"""Backend-independent section of the cofolding protocols.

    cofoldUtils.py receptor <receptor|SEQ1:SEQ2...> <ligands.json> <workDir> <chains|all>
        Writes <workDir>/receptor.json: the receptor chains (with their CA coordinates
        when a structure was given), the ligand SMILES and the free ligand chain id.

    cofoldUtils.py extract <backend> <workDir> <predDir> <outDir>
        Superposes every predicted complex onto the input receptor, writes each ligand
        pose (SDF), its predicted receptor (PDB) and the scores (<outDir>/results.json).

Both run in pwchem's RDKit env rather than in the prediction tool's own env: they need
only gemmi, numpy and rdkit, not every tool env ships gemmi, and this way one known rdkit
canonicalises the SMILES. Each backend keeps only a small prepare script of its own, for
its MSA client and input format.

Nothing here may import a prediction tool, and the ligand pose is rebuilt by molecular
graph rather than from any backend's atom naming.
"""
import glob
import json
import os
import shutil
import string
import sys
from pathlib import Path

import gemmi
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdDetermineBonds

RDLogger.DisableLog('rdApp.*')

# --------------------------------------------------------------------------- receptor
def readReceptorChains(recFile, chains):
    """[{id, seq, ca}] for every protein chain (or only `chains`) of the receptor.
    Only residues with a CA are kept, so seq and ca stay index-aligned."""
    st = gemmi.read_structure(recFile)
    out = []
    for ch in st[0]:
        if chains and ch.name not in chains:
            continue
        seq, ca = '', []
        for res in ch.first_conformer():
            info = gemmi.find_tabulated_residue(res.name)
            code = info.one_letter_code.upper() if info and info.is_amino_acid() else ''
            atom = res.find_atom('CA', '*')
            if code.isalpha() and atom:
                seq += code
                ca.append([atom.pos.x, atom.pos.y, atom.pos.z])
        if seq:
            out.append({'id': ch.name, 'seq': seq, 'ca': ca})
    return out


def sequenceChains(sequences):
    """[{id, seq, ca=None}] for ':'-separated protein sequences, named A, B, ..."""
    return [{'id': cid, 'seq': seq.strip().upper(), 'ca': None}
            for cid, seq in zip(string.ascii_uppercase, sequences.split(':')) if seq.strip()]


def readReceptor(recFile, chains):
    """Receptor chains from a structure file, or from ':'-separated sequences"""
    if os.path.isfile(recFile):
        return readReceptorChains(recFile, [] if chains == 'all' else
                                  [c.strip() for c in chains.split(',')])
    return sequenceChains(recFile)


def freeChainId(usedIds):
    """First uppercase letter not taken by a receptor chain, for the ligand"""
    return next(c for c in string.ascii_uppercase if c not in set(usedIds))


# --------------------------------------------------------------------------- ligands
def ligandSmiles(lig):
    """Canonical SMILES of the largest fragment of a ligand given as a SMILES string or
    as an SDF file"""
    mol = Chem.MolFromSmiles(lig['smiles']) if 'smiles' in lig else Chem.MolFromMolFile(lig['file'])
    if mol is None:
        return None
    mol = max(Chem.GetMolFrags(mol, asMols=True), key=lambda m: m.GetNumHeavyAtoms())
    return Chem.MolToSmiles(Chem.RemoveHs(mol))


def rebuildLigand(smiles, ligChain):
    """The predicted ligand as an RDKit mol with correct bond orders.

    Matched to the SMILES template by molecular graph, not by atom name: every backend
    names ligand atoms differently (canonical rank, per-element counter, element+index),
    and those names are an implementation detail we must not depend on. Connectivity is
    perceived from the predicted coordinates and the bond orders come from the template.
    """
    tmpl = Chem.MolFromSmiles(smiles)
    pdbBlock = ''
    for i, atom in enumerate((a for res in ligChain for a in res), start=1):
        element = atom.element.name
        pdbBlock += ('HETATM%5d %-4s LIG A   1    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n'
                     % (i, atom.name[:4], atom.pos.x, atom.pos.y, atom.pos.z, element.upper()))

    raw = Chem.MolFromPDBBlock(pdbBlock, sanitize=False, removeHs=True, proximityBonding=False)
    if raw is None:
        raise ValueError('could not read the predicted ligand atoms')
    rdDetermineBonds.DetermineConnectivity(raw)
    # the PDB reader marks every atom "no implicit Hs"; clear it or the valences stay fixed
    for a in raw.GetAtoms():
        a.SetNoImplicit(False)
        a.SetNumExplicitHs(0)
    raw.UpdatePropertyCache(strict=False)

    mol = AllChem.AssignBondOrdersFromTemplate(tmpl, raw)
    Chem.AssignStereochemistryFrom3D(mol)
    return mol


# --------------------------------------------------------------------------- geometry
def kabsch(P, Q):
    """R, t minimising |R·P + t - Q|"""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def chainCAs(model, chainIds):
    """CA coordinates of the given chains, in chain then residue order"""
    return np.array([[a.pos.x, a.pos.y, a.pos.z] for cid in chainIds
                     for res in model[cid] for a in [res.find_atom('CA', '*')] if a])


def superpose(model, predCA, refCA):
    """Move every atom of the model onto the reference frame. Returns the CA RMSD."""
    R, t = kabsch(predCA, refCA)
    for ch in model:
        for res in ch:
            for a in res:
                a.pos = gemmi.Position(*(R @ [a.pos.x, a.pos.y, a.pos.z] + t))
    return float(np.sqrt((((predCA @ R.T + t) - refCA) ** 2).sum(1).mean()))


# --------------------------------------------------------------------------- scores
def ligandPLDDT(ligChain):
    """Mean pLDDT over the ligand atoms only, on a 0-1 scale.

    Every backend writes pLDDT into the B-factor column of its predicted structure
    (Boltz and Protenix as plddt*100, Chai already on 0-100), so this is the one
    measure of pose confidence that can be read the same way for all of them.
    """
    return meanBFactor([a for res in ligChain for a in res])


def complexPLDDT(model):
    """Mean pLDDT over every atom of the model, 0-1. For backends that do not report a
    complex-wide pLDDT of their own."""
    return meanBFactor([a for ch in model for res in ch for a in res])


def meanBFactor(atoms):
    """Mean B-factor of `atoms` read as a 0-1 pLDDT, or None when the backend wrote none.
    """
    bs = [a.b_iso for a in atoms]
    if not bs or not any(bs):
        return None
    mean = float(np.mean(bs))
    return mean / 100 if mean > 1.5 else mean


def ligandIpTM(pairMatrix, ligChain, protChains):
    """Protein-ligand interface ipTM from a per-chain-pair matrix.

    Taken as the best interface the ligand makes with any single protein chain. A mean
    would punish a ligand bound to one chain of a dimer for ignoring the other, while a
    ligand sitting at a chain-chain interface still scores high with both.
    """
    def cell(a, b):
        try:
            row = pairMatrix[a] if isinstance(pairMatrix, (list, tuple)) else (
                pairMatrix.get(a, pairMatrix.get(str(a))))
            if row is None:
                return None
            return row[b] if isinstance(row, (list, tuple)) else row.get(b, row.get(str(b)))
        except (IndexError, KeyError, TypeError):
            return None

    vals = [float(v) for p in protChains for v in (cell(ligChain, p), cell(p, ligChain))
            if v is not None]
    return max(vals) if vals else None


# --------------------------------------------------------------------- receptor.json
def describe(recFile, ligandsJson, workDir, chains):
    workDir = Path(workDir)
    workDir.mkdir(parents=True, exist_ok=True)

    recChains = readReceptor(recFile, chains)
    if not recChains:
        sys.exit(f'No protein chains found in {recFile}')

    ligands = {}
    for lig in json.load(open(ligandsJson)):
        smi = ligandSmiles(lig)
        if not smi:
            print(f"Warning: could not read ligand {lig['id']} "
                  f"({lig.get('smiles') or lig['file']}); skipped")
            continue
        ligands[lig['id']] = smi
    if not ligands:
        sys.exit('No ligand could be read')

    with open(workDir / 'receptor.json', 'w') as f:
        json.dump({'ligandChain': freeChainId(ch['id'] for ch in recChains),
                   'chains': recChains, 'ligands': ligands}, f)


# --------------------------------------------------------------------- backends
def boltzModels(predDir, ligId, rec):
    """Boltz-2: <ligId>/<ligId>_model_<k>.cif, one confidence json each, one affinity json"""
    ligDir = predDir / ligId
    affFile = ligDir / f'affinity_{ligId}.json'
    aff = json.load(open(affFile)) if affFile.exists() else {}

    for k in range(len(list(ligDir.glob(f'{ligId}_model_*.cif')))):
        conf = json.load(open(ligDir / f'confidence_{ligId}_model_{k}.json'))
        yield ligDir / f'{ligId}_model_{k}.cif', {
            'confidence': conf.get('confidence_score'),
            # Boltz computes the protein-ligand ipTM itself, pooled over all protein chains
            'ligIpTM': conf.get('ligand_iptm'),
            'plddt': conf.get('complex_plddt'),
            'affinity': aff.get('affinity_pred_value'),
            'bindingProbability': aff.get('affinity_probability_binary')}


def intellifoldModels(predDir, ligId, rec):
    """IntelliFold: <ligId>/<ligId>_seed-<s>_sample-<i>.cif + _summary_confidences.json"""
    ligDir = predDir / ligId
    # its chain_pair_iptm is a square matrix indexed by chain order, receptor chains first
    nProt = len(rec['chains'])
    for cif in sorted(ligDir.glob(f'{ligId}_seed-*_sample-*.cif')):
        conf = json.load(open(str(cif)[:-4] + '_summary_confidences.json'))
        yield cif, {
            'confidence': conf.get('ranking_score'),
            'ligIpTM': ligandIpTM(conf.get('chain_pair_iptm') or {}, nProt, range(nProt)),
            'plddt': conf.get('plddt')}


def squeezed(value, ndim):
    """`value` as a nested list of exactly `ndim` dimensions, dropping a leading sample
    axis if the backend left one in, or None when it is missing or has another shape"""
    if value is None:
        return None
    arr = np.squeeze(np.array(value, dtype=float))
    return arr.tolist() if arr.ndim == ndim else None


def protenixModels(predDir, ligId, rec):
    """Protenix: <ligId>/seed_<s>/predictions/<name>_sample_<rank>.cif + summary json"""
    # chain_pair_iptm is a dense, symmetric matrix indexed by chain order, receptor chains
    # first and the ligand last: not by chain letter
    nProt = len(rec['chains'])
    for cif in sorted(Path(p) for p in glob.glob(str(predDir / ligId / 'seed_*' /
                                                     'predictions' / '*_sample_*.cif'))):
        name, rank = cif.stem.rsplit('_sample_', 1)
        confFile = cif.parent / f'{name}_summary_confidence_sample_{rank}.json'
        conf = json.load(open(confFile)) if confFile.exists() else {}
        pairs = squeezed(conf.get('chain_pair_iptm'), 2)
        chainPlddt = squeezed(conf.get('chain_plddt'), 1)
        plddt = conf.get('plddt')
        yield cif, {
            'confidence': conf.get('ranking_score'),
            'ligIpTM': ligandIpTM(pairs, nProt, range(nProt)) if pairs is not None else None,
            'ligPLDDT': float(chainPlddt[nProt]) if chainPlddt is not None and len(chainPlddt) > nProt else None,
            'plddt': plddt / 100 if isinstance(plddt, (int, float)) else None}


def chainIds(backend, rec):
    """(receptor chain ids, ligand chain id) as the backend named them in its output.

    Boltz, Protenix and IntelliFold are told the ids to use and honour them. Chai names
    chains A, B, C... by FASTA record order and ignores whatever we asked for, so with a
    receptor whose chains are not already named from A the ids would not line up.
    """
    if backend == 'chai':
        letters = [chr(ord('A') + i) for i in range(len(rec['chains']) + 1)]
        return letters[:-1], letters[-1]
    return [ch['id'] for ch in rec['chains']], rec['ligandChain']


def chaiModels(predDir, ligId, rec):
    """Chai-1: <ligId>/pred.model_idx_<i>.cif + scores.model_idx_<i>.npz"""
    ligDir = predDir / ligId
    nProt = len(rec['chains'])
    for cif in sorted(ligDir.glob('pred.model_idx_*.cif')):
        idx = cif.stem.split('_idx_')[1]
        scoreFile = ligDir / f'scores.model_idx_{idx}.npz'
        scores = dict(np.load(scoreFile)) if scoreFile.exists() else {}
        # chai indexes its per-chain matrices by chain order, receptor chains first
        pairs = scores.get('per_chain_pair_iptm')
        pairDic = {i: {j: v for j, v in enumerate(row)}
                   for i, row in enumerate(np.atleast_2d(np.squeeze(pairs)))} if pairs is not None else {}
        yield cif, {
            'confidence': float(np.squeeze(scores['aggregate_score'])) if 'aggregate_score' in scores else None,
            'ligIpTM': ligandIpTM(pairDic, nProt, range(nProt)),
            # chai writes no complex pLDDT; the ligand one is read from the B-factors
            'plddt': None}


BACKENDS = {'boltz': boltzModels, 'chai': chaiModels,
            'protenix': protenixModels, 'intellifold': intellifoldModels}


# --------------------------------------------------------------------- extraction
def extract(backend, workDir, predDir, outDir):
    workDir, predDir, outDir = Path(workDir), Path(predDir), Path(outDir)
    outDir.mkdir(parents=True, exist_ok=True)
    rec = json.load(open(workDir / 'receptor.json'))
    protIds, ligChainId = chainIds(backend, rec)
    # Sequence input has no structure: the first prediction sets the common frame
    seqInput = rec['chains'][0]['ca'] is None
    refCA = None if seqInput else np.array([xyz for ch in rec['chains'] for xyz in ch['ca']])

    results = []
    for ligId, smiles in rec['ligands'].items():
        if not (predDir / ligId).is_dir():
            print(f'Warning: no {backend} prediction for {ligId}; skipped')
            continue

        for poseId, (cif, scores) in enumerate(BACKENDS[backend](predDir, ligId, rec), start=1):
            st = gemmi.read_structure(str(cif))
            st.setup_entities()
            model = st[0]
            predCA = chainCAs(model, protIds)
            if refCA is None:
                refCA = predCA
            if len(predCA) != len(refCA):
                print(f'Warning: {cif.name} has {len(predCA)} CA against {len(refCA)} in the '
                      f'receptor; pose skipped')
                continue
            rmsd = superpose(model, predCA, refCA)

            ligChain = model[ligChainId]
            try:
                mol = rebuildLigand(smiles, ligChain)
            except Exception as e:
                print(f'Warning: could not rebuild the ligand of {cif.name} ({e}); pose skipped')
                continue
            mol.SetProp('_Name', ligId)
            poseFile = outDir / f'{ligId}_{poseId}.sdf'
            with Chem.SDWriter(str(poseFile)) as w:
                w.write(mol)

            # Taken from the backend's own summary when found  otherwise from the B-factors
            if scores.get('ligPLDDT') is None:
                scores['ligPLDDT'] = ligandPLDDT(ligChain)
            if scores.get('plddt') is None:
                scores['plddt'] = complexPLDDT(model)
            model.remove_chain(ligChainId)
            recFile = outDir / f'{ligId}_{poseId}_receptor.pdb'
            st.write_pdb(str(recFile))
            if seqInput and not (outDir / 'reference_receptor.pdb').exists():
                shutil.copy(recFile, outDir / 'reference_receptor.pdb')

            scores.update(id=ligId, poseId=poseId, receptorRMSD=rmsd,
                          poseFile=str(poseFile.resolve()), receptorFile=str(recFile.resolve()))
            results.append(scores)

    with open(outDir / 'results.json', 'w') as f:
        json.dump(results, f, indent=1)


if __name__ == '__main__':
    mode, args = sys.argv[1], sys.argv[2:]
    {'receptor': describe, 'extract': extract}[mode](*args)
