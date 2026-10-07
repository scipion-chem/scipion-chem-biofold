#!/usr/bin/env python3
"""Boltz-2 cofolding helper, run inside the Boltz conda env.

prepare <receptor|SEQ1:SEQ2...> <ligands.json> <workDir> <chains|all> <msaServerUrl|none>
    Reads the receptor chains (from a structure file, or ':'-separated protein
    sequences named A, B, ...), computes their MSA once and writes one Boltz
    YAML per ligand (all sharing that MSA) into <workDir>/inputs.
extract <workDir> <predDir> <outDir>
    Superposes every predicted complex onto the input receptor (CA Kabsch), or
    onto the first prediction when the receptor was given as sequences, and
    writes the ligand pose (SDF), the predicted receptor (PDB) and the scores
    (<outDir>/results.json). With sequence input, that first predicted receptor
    is also written as <outDir>/reference_receptor.pdb.
"""
import json
import os
import shutil
import string
import sys
from pathlib import Path

import gemmi
import numpy as np
import yaml
from rdkit import Chem
from rdkit.Chem import AllChem

from boltz.data.parse.schema import standardize


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


def ligandSmiles(molFile):
    """Canonical SMILES of the largest fragment of the first molecule in an SDF"""
    mol = Chem.MolFromMolFile(molFile)
    if mol is None:
        return None
    mol = max(Chem.GetMolFrags(mol, asMols=True), key=lambda m: m.GetNumHeavyAtoms())
    return Chem.MolToSmiles(Chem.RemoveHs(mol))


def prepare(recFile, ligandsJson, workDir, chains, msaServer):
    workDir = Path(workDir)
    inputsDir, msaDir = workDir / 'inputs', workDir / 'msa'
    inputsDir.mkdir(parents=True, exist_ok=True)
    msaDir.mkdir(parents=True, exist_ok=True)

    if os.path.isfile(recFile):
        recChains = readReceptorChains(recFile, [] if chains == 'all' else [c.strip() for c in chains.split(',')])
    else:
        recChains = sequenceChains(recFile)
    if not recChains:
        sys.exit(f'No protein chains found in {recFile}')

    # Identical chains share one entity (and one MSA)
    entities = {}
    for ch in recChains:
        entities.setdefault(ch['seq'], []).append(ch['id'])

    # The MSA is computed once here and reused by every ligand job
    msaFiles = {}
    if msaServer != 'none':
        from boltz.main import compute_msa
        names = {seq: f'receptor_{i}' for i, seq in enumerate(entities)}
        compute_msa(data={names[s]: s for s in entities}, target_id='receptor', msa_dir=msaDir,
                    msa_server_url=msaServer, msa_pairing_strategy='greedy')
        msaFiles = {s: str((msaDir / f'{names[s]}.csv').resolve()) for s in entities}

    proteins = [{'protein': {'id': ids if len(ids) > 1 else ids[0], 'sequence': seq,
                             'msa': msaFiles.get(seq, 'empty')}}
                for seq, ids in entities.items()]
    ligChain = next(c for c in string.ascii_uppercase if c not in {ch['id'] for ch in recChains})

    ligands = {}
    for lig in json.load(open(ligandsJson)):
        smi = ligandSmiles(lig['file'])
        if not smi:
            print(f"Warning: could not read {lig['file']}; ligand {lig['id']} skipped")
            continue
        ligands[lig['id']] = smi
        job = {'version': 1,
               'sequences': proteins + [{'ligand': {'id': ligChain, 'smiles': smi}}],
               'properties': [{'affinity': {'binder': ligChain}}]}
        with open(inputsDir / f"{lig['id']}.yaml", 'w') as f:
            yaml.safe_dump(job, f, sort_keys=False)

    with open(workDir / 'receptor.json', 'w') as f:
        json.dump({'ligandChain': ligChain, 'chains': recChains, 'ligands': ligands}, f)


def kabsch(P, Q):
    """R, t minimising |R·P + t - Q|"""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def templateMol(smiles):
    """Heavy-atom ligand mol named exactly as Boltz names it (boltz/data/parse/schema.py)"""
    mol = AllChem.AddHs(AllChem.MolFromSmiles(standardize(smiles)))
    for atom, rank in zip(mol.GetAtoms(), AllChem.CanonicalRankAtoms(mol)):
        atom.SetProp('name', atom.GetSymbol().upper() + str(rank + 1))
    return AllChem.RemoveHs(mol)


def posedLigand(smiles, ligChain):
    """The template mol with the predicted coordinates of the ligand chain"""
    coords = {a.name: a.pos for res in ligChain for a in res}
    mol = templateMol(smiles)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for atom in mol.GetAtoms():
        p = coords[atom.GetProp('name')]
        conf.SetAtomPosition(atom.GetIdx(), (p.x, p.y, p.z))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    Chem.AssignStereochemistryFrom3D(mol)
    return mol


def extract(workDir, predDir, outDir):
    workDir, predDir, outDir = Path(workDir), Path(predDir), Path(outDir)
    outDir.mkdir(parents=True, exist_ok=True)
    rec = json.load(open(workDir / 'receptor.json'))
    # Sequence input has no structure: the first prediction sets the common frame
    seqInput = rec['chains'][0]['ca'] is None
    refCA = None if seqInput else np.array([xyz for ch in rec['chains'] for xyz in ch['ca']])

    results = []
    for ligId, smiles in rec['ligands'].items():
        ligDir = predDir / ligId
        if not ligDir.is_dir():
            print(f'Warning: no Boltz prediction for {ligId}; skipped')
            continue
        affFile = ligDir / f'affinity_{ligId}.json'
        aff = json.load(open(affFile)) if affFile.exists() else {}

        for k in range(len(list(ligDir.glob(f'{ligId}_model_*.cif')))):
            st = gemmi.read_structure(str(ligDir / f'{ligId}_model_{k}.cif'))
            model = st[0]
            predCA = np.array([[a.pos.x, a.pos.y, a.pos.z] for ch in rec['chains']
                               for res in model[ch['id']] for a in [res.find_atom('CA', '*')]])
            if refCA is None:
                refCA = predCA
            R, t = kabsch(predCA, refCA)
            rmsd = float(np.sqrt((((predCA @ R.T + t) - refCA) ** 2).sum(1).mean()))
            for ch in model:
                for res in ch:
                    for a in res:
                        a.pos = gemmi.Position(*(R @ [a.pos.x, a.pos.y, a.pos.z] + t))

            mol = posedLigand(smiles, model[rec['ligandChain']])
            mol.SetProp('_Name', ligId)
            poseFile = outDir / f'{ligId}_{k + 1}.sdf'
            with Chem.SDWriter(str(poseFile)) as w:
                w.write(mol)

            model.remove_chain(rec['ligandChain'])
            recFile = outDir / f'{ligId}_{k + 1}_receptor.pdb'
            st.write_pdb(str(recFile))
            if seqInput and not (outDir / 'reference_receptor.pdb').exists():
                shutil.copy(recFile, outDir / 'reference_receptor.pdb')

            conf = json.load(open(ligDir / f'confidence_{ligId}_model_{k}.json'))
            results.append({'id': ligId, 'poseId': k + 1,
                            'poseFile': str(poseFile.resolve()), 'receptorFile': str(recFile.resolve()),
                            'affinity': aff.get('affinity_pred_value'),
                            'bindingProbability': aff.get('affinity_probability_binary'),
                            'confidence': conf.get('confidence_score'),
                            'ligandIpTM': conf.get('ligand_iptm'),
                            'plddt': conf.get('complex_plddt'),
                            'iplddt': conf.get('complex_iplddt'),
                            'receptorRMSD': rmsd})

    with open(outDir / 'results.json', 'w') as f:
        json.dump(results, f, indent=1)


if __name__ == '__main__':
    mode, args = sys.argv[1], sys.argv[2:]
    {'prepare': prepare, 'extract': extract}[mode](*args)
