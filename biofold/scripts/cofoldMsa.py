#!/usr/bin/env python3
"""Compute the receptor MSA once, with the prediction tool's own client.

    cofoldMsa.py <backend> <workDir> <backend args...>
        boltz       <msaServerUrl>
        intellifold <msaServerUrl>
        protenix    <msaServerMode> <protenixRoot>

This is the only part of preparing a cofolding job that needs the prediction tool's conda
env, because the MSA client ships with the tool. Everything else (the job files) is plain
text and is written by the protocol itself.

Reads <workDir>/receptor.json, written by cofoldUtils.py, and writes <workDir>/msa.json:
{receptor sequence: {job-file keys that point at its MSA}}, so a protocol can paste each
entry straight into every ligand job. Identical chains share one sequence, so one MSA.

Only the standard library is imported at the top level: each backend imports its own tool
inside its function, so this file loads in any of the tool envs.
"""
import json
import os
import sys
from pathlib import Path


def uniqueSequences(rec):
    """The distinct receptor sequences, in order of appearance"""
    return list(dict.fromkeys(ch['seq'] for ch in rec['chains']))


def boltz(workDir, rec, server):
    from boltz.main import compute_msa

    names = {seq: f'receptor_{i}' for i, seq in enumerate(uniqueSequences(rec))}
    compute_msa(data={name: seq for seq, name in names.items()}, target_id='receptor',
                msa_dir=workDir / 'msa', msa_server_url=server, msa_pairing_strategy='greedy')
    return {seq: {'msa': str((workDir / 'msa' / f'{name}.csv').resolve())}
            for seq, name in names.items()}


def intellifold(workDir, rec, server):
    from intellifold.data.inference.data_tools import compute_msa

    names = {seq: f'receptor_{i}' for i, seq in enumerate(uniqueSequences(rec))}
    compute_msa(data={name: seq for seq, name in names.items()},
                generate_template={name: False for name in names.values()},
                target_id='receptor', msa_dir=workDir / 'msa', msa_server_url=server,
                msa_pairing_strategy='greedy')
    return {seq: {'msa': str((workDir / 'msa' / f'{name}.csv').resolve())}
            for seq, name in names.items()}


def protenix(workDir, rec, mode, root):
    # Protenix keeps its weights and CCD cache here; set before importing it
    os.environ.setdefault('PROTENIX_ROOT_DIR', root)
    from runner.msa_search import update_infer_json

    msaDir = workDir / 'msa'
    msaDir.mkdir(parents=True, exist_ok=True)
    # Protenix searches from a job file, so describe a protein-only job. What it returns is
    # the same job with pairedMsaPath/unpairedMsaPath filled in, which is what we keep.
    chains = {}
    for ch in rec['chains']:
        chains.setdefault(ch['seq'], []).append(ch['id'])
    job = [{'name': 'receptor', 'sequences': [
        {'proteinChain': {'sequence': seq, 'count': len(ids), 'id': ids}}
        for seq, ids in chains.items()]}]
    recJson = msaDir / 'receptor_only.json'
    recJson.write_text(json.dumps(job, indent=1))

    updated, _ = update_infer_json(str(recJson), str(msaDir), use_msa=True, mode=mode)
    paths = {}
    for entity in json.load(open(updated))[0]['sequences']:
        chain = entity.get('proteinChain')
        if chain:
            paths[chain['sequence']] = {k: chain[k] for k in ('pairedMsaPath', 'unpairedMsaPath')
                                        if chain.get(k)}
    return paths


BACKENDS = {'boltz': boltz, 'intellifold': intellifold, 'protenix': protenix}


def main(backend, workDir, *args):
    workDir = Path(workDir)
    (workDir / 'msa').mkdir(parents=True, exist_ok=True)
    rec = json.load(open(workDir / 'receptor.json'))
    with open(workDir / 'msa.json', 'w') as f:
        json.dump(BACKENDS[backend](workDir, rec, *args), f)


if __name__ == '__main__':
    main(*sys.argv[1:])
