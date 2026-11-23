#!/usr/bin/env python3
"""Compute the receptor MSA once, with the prediction tool's own client.

    cofoldMsa.py <backend> <workDir> <backend args...>
        boltz       <msaServerUrl>
        intellifold <msaServerUrl>
        protenix    <msaServerMode> <protenixRoot>
        alphafold3  <msaServerUrl>

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


def colabfoldA3m(seq, server):
    """A3M of `seq` from a ColabFold MMseqs2 server: the UniRef hits, then the environmental ones.

    AlphaFold 3 takes an A3M and ships no client for the server, so this speaks the server's
    ticket API directly (submit, poll, download), with the standard library only.
    """
    import io
    import random
    import tarfile
    import time
    import urllib.parse
    import urllib.request

    server = server.rstrip('/')

    def call(url, data=None):
        body = urllib.parse.urlencode(data).encode() if data else None
        req = urllib.request.Request(url, data=body, headers={'User-Agent': 'scipion-chem-biofold'})
        with urllib.request.urlopen(req, timeout=120) as resp:
            return resp.read()

    query = {'q': f'>101\n{seq}\n', 'mode': 'env'}
    for _ in range(30):
        ticket = json.loads(call(f'{server}/ticket/msa', query))
        if ticket['status'] not in ('UNKNOWN', 'RATELIMIT'):
            break
        time.sleep(60)
    while ticket['status'] in ('PENDING', 'RUNNING'):
        time.sleep(5 + random.random() * 5)
        ticket = json.loads(call(f"{server}/ticket/{ticket['id']}"))
    if ticket['status'] != 'COMPLETE':
        raise RuntimeError(f"MSA server {server} ended with status {ticket['status']}")

    tar = tarfile.open(fileobj=io.BytesIO(call(f"{server}/result/download/{ticket['id']}")), mode='r:gz')
    a3ms = sorted((m for m in tar.getmembers() if m.name.endswith('.a3m')),
                  key=lambda m: 'uniref' not in m.name)  # UniRef first

    records = []
    for i, member in enumerate(a3ms):
        lines = [l for l in tar.extractfile(member).read().decode().replace('\x00', '').splitlines()
                 if l and not l.startswith('#')]
        found = [(lines[j], lines[j + 1]) for j in range(0, len(lines) - 1, 2)]
        # every file repeats the query as its first record: keep it once
        records += found if i == 0 else found[1:]
    if not records or records[0][1] != seq:
        raise RuntimeError('The MSA server did not return the query sequence first')
    return ''.join(f'{h}\n{s}\n' for h, s in records)


def alphafold3(workDir, rec, server):
    names = {seq: f'receptor_{i}' for i, seq in enumerate(uniqueSequences(rec))}
    out = {}
    for seq, name in names.items():
        a3m = workDir / 'msa' / f'{name}.a3m'
        a3m.write_text(colabfoldA3m(seq, server))
        # Only the unpaired MSA: it is the one AlphaFold 3 also uses to pair chains
        out[seq] = {'unpairedMsaPath': str(a3m.resolve()), 'pairedMsa': ''}
    return out


BACKENDS = {'boltz': boltz, 'intellifold': intellifold, 'protenix': protenix, 'alphafold3': alphafold3}


def main(backend, workDir, *args):
    workDir = Path(workDir)
    (workDir / 'msa').mkdir(parents=True, exist_ok=True)
    rec = json.load(open(workDir / 'receptor.json'))
    with open(workDir / 'msa.json', 'w') as f:
        json.dump(BACKENDS[backend](workDir, rec, *args), f)


if __name__ == '__main__':
    main(*sys.argv[1:])
