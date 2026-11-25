#!/usr/bin/env python3
"""Check the stereochemistry of each pose against the ligand it was asked to fold.

    cofoldStereo.py <workDir> <outDir>

This runs PoseBusters' own identity check on each pose written by cofoldUtils.py extract,
against a reference built straight from the input SMILES (the check compares InChI layers,
so the reference needs no coordinates). Its `stereo_tetrahedral` result is the one used by
the OpenBind blind-challenge scorer.

A ligand whose input leaves a stereocentre open is not checked at all. The check compares
whole InChI layers, so a pose always carries a handedness the open reference lacks, and
every pose of such a ligand would be reported as wrong.

Reads <outDir>/results.json and <workDir>/receptor.json; writes <outDir>/stereo.json:
{pose file: {"stereoOK": bool or null, "note": str or null}}. Dropping the failing poses
is left to the protocol.
"""
import json
import sys
from pathlib import Path

from rdkit import Chem, RDLogger
from posebusters import check_identity

RDLogger.DisableLog('rdApp.*')

#: recorded when the input itself leaves the stereochemistry open
NOT_SPECIFIED = 'the input leaves a stereocentre open'


def openCentres(mol):
    """The stereocentres the molecule leaves unspecified"""
    return [i for i, label in Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                                        useLegacyImplementation=False)
            if label == '?']


def poseStereoOK(poseFile, smiles):
    """(stereoOK, note). stereoOK is None when the pose was not checked."""
    ref = Chem.MolFromSmiles(smiles)
    if ref is None:
        return None, f'unreadable reference SMILES: {smiles}'
    # Nothing to contradict, and an InChI comparison would call every pose wrong
    if openCentres(ref):
        return None, NOT_SPECIFIED
    pose = Chem.MolFromMolFile(str(poseFile))
    if pose is None:
        return None, 'unreadable pose'
    try:
        results = check_identity(pose, ref)['results']
    except Exception as e:
        return None, f'{type(e).__name__}: {e}'
    return bool(results.get('stereo_tetrahedral', True)), None


def gate(workDir, outDir):
    workDir, outDir = Path(workDir), Path(outDir)
    rec = json.load(open(workDir / 'receptor.json'))
    results = json.load(open(outDir / 'results.json'))

    verdicts, failed, open_, unchecked = {}, [], set(), []
    for res in results:
        stereoOK, note = poseStereoOK(res['poseFile'], rec['ligands'][res['id']])
        verdicts[res['poseFile']] = {'stereoOK': stereoOK, 'note': note}
        if stereoOK is False:
            failed.append(f"{res['id']} pose {res['poseId']}")
        elif note == NOT_SPECIFIED:
            open_.add(res['id'])
        elif stereoOK is None:
            unchecked.append(f"{res['id']} pose {res['poseId']}: {note}")

    with open(outDir / 'stereo.json', 'w') as f:
        json.dump(verdicts, f, indent=1)

    checked = sum(1 for v in verdicts.values() if v['stereoOK'] is not None)
    print(f'Stereochemistry: {checked - len(failed)} of {checked} checked poses match the input '
          f'stereochemistry ({len(results)} poses in all)')
    if failed:
        print(f'Wrong stereochemistry ({len(failed)}): {", ".join(failed)}')
    if open_:
        print(f'Not checked, {NOT_SPECIFIED}, so any pose is acceptable: '
              f'{", ".join(sorted(open_))}')
    for line in unchecked:
        print(f'Warning: could not check {line}')


if __name__ == '__main__':
    gate(*sys.argv[1:])
