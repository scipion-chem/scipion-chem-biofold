import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


class TestMsaSubsamplingConfiguration(unittest.TestCase):
    def test_boltz_subsampling_controls_are_present(self):
        text = (REPO_ROOT / 'biofold' / 'protocols' / 'protocol_boltz.py').read_text()

        self.assertIn("group.addParam('subsampleMsa'", text)
        self.assertIn("group.addParam('numSubsampledMsa'", text)
        self.assertIn("--subsample_msa", text)
        self.assertIn("--num_subsampled_msa", text)

    def test_other_models_expose_msa_depth_flags(self):
        expectations = {
            REPO_ROOT / 'biofold' / 'protocols' / 'protocol_chai.py': (
                "form.addParam('maxMsaSeqs'",
                "--max_msa_seqs",
            ),
            REPO_ROOT / 'biofold' / 'protocols' / 'protocol_intellifold.py': (
                "group.addParam('msaMaxDepth'",
                "--msa_max_depth",
            ),
            REPO_ROOT / 'biofold' / 'protocols' / 'protocol_protenix.py': (
                "group.addParam('maxMsaDepth'",
                "--max_msa_depth",
            ),
        }

        for path, snippets in expectations.items():
            with self.subTest(path=path.name):
                text = path.read_text()
                for snippet in snippets:
                    self.assertIn(snippet, text)


if __name__ == '__main__':
    unittest.main()
