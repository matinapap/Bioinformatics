import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from bioseq.algorithms import align, build_profile, generate_sequences, progressive_align


class BioseqTests(unittest.TestCase):
    def test_global_alignment_handles_terminal_gaps(self):
        left, right, score = align("ACGT", "ACG")
        self.assertEqual((left, right, score), ("ACGT", "ACG-", 1))
        self.assertEqual(align("", "AG"), ("--", "AG", -4))

    def test_progressive_alignment_preserves_input(self):
        sequences = ["ACGT", "AGT", "ACGTT", "ACT"]
        alignment, tree = progressive_align(sequences)
        self.assertEqual([row.replace("-", "") for row in alignment], sequences)
        self.assertEqual(len({len(row) for row in alignment}), 1)
        self.assertTrue(tree.newick().startswith("("))

    def test_profile_probabilities_and_deletion_path(self):
        profile = build_profile(["ACG", "A-G", "ATG", "ACG"])
        self.assertEqual(profile["match_columns_1_based"], [1, 2, 3])
        self.assertIn("D2", profile["transitions"])
        for probabilities in list(profile["emissions"].values()) + list(profile["transitions"].values()):
            self.assertAlmostEqual(sum(probabilities.values()), 1)

    def test_seeded_demo_writes_readable_outputs(self):
        self.assertEqual(generate_sequences(seed=7), generate_sequences(seed=7))
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            subprocess.run([sys.executable, "-m", "bioseq", "demo", "--output", str(output)],
                           check=True, capture_output=True, text=True)
            self.assertEqual(len((output / "dataset_a.txt").read_text().splitlines()), 15)
            self.assertEqual(len((output / "dataset_b.txt").read_text().splitlines()), 35)
            self.assertEqual(len((output / "alignment.txt").read_text().splitlines()), 15)
            self.assertIn("emissions", json.loads((output / "profile.json").read_text()))
            self.assertTrue((output / "guide_tree.nwk").read_text().endswith(";\n"))

    def test_invalid_input_reports_error(self):
        with self.assertRaises(ValueError):
            progressive_align(["ACGT", "ACXT"])
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "empty.txt"
            source.write_text("\n", encoding="utf-8")
            result = subprocess.run([sys.executable, "-m", "bioseq", "analyze", str(source)],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn("contains no sequences", result.stderr)
            output = Path(directory) / "unused"
            result = subprocess.run([sys.executable, "-m", "bioseq", "demo", "--threshold", "1.5",
                                     "--output", str(output)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 2)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
