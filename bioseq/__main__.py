"""Command-line entry point for the reproducible portfolio demonstration."""

import argparse
import json
from pathlib import Path

from .algorithms import build_profile, generate_sequences, progressive_align


def _write_lines(path: Path, lines):
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _read_sequences(path: Path):
    sequences = [line.strip().upper() for line in path.read_text(encoding="utf-8").splitlines()
                 if line.strip()]
    if not sequences:
        raise ValueError("{} contains no sequences".format(path))
    return sequences


def _threshold(value):
    number = float(value)
    if not 0 < number <= 1:
        raise argparse.ArgumentTypeError("must be greater than 0 and at most 1")
    return number


def analyze(sequences, output: Path, threshold: float):
    output.mkdir(parents=True, exist_ok=True)
    alignment, tree = progressive_align(sequences)
    profile = build_profile(alignment, threshold)
    _write_lines(output / "alignment.txt", alignment)
    (output / "guide_tree.nwk").write_text(tree.newick() + ";\n", encoding="utf-8")
    (output / "profile.json").write_text(json.dumps(profile, indent=2) + "\n", encoding="utf-8")
    print("Aligned {} sequence{} across {} columns".format(
        len(alignment), "" if len(alignment) == 1 else "s", len(alignment[0])))
    print("Profile contains {} match states".format(len(profile["match_columns_1_based"])))
    print("Results: {}".format(output.resolve()))


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="bioseq",
        description="DNA alignment and profile-HMM demonstration",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    demo = subparsers.add_parser("demo", help="generate a seeded dataset and analyze its first 15 sequences")
    demo.add_argument("--seed", type=int, default=2024)
    demo.add_argument("--output", type=Path, default=Path("results"))
    demo.add_argument("--threshold", type=_threshold, default=0.7)
    existing = subparsers.add_parser("analyze", help="analyze one DNA sequence per input line")
    existing.add_argument("input", type=Path)
    existing.add_argument("--output", type=Path, default=Path("results"))
    existing.add_argument("--threshold", type=_threshold, default=0.7)
    args = parser.parse_args(argv)
    try:
        if args.command == "demo":
            sequences = generate_sequences(seed=args.seed)
            args.output.mkdir(parents=True, exist_ok=True)
            _write_lines(args.output / "all_sequences.txt", sequences)
            _write_lines(args.output / "dataset_a.txt", sequences[:15])
            _write_lines(args.output / "dataset_b.txt", sequences[15:])
            analyze(sequences[:15], args.output, args.threshold)
        else:
            analyze(_read_sequences(args.input), args.output, args.threshold)
    except (OSError, ValueError) as error:
        parser.exit(2, "error: {}\n".format(error))


if __name__ == "__main__":
    main()
