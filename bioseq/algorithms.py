"""Small, dependency-free implementations of the project's core algorithms."""

from collections import Counter, defaultdict
from dataclasses import dataclass
from random import Random
from typing import Dict, List, Optional, Sequence, Tuple

DNA = "ACGT"
MOTIFS = ("AATTGA", "CGCTTAT", "GGACTCAT", "TTATTCGTA")


def generate_sequences(count: int = 50, seed: int = 2024) -> List[str]:
    """Create shuffled DNA strings with mutated copies of four shared motifs."""
    if count < 1:
        raise ValueError("count must be positive")
    rng = Random(seed)
    sequences = []
    for _ in range(count):
        parts = ["".join(rng.choices(DNA, k=rng.randint(1, 3)))]
        for motif in MOTIFS:
            letters = list(motif)
            for _ in range(rng.randint(1, 2)):
                position = rng.randrange(len(letters))
                replacement = rng.choice([base for base in DNA if base != letters[position]] + [""])
                letters[position:position + 1] = list(replacement)
            parts.append("".join(letters))
            parts.append("".join(rng.choices(DNA, k=rng.randint(1, 2))))
        sequences.append("".join(parts))
    rng.shuffle(sequences)
    return sequences


def align(a: str, b: str, match: int = 1, mismatch: int = -1,
          gap: int = -2) -> Tuple[str, str, int]:
    """Needleman-Wunsch global alignment with a linear gap penalty."""
    rows, cols = len(a) + 1, len(b) + 1
    scores = [[0] * cols for _ in range(rows)]
    moves = [[""] * cols for _ in range(rows)]
    for i in range(1, rows):
        scores[i][0], moves[i][0] = i * gap, "U"
    for j in range(1, cols):
        scores[0][j], moves[0][j] = j * gap, "L"
    for i in range(1, rows):
        for j in range(1, cols):
            choices = (
                (scores[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch), "D"),
                (scores[i - 1][j] + gap, "U"),
                (scores[i][j - 1] + gap, "L"),
            )
            scores[i][j], moves[i][j] = max(choices, key=lambda choice: choice[0])
    left, right = [], []
    i, j = len(a), len(b)
    while i or j:
        move = moves[i][j]
        if move == "D":
            left.append(a[i - 1])
            right.append(b[j - 1])
            i, j = i - 1, j - 1
        elif move == "U":
            left.append(a[i - 1])
            right.append("-")
            i -= 1
        else:
            left.append("-")
            right.append(b[j - 1])
            j -= 1
    return "".join(reversed(left)), "".join(reversed(right)), scores[-1][-1]


def _distance(a: str, b: str) -> float:
    aligned_a, aligned_b, _ = align(a, b)
    return sum(x != y for x, y in zip(aligned_a, aligned_b)) / len(aligned_a)


@dataclass
class Tree:
    members: Tuple[int, ...]
    height: float = 0.0
    left: Optional["Tree"] = None
    right: Optional["Tree"] = None

    def newick(self, parent_height: Optional[float] = None) -> str:
        if self.left is None:
            label = "seq{:02d}".format(self.members[0] + 1)
        else:
            label = "({},{})".format(self.left.newick(self.height), self.right.newick(self.height))
        if parent_height is not None:
            label += ":{:.5f}".format(max(0.0, parent_height - self.height))
        return label


def guide_tree(sequences: Sequence[str]) -> Tree:
    """Build an average-linkage (UPGMA) guide tree from alignment distances."""
    if not sequences:
        raise ValueError("at least one sequence is required")
    pair = {(i, j): _distance(sequences[i], sequences[j])
            for i in range(len(sequences)) for j in range(i + 1, len(sequences))}
    clusters = [Tree((i,)) for i in range(len(sequences))]
    while len(clusters) > 1:
        def average(first: Tree, second: Tree) -> float:
            values = [pair[min(i, j), max(i, j)] for i in first.members for j in second.members]
            return sum(values) / len(values)

        distance, i, j = min(
            (average(clusters[i], clusters[j]), i, j)
            for i in range(len(clusters)) for j in range(i + 1, len(clusters))
        )
        left, right = clusters[i], clusters[j]
        merged = Tree(left.members + right.members, distance / 2, left, right)
        clusters = [cluster for k, cluster in enumerate(clusters) if k not in (i, j)] + [merged]
    return clusters[0]


def _consensus(profile: Sequence[str]) -> str:
    return "".join(Counter(base for base in column if base != "-").most_common(1)[0][0]
                   for column in zip(*profile))


def _merge_profiles(first: Sequence[str], second: Sequence[str]) -> List[str]:
    first_consensus, second_consensus = align(_consensus(first), _consensus(second))[:2]
    merged = [[] for _ in range(len(first) + len(second))]
    first_column = second_column = 0
    for base_a, base_b in zip(first_consensus, second_consensus):
        for index, sequence in enumerate(first):
            merged[index].append(sequence[first_column] if base_a != "-" else "-")
        for index, sequence in enumerate(second, len(first)):
            merged[index].append(sequence[second_column] if base_b != "-" else "-")
        first_column += base_a != "-"
        second_column += base_b != "-"
    return ["".join(row) for row in merged]


def progressive_align(sequences: Sequence[str]) -> Tuple[List[str], Tree]:
    """Align profiles along the UPGMA tree, preserving every original base."""
    if not sequences or any(not sequence or set(sequence) - set(DNA) for sequence in sequences):
        raise ValueError("sequences must be nonempty A/C/G/T strings")
    tree = guide_tree(sequences)

    def visit(node: Tree) -> List[str]:
        if node.left is None:
            return [sequences[node.members[0]]]
        return _merge_profiles(visit(node.left), visit(node.right))

    rows = visit(tree)
    by_index = dict(zip(tree.members, rows))
    return [by_index[index] for index in range(len(sequences))], tree


def _probabilities(counts: Counter, alphabet: Sequence[str], pseudocount: float = 0) -> Dict[str, float]:
    total = sum(counts.values()) + pseudocount * len(alphabet)
    return {item: (counts[item] + pseudocount) / total for item in alphabet}


def build_profile(alignment: Sequence[str], threshold: float = 0.7) -> Dict:
    """Estimate match/insert emissions and observed HMM state transitions."""
    if not 0 < threshold <= 1:
        raise ValueError("threshold must be in (0, 1]")
    if not alignment or not alignment[0] or len({len(row) for row in alignment}) != 1:
        raise ValueError("alignment must contain equally sized, nonempty rows")
    if any(set(row) - set(DNA + "-") for row in alignment):
        raise ValueError("alignment contains invalid symbols")
    columns = list(zip(*alignment))
    match_columns = [index for index, column in enumerate(columns)
                     if sum(base != "-" for base in column) / len(alignment) >= threshold]
    if not match_columns:
        raise ValueError("alignment has no match columns at this threshold")
    match_index = {column: number for number, column in enumerate(match_columns, 1)}
    emissions = {"M{}".format(i): Counter() for i in range(1, len(match_columns) + 1)}
    emissions.update({"I{}".format(i): Counter() for i in range(len(match_columns) + 1)})
    transitions = defaultdict(Counter)
    for row in alignment:
        previous, position = "B", 0
        for index, base in enumerate(row):
            if index in match_index:
                position = match_index[index]
                state = "M{}".format(position) if base != "-" else "D{}".format(position)
            elif base != "-":
                state = "I{}".format(position)
            else:
                continue
            transitions[previous][state] += 1
            if base != "-":
                emissions[state][base] += 1
            previous = state
        transitions[previous]["E"] += 1
    return {
        "match_columns_1_based": [index + 1 for index in match_columns],
        "match_occupancy": [sum(base != "-" for base in columns[index]) / len(alignment)
                            for index in match_columns],
        "emissions": {state: _probabilities(counts, DNA, 0.5)
                      for state, counts in emissions.items()},
        "transitions": {state: _probabilities(counts, sorted(counts))
                        for state, counts in sorted(transitions.items())},
    }
