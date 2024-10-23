import numpy as np

file_name = 'auxiliary2024/datasetA.txt'

try:
    with open(file_name, 'r') as file:
        # Read all lines in the file
        lines = file.readlines()
        
        # Strip whitespace characters and newline characters from each line, then store in datasetA
        datasetA = [line.strip() for line in lines]
        
except FileNotFoundError:
    print(f"The file {file_name} does not exist.")
except Exception as e:
    print(f"An error occurred: {e}")


def global_alignment(A, B, alpha=2):
    m, n = len(A), len(B)
    table = np.zeros((m + 1, n + 1))
    gap_penalty = -alpha
    match_score = 1
    mismatch_penalty = -alpha / 2

    # Initialization
    for i in range(1, m + 1):
        table[i][0] = i * gap_penalty
    for j in range(1, n + 1):
        table[0][j] = j * gap_penalty

    # Scoring
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = table[i - 1][j - 1] + (match_score if A[i - 1] == B[j - 1] else mismatch_penalty)
            delete = table[i - 1][j] + gap_penalty
            insert = table[i][j - 1] + gap_penalty
            table[i][j] = max(match, delete, insert)

    # Traceback
    align1, align2 = [], []
    i, j = m, n

    while i > 0 and j > 0:
        score_current = table[i][j]
        score_diagonal = table[i-1][j-1]
        score_left = table[i][j-1]
        score_up = table[i-1][j]

        max_score = max(score_diagonal, score_up, score_left)
        if max_score == score_diagonal:
            align1.append(A[i-1])
            align2.append(B[j-1])
            i -= 1
            j -= 1
        elif max_score == score_left:
            align1.append("-")
            align2.append(B[j-1])
            j -= 1
        elif max_score == score_up:
            align1.append(A[i-1])
            align2.append("-")
            i -= 1

    align1.reverse()
    align2.reverse()
    return ''.join(align1), ''.join(align2), int(table[m][n])


def pairwise_distance_matrix(sequences, alpha=2):
    n = len(sequences)
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            _, _, score = global_alignment(sequences[i], sequences[j], alpha)
            dist_matrix[i][j] = dist_matrix[j][i] = score 
    return dist_matrix


def neighbor_joining(dist_matrix):
    n = len(dist_matrix)
    clusters = [[i] for i in range(n)]
    while len(clusters) > 1:
        min_dist = float('inf')
        a, b = -1, -1
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                dist = -sum(dist_matrix[p][q] for p in clusters[i] for q in clusters[j]) / (len(clusters[i]) * len(clusters[j]))
                if dist < min_dist:
                    min_dist = dist
                    a, b = i, j
        new_cluster = clusters[a] + clusters[b]
        clusters.append(new_cluster)
        clusters.pop(max(a, b))
        clusters.pop(min(a, b))
    if len(clusters[0]) % 2 == 1:
        temp = clusters[0][0]
        clusters[0].pop(0)
        clusters[0].append(temp)

    return clusters[0]


def progressive_alignment(sequences):
    dist_matrix = pairwise_distance_matrix(sequences)
    guide_tree = neighbor_joining(dist_matrix)

    aligned_sequences = [sequences[guide_tree[0]], sequences[guide_tree[1]]]
    align1, align2, _ = global_alignment(aligned_sequences[0], aligned_sequences[1])
    aligned_sequences = [align1, align2]

    for idx in guide_tree[2:]:
        new_alignments = []
        for aligned_seq in aligned_sequences:
            align1, align2, _ = global_alignment(sequences[idx], aligned_seq)
            new_alignments.append(align2)
        new_alignments.insert(0, align1)
        aligned_sequences = new_alignments

    max_length = max(len(seq) for seq in aligned_sequences)
    aligned_sequences = [seq.ljust(max_length, '-') for seq in aligned_sequences]
    
    return aligned_sequences

result = progressive_alignment(datasetA)
print("Multiple Sequence Alignment Result:")
for aligned_seq in result:
    print(aligned_seq)

print('')
try:
    output_file = 'auxiliary2024/multiple_alignment_result.txt'
    with open(output_file, 'w') as f:
        for aligned_seq in result:
            f.write(aligned_seq + '\n')

    print(f"Results written to {output_file}")
except FileNotFoundError:
    print(f"The file {file_name} does not exist.")
except Exception as e:
    print(f"An error occurred: {e}")
