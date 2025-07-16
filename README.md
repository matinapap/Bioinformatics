# Bioinformatics Project – Sequence Alignment and HMM Profile Construction

## 📌 Description

This project implements various bioinformatics algorithms including synthetic DNA sequence generation, pairwise and multiple sequence alignment, phylogenetic tree construction using the Neighbor-Joining algorithm, and the construction of a Profile Hidden Markov Model (HMM). It was developed as part of a university course assignment in Spring 2023-2024.

## 🧬 Project Structure

### 1. **Synthetic DNA Sequence Generation**

- Uses the `random` and `sys` libraries.
- Random DNA strings are created using the nucleotides A, C, G, and T.
- Patterns are randomly modified and inserted into the string using substitution or insertion.
- A total of 50 sequences are generated and split into two datasets:
  - `datasetA.txt` (15 sequences)
  - `datasetB.txt` (35 sequences)
  - Full dataset saved in `FullDataset.txt`

### 2. **Global Alignment**

- Implements a global alignment algorithm (`global_alignment(A, B, alpha=2)`) that uses dynamic programming.
- Scores:
  - +1 for match
  - -α/2 for mismatch
  - Gap penalties included
- Returns aligned sequences and alignment score.

### 3. **Pairwise Distance Matrix**

- Calculates similarity scores between all sequence pairs in a dataset.
- Builds an n×n matrix storing alignment scores for each sequence pair.
- Symmetric matrix.

### 4. **Neighbor Joining Algorithm**

- Constructs a phylogenetic tree based on the distance matrix.
- Iteratively merges the most similar clusters until all sequences are grouped.

### 5. **Progressive Alignment**

- Builds a multiple sequence alignment using a guide tree from the Neighbor-Joining algorithm.
- Aligns the most similar sequences first, then adds others progressively.
- Results saved in `multiple_alignment_result.txt`.

### 6. **Profile HMM Construction**

#### Conserved Region Detection

- Columns are categorized into:
  - `match`: ≥70% consensus
  - `delete`: ≥70% consensus with gaps
  - `insert`: other cases

#### Emission Probabilities

- For each conserved (match) column, calculates the probability distribution of nucleotides (A, C, G, T).
- Stores in an emission probability table.

#### Transition Probabilities

- Builds a transition matrix that describes the probabilities of moving between:
  - Match → Match, Insert, Delete
  - Insert → Match, Insert, Delete
  - Delete → Match, Insert, Delete

### 7. **Viterbi Algorithm Preparation**

- Calculates emission probabilities for Match, Insert, and Delete states.
- These are used to determine the most probable path through the HMM for a given observed sequence.

---

## 🗂️ Output Files

- `FullDataset.txt`: All 50 generated sequences.
- `datasetA.txt` and `datasetB.txt`: Split sequence sets.
- `multiple_alignment_result.txt`: Multiple sequence alignment result.
- Emission and transition probability tables printed or optionally saved to files.

---

## 📄 Documentation

📘 [Project Documentation (PDF)]([[https://github.com/matinapap/Learning-Python-Android-App/blob/main/%CE%95%CE%BA%CF%80%CE%B1%CE%B9%CE%B4%CE%B5%CF%85%CF%84%CE%B9%CE%BA%CF%8C%20%CE%9B%CE%BF%CE%B3%CE%B9%CF%83%CE%BC%CE%B9%CE%BA%CF%8C_%CE%91%CE%BD%CE%B1%CF%86%CE%BF%CF%81%CE%AC.pdf](https://github.com/matinapap/Bioinformatics/blob/main/bioinformatics_doc.pdf)](https://github.com/matinapap/Bioinformatics/blob/main/bioinformatics_doc.pdf))

---

## 👥 Contributors

<table>
  <tr>
    <td align="center">
      <a href="https://github.com/matinapap">
        <img src="https://github.com/matinapap.png" width="100px;" alt="Matina Papadakou"/><br />
        <sub><b>Matina Papadakou</b></sub>
      </a>
      <p>P21127</p>
    </td>
    <td align="center">
      <a href="https://github.com/Georgechrp">
        <img src="https://github.com/Georgechrp.png" width="100px;" alt="George Christopoulos"/><br />
        <sub><b>Dimitrios Fotiadis</b></sub>
      </a>
      <p>P21183</p>
    </td>
  </tr>
</table>

## 🛠️ Requirements

- Python 3.x
- NumPy

Install dependencies (if needed):

```bash
pip install numpy
