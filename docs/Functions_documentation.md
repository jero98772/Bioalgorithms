# Bioinformatics-toolspy Code Documentation

## Introduction
This document provides documentation for the bioinformatics code snippets provided. These snippets cover various bioinformatics algorithms and functions for sequence analysis, alignment, and motif finding.

## Functions

### `levenshtein_distance(string1, string2)`
Calculates the Levenshtein distance between two input strings.

**Example:**
```python
import bioinformatics
print(bioinformatics.levenshtein_distance("kitten","sitting"))  # Output: 3
```

### `sequence_alignment(sequence1, sequence2, gap_penalty, mismatch_penalty)`
Performs sequence alignment between two sequences with specified gap and mismatch penalties.

**Example:**
```python
import bioinformatics
print(bioinformatics.sequence_alignment("AGGGCT", "AGGCA", 3, 2))  # Output: (3, 'AGGGCT-', '-AGGCA')
```

### `longest_common_subsequences(sequence_list)`
Finds the longest common subsequence among a list of input sequences.

**Example:**
```python
import bioinformatics
print(bioinformatics.longest_common_subsequences(["ACCGAAGG","ACCGAACC","CCACCGAAGG","GGACCGAACC"]))  # Output: 'ACCGA'
```

### `longest_common_subsequence(sequence1, sequence2)`
Finds the longest common subsequence between two input sequences.

**Example:**
```python
import bioinformatics
print(bioinformatics.longest_common_subsequence("ACCGAAGG", "ACCGAACC"))  # Output: 'ACCGAA'
```

### `commun_patterns(pattern_list)`
Finds common patterns among a list of input patterns.

**Example:**
```python
import bioinformatics
print(bioinformatics.commun_patterns(["XAaXV","XAsXV","XAcXV"]))  # Output: 'XAXV'
```

### `reconstruct_from_kmers(k, kmers)`
Reconstructs a string from a collection of k-mers.

**Example:**
```python
import bioinformatics
print(bioinformatics.reconstruct_from_kmers(3,["AAT","ATG", "TGC", "GCT", "CTA"]))  # Output: 'AATGCTA'
```

### `translate_rna_to_aminoacid(rna_sequence)`
Translates an RNA sequence into an amino acid sequence.

**Example:**
```python
import bioinformatics
print(bioinformatics.translate_rna_to_aminoacid("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"))  # Output: 'MAAMGSST*'
```

### `de_bruijn_collection(kmers, prefix_func, suffix_func)`
Constructs a de Bruijn graph from a collection of k-mers.

**Example:**
```python
import bioinformatics
print(bioinformatics.de_bruijn_collection(["ATG", "ATG", "TGT", "TGG", "CAT", "GGA", "GAT", "AGA"], lambda kmer: kmer[:-1], lambda kmer: kmer[1:]))  # Output: {'AT': ['TG', 'TG', 'TG'], 'TG': ['GT', 'GG'], 'GT': ['TG'], 'GG': ['GA'], 'GA': ['AT'], 'CA': ['AT']}
```

### `find_eulerian_cycle(graph)`
Finds an Eulerian cycle in a graph.

**Example:**
```python
import bioinformatics
print(bioinformatics.find_eulerian_cycle({"AAT": ["ATG"],"ATG": ["TGC"],"TGC": ["GCT"],"GCT": ["CTA"],"CTA": ["TAC"],"TAC": ["ACG"],"ACG": ["CGA"],"CGA": ["GAT"],"GAT": ["ATG"]}))  # Output: ['AAT', 'ATG', 'TGC', 'GCT', 'CTA', 'TAC', 'ACG', 'CGA', 'GAT', 'ATG']
```

### `de_bruijn(k, sequence)`
Constructs a de Bruijn graph from a sequence.

**Example:**
```python
import bioinformatics
print(bioinformatics.de_bruijn(3, "ATGATCAAG"))  # Output: {'ATG': ['TGA'], 'TGA': ['GAT'], 'GAT': ['ATC'], 'ATC': ['TCA'], 'TCA': ['CAA', 'CAG'], 'CAA': ['AAG'], 'AAG': ['AG']}
```

### `grph_kmers(kmers)`
Constructs a graph from a collection of k-mers.

**Example:**
```python
import bioinformatics
print(bioinformatics.grph_kmers(["ACCGA", "CCGAA", "CGAAG", "GAAGC", "AAGCT"]))  # Output: {'ACCGA': ['CCGAA'], 'CCGAA': ['CGAAG'], 'CGAAG': ['GAAGC'], 'GAAGC': ['AAGCT']}
```

### `reconstruct(kmers)`
Reconstructs a string from a collection of overlapping k-mers.

**Example:**
```python
import bioinformatics
print(bioinformatics.reconstruct(["ACCGA", "CCGAA", "CGAAG", "GAAGC", "AAGCT"]))  # Output: 'ACCGAAGCT'
```

### `kmer_composition(k, sequence)`
Finds the k-mer composition of a sequence.

**Example:**
```python
import bioinformatics
print(bioinformatics.kmer_composition(2,'CAATCCAAC'))  # Output: ['CA', 'AA', 'AT', 'TT', 'TC', 'CC', 'CA', 'AA', 'AC']
```

### `distance_between_pattern_and_strings(pattern, string_list)`
Calculates the total Hamming distance between a pattern and a list of strings.

**Example:**
```python
import bioinformatics
print(bioinformatics.distance_between_pattern_and_strings("AA",['AAATTGACGCAT','GACGAAAAACGTT','CGTCAGCGCCTG''GCTGAGCAAAGG','AGTACGGGACAG']))  # Output: 14
```

### `gibbs(k, t, n, dna, iterations)`
Performs Gibbs sampling for motif finding.

**Example:**
```python
import bioinformatics
print(bioinformatics.gibbs(4, 5, 10,["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"],1000))  # Output: ['TCAG', 'TCAG', 'TCAG', 'TCAG', 'TCAG']
```

### `randomized_motif_search(k, t, dna, iterations)`
Performs randomized motif search for motif finding.

**Example:**
```python
import bioinformatics
print(bioinformatics.randomized

_motif_search(3,5,["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"],1))  # Output: ['TTC', 'TTC', 'TTC', 'TTC', 'TTC']
```

### `greedy_motif_search(k, t, dna)`
Performs greedy motif search for motif finding.

**Example:**
```python
import bioinformatics
print(bioinformatics.greedy_motif_search(3,5,["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]))  # Output: ['TTC', 'TTC', 'TTC', 'TTC', 'TTC']
```

### `most_probable(dna_string, k, profile_matrix)`
Finds the most probable k-mer in a DNA string given a profile matrix.

**Example:**
```python
import bioinformatics
print(bioinformatics.most_probable("ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT",5,[[0.2, 0.2, 0.3, 0.2, 0.3],[0.4, 0.3, 0.1, 0.5, 0.1],[0.3, 0.3, 0.5, 0.2, 0.4],[0.1, 0.2, 0.1, 0.1, 0.2]]))  # Output: 'CCGAG'
```

### `median_string(dna, k)`
Finds the median string in a collection of DNA strings.

**Example:**
```python
import bioinformatics
print(bioinformatics.median_string(['AAATTGACGCAT','GACGACCACGTT','CGTCAGCGCCTG''GCTGAGCACCGG','AGTACGGGACAG'],6))  # Output: 'GACGAC'
```

### `enumerate_motifs(dna, k, d)`
Enumerates all motifs of length k with at most d mismatches in a collection of DNA strings.

**Example:**
```python
import bioinformatics
print(bioinformatics.enumerate_motifs(["ATTTGGC","TGCCTTA","CGGTATC","GAAAATT"],3,1))  # Output: ['ATA', 'ATT', 'GTT', 'TTA', 'TTG']
```

### `pattern_to_number(pattern)`
Converts a DNA pattern to its corresponding integer.

**Example:**
```python
import bioinformatics
print(bioinformatics.pattern_to_number("CC"))  # Output: 13
```

### `number_to_pattern(number, k)`
Converts an integer to its corresponding DNA pattern of length k.

**Example:**
```python
import bioinformatics
print(bioinformatics.number_to_pattern(5,2))  # Output: 'GG'
```

### `generate_frequency_array(dna_string, k)`
Generates the frequency array of k-mers in a DNA string.

**Example:**
```python
import bioinformatics
print(bioinformatics.generate_frequency_array("AAACAGATCACCCGCTGAGCGGGTTATCTGTT",1))  # Output: [5, 4, 2, 1, 0, 0, 2, 2, 1, 2]
```

### `reverse_complement(dna_string)`
Finds the reverse complement of a DNA string.

**Example:**
```python
import bioinformatics
print(bioinformatics.reverse_complement("AAAAAGCATAAACATTAAAGAG"))  # Output: 'CTCTTTAATGTTTATGCTTTTT'
```

### `frequent_words_mismatch(dna_string, k, d)`
Finds the most frequent k-mers with at most d mismatches in a DNA string.

**Example:**
```python
import bioinformatics
print(bioinformatics.frequent_words_mismatch("ACGTTGCATGTCGCATGATGCATGAGAGCT",4,1))  # Output: ['ATGT', 'GATG', 'ATGC']
```

### `approximate_pattern_matching(pattern, text, d)`
Finds all approximate occurrences of a pattern in a text with at most d mismatches.

**Example:**
```python
import bioinformatics
print(bioinformatics.approximate_pattern_matching("AAAAAGCATAAACATTAAAGAG","AAAAA",0))  # Output: [0, 1, 2, 3, 4, 5, 6, 21]
```

### `approximate_pattern_count(pattern, text, d)`
Counts the number of approximate occurrences of a pattern in a text with at most d mismatches.

**Example:**
```python
import bioinformatics
print(bioinformatics.approximate_pattern_count("AAAAAGCATAAACATTAAAGAG","AAAAA",0))  # Output: 8
```

### `min_skew(genome)`
Finds the positions in a genome where the skew diagram attains its minimum value.

**Example:**
```python
import bioinformatics
print(bioinformatics.min_skew("CATGGGCATCGGCCATACGCC"))  # Output: [1, 6]
```

### `clump_finding(genome, k, L, t)`
Finds patterns forming clumps in a genome.

**Example:**
```python
import bioinformatics
print(bioinformatics.clump_finding("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAG",5,50,4))  # Output: ['CGACA']
```

### `pattern_count(text, pattern)`
Counts the occurrences of a pattern in a text.

**Example:**
```python
import bioinformatics
print(bioinformatics.pattern_count("cgatatatccatag","ata"))  # Output: 3
```

### `frequent_words(text, k)`
Finds the most frequent k-mers in a text.

**Example:**
```python
import bioinformatics
print(bioinformatics.frequent_words("actgactcccaccccc",3))  # Output: ['ccc']
```

### `pattern_count_positions(text, pattern)`
Finds the positions of all occurrences of a pattern in a text.

**Example:**
```python
import bioinformatics
print(bioinformatics.pattern_count_positions("cgatatatccatag","ata"))  # Output: [2, 4, 10]
```

### `hamming_distances(pattern, sequence)`
Calculates the Hamming distances between a pattern and a sequence.

**Example:**
```python
import bioinformatics
print(bioinformatics.hamming_distances("cgatatatccatag","ata"))  # Output: [3, 2, 1, 2, 2, 1, 2, 1, 0, 1, 2, 3, 1, 2]
```

