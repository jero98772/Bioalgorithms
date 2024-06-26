# Bioalgorithms
<p align="center"><img src="https://raw.githubusercontent.com/jero98772/Bioalgorithms/main/docs/logo.jpeg" width="500" height="500"></p>

## Overview
Bioalgorithms is a powerful library tailored for Python users but built in Rust, harnessing the efficiency and speed of compiled languages. With over 30 functions crafted for both bioinformatics and general programming tasks, this library integrates dynamic programming, greedy techniques, and optimized data structures such as hashmaps.

Beyond Python, Bioalgorithms offers compatibility with Rust, empowering developers to leverage its capabilities across different programming paradigms.

## Features

- **Efficiency**: Utilizes Rust's compiled nature for enhanced performance compared to interpreted languages like Python and JavaScript.
- **Comprehensive Functionality**: Offers a diverse set of functions optimized for various bioinformatics and programming tasks.
- **Algorithmic Techniques**: Implements advanced algorithms including dynamic programming and greedy strategies.
- **Data Structures**: Leverages efficient data structures such as hashmaps to streamline problem-solving.

## Installation

### Prerequisites

- Git
- Python 3.x
- Rust (for compiling Rust code)
- `maturin` (for Rust bindings)

### Installation Steps

      pip install bioalgorithms

or 

1. **Clone Repository:**
   ```
   git clone https://github.com/jero98772/Bioalgorithms.git
   ```

2. **Set Up Virtual Environment:**
   ```
   cd Bioinformatics-toolspy
   python -m venv env
   source env/bin/activate
   pip install maturin
   ```

3. **Compile Rust Code:**
   ```
   maturin develop
   ```

4. **Integration with Your Code:**
   Place your Python file temporarily inside the folder and execute:
   ```
   python <your_file>
   ```

## Code Documentation

For comprehensive documentation and usage guidelines, please refer to the [Bioalgorithms Documentation](https://github.com/jero98772/Bioalgorithms/blob/main/docs/Functions_documentation.md). This documentation covers installation instructions, detailed explanations of available functions, code examples, and more. Whether you're a beginner exploring bioinformatics algorithms or an experienced developer seeking to optimize your workflows, the documentation serves as a valuable resource to unleash the full potential of Bioalgorithms.

Example of a code in python using all functions:

      import bioalgorithms
      print(bioalgorithms.levenshtein_distance("kitten","sitting"))
      print(bioalgorithms.sequence_aligment("AGGGCT","AGGCA",3,2))
      print(bioalgorithms.longest_commons_subsequences(["ACCGAAGG","ACCGAACC","CCACCGAAGG","GGACCGAACC"]))
      print(bioalgorithms.longest_common_subsequence("ACCGAAGG","ACCGAACC"))
      print(bioalgorithms.commun_patters(["XAaXV","XAsXV","XAcXV"]))
      print(bioalgorithms.reconstruct_from_kmers(3,["AAT","ATG", "TGC", "GCT", "CTA"]))#check this output
      print(bioalgorithms.translate_rna_to_aminoacid("AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"))
      print(bioalgorithms.de_bruijn_collection(["ATG", "ATG", "TGT", "TGG", "CAT", "GGA", "GAT", "AGA"], lambda kmer: kmer[:-1], lambda kmer: kmer[1:]))
      print(bioalgorithms.find_eulerian_cycle({"AAT": ["ATG"],"ATG": ["TGC"],"TGC": ["GCT"],"GCT": ["CTA"],"CTA": ["TAC"],"TAC": ["ACG"],"ACG": ["CGA"],"CGA": ["GAT"],"GAT": ["ATG"]}))
      print(bioalgorithms.de_bruijn(3, "ATGATCAAG"))
      print(bioalgorithms.grph_kmers(["ACCGA", "CCGAA", "CGAAG", "GAAGC", "AAGCT"]))
      print(bioalgorithms.reconstruct(["ACCGA", "CCGAA", "CGAAG", "GAAGC", "AAGCT"]))
      print(bioalgorithms.kmer_composition(2,'CAATCCAAC'))
      print(bioalgorithms.distance_between_pattern_and_strings("AA",['AAATTGACGCAT','GACGAAAAACGTT','CGTCAGCGCCTG''GCTGAGCAAAGG','AGTACGGGACAG']))
      print(bioalgorithms.gibbs(4, 5, 10,["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"],1000))
      print(bioalgorithms.randomized_motif_search(3,5,["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"],1))
      print(bioalgorithms.greedy_motif_search(3,5,["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]))
      print(bioalgorithms.most_probable("ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT",5,[[0.2, 0.2, 0.3, 0.2, 0.3],[0.4, 0.3, 0.1, 0.5, 0.1],[0.3, 0.3, 0.5, 0.2, 0.4],[0.1, 0.2, 0.1, 0.1, 0.2]]))
      print(bioalgorithms.median_string(['AAATTGACGCAT','GACGACCACGTT','CGTCAGCGCCTG''GCTGAGCACCGG','AGTACGGGACAG'],6))
      print(bioalgorithms.enumerate_motifs(["ATTTGGC","TGCCTTA","CGGTATC","GAAAATT"],3,1))
      print(bioalgorithms.pattern_to_number("CC"))
      print(bioalgorithms.number_to_pattern(5,2))
      print(bioalgorithms.generate_frequency_array("AAACAGATCACCCGCTGAGCGGGTTATCTGTT",1))
      print(bioalgorithms.reverse_complement("AAAAAGCATAAACATTAAAGAG"))
      print(bioalgorithms.frequent_words_mismatch("ACGTTGCATGTCGCATGATGCATGAGAGCT",4,1))
      print(bioalgorithms.approximate_pattern_matching("AAAAAGCATAAACATTAAAGAG","AAAAA",0))
      print(bioalgorithms.approximate_pattern_count("AAAAAGCATAAACATTAAAGAG","AAAAA",0))
      print(bioalgorithms.min_skew("CATGGGCATCGGCCATACGCC"))
      print(bioalgorithms.clump_finding("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAG",5,50,4))
      print(bioalgorithms.pattern_count("cgatatatccatag","ata"))
      print(bioalgorithms.frequent_words("actgactcccaccccc",3))
      print(bioalgorithms.pattern_count_positions("cgatatatccatag","ata"))
      print(bioalgorithms.hamming_distances("cgatatatccatag","ata"))



## Purpose

Bioalgorithms is a project aimed at enhancing bioinformatics algorithm implementations, focusing on optimizing Python with Rust and PyO3 bindings. It aims to provide a complement to existing libraries like [Biopython](https://github.com/biopython/biopython), featuring functions both in pure Python and Rust for improved performance.

## Contribution

Contributions to this project are welcome. If you encounter any issues or have suggestions for improvement, please open an issue. For inquiries or assistance, contact [jero98772@protonmail.com](mailto:jero98772@protonmail.com).

While the library shows promise, its future development may be uncertain due to other commitments.


### Training material and documentation

book bioinformtics algorithms an active learning aproch solutions by Phillip Compeau and Pavel Pevzner

[geeksforgeeks](www.geeksforgeeks.org)

### References

[Biocomp BookClub ](https://github.com/juanjo255/Biocomp-BookClub)

[https://github.com/weka511/bioinformatics](https://github.com/weka511/bioinformatics)

### License

Bioalgorithms y is licensed under the GNU General Public License v3.0 (GPL-3.0). This license grants users the freedom to copy, modify, and distribute the software according to the terms outlined in the license agreement. If you choose to copy or redistribute the codebase, all you need to do is ensure that my name, jero98772, is included in the credits or acknowledgments. This requirement ensures proper attribution and helps maintain transparency in the open-source community.
