import bioinformatics

def custom_head(kmer):
    return kmer[:-1]
def custom_tail(kmer):
    return kmer[1:]
print(bioinformatics.de_bruijn_collection(["ATG", "ATG", "TGT", "TGG", "CAT", "GGA", "GAT", "AGA"], custom_head, custom_tail))
print(bioinformatics.find_eulerian_cycle({"AAT": ["ATG"],"ATG": ["TGC"],"TGC": ["GCT"],"GCT": ["CTA"],"CTA": ["TAC"],"TAC": ["ACG"],"ACG": ["CGA"],"CGA": ["GAT"],"GAT": ["ATG"]}))
print(bioinformatics.de_bruijn(3, "ATGATCAAG"))
print(bioinformatics.grph_kmers(["ACCGA", "CCGAA", "CGAAG", "GAAGC", "AAGCT"]))
print(bioinformatics.reconstruct(["ACCGA", "CCGAA", "CGAAG", "GAAGC", "AAGCT"]))
print(bioinformatics.kmer_composition(2,'CAATCCAAC'))
print(bioinformatics.distance_between_pattern_and_strings("AA",['AAATTGACGCAT','GACGAAAAACGTT','CGTCAGCGCCTG''GCTGAGCAAAGG','AGTACGGGACAG']))
print(bioinformatics.gibbs(4, 5, 10,["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"],1000))
print(bioinformatics.randomized_motif_search(3,5,["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"],1))
print(bioinformatics.greedy_motif_search(3,5,["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]))
print(bioinformatics.most_probable("ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT",5,[[0.2, 0.2, 0.3, 0.2, 0.3],[0.4, 0.3, 0.1, 0.5, 0.1],[0.3, 0.3, 0.5, 0.2, 0.4],[0.1, 0.2, 0.1, 0.1, 0.2]]))
print(bioinformatics.median_string(['AAATTGACGCAT','GACGACCACGTT','CGTCAGCGCCTG''GCTGAGCACCGG','AGTACGGGACAG'],6))
print(bioinformatics.enumerate_motifs(["ATTTGGC","TGCCTTA","CGGTATC","GAAAATT"],3,1))
print(bioinformatics.pattern_to_number("CC"))
print(bioinformatics.number_to_pattern(5,2))
print(bioinformatics.generate_frequency_array("AAACAGATCACCCGCTGAGCGGGTTATCTGTT",1))
print(bioinformatics.reverse_complement("AAAAAGCATAAACATTAAAGAG"))
print(bioinformatics.frequent_words_mismatch("ACGTTGCATGTCGCATGATGCATGAGAGCT",4,1))
print(bioinformatics.approximate_pattern_matching("AAAAAGCATAAACATTAAAGAG","AAAAA",0))
print(bioinformatics.approximate_pattern_count("AAAAAGCATAAACATTAAAGAG","AAAAA",0))
print(bioinformatics.min_skew("CATGGGCATCGGCCATACGCC"))
print(bioinformatics.clump_finding("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAG",5,50,4))
print(bioinformatics.pattern_count("cgatatatccatag","ata"))
print(bioinformatics.frequent_words("actgactcccaccccc",3))
print(bioinformatics.pattern_count_positions("cgatatatccatag","ata"))
print(bioinformatics.hamming_distances("cgatatatccatag","ata"))

