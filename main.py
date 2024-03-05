import bioinformatics

print(bioinformatics.approximate_pattern_matching("AAAAAGCATAAACATTAAAGAG","AAAAA",0))
print(bioinformatics.approximate_pattern_count("AAAAAGCATAAACATTAAAGAG","AAAAA",0))
print(bioinformatics.min_skew("CATGGGCATCGGCCATACGCC"))
print(bioinformatics.clump_finding("CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAG",5,50,4))
print(bioinformatics.pattern_count("cgatatatccatag","ata"))
print(bioinformatics.frequent_words("actgactcccaccccc",3))
print(bioinformatics.pattern_count_positions("cgatatatccatag","ata"))
print(bioinformatics.hamming_distances("cgatatatccatag","ata"))

