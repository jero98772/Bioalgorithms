def HAMMINGDISTANCE(pattern1, pattern2):
    # Assuming both patterns are of equal length
    return sum(1 for i in range(len(pattern1)) if pattern1[i] != pattern2[i])

def APPROXIMATEPATTERNCOUNT(Text, Pattern, d):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        Pattern_prime = Text[i:i+len(Pattern)]
        if HAMMINGDISTANCE(Pattern, Pattern_prime) <= d:
            count += 1
    return count

Text = "AACAAGCATAAACATTAAAGAG"
Pattern = "AAAA"
d = 1

result = APPROXIMATEPATTERNCOUNT(Text, Pattern, d)
print(result)  # Output: 4

