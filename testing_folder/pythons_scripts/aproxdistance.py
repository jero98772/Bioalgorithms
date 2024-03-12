import numpy as np

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
#print(result)  # Output: 4



def approx(kmer, substring, k, n):
    count = 0
    for i in range(k):
        if substring[i] != kmer[i]:
            count += 1
        if count > n:
            return False
    return True

def Frequent_Words_Mismatch(DNA,k,n):
    counts = {}
    for i in range(0, len(DNA) - k + 1):
        kmer = DNA[i:i + k]
        if kmer in counts:
            counts[kmer] += 1
        else:
            counts[kmer] = 1
    updateCounts = {}
    for a in counts:
        c = 0
        for b in counts:
            if approx(a, b, k, n):
                c += counts.get(b)
        updateCounts[a] = c
    #print(updateCounts)
    frequent = max(updateCounts.values())
    ans = []
    for k in updateCounts:
        if updateCounts[k] == frequent:
            ans.append(k)

    output = ' '.join(ans)
    #print(" ")
    #print(output)

#print(Frequent_Words_Mismatch("ACGTTGCATGTCGCATGATGCATGAGAGCT",4,1))





from collections import defaultdict

def HammingDistance(s1, s2):
    d = sum([1 for i in range(len(s1)) if s1[i]!=s2[i]])
    return d

def ReversePattern(pattern):
    chain="ATGATCAAG"
    newchain=[]
    for i in range(len(chain)):
        char=chain[len(chain)-i-1]
        if char.lower()=="a":newchain.append("t")
        elif char.lower()=="t":newchain.append("a")
        elif char.lower()=="g":newchain.append("c")
        elif char.lower()=="c":newchain.append("g")
    return ' '.join(str(i) for i in newchain)

def neighbour(pattern, mismatch, words):
    if mismatch == 0:
        words.add(pattern)
    else:
        bases = ['A', 'T', 'C', 'G']
        for i in range(len(pattern)):
            for j in range(len(bases)):
                new_pattern = pattern[:i] + bases[j] + pattern[i+1:]
                if mismatch <= 1:
                    words.add(new_pattern)
                else:
                    neighbour(new_pattern, mismatch-1, words)

def FindMostFrequentPattern(text, k, d):
    allfrequentwords = defaultdict(int)
    for i in range(len(text) - k + 1):
        frequentwords = set()
        neighbour(text[i:i + k], d, frequentwords)

        for words in frequentwords:
            allfrequentwords[words] += 1

    for t in allfrequentwords.keys():
        reverse_k = ReversePattern(t)
        for i in range(len(text) - k + 1):
            if HammingDistance(text[i:i + k], reverse_k) <= d:
                allfrequentwords[t] += 1

    result = set()
    for t in allfrequentwords.keys():
        if allfrequentwords[t] == max(allfrequentwords.values()):
            result.add(t)
            result.add(ReversePattern(t))
    for i in result:
        print(i)
        

def patternToNumber(kmer):
    '''BA1L Implement PatternToNumber'''
    n=0
    for letter in kmer:
        n*=4
        n+=bases.find(letter)
    return n
def generateFrequencyArray(text,k):

    frequencies=np.zeros((4**k))
    for i in range(len(text)-k+1):
        frequencies[patternToNumber(text[i:i+k])] += 1
    return frequencies

#faster frecuency array

if __name__ == "__main__":
    # text = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
    k, d = 4, 1
    text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
    FindMostFrequentPattern(text, k, d)

