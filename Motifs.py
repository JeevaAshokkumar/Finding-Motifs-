# Copy your Score(Motifs), Count(Motifs), Profile(Motifs), and Consensus(Motifs) functions here.
def Count(Motifs):
    count = {} # initializing the count dictionary
    # your code here
    # Set k equal to the length of Motifs[0]
    k = len(Motifs[0])
    # Iterate over all nucleotides symbol and create a list of zeroes corresponding to count[symbol]
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    # Iterate over all elements symbol = Motifs[i][j] of the count matrix and add 1 to count[symbol][j]
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count
# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile={}
    for symbol in "ACGT":
        profile[symbol] = []
        for j in range(k):
             profile[symbol].append(0)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            profile[symbol][j] += 1
    for symbol in "ACGT":
        for j in range(k):
            profile[symbol][j]=profile[symbol][j]/t
    return profile
def Consensus(Motifs):
    # insert your code here
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus
# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score(Motifs):
    # Insert code here
    consensus = Consensus(Motifs)  # Get the consensus string using Consensus function
    k = len(consensus)
    t = len(Motifs)
    score = 0
    # Sum the number of symbols in the j-th column of Motifs that do not match the symbol in position j of the consensus string
    for j in range(k):
        for i in range(t):
            if Motifs[i][j] != consensus[j]:
                score += 1
    return score
# Then copy your ProfileMostProbableKmer(Text, k, Profile) and Pr(Text, Profile) functions here.
def Pr(Text, Profile):
    # insert your code here
    p = 1.0
    for i in range(len(Text)):
        symbol = Text[i]
        column = i
        # Multiply p by the value of Profile corresponding to symbol Text[i] and column i
        p *= Profile[symbol][column]
    return p
# Write your ProfileMostProbableKmer() function here.
# The profile matrix assumes that the first row corresponds to A, the second corresponds to C,
# the third corresponds to G, and the fourth corresponds to T.
# You should represent the profile matrix as a dictionary whose keys are 'A', 'C', 'G', and 'T' and whose values are lists of floats
def ProfileMostProbablePattern(text, k, profile):
    max_probability = -1.0  # Initialize max probability as a negative value
    most_probable_kmer = ""
    for i in range(len(text) - k + 1):
        kmer = text[i:i+k]
        probability = Pr(kmer, profile)  # Use the Pr function to compute probability
        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer
    return most_probable_kmer
# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearch(Dna, k, t):
    # type your GreedyMotifSearch code here.
    n = len(Dna[0])
    best_motifs = [Dna[i][0:k] for i in range(t)]
    for i in range(n - k + 1):
        motifs = [Dna[0][i:i+k]]
        for j in range(1, t):
            profile = Profile(motifs[0:j])
            motifs.append(ProfileMostProbablePattern(Dna[j], k, profile))
        if Score(motifs) < Score(best_motifs):
            best_motifs = motifs
    return best_motifs 
def CountWithPseudocounts(Motifs):
    count = {symbol: [1] * len(Motifs[0]) for symbol in "ACGT"}
    for i in range(len(Motifs)):
        for j in range(len(Motifs[i])):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs) + 4
    count = CountWithPseudocounts(Motifs)
    profile = {symbol: [float(c) / t for c in count[symbol]] for symbol in "ACGT"}
    return profile
def ScoreWithPseudocounts(Motifs):
    count = CountWithPseudocounts(Motifs)
    score = 0
    for j in range(len(Motifs[0])):
        max_count = max(count[symbol][j] for symbol in "ACGT")
        score += sum(count[symbol][j] for symbol in "ACGT") - max_count
    return score
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = [string[:k] for string in Dna]    
    for i in range(len(Dna[0]) - k + 1):
        Motifs = [Dna[0][i:i+k]]        
        for j in range(1, t):
            profile = ProfileWithPseudocounts(Motifs)
            motif = ProfileMostProbablePattern(Dna[j], k, profile)
            Motifs.append(motif)   
        if ScoreWithPseudocounts(Motifs) < ScoreWithPseudocounts(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
def Motifs(Profile, Dna):
    k = len(Profile["A"])
    most_probable_kmers = []
    for text in Dna:
        most_probable = ProfileMostProbablePattern(text, k, Profile)
        most_probable_kmers.append(most_probable)
    return most_probable_kmers
# For example data set
Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC",
       "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG",
       "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC",
       "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC",
       "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG",
       "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA",
       "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA",
       "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG",
       "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG",
       "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]
# set t equal to the number of strings in Dna and k equal to 15
t = len(Dna)
k = 15
# Call GreedyMotifSearchWithPseudocounts(Dna, k, t) and store the output in a variable called Motifs
motifs = GreedyMotifSearchWithPseudocounts(Dna, k, t)
# Print the Motifs variable
#print(motifs)
# Print Score(Motifs)
#print(Score(motifs))
profile=Profile(motifs)
#print(Motifs(profile, Dna))
import random
# Input:  A list of strings Dna, and integers k and t
# Output: RandomMotifs(Dna, k, t)
# HINT:   You might not actually need to use t since t = len(Dna), but you may find it convenient
def RandomMotifs(Dna, k, t):
    random_motifs = []
    t = len(Dna)
    l = len(Dna[0])
    for i in range(t):
        ran_num = random.randint(0,l-k)
        random_motifs.append(Dna[i][ran_num:ran_num+k])
    return random_motifs
#print(RandomMotifs(Dna, k, t))
def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs
#print(RandomizedMotifSearch(Dna, k, t))
# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
def Normalize(Probabilities):
    # your code here
    total_probability = sum(Probabilities.values())  # Calculate the sum of all probabilities
    # Normalize each probability by dividing by the total probability
    normalized_probabilities = {key: value / total_probability for key, value in Probabilities.items()}        
    return normalized_probabilities
# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    import random
    # Generate a random number between 0 and 1
    p = random.uniform(0, 1)
    # Iterate through k-mers and their probabilities
    for kmer, probability in Probabilities.items():
        p -= probability
        if p <= 0:
            return kmer
# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    # your code here
    n=len(Text)
    probabilities={}
    for i in range(n-k+1):
        probabilities[Text[i:i+k]]=Pr(Text[i:i+k],profile)
    probabilities=Normalize(probabilities)
    return WeightedDie(probabilities)
def GibbsSampler(Dna, k, t, N):
    Motifs=RandomMotifs(Dna, k, t)
    BestMotifs=Motifs
    for j in range(N):
        i=random.randint(1,t)
        text=Motifs.pop(i-1)
        profile=ProfileWithPseudocounts(Motifs)
        kmer=ProfileMostProbableKmer(text, k, profile)
        Motifs.insert(i-1,kmer)
        if Score(Motifs)<Score(BestMotifs):
            BestMotifs=Motifs
    return BestMotifs

