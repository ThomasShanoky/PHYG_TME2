

def kmer2str(val, k): #fonction test
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = [] #stocker les caractère du k-mer
    for _ in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2 

    str_val.reverse()
    return "".join(str_val)


def create_index(list_of_seq, k):
    """Create a dictionary index of k-mers from a sequence
    :param list list_of_seq: The list of sequence(s) to index
    :param int k: The number of nucleotides involved into the kmer.
    :return dict: A dictionary index of k-mers"""
    index = {}
    mask = (1 << (2 * k)) - 1  
    kmer = encode_kmer(list_of_seq[0][:k], k)  
    kmer_rc = reverse_complement(kmer, k)  
    kmer_to_add = min(kmer, kmer_rc)  
    index[kmer_to_add] = 1
    for seq in list_of_seq:
        for i in range(1, len(seq) - k + 1):
            kmer = ((kmer << 2) & mask) + encode_nuc(seq[i + k - 1])
            kmer_rc = reverse_complement(kmer, k) 
            kmer_to_add = min(kmer, kmer_rc)
            if kmer_to_add in index:
                index[kmer_to_add]+=1
            else:
                index[kmer_to_add] = 1
    return index


# Fonction pour encoder un nucléotide en entier selon l'ordre 'A', 'C', 'T', 'G'
def encode_nuc(letter):
    """Encode un nucléotide en entier"""
    encoding = {'A': 0, 'C': 1, 'T': 2, 'G': 3}  # Ordre ACTG
    return encoding[letter]


# Fonction pour encoder un k-mer sous forme d'entier
def encode_kmer(seq, k):
    """Encode un k-mer en entier."""
    kmer = 0
    for letter in seq[:k]:
        kmer <<= 2
        kmer += encode_nuc(letter)
    return kmer


#Fonction pour calculer le complément inverse d'un k-mer encodé
def reverse_complement(kmer, k):
    complement = 0
    for _ in range(k):
        # Extraire les 2 bits du nucléotide actuel
        nuc = kmer & 0b11
        # Appliquer le complément dans l'ordre 'A', 'C', 'T', 'G': 
        # 00 (A) <-> 10 (T), 01 (C) <-> 11 (G)
        if nuc == 0:  # A
            complement = (complement << 2) | encode_nuc('T')  # A -> T
        elif nuc == 1:  # C
            complement = (complement << 2) | encode_nuc('G')  # C -> G
        elif nuc == 2:  # T
            complement = (complement << 2) | encode_nuc('A')  # T -> A
        elif nuc == 3:  # G
            complement = (complement << 2) | encode_nuc('C')  # G -> C
        kmer >>= 2  # Décaler le k-mer pour le prochain nucléotide
    return complement


# Fonction pour générer les k-mers canoniques d'une séquence
def stream_kmers(list_of_seq, k):
    """Génère les k-mers canoniques d'une séquence."""
    mask = (1 << (2 * k)) - 1  # Masque pour limiter les bits du k-mer à k
    kmer = encode_kmer(list_of_seq[0][:k], k)  # Encodage du premier k-mer
    kmer_rc = reverse_complement(kmer, k)  # Complément inverse du premier k-mer
    yield min(kmer, kmer_rc)  # Générer le plus petit k-mer canonique

    for seq in list_of_seq:
        for i in range(1, len(seq) - k + 1):
            # Décale à gauche et ajoute le nouveau nucléotide encodé
            kmer = ((kmer << 2) & mask) + encode_nuc(seq[i + k - 1])
            kmer_rc = reverse_complement(kmer, k)  # Complément inverse du k-mer
            # Générer le plus petit entre le k-mer et son complément inverse
            yield min(kmer, kmer_rc)