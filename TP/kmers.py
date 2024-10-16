def kmer2str(val, k):
    """ Transform a kmer integer into a its string representation
    :param int val: An integer representation of a kmer
    :param int k: The number of nucleotides involved into the kmer.
    :return str: The kmer string formatted
    """
    letters = ['A', 'C', 'T', 'G']
    str_val = []
    for _ in range(k):
        str_val.append(letters[val & 0b11])
        val >>= 2

    str_val.reverse()
    return "".join(str_val)


def encode_nucl(nucl):
    """ Encode a nucleotide into a 2-bit integer
    :param str nucl: The nucleotide to encode
    :return (int, int): The encoded nucleotide and its reverse complement
    """
    encoded = (ord(nucl) >> 1) & 0b11 # Extract the two bits of the ascii code that represent the nucleotide
    rencoded = (encoded + 2) & 0b11 # Complement encoding with bit tricks. Avoid slow if statement.

    return encoded, rencoded


def xorshift(val): #fonction hash : permet de mélanger la distribution
    val ^= val << 13
    val &= 0xFFFFFFFFFFFFFFFF
    val ^= val >> 7
    val ^= val << 17
    val &= 0xFFFFFFFFFFFFFFFF
    return val


def stream_kmers(seq, k):
    # Initialize the kmer and its reverse complement
    kmer = 0
    rkmer = 0

    # Add the first k-1 nucleotides to the first kmer and its reverse complement
    for i in range(k-1):
        nucl, rnucl = encode_nucl(seq[i])
        kmer |= nucl << (2*(k-2-i))
        rkmer |= rnucl << (2*(i+1))

    mask = (1 << (2*(k-1))) - 1

    # Yield the kmers
    for i in range(k-1, len(seq)):
        nucl, rnucl = encode_nucl(seq[i])
        # Remove the leftmost nucleotide from the kmer 
        kmer &= mask
        # Shift the kmer to make space for the new nucleotide
        kmer <<= 2
        # Add the new nucleotide to the kmer
        kmer |= nucl
        # Make space for the new nucleotide in the reverse kmer (remove the rightmost nucleotide by side effect)
        rkmer >>= 2
        # Add the new nucleotide to the reverse kmer
        rkmer |= rnucl << (2*(k-1))

        yield min(xorshift(kmer), xorshift(rkmer))


def filter_smallest(seq, k, s): #s : sketch size
    smallests = [float('inf') for _ in range(s)]
    max_el = float('inf')
    for kmer in stream_kmers(seq, k):
        if kmer < max_el:
            index = smallests.index(max_el)
            smallests[index] = kmer
            max_el = max(smallests)
    return smallests
