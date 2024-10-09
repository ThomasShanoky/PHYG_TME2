from TP.loading import load_directory
from TP.kmers import stream_kmers, kmer2str, create_index
import time


def jaccard(fileA, fileB, k):
    j = 0
    index = create_index(fileA, k) #dictionnaire
    intersect = 0
    union = sum(index.values()) #le nombre de k-mer dans seq1
    for kmer in stream_kmers(fileB, k): 
        if (kmer in index) and (index[kmer] > 0):
            intersect += 1
            index[kmer] -= 1 #on ne compte pas plusieurs fois les mêmes kmers, on épuise le stock au fur et à mesure
        else:
            union += 1
    j = intersect / union
    return j


if __name__ == "__main__":
    print("Computation of Jaccard similarity between files")

    # Load all the files in a dictionary
    files = load_directory("data")
    k = 21
    
    print("Computing Jaccard similarity for all pairs of samples")
    filenames = list(files.keys())
    n = len(filenames)
    jaccard_matrix = [[1.0 for _ in range(n)] for _ in range(n)]

    init = time.time()
    for i in range(len(files)):
        for j in range(i+1, len(files)):
            dist_j = jaccard(files[filenames[i]], files[filenames[j]], k)
            print(filenames[i], filenames[j], dist_j, time.time() - init)
            jaccard_matrix[i][j] = dist_j
            jaccard_matrix[j][i] = dist_j

    print("Jaccard Distance Matrix:")
    print(" " * 15 + " ".join(f"{name:15}" for name in filenames))
    for i in range(n):
        row = f"{filenames[i]:15}" + " ".join(f"{jaccard_matrix[i][j]:15.6f}" for j in range(n))
        print(row)
