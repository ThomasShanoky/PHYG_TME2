# Rendu de Thomas LOUVET et Antoine GRISLAIN

## Jaccard Distance Matrix:

|               |GCA_000008865.2| GCA_000069965.1| GCA_000013265.1| GCA_030271835.1| GCA_000005845.2|
|GCA_000008865.2|       1.000000|        0.002314|        0.307048|        0.002318|        0.436482|
|GCA_000069965.1|       0.002314|        1.000000|        0.002437|        0.031134|        0.002567|
|GCA_000013265.1|       0.307048|        0.002437|        1.000000|        0.002434|        0.341008|
|GCA_030271835.1|       0.002318|        0.031134|        0.002434|        1.000000|        0.002577|
|GCA_000005845.2|       0.436482|        0.002567|        0.341008|        0.002577|        1.000000|

On peut voir grâce à la matrice des distances de Jaccard que les deux séquences les plus proches sont GCA_000008865.2 et GCA_000005845.2. 
De manière générale les valeurs sont assez proches de 0 ce qui signifie que ces ensembles partagent très peu de k-mers et sont donc très différent. 
On remarque tout de même quelques valeurs autour de 0.3 0.4 qui indiquent une similarité modéré.

## Description des méthodes implémentées 

### Fichier main :

jaccard(fileA, fileB, k) : Cette fonction calcule la similarité de Jaccard entre deux séquences en utilisant des k-mers de taille k. Elle génère un dictionnaire d'index à partir du premier fichier et compare les k-mers du second fichier pour déterminer l'intersection et l'union.

Main : Charge les fichiers de séquences dans un dictionnaire depuis le répertoire "data" et concatène les séquences par fichier. Calcule ensuite la similarité de Jaccard pour chaque paire de fichiers et affiche une matrice des distances.

### Fichier kmers :

create_index(seq, k) : Crée un index sous forme de dictionnaire contenant les k-mers de taille k extraits d'une séquence donnée.

encode_nuc(letter) : Encode un nucléotide (A, C, T, G) en un entier, selon l'ordre ACTG : A = 0, C = 1, T = 2, G = 3.

encode_kmer(seq, k) : Encode une séquence de k-mer en entier en concaténant les valeurs entières de chaque nucléotide de la séquence.

reverse_complement(kmer, k) : Calcule le complément inverse d'un k-mer encodé sous forme d'entier en inversant et en prenant le complément de chaque nucléotide.

stream_kmers(seq, k) : Génère un flux de k-mers à partir d'une séquence, permettant une itération sur les k-mers d'une séquence sans charger tous les k-mers en mémoire.

kmer2str(kmer) : Convertit un k-mer représenté en format numérique en chaîne de caractères, facilitant la manipulation et l'affichage.

### Fichier loading :

load_fasta(file_pointer) : Charge un fichier FASTA (ou une archive gzip) et extrait les séquences sous forme de liste de chaînes.

load_directory(directory) : Parcourt un répertoire contenant des fichiers FASTA et charge tous les fichiers trouvés dans un dictionnaire où chaque clé est un nom de fichier et la valeur est une liste de séquences.