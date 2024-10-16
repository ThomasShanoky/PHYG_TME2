# Rendu du deuxième TP de Thomas LOUVET et Antoine GRISLAIN



## Résultats 

Nous avons executés le nouvel algorithme (avec les nouvelles fonctions xorshift et filter_smallest) sur les 5 génomes précédemment utilisés afin de vérifier la sortie de l'algorithme en comparaison avec les vraies valeurs des distances de Jacard calculées. La matrice obtenue est : 

|                   | GCA_000013265.1 | GCA_000005845.2 | GCA_000069965.1 | GCA_000008865.2 | GCA_030271835.1 |
|-------------------|-----------------|-----------------|-----------------|-----------------|-----------------|
| GCA_000013265.1   | 1.000000        | 0.192511        | 0.001201        | 0.127192        | 0.001569        |
| GCA_000005845.2   | 0.192511        | 1.000000        | 0.001937        | 0.162825        | 0.003261        |
| GCA_000069965.1   | 0.001201        | 0.001937        | 1.000000        | 0.000861        | 0.020721        |
| GCA_000008865.2   | 0.127192        | 0.162825        | 0.000861        | 1.000000        | 0.001151        |
| GCA_030271835.1   | 0.001569        | 0.003261        | 0.020721        | 0.001151        | 1.000000        |

On remarque que les valeurs approximées semblent être du même ordre que les vraies distances. Cependant, l'ordre de certaines distances ne sont pas respectées (c'est le cas par exemple entre d(GCA_000005845.2, GCA_000013265.1) et d(GCA_000005845.2, GCA_000008865.2)).

Nous n'avons pas pu executer le nouvel algorithme sur les séquences du singe, de la souris et de l'homme puisque les séquences n'étaient pas téléchargeables (données conrompues ou bien les séquences avaient des nucléotides manquants remplacés par des N ou M). 

## Nouvelles fonction et différences avec les anciennes 

### kmers.py

#### encode_nucl(nucl)
Fonction : Encode un nucléotide en un entier de 2 bits et calcule simultanément son complément inverse (reverse complement). Cette méthode utilise des opérations sur des bits pour éviter les instructions conditionnelles, ce qui améliore la performance.
Différence avec encode_nuc : Dans la version précédente, la fonction encode_nuc utilisait un dictionnaire pour associer des lettres (A, C, T, G) à des entiers. Cette nouvelle version utilise des manipulations directes sur les bits ASCII du nucléotide, ce qui est plus rapide et efficace.

#### xorshift(val)
Fonction : Une fonction de hachage (hashing) qui applique l'algorithme de XOR Shift pour mélanger les bits d'un entier, afin de créer une distribution plus uniforme des valeurs des k-mers.
Nouvelle fonction : Cette fonction n'existait pas dans la version précédente. Elle est ajoutée pour améliorer la distribution des valeurs de k-mers, ce qui peut être important lors du calcul de sketchs (empreintes compactes des k-mers).

#### filter_smallest(seq, k, s)
Fonction : Cette nouvelle fonction identifie les s plus petits k-mers dans une séquence donnée, parmi tous les k-mers générés par stream_kmers. Elle maintient une liste des plus petits k-mers et retourne cette liste.
Nouvelle fonction : Elle permet de sélectionner un échantillon réduit des plus petits k-mers, ce qui est utile dans le contexte de la minimisation de l'espace mémoire (sketching) ou pour des méthodes rapides de comparaison de séquences génomiques.

#### stream_kmers(seq, k)
Fonction : Génère les k-mers et leurs compléments inverses à partir d'une séquence, en utilisant des opérations bitwise pour encoder et manipuler les k-mers. La version actuelle utilise la nouvelle fonction xorshift pour hacher les k-mers avant de les retourner.
Différence avec l'ancienne version : L'ancienne version générait les k-mers en encodant directement chaque nucléotide et calculait les k-mers canoniques. La nouvelle version intègre une phase de hachage (XOR Shift), rendant potentiellement les résultats plus adaptés à des applications comme le sketching ou des comparaisons rapides.

### __main__.py

#### jaccard_from_sorted_lists(lstA, lstB)

Fonction : Cette fonction calcule la similarité de Jaccard entre deux listes de k-mers déjà triées (lstA et lstB). Elle procède par une approche à deux pointeurs pour trouver les k-mers communs (intersection) et les k-mers totaux (union).
Différence avec la version précédente : L'ancienne version utilisait directement des ensembles de k-mers, ce qui consommait plus de mémoire et de temps. Cette version optimise la comparaison en utilisant des listes triées, réduisant ainsi les coûts liés à l'intersection et à l'union des k-mers.

#### Paramètres k et s :

Fonction : Définir les valeurs de k = 21 pour la taille des k-mers et s = 10,000 pour la taille de l'échantillon (sketch) permet de réduire la quantité de k-mers considérée dans chaque fichier.
Différence : Cette version utilise le paramètre s pour limiter la quantité de k-mers à comparer, grâce à la nouvelle fonction filter_smallest du fichier kmers.py. Cette approche permet de réduire la consommation mémoire et d'accélérer les comparaisons

#### Génération et tri des k-mers :

Fonction : Les k-mers sont générés à partir des séquences des fichiers avec la fonction filter_smallest, puis triés pour faciliter la comparaison.
Différence : Contrairement à la version précédente qui stockait tous les k-mers sans optimisation, cette version ne conserve que les s plus petits k-mers, ce qui optimise l'espace mémoire et le temps de calcul.

#### Calcul de la similarité de Jaccard :

Fonction : La similarité de Jaccard est calculée pour chaque paire de fichiers à l'aide de la fonction jaccard_from_sorted_lists. Les résultats sont ensuite stockés dans une matrice symétrique.
Différence : La méthode précédente utilisait une approche moins efficace pour les comparaisons, sans le tri préalable des k-mers. Désormais, les listes triées et la fonction de hachage améliorent la performance des comparaisons.