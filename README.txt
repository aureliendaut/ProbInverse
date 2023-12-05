Résolution numérique de l'équation de la chaleur en 2D

Ce programme résout numériquement l'équation de la chaleur en 2D avec différentes conditions de bord et sources. Les paramètres sont spécifiés dans le fichier 'parametres.txt', qui doit être fourni en argument lors de l'exécution du programme. Pour l'instant la résolution s'effectue seulement de manière séquentielle.

Organisation du code

Le fichier 'main.c' contient la lecture des paramètres de résolution présente dans le fichier 'parametres.txt'. La boucle principale fait appel à la fonction 'BuildSecondMember' et 'GC' qui permettent, respectivement, de construire le second membre lors de la résolution de "Au=b" et le gradient conjugué résout ce système. Pour éviter de stocker la matrice, le fichier 'matrice.c' permet de calculer le produit Ax pour la matrice A de notre problème. La fonction 'BuildSecondMember' se trouve dans le fichier 'matrice.c'. Les fonctions pour construire les bords et le terme source se trouvent dans le fichier 'function.c'. Les fichiers '.txt' sont les fichiers d'instructions pour gnuplot.

Compilation

Pour compiler le programme, il suffit d'exécuter la commande 'make' à partir du terminal. Cela va créer un exécutable nommé 'chp'.
2 options de compilation sont possibles, 1 pour le débuggage : 'DEBUG_FLAG' ; 1 pour l'optimisation des calculs : 'OPTIM_FLAG'. Pour changer celle-ci il faut changer le nom dans le Makefile de 'CXX_FLAGS'.

Exécution

Pour exécuter le programme, utilisez la commande './chp "parametres.txt" ', où 'parametres.txt' contient les paramètres suivants :
- 'Nx' : nombre de points de discrétisation dans la direction x.
- 'Ny' : nombre de points de discrétisation dans la direction y.
- 'Lx' : longueur du domaine dans la direction x.
- 'Ly' : longueur du domaine dans la direction y.
- 'D' : diffusivité thermique.
- 'dt' : pas de temps.
- 'tfinal' : temps final de la simulation.
- 'cas' : spécifie les conditions de bord et sources.
Les différents choix possibles pour 'cas' sont les suivants :
- 'cas=1' : terme source stationnaire et conditions de bord de Dirichlet homogènes.
- 'cas=2' : terme source stationnaire et conditions de bord de Dirichlet non homogènes.
- 'cas=3' : terme source instationnaire et conditions de bord de Dirichlet non homogènes.

Visualisation des résultats

Pour visualiser les solutions numériques, utilisez la commande 'make images'. Cela va générer des fichiers PNG pour chaque instant de temps de la simulation.
Pour visualiser les solutions exactes (si disponibles), utilisez la commande 'make images_exactes'. Cela va générer des fichiers PNG pour chaque instant de temps de la simulation.

Attention
Les temps indiqués en légende des images ne sont pas automatiquement modifiés en fonction du pas de temps, du temps final de résolution et le nombre d'images crée aussi. Il faut pour cela modifier les titres, le pas de temps et le nombre d'images à la main dans les fichiers 'gnuplot.txt' et 'gnuplot_exacte.txt'. De plus, pour le cas 3, n'ayant pas de solution exacte, il n'est pas possible de faire la commande 'make images_exactes' car les fichiers sont vides.

Validation des résultats

Les résultats ont été validés à l'aide des cas 1 et 2. Pour ces 2 cas, il est possible de calculer la solution exacte et donc, il est possible de vérifier l'ordre des schémas numériques en traçant log(erreur) en fonction de log(dt) (ou log du pas d'espace pour vérifier l'ordre d'espace) où l'erreur est définie comme la norme de la solution exacte - la solution numérique. Ce code obtient un ordre de 0.99 en temps.

Résultats obtenus

Les images suivantes présentent quelques résultats possibles d'obtenir avec le code pour les 3 cas présentés précédemment (Lire le fichier 'README.md' pour les visualiser ou bien elles sont visible dans le dossier 'image).

Nettoyage

Pour supprimer l'exécutable et les images générées, utilisez la commande 'make clean'.
