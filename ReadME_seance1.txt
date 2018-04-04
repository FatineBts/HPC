######################################
Prémices : 

1) On copie le dossier : mkdir /tmp/3302011
2) On lance le makefile : make
3) On crée un dossier avec le numéro étudiant dans tmp sur notre racine : mkdir /tmp/3302011
4) On lance shalw avec export dans le fichier num_etudiant créé : ./bin/shalw --export --export-path /tmp/3302011 (simulation complète)
On obtient shalw_256x256_T1000.sav dans /tmp/3302011
5) On utilise visu avec l'image enregistrée dans le fichier du num_etudiant : ./visu.py /tmp/3302011/shalw_256x256_T1000.sav

Ainsi, le fichier visu nous permet de visualiser ce qu'on devrait obtenir si on réorientait le calcul de notre executable shawl vers un fichier .ras. 

#################
Dossiers : 
- bin  : fichiers executables
- inc : fichiers headers 
- obj : fichiers objets 
- src : fichiers .c et .cpp

Fichiers : 
- Makefile 
- visu.py 

##########################
Explication du but des fichiers : 

- export.c : crée le fichier .ras 
- init.c : initialise la grille avec la méthode de Gauss 
- memory. c : fait une allocation et libération de mémoire 
- parse_args.cpp : permet de récupérer les éléments tapés au clavier et si on tape 'help', nous affiche les éléments et leur explication
- shalw.cpp : fichier main  

Explication du code : 

- Fichier forward.c : HFIL, VFIL : macro accesseurs. HFIL est une macro qui prend t, i et j et qui associe un calcul donnant le pointeur vers le i, j au temps de pas t vers la structure associée. 

############################
Partie 1 : parallélisation MPI

1. On commence par le découpage par bande : on utilise le cas 2 

—  Cas 2, pour les mesures de performances en parallèle :
$ ./bin/shalw -x 8192 -y 8192 -t 20
Il est strictement interdit de générer un export pour ce cas. Si nécessaire ou pertinent, on pourra aussi s’intéresser
à des grilles de dimensions plus petites (mais toujours sans export).


Modifications faites : 

- Création des dossiers
- Ajout de la bu mpi dans le fichier header shalx.hpp
- Ajout de la variable p int dans shalw.cpp

- Modification du main dans shalw.cpp : 

* passage de w à size_y 
* passage de local_h à size_x 

Comme size_x et size_y ne sont pas modifiés dans le reste du code, on fait des modifications pour dire qu'il s'agit des variables x et y d'une bande (donc dimensions locales). 

La méthode utilisée par bande repose sur la méthode du TP4.c 