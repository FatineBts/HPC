# HPC

EN: the “shallow water” model represents the vertically homogenous flow of a fluid. The equation system is easy to understand but more difficult to solve. Moreover, it’s still a subject of research. The project can be decomposed in two parts. The first one corresponds to the MPI parallelization: band and block decomposition, blocking and non blocking communications. The second one releases on the utilisation of OpenMP and SIMD. 

FR : le modèle « shallow water » permet de représenter l’écoulement d’un fluide homogène à la verticale. Le système d’équation est relativement simple à poser mais leur résolution est toujours aujourd’hui un sujet de recherche. Le projet est décomposé en deux parties. La première partie correspond à la parallélisation MPI : décompositions par bande et par bloc, bloquante et non bloquante. La seconde partie correspond à l’utilisation d’OpenMP et SIMD : prise en compte de la hiérarchie mémoire, parallélisation MPI+OpenMP et SIMD.
