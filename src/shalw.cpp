#include <stdlib.h>
#include <shalw.h>
#include <parse_args.hpp>
#include <memory.h>
#include <init.h>
#include <forward.h>
#include <export.h>
#include <mpi.h>
#include <math.h>   /* pour le rint */
#include <time.h>   /* chronometrage */
#include "sys/time.h"

double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy;
double *g_hFil;//, *g_uFil, *g_vFil, *g_hPhy, *g_uPhy, *g_vPhy;
int size_x, size_y;
int g_size_x, g_size_y, nb_steps;
double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
bool file_export;
std::string export_path;

double my_gettimeofday(){
  struct timeval tmp_time;
  gettimeofday(&tmp_time, NULL);
  return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

int main(int argc, char **argv) {
 	/* Par processeur */ 
  int rang; // rang
  int NP; // NP = nombre de processus



 /* Variables liees au chronometrage */
  double debut=0, fin=0;

   /* debut du chronometrage */
   
 
  parse_args(argc, argv);
  printf("Command line options parsed\n");

  gauss_init();
  printf("State initialised\n"); 

    /* Initialisation MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &NP);
  MPI_Comm_rank(MPI_COMM_WORLD, &rang); 
	
  if (rang ==0)
  {
    debut = my_gettimeofday(); 
		alloc();
		printf("Memory allocated by rank 0\n");  
		/* Initialisations / calcul */		

	}

  /* Permet de donner la taille de sixe_x et size_y*/
  size_y = g_size_y ;
  size_x = (rang==0 || rang==NP-1)?(g_size_x/NP +1):(g_size_x/NP +2);
  
  loc_alloc();
  printf("Local memory allocated. Rang = %d \n", rang); 
   
//  printf("Avant Scatter \n");
	//MPI_Scatter(g_hFil /*sbuf*/, size_x/NP*size_y /*scount*/, MPI_DOUBLE /*sdtype*/, hFil+size_y*(rang!=0) /*rbuf*/, size_x/NP*size_y /*rcount*/, MPI_DOUBLE /*rdtype*/, 0 /*root*/, MPI_COMM_WORLD /*comm*/);
  //printf("Après Scatter \n");

	forward(NP, rang); // MPI send and receive 
	printf("State computed\n");



  if (rang ==0)
  {
    /* fin du chronometrage */
    fin = my_gettimeofday();
    printf("Temps total de calcul : %g seconde(s) \n", fin - debut);
    dealloc(); // MÉMOIRE 

  }
  dealloc_2(); 
  printf("Memory freed\n");
  MPI_Finalize();
  
  return EXIT_SUCCESS;
}
