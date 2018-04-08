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
  int test = 0; // Vaut true si c'est une ligne et false si c'est une colonne 
  //int test2=1; 

    /* Initialisation MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &NP);
  MPI_Comm_rank(MPI_COMM_WORLD, &rang); 
  
 /* Variables liees au chronometrage */
  double debut=0, fin=0;

   /* debut du chronometrage */
   
   debut = my_gettimeofday();

  parse_args(argc, argv);
  printf("Command line options parsed\n");
  int i=NP; 

  /* Permet de donner la taille de sixe_x et size_y*/
  size_y = g_size_y;
  size_x = g_size_x; 

while(i>1)
{

    if(test == 0)
    {
      size_x/=2; 
      test=1;
    }

    if(test ==1)
    {
      size_y/=2; 
      test=0; 
    }

  i=i/2;  
}
     
  if (rang==0)
  {
    alloc();
    printf("Memory allocated by rank 0\n"); 
    gauss_init();
    printf("State initialised\n");  
   
    /* Initialisations / calcul */    
  }

  alloc_2();
  printf("Local memory allocated. Rang = %d \n", rang); 

  MPI_Scatter(g_hFil /*sbuf*/, size_x*size_y /*scount*/, MPI_DOUBLE /*sdtype*/, hFil+size_y*(rang!=0) /*rbuf*/, size_x*size_y /*rcount*/, MPI_DOUBLE /*rdtype*/, 0 /*root*/, MPI_COMM_WORLD /*comm*/);
 
  forward(NP, rang); // MPI send and receive 
  printf("State computed\n");

  /* fin du chronometrage */
  fin = my_gettimeofday();
  printf("Temps total de calcul : %g seconde(s) \n", fin - debut);

  if(rang==0)
    {
	   dealloc(); // MÉMOIRE 
	   printf("Memory freed 1\n");
    }
   
  dealloc_2(); 
  printf("Memory freed 2\n");
   
  MPI_Finalize(); 

  return EXIT_SUCCESS;
}
