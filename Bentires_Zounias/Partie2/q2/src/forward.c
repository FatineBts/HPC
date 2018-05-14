#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <stdlib.h>
#include <export.h>
#include <mpi.h>
#include <omp.h>

double hFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //HPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return HPHY(t, i, j);
  return HPHY(t - 1, i, j) +
    alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
}

double uFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //UPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return UPHY(t, i, j);
  return UPHY(t - 1, i, j) +
    alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
}

double vFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //VPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return VPHY(t, i, j);
  return VPHY(t - 1, i, j) +
    alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
}

double hPhy_forward(int t, int i, int j) {
  double c, d;
  
  c = 0.;
  if (i > 0)
    c = UPHY(t - 1, i - 1, j);

  d = 0.;
  if (j < size_y - 1)
    d = VPHY(t - 1, i, j + 1);

  return HFIL(t - 1, i, j) -
    dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx +
		 (d - VPHY(t - 1, i, j)) / dy);
}

double uPhy_forward(int t, int i, int j) {
  double b, e, f, g;
  
  if (i == size_x - 1)
    return 0.;

  b = 0.;
  if (i < size_x - 1)
    b = HPHY(t - 1, i + 1, j);

  e = 0.;
  if (j < size_y - 1)
    e = VPHY(t - 1, i, j + 1);

  f = 0.;
  if (i < size_x - 1)
    f = VPHY(t - 1, i + 1, j);

  g = 0.;
  if (i < size_x - 1 && j < size_y - 1)
    g = VPHY(t - 1, i + 1, j + 1);

  return UFIL(t - 1, i, j) +
    dt * ((-grav / dx) * (b - HPHY(t - 1, i, j)) +
	  (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) -
	  (dissip * UFIL(t - 1, i, j)));
}

double vPhy_forward(int t, int i, int j) {
  double c, d, e, f;

  if (j == 0)
    return 0.;

  c = 0.;
  if (j > 0)
    c = HPHY(t - 1, i, j - 1);

  d = 0.;
  if (i > 0 && j > 0)
    d = UPHY(t - 1, i -1, j -1);

  e = 0.;
  if (i > 0)
    e = UPHY(t - 1, i - 1, j);

  f = 0.;
  if (j > 0)
    f = UPHY(t - 1, i, j - 1);

  return VFIL(t - 1, i, j) +
    dt * ((-grav / dy) * (HPHY(t - 1, i, j) - c) -
	  (pcor / 4.) * (d + e + f + UPHY(t - 1, i, j)) -
	  (dissip * VFIL(t - 1, i, j)));
}

void forward(int NP, int rang) {
  MPI_Status status[24];
  MPI_Request request[24];
  int TAG_FIRST_ROW = 0;
  int TAG_LAST_ROW = 1;
  FILE *file = NULL;
  double svdt = 0.;
  int t = 0;
  
  if(rang==0)
  {
    printf("Début if\n");
    if (file_export) {
      file = create_file();
      export_step(file, t);
  }
  }
  printf("Avant tous les send\n");
  
  for (t = 1; t < nb_steps; t++)
  { //&HFIL(t,sixe_x, 0) + utiliser MPI_SendRecv

    if (t == 1) {
      svdt = dt;
      dt = 0;
    }
    if (t == 2){
      dt = svdt / 2.;
    }

  #pragma omp parallel for

    for (int j = 0; j < size_y; j++) {

        #pragma omp parallel for

      for (int i = 0; i <size_x; i++) {
	   HPHY(t, i, j) = hPhy_forward(t, i, j);
	   UPHY(t, i, j) = uPhy_forward(t, i, j);
	   VPHY(t, i, j) = vPhy_forward(t, i, j);
	   HFIL(t, i, j) = hFil_forward(t, i, j);
	   UFIL(t, i, j) = uFil_forward(t, i, j);
	   VFIL(t, i, j) = vFil_forward(t, i, j);
          }
        }
                                       

  /* int MPI_Isend(
  void *buf,
  int count,
  MPI_Datatype datatype,
  int dest,
  int tag,
  MPI_Comm comm,
  MPI_Request *request

  int MPI_Irecv(
  void *buf,
  int count,
  MPI_Datatype datatype,
  int source,
  int tag,
  MPI_Comm comm,
  MPI_Request *request
);
);*/



  if (rang!=0)
      {
        MPI_Isend(&HFIL(t,1,0),size_y, MPI_DOUBLE, rang-1, TAG_FIRST_ROW, MPI_COMM_WORLD,&request[0]);
        MPI_Irecv(&HFIL(t,0,0),size_y, MPI_DOUBLE, rang-1, TAG_LAST_ROW, MPI_COMM_WORLD,&request[1]);

        MPI_Isend(&UFIL(t,1,0),size_y, MPI_DOUBLE, rang-1, TAG_FIRST_ROW, MPI_COMM_WORLD,&request[2]);
        MPI_Irecv(&UFIL(t,0,0),size_y, MPI_DOUBLE, rang-1, TAG_LAST_ROW, MPI_COMM_WORLD,&request[3]);

        MPI_Isend(&VFIL(t,1,0),size_y, MPI_DOUBLE, rang-1, TAG_FIRST_ROW, MPI_COMM_WORLD,&request[4]);
        MPI_Irecv(&VFIL(t,0,0),size_y, MPI_DOUBLE, rang-1, TAG_LAST_ROW, MPI_COMM_WORLD,&request[5]);

        MPI_Isend(&HPHY(t,1,0),size_y, MPI_DOUBLE, rang-1, TAG_FIRST_ROW, MPI_COMM_WORLD,&request[6]);
        MPI_Irecv(&HPHY(t,0,0),size_y, MPI_DOUBLE, rang-1, TAG_LAST_ROW, MPI_COMM_WORLD,&request[7]);

        MPI_Isend(&UPHY(t,1,0),size_y, MPI_DOUBLE, rang-1, TAG_FIRST_ROW, MPI_COMM_WORLD,&request[8]);
        MPI_Irecv(&UPHY(t,0,0),size_y, MPI_DOUBLE, rang-1, TAG_LAST_ROW, MPI_COMM_WORLD,&request[9]);

        MPI_Isend(&VPHY(t,1,0),size_y, MPI_DOUBLE, rang-1, TAG_FIRST_ROW, MPI_COMM_WORLD,&request[10]);
        MPI_Irecv(&VPHY(t,0,0),size_y, MPI_DOUBLE, rang-1, TAG_LAST_ROW, MPI_COMM_WORLD,&request[11]);
      }
 

  if (rang!=NP-1) 
      {
        MPI_Isend(&HFIL(t,(size_x-2),0),size_y, MPI_DOUBLE, rang+1, TAG_LAST_ROW, MPI_COMM_WORLD,&request[12-12*(rang==0)]);
        MPI_Irecv(&HFIL(t,(size_x-1),0),size_y, MPI_DOUBLE, rang+1, TAG_FIRST_ROW, MPI_COMM_WORLD,&request[13-12*(rang==0)]);

        MPI_Isend(&UFIL(t,(size_x-2),0),size_y, MPI_DOUBLE, rang+1, TAG_LAST_ROW, MPI_COMM_WORLD,&request[14-12*(rang==0)]);
        MPI_Irecv(&UFIL(t,(size_x-1),0),size_y, MPI_DOUBLE, rang+1, TAG_FIRST_ROW, MPI_COMM_WORLD,&request[15-12*(rang==0)]);

        MPI_Isend(&VFIL(t,(size_x-2),0),size_y, MPI_DOUBLE, rang+1, TAG_LAST_ROW, MPI_COMM_WORLD,&request[16-12*(rang==0)]);
        MPI_Irecv(&VFIL(t,(size_x-1),0),size_y, MPI_DOUBLE, rang+1, TAG_FIRST_ROW, MPI_COMM_WORLD,&request[17-12*(rang==0)]);

        MPI_Isend(&HPHY(t,(size_x-2),0),size_y, MPI_DOUBLE, rang+1, TAG_LAST_ROW, MPI_COMM_WORLD,&request[18-12*(rang==0)]);
        MPI_Irecv(&HPHY(t,(size_x-1),0),size_y, MPI_DOUBLE, rang+1, TAG_FIRST_ROW, MPI_COMM_WORLD,&request[19-12*(rang==0)]);

        MPI_Isend(&UPHY(t,(size_x-2),0),size_y, MPI_DOUBLE, rang+1, TAG_LAST_ROW, MPI_COMM_WORLD,&request[20-12*(rang==0)]);
        MPI_Irecv(&UPHY(t,(size_x-1),0),size_y, MPI_DOUBLE, rang+1, TAG_FIRST_ROW, MPI_COMM_WORLD,&request[21-12*(rang==0)]);

        MPI_Isend(&VPHY(t,(size_x-2),0),size_y, MPI_DOUBLE, rang+1, TAG_LAST_ROW, MPI_COMM_WORLD,&request[22-12*(rang==0)]);
        MPI_Irecv(&VPHY(t,(size_x-1),0),size_y, MPI_DOUBLE, rang+1, TAG_FIRST_ROW, MPI_COMM_WORLD,&request[23-12*(rang==0)]);
	 }

    MPI_Waitall(24-12*((rang!=0)||(rang!=NP-1)),request,status);


      if(file_export) {
      MPI_Gather(&HFIL(t,(rang!=0), 0),(gsize_x/NP)*gsize_y, MPI_DOUBLE, &ghfil(t, 0, 0), (gsize_x/NP)*gsize_y, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }


    if(rang==0)
    {
      if (file_export) {
      export_step(file, t);
        }
    }

    
    if (t == 2) {
      dt = svdt;
    }
  }

    if(rang==0)
    {
      if (file_export) 
      {
      finalize_export(file);
      printf("\n\n");    
      }
    }

}
