#include <stdlib.h>
#include <shalw.h>

void alloc(void) {
  g_hFil = (double *) calloc(2*g_size_x*g_size_y, sizeof(double)); // on utilise deux grilles. En fonction de t, on accède soit à la première, soit à la deuxième.
  // g_uFil = (double *) calloc(2*g_size_x*g_size_y, sizeof(double));
  // g_vFil = (double *) calloc(2*g_size_x*g_size_y, sizeof(double));
  // g_hPhy = (double *) calloc(2*g_size_x*g_size_y, sizeof(double));
  // g_uPhy = (double *) calloc(2*g_size_x*g_size_y, sizeof(double));
  // g_vPhy = (double *) calloc(2*g_size_x*g_size_y, sizeof(double));
}

void loc_alloc(void) {
  hFil = (double *) calloc(2*size_x*size_y, sizeof(double)); 
  uFil = (double *) calloc(2*size_x*size_y, sizeof(double));
  vFil = (double *) calloc(2*size_x*size_y, sizeof(double));
  hPhy = (double *) calloc(2*size_x*size_y, sizeof(double));
  uPhy = (double *) calloc(2*size_x*size_y, sizeof(double));
  vPhy = (double *) calloc(2*size_x*size_y, sizeof(double));
}

void dealloc(void) {
  // free(g_uFil);
  free(g_hFil);
  // free(g_vFil);
  // free(g_hPhy);
  // free(g_uPhy);
  // free(g_vPhy);

}

void dealloc_2(void)
{
  free(hFil);
  free(uFil);
  free(vFil);
  free(hPhy);
  free(uPhy);
  free(vPhy);
}
