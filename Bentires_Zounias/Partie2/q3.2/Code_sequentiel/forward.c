#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <export.h>
#include <immintrin.h> 

_m256d hFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //HPHY(t - 1, i, j) est encore nul
  _m256d h_hphy,h_hphy1,h_hfil1,m,m1,deux,alpha2,add,add2,sous; 

  h_hphy =  _mm256_load_pd(&HPHY(t,i,j*4));

  if (t <= 2)
    return h_hphy;

  h_hphy1 =_mm256_load_pd(&HPHY(t-1,i,j*4));
  h_hfil1 =_mm256_load_pd(&HFIL(t-1,i,j*4));
  alpha2 = _mm256_set1_pd(alpha);
  deux=_mm256_set1_pd(2);

  m=_mm256_mul_pd(alpha2,h_hfil1); 
  m1=_mm256_mul_pd(deux,h_hphy1);
  sous= _mm256_sub_pd(m,m1);
  add=_mm256_add_pd(h_hphy1,sous); 
  add2=_mm256_add_pd(add,h_hphy); 


  return (add2);
  //return HPHY(t - 1, i, j) +
   // alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
}

_m256d uFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //UPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return UPHY(t, i, j);
  return UPHY(t - 1, i, j) +
    alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
}

_m256d vFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //VPHY(t - 1, i, j) est encore nul
  if (t <= 2)
    return VPHY(t, i, j);
  return VPHY(t - 1, i, j) +
    alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
}

_m256d hPhy_forward(int t, int i, int j) {
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

_m256d uPhy_forward(int t, int i, int j) {
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

_m256d vPhy_forward(int t, int i, int j) {
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

void forward(void) {
  __m256d  h_hfil, h_hphy, h_uphy, h_vphy,h_ufil, h_vfil; 
  FILE *file = NULL;
  double svdt = 0.;
  int t = 0;
  
  if (file_export) {
    file = create_file();
    export_step(file, t);
  }
    
  for (t = 1; t < nb_steps; t++) {
    if (t == 1) {
      svdt = dt;
      dt = 0;
    }
    if (t == 2){
      dt = svdt / 2.;
    }

    for (int i = 0; i < size_x; i++) {
      for (int j = 0; j < size_y/4; j++) {

        h_hphy = hPhy_forward(t, i, j);
        h_uphy = uPhy_forward(t, i, j);
        h_vphy = vPhy_forward(t, i, j);
        h_hfil = hFil_forward(t, i, j);
        h_ufil = uFil_forward(t, i, j);
        h_vfil = vFil_forward(t, i, j);


      _mm256_store_pd(&HPHY(t,i,j*4),h_hphy); // il faut faire des sauts de 4 
      _mm256_store_pd(&UPHY(t,i,j*4),h_uphy); 
      _mm256_store_pd(&VPHY(t,i,j*4),h_vphy); 
      _mm256_store_pd(&HFIL(t,i,j*4),h_hfil); 
      _mm256_store_pd(&UFIL(t,i,j*4),h_ufil); 
      _mm256_store_pd(&VFIL(t,i,j*4),h_vfil); 
      }
    }

    if (file_export) {
      export_step(file, t);
    }
    
    if (t == 2) {
      dt = svdt;
    }
  }

  if (file_export) {
    finalize_export(file);
  }
}
