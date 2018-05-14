#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <export.h>
#include <immintrin.h> 


int execution(int t, int i, int j)
{

  __m256d hf1,hf2, h_alpha,h_hfil, h_hphy,h2_hphy, tmp; 

  h_alpha=_mm256_set1_pd(alpha); 
  h_hfil=_mm256_set1_pd(HFIL(t - 1, i, j));
  h_hphy=_mm256_set1_pd(HPHY(t - 1, i, j));
  h2_hphy=_mm256_set1_pd(HPHY(t, i, j));

  if (t <= 2)
  {
    HPHY(t, i, j)=HPHY(t, i, j);
    UPHY(t, i, j)=UPHY(t, i, j); 
    VPHY(t, i, j)=VPHY(t, i, j); 
  }

  hf1 = _mm256_mul_pd(h_alpha ,h_hfil);
  hf2 = hf1 - 2*h_hphy;
  tmp=_mm256_add_pd(h_hphy,h2_hphy);

  HPHY(t, i, j)=_mm256_add_pd(tmp,hf2);
  printf("%p",_mm256_add_pd(tmp,hf2));
  UPHY(t, i, j)=UPHY(t - 1, i, j) + alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
  VPHY(t, i, j)=VPHY(t - 1, i, j) + alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));

  double c, d;
  double b, e, f, g;
  
  c = 0.;
  if (i > 0)
    c = UPHY(t - 1, i - 1, j);

  d = 0.;
  if (j < size_y - 1)
    d = VPHY(t - 1, i, j + 1);

  HPHY(t,i,j)=HFIL(t - 1, i, j) - dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx +(d - VPHY(t - 1, i, j)) / dy);


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

  UPHY(t,i,j)=UFIL(t - 1, i, j) +dt * ((-grav / dx) * (b - HPHY(t - 1, i, j)) + (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) -(dissip * UFIL(t - 1, i, j)));


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

  VPHY(t,i,j)=VFIL(t - 1, i, j) +dt * ((-grav / dy) * (HPHY(t - 1, i, j) - c) -(pcor / 4.) * (d + e + f + UPHY(t - 1, i, j)) -(dissip * VFIL(t - 1, i, j)));

}

void forward(void) {
  FILE *file = NULL;
  double svdt = 0.;
  int t = 0;
  __m256d HP, UP, VP, HF, UF, VF; 
  
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

    for (int j = 0; j < size_y; j++) {
      for (int i = 0; i < size_x; i++) {

        execution(t,i,j); 

        HP = _mm256_load_pd(&HPHY(t,i,j));
        UP = _mm256_load_pd(&UPHY(t,i,j));
        VP = _mm256_load_pd(&VPHY(t,i,j));
        HF = _mm256_load_pd(&HFIL(t,i,j));
        UF = _mm256_load_pd(&UFIL(t,i,j));
        VF = _mm256_load_pd(&VFIL(t,i,j));
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
