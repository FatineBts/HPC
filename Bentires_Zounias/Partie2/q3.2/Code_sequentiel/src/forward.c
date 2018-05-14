#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <export.h>
#include <immintrin.h> 

void forward(void) {
  FILE *file = NULL;
  double svdt = 0.;
    double c, d, b, e, f, g;
    int toto=1; 
  
  int t = 0;

 __m256d  h_hfil, h_hphy, h_ufil, h_vfil, h2_vfil, h2_hphy, h2_uphy, h2_vphy;
 __m256d h2_hfil, h2_ufil, h_dt, h_i_dx, h_d, h_uphy, h_vphy, h_hmoy, h_c, h_i_dy, t1, t2, t3, t4, t5; 


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


/*On a inversé les boucles*/
//////////////////////////////////////////////////////////////////////////////////////////////////////
    for (int i = 0; i < size_x/4; i++) { // /4 car vecteur de 32 et chaque individu de 8 octets donc 4 éléments 
      for (int j = 0; j < size_y/4; j++) {

      // on charge les macros 
      HP = _mm256_load_pd(&HPHY(t,i,j)); // il faudra multiplier par 4 pour faire des sauts de 4 
      UP = _mm256_load_pd(&UPHY(t,i,j));
      VP = _mm256_load_pd(&VPHY(t,i,j));
      HF = _mm256_load_pd(&HFIL(t,i,j));
      UF = _mm256_load_pd(&UFIL(t,i,j));
      VF = _mm256_load_pd(&VFIL(t,i,j));

  // variables 

        // set1 = initialise un vector avec la valeur d'un double 
  h_hfil=_mm256_set1_pd(HFIL(t - 1, i, j));
  h_ufil=_mm256_set1_pd(UFIL(t - 1, i, j));
  h_vfil=_mm256_set1_pd(VFIL(t - 1, i, j));
  h_hphy=_mm256_set1_pd(HPHY(t - 1, i, j));
  h_uphy=_mm256_set1_pd(UPHY(t - 1, i, j));
  h_vphy=_mm256_set1_pd(VPHY(t - 1, i, j));
   

  h2_hphy=_mm256_set1_pd(HPHY(t, i, j));
  h2_uphy=_mm256_set1_pd(UPHY(t, i, j));
  h2_vphy=_mm256_set1_pd(VPHY(t, i, j));
  h2_hfil=_mm256_set1_pd(HFIL(t, i, j));
  h2_ufil=_mm256_set1_pd(UFIL(t, i, j));
  h2_vfil=_mm256_set1_pd(VFIL(t, i, j));
  h_dt=_mm256_set1_pd(dt);
  h_hmoy=_mm256_set1_pd(hmoy);
  h_c=_mm256_set1_pd(c);
  h_i_dx=_mm256_set1_pd(1./dx);
  h_d=_mm256_set1_pd(d);
  h_i_dy=_mm256_set1_pd(dy);


  // Partie 1 : hPhy_forward


  c = 0.;
  if (i > 0)
    c = UPHY(t - 1, i - 1, j);

  d = 0.;
  if (j < size_y - 1)
    d = VPHY(t - 1, i, j + 1);

  t1=_mm256_mul_pd(h_dt ,h_hmoy); 
  t2=_mm256_mul_pd((h_uphy - h_c),h_i_dx);
  t3=_mm256_mul_pd(t1,t2);
  t4=_mm256_mul_pd((h_d - h_vphy),h_i_dy);
  t5=_mm256_add_pd(t3,t4);

  h2_hphy= h_hphy-t5;


  //Partie 2 : uPhy_forward

  if (i == size_x - 1)
    UPHY(t,i,j)=0.;

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
 
    h2_uphy=h_ufil +h_dt* ((-grav * h_i_dx) * (b - h_hphy) + (pcor / 4.) * (h_vphy + e + f + g) -(dissip * h_ufil));
 
   if(toto==1)
  {
    printf("toto\n");
    toto=2; 
  }

  //Partie 3 : vPhy_forward

  if (j == 0)
    VPHY(t,i,j)=0.;

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
 
    h2_vphy=h_vfil +dt * ((-grav / dy) * (h_hphy - c) -(pcor / 4.) * (d + e + f + h_uphy) -(dissip * h_vfil));


  //Partie 4 : hFil_forward

  if (t <= 2)
    h2_hfil=h2_hphy;

  else 
    h2_hfil=h_hphy +alpha * (h_hfil - 2 * h_hphy + h2_hphy);//_mm256_add_pd(tmp,hf2);

  //Partie 5 : uFil_forward 

  if(t<=2)
    h2_ufil=h2_uphy;

  else 
    h2_ufil=h_uphy + alpha * (h_ufil - 2 * h_uphy + h2_uphy);


  //Partie 6 : vFil_forward

  if(t<=2)
    h2_vfil=h2_vphy;

  else 
    h2_vfil=h_vphy + alpha * (h_vfil - 2 * h_vphy + h2_vphy);


 // hf1 = _mm256_mul_pd(h_alpha ,h_hfil);
  //hf2 = hf1 - 2*h_hphy;
  //tmp=_mm256_add_pd(h_hphy,h2_hphy);



      //pour récupérer les infos de la mémoire 
  // déplace ce qui est contenu dans le vector vers la variable double 
      _mm256_store_pd(&HPHY(t,i,j),h2_hphy); // il faut faire des sauts de 4 
      _mm256_store_pd(&UPHY(t,i,j),h2_uphy); 
      _mm256_store_pd(&VPHY(t,i,j),h2_vphy); 
      _mm256_store_pd(&HFIL(t,i,j),h2_hfil); 
      _mm256_store_pd(&UFIL(t,i,j),h2_ufil); 
      _mm256_store_pd(&VFIL(t,i,j),h2_vfil); 
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
