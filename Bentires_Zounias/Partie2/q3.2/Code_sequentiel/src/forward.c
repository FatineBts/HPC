#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <export.h>
#include <immintrin.h> 


inline __m256d hFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //HPHY(t - 1, i, j) est encore nul

  __m256d h_hphy,h_hphy1,h_hfil1,m,m1,deux,alpha2,add,res,sous; 

  h_hphy =  _mm256_loadu_pd(&HPHY(t,i,j));

  if (t <= 2)
  {
    res =  h_hphy;
  }
  else
  {
    h_hphy1 =_mm256_loadu_pd(&HPHY(t-1,i,j));
    h_hfil1 =_mm256_loadu_pd(&HFIL(t-1,i,j));
    alpha2 = _mm256_set1_pd(alpha);
    deux=_mm256_set1_pd(2.);

    m=_mm256_mul_pd(deux,h_hphy1); // 2 * HPHY(t - 1, i, j)
    sous= _mm256_sub_pd(h_hfil1,m); //HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j)
    add=_mm256_add_pd(h_hphy,sous);  //HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j)
    m1=_mm256_mul_pd(alpha2,add);
    res=_mm256_add_pd(m1,h_hphy1); 
  }


  return (res);
  //return HPHY(t - 1, i, j) +
   // alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
}

inline __m256d uFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //UPHY(t - 1, i, j) est encore nul

  __m256d h_uphy, h_uphy1, h_ufil1, m,m1,deux,alpha2,add,res,sous;
  h_uphy =  _mm256_loadu_pd(&UPHY(t, i, j));

  if (t <= 2)
  {
    res = h_uphy;
    //return UPHY(t, i, j);
  }
  else
  {
    h_uphy1 =  _mm256_loadu_pd(&UPHY(t-1, i, j));
    h_ufil1 = _mm256_loadu_pd(&UFIL(t-1, i, j));
    alpha2 = _mm256_set1_pd(alpha);
    deux=_mm256_set1_pd(2.);

 
    m=_mm256_mul_pd(deux,h_uphy1); // 2 * UPHY(t - 1, i, j)
    sous= _mm256_sub_pd(h_ufil1,m); // UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j)
    add=_mm256_add_pd(h_uphy,sous); // UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j)
    m1 = _mm256_mul_pd(alpha2,add); //  alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
    res=_mm256_add_pd(h_uphy1,m1); 
  }
  return(res);
  //return UPHY(t - 1, i, j) +
  //  alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
}

inline __m256d vFil_forward(int t, int i, int j) {
  //Phase d'initialisation du filtre
  //VPHY(t - 1, i, j) est encore nul
  __m256d h_vphy, h_vphy1, h_vfil1, m,m1,deux,alpha2,add,res,sous;
  h_vphy =  _mm256_loadu_pd(&VPHY(t, i, j));

  if (t <= 2)
    res = h_vphy;
    //return VPHY(t, i, j);
  else
  {
    h_vphy1 = _mm256_loadu_pd(&VPHY(t-1, i, j));
    h_vfil1 = _mm256_loadu_pd(&VFIL(t-1, i, j));
    alpha2 = _mm256_set1_pd(alpha);
    deux=_mm256_set1_pd(2.);

    m=_mm256_mul_pd(deux,h_vphy1); 
    sous= _mm256_sub_pd(h_vfil1,m); 
    add=_mm256_add_pd(h_vphy,sous); 
    m1 = _mm256_mul_pd(alpha2,add); 
    res=_mm256_add_pd(h_vphy1,m1); 
  }
  return(res);
    
  //return VPHY(t - 1, i, j) +
    //alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
}

inline __m256d hPhy_forward(int t, int i, int j) {


  __m256d dx_i, dy_i, hmoy2, h_uphy1, res, h_vphy1, h_hfil1, m, sous, m2, m3, sous2, m4, add,c,d, dt1;

  dt1 = _mm256_set1_pd(dt);
  dx_i = _mm256_set1_pd(1./dx);
  dy_i = _mm256_set1_pd(1./dy);
  hmoy2 = _mm256_set1_pd(hmoy);
  h_uphy1=_mm256_loadu_pd(&UPHY(t - 1, i, j));
  h_vphy1=_mm256_loadu_pd(&VPHY(t - 1, i, j));
  h_hfil1 =_mm256_loadu_pd(&HFIL(t - 1, i, j));

  c= _mm256_set1_pd(0.);

  if (i > 0)
    c =  _mm256_loadu_pd(&UPHY(t - 1, i - 1, j));

  d = _mm256_set1_pd(0.);

  if (j < size_y - 1)
    d = _mm256_loadu_pd(&VPHY(t - 1, i, j + 1));

  m=_mm256_mul_pd(dt1,hmoy2); //dt * hmoy
  sous=_mm256_sub_pd(h_uphy1,c); //(UPHY(t - 1, i, j) - c) 
  m2=_mm256_mul_pd(sous,dx_i); //((UPHY(t - 1, i, j) - c) / dx
  sous2=_mm256_sub_pd(d,h_vphy1); //(d - VPHY(t - 1, i, j))
  m4=_mm256_mul_pd(sous2,dy_i);  //(d - VPHY(t - 1, i, j)) / dy);
  add=_mm256_add_pd(m2,m4);  ////(d - VPHY(t - 1, i, j))+(d - VPHY(t - 1, i, j)) / dy)
  m3 = _mm256_mul_pd(m,add);
  
  
  //dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx + (d - VPHY(t - 1, i, j)) / dy);
  res=_mm256_sub_pd(h_hfil1,m3);


  return(res);
  //return HFIL(t - 1, i, j) -
    //dt * hmoy * ((UPHY(t - 1, i, j) - c) / dx +
     //(d - VPHY(t - 1, i, j)) / dy);
}

inline __m256d uPhy_forward(int t, int i, int j) {

  __m256d h_vphy, h_ufil, dissip2, h_hphy1, dx_i, m_grav, addo2, pcor1, quatre, m, m2, m3, mm, mm3, dt1,s, add, add1, add2, sous, amo, res, g, b, e, f;

  h_vphy=_mm256_loadu_pd(&VPHY(t - 1, i, j));
  h_ufil = _mm256_loadu_pd(&UFIL(t - 1, i, j));
  h_hphy1=_mm256_loadu_pd(&HPHY(t - 1, i, j));

  dt1 = _mm256_set1_pd(dt);
  dx_i= _mm256_set1_pd(1./dx);
  m_grav= _mm256_set1_pd(-grav);
  dissip2= _mm256_set1_pd(dissip);
  pcor1= _mm256_set1_pd(pcor);
  quatre=_mm256_set1_pd(1./4.);

  if (i == size_x - 1)
    return (_mm256_set1_pd(0.));

  b = _mm256_set1_pd(0.);
  if (i < size_x - 1)
    b = _mm256_loadu_pd(&HPHY(t - 1, i + 1, j));

  e = _mm256_set1_pd(0.);
  if (j < size_y - 1)
    e = _mm256_loadu_pd(&VPHY(t - 1, i, j + 1));

  f = _mm256_set1_pd(0.);
  if (i < size_x - 1)
    f = _mm256_loadu_pd(&VPHY(t - 1, i + 1, j));

  g = _mm256_set1_pd(0.);
  if (i < size_x - 1 && j < size_y - 1)
    g = _mm256_loadu_pd(&VPHY(t - 1, i + 1, j + 1));


  m=_mm256_mul_pd(m_grav,dx_i); //((-grav / dx)
  s=_mm256_sub_pd(b,h_hphy1); //(b - HPHY(t - 1, i, j))
  m3=_mm256_mul_pd(s,m);// (grav / dx) * (b - HPHY(t - 1, i, j))


  mm= _mm256_mul_pd(pcor1,quatre);//(pcor / 4.) 
  add=_mm256_add_pd(h_vphy,e); //VPHY(t - 1, i, j) + e 
  add1=_mm256_add_pd(add,f); //VPHY(t - 1, i, j) + e + f
  add2=_mm256_add_pd(add1,g); //(VPHY(t - 1, i, j) + e + f + g)
  mm3=_mm256_mul_pd(mm,add2);//(pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g)

  m2=_mm256_add_pd(mm3,m3);//  (grav / dx) * (b - HPHY(t - 1, i, j)) + (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g)

  res=_mm256_mul_pd(dissip2,h_ufil); //(dissip * UFIL(t - 1, i, j)));
  sous = _mm256_sub_pd(m2,res); // (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) - (dissip * UFIL(t - 1, i, j)));


  amo=_mm256_mul_pd(dt1,sous);   
 

  addo2=_mm256_add_pd(amo,h_ufil); //UFIL(t - 1, i, j) +dt * ((-grav / dx) * (b - HPHY(t - 1, i, j)) + (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) - (dissip * UFIL(t - 1, i, j)));


  return (addo2);
  //return UFIL(t - 1, i, j) +
    //dt * ((-grav / dx) * (b - HPHY(t - 1, i, j)) +
    //(pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) -
    //(dissip * UFIL(t - 1, i, j)));
}

inline __m256d vPhy_forward(int t, int i, int j) {
  

  __m256d h_hphy1,c, d, e, f, res, alpha1, alpha2, sous1, mult1, add1, add2, add, h_uphy1, mult2, h_vfil1, mult3;
 

  if (j == 0)
    return _mm256_set1_pd(0.);

  c = _mm256_set1_pd(0.);
  if (j > 0)
    c = _mm256_loadu_pd(&HPHY(t-1, i, j-1));

  d = _mm256_set1_pd(0.);
  if (i > 0 && j > 0)
    d = _mm256_loadu_pd(&UPHY(t-1, i-1, j-1));

  d = _mm256_set1_pd(0.);
  if (i > 0)
    e = _mm256_loadu_pd(&UPHY(t-1, i-1, j));

  d = _mm256_set1_pd(0.);
  if (j > 0)
    f = _mm256_loadu_pd(&UPHY(t-1, i, j-1));

  alpha1 = _mm256_set1_pd((-grav / dy)); 
  h_hphy1 = _mm256_loadu_pd(&HPHY(t-1, i, j));
  sous1 =  _mm256_sub_pd(h_hphy1,c);
  mult1 = _mm256_mul_pd(alpha1,sous1); // ((-grav / dy) * (HPHY(t - 1, i, j) - c)

  alpha2 = _mm256_set1_pd(-pcor / 4.); //(pcor / 4.)
  h_uphy1 = _mm256_loadu_pd(&UPHY(t-1, i, j));
  add1=_mm256_add_pd(d,e); //d + e 
  add2=_mm256_add_pd(f,h_uphy1); //f + UPHY(t - 1, i, j))
  add = _mm256_add_pd(add1,add2); //d+e+f + UPHY(t - 1, i, j))
  mult2 = _mm256_mul_pd(alpha2,add);//(pcor / 4.)*(d+e+f + UPHY(t - 1, i, j)))

  alpha1 = _mm256_set1_pd(-dissip);
  h_vfil1 = _mm256_loadu_pd(&VFIL(t-1, i, j));
  mult3 = _mm256_mul_pd(alpha1,h_vfil1);//(dissip * VFIL(t - 1, i, j)));

  add1 = _mm256_add_pd(mult2,mult3);
  add2 = _mm256_add_pd(mult1,add1);
  alpha1 = _mm256_set1_pd(dt); 
  mult1 = _mm256_mul_pd(alpha1,add2);
  res = _mm256_add_pd(mult1,h_vfil1);
  return(res);


  //return VFIL(t - 1, i, j) +
    //dt * ((-grav / dy) * (HPHY(t - 1, i, j) - c) -
	  //(pcor / 4.) * (d + e + f + UPHY(t - 1, i, j)) -
	  //(dissip * VFIL(t - 1, i, j)));
}


void forward(void) {
  __m256d  h_hfil, h_hphy, h_uphy, h_vphy,h_ufil, h_vfil; 
  FILE *file = NULL;
  double svdt = 0.;
  int t = 0;
  //h_hphy = hPhy_forward(1, 0, 0);
  //print256_num(h_hphy);
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
      for (int j = 0; j < size_y; j++) {
        h_hphy = hPhy_forward(t, i, j);
        _mm256_storeu_pd(&HPHY(t,i,j),h_hphy); 
        h_uphy = uPhy_forward(t, i, j);
        _mm256_storeu_pd(&UPHY(t,i,j),h_uphy); 
        h_vphy = vPhy_forward(t, i, j);
        _mm256_storeu_pd(&VPHY(t,i,j),h_vphy); 
        h_hfil = hFil_forward(t, i, j);
        _mm256_storeu_pd(&HFIL(t,i,j),h_hfil);
        h_ufil = uFil_forward(t, i, j);
        _mm256_storeu_pd(&UFIL(t,i,j),h_ufil);
        h_vfil = vFil_forward(t, i, j);
        _mm256_storeu_pd(&VFIL(t,i,j),h_vfil);
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
