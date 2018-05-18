#include <stdio.h>
#include <math.h>
#include <shalw.h>
#include <export.h>
#include <mpi.h>
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

void forward(int NP, int rang) {
  __m256d  h_hfil, h_hphy, h_uphy, h_vphy,h_ufil, h_vfil; 
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

    for (int j = 0; j < size_y; j++) {
      for (int i = 0; i <size_x; i++) {
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
