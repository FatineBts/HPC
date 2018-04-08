#include <string>
extern double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy; // variables locales 
extern int size_x, size_y, nb_steps;
extern int g_size_x, g_size_y;
extern double *g_hFil;
extern double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
extern bool file_export;
extern std::string export_path;
extern int NP, rang; 

#define HFIL(t, i, j) hFil[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define UFIL(t, i, j) uFil[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define VFIL(t, i, j) vFil[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define HPHY(t, i, j) hPhy[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define UPHY(t, i, j) uPhy[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]
#define VPHY(t, i, j) vPhy[ (j) +			\
			    (i) * size_y +		\
			    ((t)%2) * size_x * size_y ]


#define G_HFIL(t, i, j) g_hFil[ (j) +			\
			    (i) * g_size_y +		\
			    ((t)%2) * g_size_x * g_size_y ]

