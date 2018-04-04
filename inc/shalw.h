#include <string>
extern double *hFil, *uFil, *vFil, *hPhy, *uPhy, *vPhy; // variables locales 
extern int size_x, size_y, nb_steps;
extern int g_size_x, g_size_y;
extern double *g_hFil;//, *g_uFil, *g_vFil, *g_hPhy, *g_uPhy, *g_vPhy; // variables globales 
extern double dx, dy, dt, pcor, grav, dissip, hmoy, alpha, height, epsilon;
extern bool file_export;
extern std::string export_path;

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
#define G_UFIL(t, i, j) g_uFil[ (j) +			\
 			    (i) * g_size_y +		\
 			    ((t)%2) * g_size_x * g_size_y ]
#define G_VFIL(t, i, j) g_vFil[ (j) +			\
 			    (i) * g_size_y +		\
 			    ((t)%2) * g_size_x * g_size_y
#define G_HPHY(t, i, j) g_hPhy[ (j) +			\
			    (i) * g_size_y +		\
 			    ((t)%2) * g_size_x * g_size_y ]
#define G_UPHY(t, i, j) g_uPhy[ (j) +			\
 			    (i) * g_size_y +		\
 			    ((t)%2) * g_size_x * g_size_y ]
 #define G_VPHY(t, i, j) g_vPhy[ (j) +			\
 			    (i) * g_size_y +		\
 			    ((t)%2) * g_size_x * g_size_y ]

