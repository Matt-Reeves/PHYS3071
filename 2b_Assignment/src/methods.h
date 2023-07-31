#ifndef _HASMETHODS
#define _HASMETHODS
#include <stdlib.h>
#include <complex.h>
#include "fftw3.h"

void fftshift(int N, double *k);
void schrodingerDynamics( bool groundState,
                         int Nt, int N, 
                         double& tnow, double dt, double dx, double g,
                         double* V, double* k, 
                         std::complex<double>* psi,
                         fftw_plan& p, fftw_plan& Pinv );

#endif