#include <stdlib.h>
#include <complex.h>
#include "fftw3.h"
#include <iostream>


void fftshift(int N, double *k)
{
    //Swap the left and right halves of k -- only works for N even. 
    for (int j = 0; j < N / 2; j++)
    {
        double temp1 = k[j];
        double temp2 = k[j + N / 2];
        k[j] = temp2;
        k[j + N / 2] = temp1;
    }
}

extern std::complex<double> i;
void schrodingerDynamics( bool groundState,
                         int Nt, int N, 
                         double& tnow, double dt, double dx, double g,
                         double* V, double* k, 
                         std::complex<double>* psi,
                         fftw_plan& p, fftw_plan& Pinv ){

    for (int z = 0; z < Nt; z++)
        {
            for (int j = 0; j < N; j++)
            {
                psi[j] *= std::exp( -i * (V[j] + g*std::norm(psi[j]) ) * dt * 0.5  );
            }
            fftw_execute(p);
            for (int j = 0; j < N; j++)
            {
                psi[j] *= std::exp( -i * k[j] * k[j] * 0.5 * dt );
            }
            fftw_execute(Pinv);
            for (int j = 0; j < N; j++){
                 psi[j] /= N; 
                 psi[j] *= std::exp( -i * (V[j] + g*std::norm(psi[j]) ) * dt * 0.5  );
             }
             tnow += dt;

             if (groundState == true)
             {
                double normalization = 0.0;
                for (int j =0; j<N; j++) 
                    normalization += std::norm(psi[j]);
                
                normalization *= dx; 
                normalization = sqrt(normalization);
                for (int j =0; j<N; j++)
                    psi[j] /= normalization;
             }
            
            
        }
}
