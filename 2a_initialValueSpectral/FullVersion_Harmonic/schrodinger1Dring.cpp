
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex.h>
#include <random>
#include <iomanip>

#include "fftw3.h"

//User-defined functions used in main()
void fftshift(int N, double* k){
    
    for (int j =0; j<N/2; j++)
    {
        double temp1 = k[j];
        double temp2 = k[j+N/2];
        k[j] = temp2;
        k[j+N/2] = temp1;
    }

}

int main()
{
    //------------------------------- PROBLEM SETUP ----------------------------------------
    // Simulation parameters
    const int N = 256;                   // Number of points / Fourier modes
    const double L = 2.0*M_PI;           // Domain size
    
    const double dt = pow(2,-7);         // Small timestep for integration
    const double DT = 20;                // Large timestep for saving outputs
    const int    Nt = (int) (DT/dt);     // Number of steps per DT
    const double tf = 100.0;             // Final integration time
    double tnow = 0.0;                   // Curent simulation time  
    
    const std::complex<double> i (0.0, 1.0); // imaginary unit
    
    //------------------------------------  GRIDS ---------------------------------------
    //x-space grids
    const double dx = L / N;
    double x[N];

    // k-space grids
    const double dk = 2 * M_PI / L;
    const double kmax = 2 * M_PI / dx;
    double k[N];

    // wavefunctions in x- and -kspace
    std::complex<double> psi0[N];  // initial wavefunction
    std::complex<double>  psi[N];  // x-space
    std::complex<double> Fpsi[N];  // k-space

    //------------------------------------- FFTW SETUP -------------------------------
    // Setup for FFTs ...
    fftw_plan F;
    fftw_plan Finv;
    fftw_plan F2;
    fftw_complex *in = reinterpret_cast<fftw_complex *>(&psi[0]); // typecast to fftw_complex pointer to call fftw
    fftw_complex *out = reinterpret_cast<fftw_complex *>(&Fpsi[0]);
    
    std::cout << "Planning FFTs ... "; 
    F = fftw_plan_dft_1d(N, in, in, FFTW_FORWARD, FFTW_PATIENT);       //in-place  FFT for simulation 
    Finv = fftw_plan_dft_1d(N, in, in, FFTW_BACKWARD, FFTW_PATIENT);   //in-place IFFT for simiulation
    F2 = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);     //out-of-place FFT for analysis. 
    std::cout << " done." << std::endl;
    //-------------------------------------- OUTPUT ------------------------------------
    
    //Savefile for output
    std::ofstream saveFile, FsaveFile, gridFile;
    saveFile.open  ("xspace_psi.dat");   saveFile.precision(16);
    FsaveFile.open ("kspace_psi.txt");   FsaveFile.precision(16);
    gridFile.open  ("simGrids.txt");     gridFile.precision(16);
    
    //----------------------------- ARRAY INITIALIZATION -------------------------------

    // Initialize the arrays
    for (int j = 0; j < N; j++)
    {
        x[j] = -L/2.0 + j*dx;
        k[j] = -kmax/2.0 + j*dk;
        psi0[j] = std::complex<double> (exp(-0.5*(x[j]+a) * (x[j]+a)/l_ho/l_ho)  , 0.0) ;
        psi0[j] /=  pow((M_PI)*l_ho*l_ho,0.25) ;  // normalization 
        psi[j] =  psi0[j];
    }
    // rearrange the k-vector to the output format of forward fft
    fftshift(N,k);  // (e.g. [1 2 3 4 5 6] becomes [4 5 6 1 2 3])  
    
    //--------------------------- SIMULATION LOOP ---------------------------------  
    while (tnow < tf)
    {
            fftw_execute(F);
            for (int j = 0; j < N; j++){
                psi[j] *= std::exp(-i * k[j]*k[j] * 0.5 * DT);
            }
            fftw_execute(Finv);
            for (int j = 0; j < N; j++){
                psi[j] /= N; // renormalize -- (fftw computes unnormalized transform; forward then backward multiples by N)
            }
            tnow += DT; 
        }
        
        //Save to output
        for (int j = 0; j < N; j++) saveFile  << std::real(psi[j])  << " " ; saveFile << std::endl; 
        for (int j = 0; j < N; j++) saveFile  << std::imag(psi[j])  << " " ; saveFile  << std::endl; 
        
        std::cout << tnow << std::endl;
        
    }
    saveFile.close(); 
    
//-------------------------------------- OUPUT AND CLEANUP ------------------------------------  
    
    for (int j = 0; j < N; j++)
    {
        gridFile <<               V[j] << " " 
                 <<               k[j] << " " 
                 <<               x[j] << " " 
                 << std::endl;
    }
    gridFile.close();

    // Cleanup
    std::cout << "Simulation finished. " << std::endl;
    fftw_destroy_plan(F);
    fftw_destroy_plan(Finv);
    fftw_destroy_plan(F2);
    //fftw_free(in);
    //fftw_free(out);

    return 0;
}
