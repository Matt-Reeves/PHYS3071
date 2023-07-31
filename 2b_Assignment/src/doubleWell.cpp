
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex.h>

#include "methods.h"
#include "fftw3.h"

std::complex<double> i (0.0,1.0);

int main(int argc, char** argv)
{
    //------------------------------- PROBLEM SETUP ----------------------------------------
    // Simulation parameters
    const int N = 512;                  // Number of points / Fourier modes
    const double L = 30.0;              // Domain size
    const double V0 = 1.5;              // Well Depth
    const double sigma = 1.0;           // Well Width (always 1)
    const double a = 2;                 // Well Displacement
    const double l_ho = 1.0;
    const double g = 0.0;
    
    double dt = pow(2, -6);             // Small timestep for integration
    const int Nt = 100;                 // Number of steps per DT
    const double DT = Nt * dt;          // Large timestep for saving outputs
    const double tf = 1000.0;           // Final integration time
    double tnow = 0.0;                  // Curent simulation time

    //------------------------------------  GRIDS ---------------------------------------
    // x-space grids
    const double dx = L / N;
    double x[N];
    double V[N];
    double V2[N];

    // k-space grids
    const double dk = 2 * M_PI / L;
    const double kmax = 2 * M_PI / dx;
    double k[N];

    // wavefunctions in x- and -kspace
    std::complex<double> psi0[N]; // initial wavefunction
    std::complex<double> psi[N];  // x-space
    std::complex<double> Fpsi[N]; // k-space

    //------------------------------------- FFTW SETUP -------------------------------
    // Setup for FFTs ...
    fftw_plan p;
    fftw_plan Pinv;
    fftw_plan F;
    fftw_complex *in = reinterpret_cast<fftw_complex *>(&psi[0]); // typecast to fftw_complex pointer to call fftw
    fftw_complex *out = reinterpret_cast<fftw_complex *>(&Fpsi[0]);

    std::cout << "Planning FFTs ... " << std::flush;
    p = fftw_plan_dft_1d(N, in, in, FFTW_FORWARD, FFTW_PATIENT);     // in-place  FFT for simulation
    Pinv = fftw_plan_dft_1d(N, in, in, FFTW_BACKWARD, FFTW_PATIENT); // in-place IFFT for simiulation
    F = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);   // out-of-place FFT for analysis.
    std::cout << " done." << std::endl;
    //-------------------------------------- OUTPUT ------------------------------------

    // Savefile for output
    std::ofstream saveFile, gridFile, gFile;
    saveFile.open("./data/dynamics.dat");
    saveFile.precision(16);
    gFile.open("./data/groundState.dat");
    gFile.precision(16);
    gridFile.open("./data/simGrids.txt");
    gridFile.precision(16);

    // eFile.open("./data/expectationValues.txt");
    // eFile.precision(16);

    //----------------------------- ARRAY INITIALIZATION -------------------------------

    // Initialize the arrays
    for (int j = 0; j < N; j++)
    {
        x[j] = -L / 2.0 + j * dx;
        k[j] = -kmax / 2.0 + j * dk;
        psi0[j] = std::complex<double>(exp(-0.5 * (x[j] + a) * (x[j] + a) / l_ho / l_ho), 0.0);
        psi0[j] /= pow(M_PI * l_ho * l_ho, 0.25); // normalization  
        psi[j] = psi0[j];

        V[j] = -V0*exp(- (x[j] + a )*(x[j] + a)/sigma/sigma );
        V2[j] = -V0*(exp(- (x[j] + a )*(x[j] + a)/sigma/sigma ) 
                  + exp(- (x[j] - a )*(x[j] - a)/sigma/sigma ) );
    }
    // rearrange the k-vector to the output format of forward fft
    fftshift(N, k); // (e.g. [1 2 3 4 5 6] becomes [4 5 6 1 2 3])

    //--------------------------- SIMULATION LOOP -- Groundstate ---------------------------------
    
    i = std::complex<double>(1.0,0.0);
    std::cout << "Solving groundstate..." << std::flush ; 
    while (tnow < tf)
        schrodingerDynamics( true,  Nt,  N, tnow, dt,  dx,  g, V,  k,  psi, p, Pinv );  
    std::cout << "done." << std::endl;

    for (int j = 0; j < N; j++) gFile << std::norm(psi[j]) << " ";
        gFile << std::endl;
    gFile.close();
    //--------------------------- SIMULATION LOOP -- Dynamics ---------------------------------
    tnow = 0.0;
    i = std::complex<double> (0,1.0);
    std::cout << "Solving dynamics..." << std::flush ; 
    while (tnow < tf)
    {
        schrodingerDynamics(false,  Nt,  N, tnow, dt,  dx,  g, V2,  k,  psi, p, Pinv );
        // Save to output
        for (int j = 0; j < N; j++) saveFile << std::norm(psi[j]) << " ";
        saveFile << std::endl;
    }
    std::cout << "done." << std::endl;
    saveFile.close();
    //-------------------------------------- OUPUT AND CLEANUP ------------------------------------

    for (int j = 0; j < N; j++)
    {
        gridFile << V[j] << " "
                 << k[j] << " "
                 << x[j] << " "
                 << std::endl;
    }
    gridFile.close();

    // Cleanup
    std::cout << "Simulation finished. " << std::endl;
    fftw_destroy_plan(p);
    fftw_destroy_plan(Pinv);
    // fftw_free(in);
    // fftw_free(out);

    return 0;
}
