
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex.h>

#include "fftw3.h"

// User-defined functions used in main()
void fftshift(int N, double *k)
{
    //fftshift flips the left and right halves of  a vector,
    // to be consistent with the outout of a forward fft computed by fftw. 
    if (N % 2 != 0)
    {
        perror("Error, N must be even");
        abort();
    }
    double temp;
    for (int j = 0; j < N / 2; j++)
    {
        temp = k[j];
        k[j] = k[j + N / 2];
        k[j + N / 2] = temp;
    }
}

int main()
{
    //------------------------------- PROBLEM SETUP ----------------------------------------
    // Simulation parameters

    const int N = 256;              // Number of points / Fourier modes
    const double L = 2.0 * M_PI;    // Domain size
    const double sigma = 0.1;       // Wavepacket Width
    const double DT = pow(2, -8);   // Timestep for saving outputs
    const double tf = 0.25;         // Final integration time
    double tnow = 0.0;              // Current simulation time

    const std::complex<double> i(0.0, 1.0); // imaginary unit

    //------------------------------------  GRIDS ---------------------------------------
    // x-space grids
    const double dx = L / N;
    double x[N];

    // k-space grids
    const double dk = 2 * M_PI / L;
    const double k_range = 2 * M_PI / dx;
    double k[N];

    // wavefunctions in x- and -kspace
    std::complex<double> psi0[N]; // initial wavefunction
    std::complex<double> psi[N];  // x-space
    std::complex<double> Fpsi[N]; // k-space

    //------------------------------------- FFTW SETUP -------------------------------
    // Setup for FFTs ...
    fftw_plan F;
    fftw_plan Finv;
    fftw_plan F2;
    fftw_complex *in = reinterpret_cast<fftw_complex *>(&psi[0]); // typecast to fftw_complex pointer to call fftw
    fftw_complex *out = reinterpret_cast<fftw_complex *>(&Fpsi[0]);

    std::cout << "Planning FFTs ... ";
    F = fftw_plan_dft_1d(N, in, in, FFTW_FORWARD, FFTW_PATIENT);     // in-place  FFT for simulation
    Finv = fftw_plan_dft_1d(N, in, in, FFTW_BACKWARD, FFTW_PATIENT); // in-place IFFT for simiulation
    F2 = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);  // out-of-place FFT for analysis.
    std::cout << " done." << std::endl;

    //----------------------------- ARRAY INITIALIZATION -------------------------------

    // Initialize the arrays
    for (int j = 0; j < N; j++)
    {
        x[j] = -L / 2.0 + j * dx;
        k[j] = -k_range / 2.0 + j * dk;
        psi0[j] = std::complex<double>(exp(-0.5 * x[j] * x[j] / sigma / sigma), 0.0);
        psi0[j] /= pow((M_PI)*sigma * sigma, 0.25); // normalization
        psi[j] = psi0[j];
    }
    // rearrange the k-vector to the output format of forward fft
    fftshift(N, k); // (e.g. [1 2 3 4 5 6] becomes [4 5 6 1 2 3])

    //-------------------------------------- OUTPUT ------------------------------------

    // Savefile for output
    std::ofstream saveFile;
    saveFile.precision(16);
    std::string saveDir = "./data/";
    std::string fileType = ".dat";
    std::string filePrefix = "output";
    std::string fileName = saveDir + filePrefix + fileType;

    // Write initial state
    saveFile.open(fileName);
    for (int j = 0; j < N; j++)
        saveFile << std::norm(psi0[j]) << " ";
    saveFile << std::endl;

    //--------------------------- SIMULATION LOOP ---------------------------------

    while (tnow < tf)
    {
        fftw_execute(F);
        for (int j = 0; j < N; j++)
        {
            psi[j] *= std::exp(-i * k[j] * k[j] * 0.5 * DT);
        }
        fftw_execute(Finv);
        for (int j = 0; j < N; j++)
        {
            psi[j] /= N; /* renormalize -- (fftw is unnormalized ; forward then backward multiples by N) */
        }
        tnow += DT;

        //--------------------------------------  OUTPUT ------------------------------------
        // Save to output

        for (int j = 0; j < N; j++)
            saveFile << std::norm(psi[j]) << " ";
        saveFile << std::endl;

        std::cout << tnow << std::endl;
    }
    std::cout << "Simulation finished. " << std::endl;

    //--------------------------------------  CLEANUP ------------------------------------
    // Cleanup
    saveFile.close();
    fftw_destroy_plan(F);
    fftw_destroy_plan(Finv);
    fftw_destroy_plan(F2);

    // fftw_free(in); //cleanup if using fftw_malloc
    // fftw_free(out);

    return 0;
}
