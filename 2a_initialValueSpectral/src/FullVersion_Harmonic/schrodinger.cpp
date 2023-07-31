
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
    //Swap the left and right halves of k -- only works for N even. 
    for (int j = 0; j < N / 2; j++)
    {
        double temp1 = k[j];
        double temp2 = k[j + N / 2];
        k[j] = temp2;
        k[j + N / 2] = temp1;
    }
}

int main()
{
    //------------------------------- PROBLEM SETUP ----------------------------------------
    // Simulation parameters
    const int N = 512;       // Number of points / Fourier modes
    const double L = 20.0;    // Domain size
    const double omega = 0.5; // Oscillator frequency
    const double l_ho = 1.0;  // Oscillator length
    const double a = 0.0;     // Well Displacement

    const double dt = pow(2, -8); // Small timestep for integration
    const int Nt = 30;            // Number of steps per DT
    const double DT = Nt * dt;    // Large timestep for saving outputs
    const double tf = 100.0;      // Final integration time
    double tnow = 0.0;            // Curent simulation time

    const std::complex<double> i(0.0, 1.0); // imaginary unit

    //------------------------------------  GRIDS ---------------------------------------
    // x-space grids
    const double dx = L / N;
    double x[N];
    double V[N];

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

    std::cout << "Planning FFTs ... ";
    p = fftw_plan_dft_1d(N, in, in, FFTW_FORWARD, FFTW_PATIENT);     // in-place  FFT for simulation
    Pinv = fftw_plan_dft_1d(N, in, in, FFTW_BACKWARD, FFTW_PATIENT); // in-place IFFT for simiulation
    F = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);   // out-of-place FFT for analysis.
    std::cout << " done." << std::endl;
    //-------------------------------------- OUTPUT ------------------------------------

    // Savefile for output
    std::ofstream saveFile, gridFile, eFile;
    saveFile.open("xspace_psi.dat");
    saveFile.precision(16);
    gridFile.open("simGrids.txt");
    saveFile.precision(16);
    eFile.open("expectationValues.txt");
    eFile.precision(16);

    //----------------------------- ARRAY INITIALIZATION -------------------------------

    // Initialize the arrays
    for (int j = 0; j < N; j++)
    {
        x[j] = -L / 2.0 + j * dx;
        k[j] = -kmax / 2.0 + j * dk;
        psi0[j] = std::complex<double>(exp(-0.5 * (x[j] + a) * (x[j] + a) / l_ho / l_ho), 0.0);
        psi0[j] /= pow((M_PI)*l_ho * l_ho, 0.25); // normalization
        psi[j] = psi0[j];
        V[j] = 0.5 * x[j] * x[j] * omega * omega;
    }
    // rearrange the k-vector to the output format of forward fft
    fftshift(N, k); // (e.g. [1 2 3 4 5 6] becomes [4 5 6 1 2 3])

    //--------------------------- SIMULATION LOOP ---------------------------------
    while (tnow < tf)
    {
        for (int z = 0; z < Nt; z++)
        {
            for (int j = 0; j < N; j++)
            {
                psi[j] *= std::exp(-i * V[j] * dt);
            }
            fftw_execute(p);
            for (int j = 0; j < N; j++)
            {
                psi[j] *= std::exp(-i * k[j] * k[j] * 0.5 * dt);
            }
            for (int j = 0; j < N; j++){
                 psi[j] /= N; 
             }
            fftw_execute(Pinv);
            tnow += dt;

            // Momentum space wavefunction for analysis;
            fftw_execute(F);
            double k_norm = 0.0;
            for (int j = 0; j < N; j++)
                k_norm += std::norm(Fpsi[j]);
            k_norm *= dk;
            for (int j = 0; j < N; j++)
                Fpsi[j] /= sqrt(k_norm);

            double X = 0.0, P = 0.0;
            double T = 0.0, V = 0.0;
            for (int j = 0; j < N; j++)
            {
                V += 0.5 * std::norm(psi[j]) * x[j] * x[j] * omega * omega * dx; //p^2/2 is kinetic energy
                T += 0.5 * std::norm(Fpsi[j]) * k[j] * k[j] * dk; //x^2/2 is potential energy
                X += std::norm(psi[j]) * x[j] * dx;
                P += std::norm(Fpsi[j]) * k[j] * dk;
            }
            eFile << T << " " << V << " " << X << " " << P << std::endl;
        }

        // Save to output
        for (int j = 0; j < N; j++) saveFile << std::norm(psi[j]) << " ";
        saveFile << std::endl;
        
        std::cout << tnow << std::endl;
    }
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
