#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>

#include "methods.h"

int main(){

//---------------------------- PROBLEM PARAMETERS -----------------------------
    int    N = 513;
    int BC_type = 0;         // type of boundary conditions --- see methods.cpp
    double L = 1.0;
    double dx = L / (N-1);
    double sigma = 0.05;
    double dt = 0.5*dx; //pow(2.0, -20);
    double DT = pow(2.0,-6.0);
    int Nt = (int) (DT/dt);  
    double t = 0.0; 
    double tfinal = 2.0; 
    int NT = (int) (tfinal/DT);
    double x0 = 0.5;
//----------------------------- ARRAY INITIALIZATION ------------------------

    double x[N];     // Spatial Grid
    double u0[N];    // Initial Condition
    double u[N];     // Wave field
    double du_dt[N]; // Wave field time derivative

    for (int j = 0; j < N; j++)
    {
        x[j] = j * dx;
        u0[j] =  + exp(-(x[j]-x0) * (x[j]-x0) / sigma / sigma);
        u[j] = u0[j];
        du_dt[j] = 2*(x[j]-x0)/sigma/sigma*u0[j]; // time derivative
        
    }
//-------------------------------------- OUTPUT ------------------------------------
std::ofstream saveFile, energyFile; 
saveFile.precision(16);
saveFile.open("waveEquation.txt");
energyFile.precision(16);
energyFile.open("energy.txt");
for (int j = 0; j<N; j++){ saveFile << u0[j] << " "; }
saveFile << std::endl;
//----------------------------- SIMULATION LOOP  ----------------------------
    double v[N];
    for (int j =0; j<N; j++) v[j] = du_dt[j];
    v[0] = v[N-1] = 0;

    for(int zz =0; zz<NT; zz++)
    {
        //forwardEuler(BC_type,N,Nt, dx,t,dt, u,du_dt); 
        leapFrog(BC_type, N, Nt,dx, t, dt, u, v);
        //rk4(BC_type,N,Nt, dx,t,dt, u,du_dt);
        
        //Save Output
        std::cout << "t = " << t << std::endl;
        for (int j = 0; j<N; j++) saveFile << u[j] << " ";
        saveFile << std::endl; 
        
        // Calculate Energy --- uses RK4 but similar idea for leapfrog. 
        double KE = 0; double PE = 0;
        double E = calculateEnergy(N,dx,u,du_dt,KE,PE);
        energyFile << E << " " << KE << " " << " " << PE << std::endl;
    }

    return 0; 
}