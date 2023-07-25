#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>

#include "methods.h"

int main(){

//---------------------------- PROBLEM PARAMETERS -----------------------------
    int    N = 129;
    int BC_type = 1;         // type of boundary conditions --- see methods.cpp
    double L = 1.0;
    double dx = L / (N-1);
    double sigma = 0.05;
    double dt = pow(2.0, -14);
    double DT = pow(2.0,-6.0);
    int Nt = (int) (DT/dt);  
    double t = 0.0; 
    double tfinal = 20.0; 
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
std::ofstream saveFile; 
saveFile.precision(16);
saveFile.open("waveEquation.txt");
for (int j = 0; j<N; j++){ saveFile << u0[j] << " "; }
saveFile << std::endl;
//----------------------------- SIMULATION LOOP  ----------------------------
    double v[N];
    for (int j =0; j<N; j++) v[j] = du_dt[j];
    v[0] = v[N-1] = 0;

    for(int zz =0; zz<NT; zz++)
    {
        forwardEuler(BC_type,N,Nt, dx,t,dt, u,du_dt); 
        
        //Save Output
        std::cout << "t = " << t << std::endl;
        for (int j = 0; j<N; j++) saveFile << u[j] << " ";
        saveFile << std::endl; 
        
    }

    return 0; 
}