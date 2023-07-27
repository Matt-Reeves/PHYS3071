#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

int main(){

    // Basic Parameters
    int N = 101;
    double L = 1.0;
    double dx = L /((double) (N-1));
    
    double x[N];
    double phi[N][N]; 
    double f[N][N];
    double phi2[N][N];
    double phi_exact[N][N];


    for (int j =0; j<N; j++){
        x[j] = j*dx;
    }
    
    for (int j =0; j<N; j++){
        for (int k =0; k<N; k++){        
            phi[j][k] = 0.0;//0.5*sin(M_PI*x[j]) +0.5*sin(M_PI*x[k]);
            phi2[j][k] = 0.0;//0.5*sin(M_PI*x[j])+ 0.5*sin(M_PI*x[k]);
            f[j][k] = (2 + M_PI * M_PI * x[k] * (1 - x[k])) * sin(M_PI * x[j]) 
                     +(2 + M_PI * M_PI * x[j] * (1 - x[j])) * sin(M_PI * x[k]);

            phi_exact[j][k] = x[k]*(1- x[k])*sin(M_PI*x[j]) + x[j]*(1- x[j])*sin(M_PI*x[k]);
    }
    }


int Niter = 150000;
double SSE[Niter]; 
for(int j = 0; j<Niter; j++)  SSE[j] = 0.0;

for (int z = 0; z<Niter; z++){
    for (int j = 0; j<N; j++){
        for (int k = 0; k<N; k++){
            
                if (j ==0 || k == 0 || j == N-1 || k == N-1 )
                {   
                    phi[j][k] = 0.0; 
                }
                else
                {
                    phi[j][k] = 0.25*( phi[j+1][k] + phi[j-1][k] 
                             + phi[j][k+1] + phi[j][k-1] 
                             + dx*dx*f[j][k] );
                }
            SSE[z] += pow(phi[j][k] - phi_exact[j][k],2.0); 
            
        }
    }
    
}

std::ofstream saveFile, errFile;
  saveFile.precision(16);
  saveFile.open("GaussSeidel.txt");
  
  for (int j = 0; j<N; j++) for (int k=0; k<N; k++ ) saveFile << phi[j][k] << std::endl;
  for (int j = 0; j<N; j++) for (int k=0; k<N; k++ ) std::cout << f[j][k] << std::endl;
  
errFile.precision(16);
errFile.open("GSErr.txt");
for (int z =0; z<Niter; z++) errFile << SSE[z] <<std::endl;

  saveFile.close();
  errFile.close();

    
    





    return 0; 
}