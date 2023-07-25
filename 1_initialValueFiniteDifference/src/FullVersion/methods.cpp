#include "methods.h"
#include <iostream>

void imposeBoundaryConditionsLeapFrog(int BC_type, int N, double dx, double alpha2,
                                      double *u, double *v, double *temp1, double *temp2)
{
    int L = 0, R = N - 1;
    switch (BC_type)
    {
    case 0: // Periodic Boundary Conditions
        temp1[L] = u[L];
        temp2[L] = -v[L] + 2 * (1 - alpha2) * u[L] + alpha2 * (u[L + 1] + u[R]); 

        temp1[R] = u[R];
        temp2[R] = -v[R] + 2 * (1 - alpha2) * u[R] + alpha2 * (u[R - 1] + u[L]); 
        break;
    case 1: // Dirichlet Boundary Conditions

        temp1[L] = temp1[R] = 0.0;
        temp2[L] = temp2[R] = 0.0;
        break;

    case 2: // "Neumann Boundary Conditions"
        
        //if (t< dt){
        // temp1[L] = u[L];
        // temp2[L] = dt * v[L] + (1 - alpha2) * u[L] +  alpha2 * (u[L + 1]);

        // temp1[R] = u[R];
        // temp2[R] = dt * v[R] + (1 - alpha2) * u[R] +  alpha2 * (u[R - 1]);
        //}
        
        temp1[L] = u[L];
        temp2[L] = -v[L] + 2 * (1 - alpha2) * u[L] + 2.0 * alpha2 * (u[L + 1]); // alpha2*(u[L+1] + uL);

        temp1[R] = u[R];
        temp2[R] = -v[R] + 2 * (1 - alpha2) * u[R] + 2.0 * alpha2 * (u[R - 1]); //+ alpha2*(uR + u[R-1]);
        break;
    case 3: // Mixed Boundary Conditions
        temp1[L] = temp2[L] = 0.0;

        temp1[R] = u[R];
        temp2[R] = -v[R] + 2 * (1 - alpha2) * u[R] + 2.0 * alpha2 * (u[R - 1]);
        break;
    }
    return;
}

void imposeBoundaryConditions(int BC_type, int N, double dx,
                              double *u, double *du_dt, double *U_t, double *U_tt)
{
    int L = 0, R = N - 1;
    switch (BC_type)
    {
    case 0: // Periodic Boundary Conditions
        U_t[L] = du_dt[L];
        U_t[R] = du_dt[R];

        U_tt[L] = (u[R] - 2 * u[L] + u[L + 1]) / dx / dx;
        U_tt[R] = (u[R - 1] - 2 * u[R] + u[L]) / dx / dx;
        break;
    case 1: // Dirichlet Boundary Conditions
        U_t[L] = 0;
        U_t[R] = 0;

        U_tt[L] = 0;
        U_tt[R] = 0;
        break;
    case 2: // "Neumann Boundary Conditions"
        U_t[L] = du_dt[L];
        U_t[R] = du_dt[R];

        U_tt[L] = (2 * u[L + 1] - 2 * u[L]) / dx / dx;
        U_tt[R] = (2 * u[R - 1] - 2 * u[R]) / dx / dx;
        break;
    case 3: // Mixed Boundary Conditions
        U_t[L] = 0;
        U_t[R] = du_dt[R];

        U_tt[L] = 0;
        U_tt[R] = (2 * u[R - 1] - 2 * u[R]) / dx / dx;
        break;
    }
    return;
}

void waveEquation(int BC, int N, double dx,
                  double *u, double *du_dt, double *U_t, double *U_tt)
{
    // Dynamics of the interior points ...
    for (int j = 1; j < N - 1; j++)
    {
        U_t[j] = du_dt[j];
        U_tt[j] = (u[j - 1] - 2 * u[j] + u[j + 1]) / dx / dx;
    }
    //..and the boundaries.
    imposeBoundaryConditions(BC, N, dx, u, du_dt, U_t, U_tt);
    return;
}

void forwardEuler(int BC, int N, int Nt,
                  double dx, double &t, double dt, double *u, double *du_dt)
{
    // time derivatives of u and du_dt for integration scheme
    double U_t[N];
    double U_tt[N];

    for (int k = 0; k < Nt; k++)
    {
        waveEquation(BC, N, dx, u, du_dt, U_t, U_tt);
        for (int j = 0; j < N; j++)
        {
            u[j] += dt * U_t[j];
            du_dt[j] += dt * U_tt[j];
        }
    }
    t += (Nt * dt);
    return;
}

void leapFrog(int BC, int N, int Nt,
              double dx, double &t, double dt, double *u, double *v)
{

    double temp1[N], temp2[N];
    double alpha2 = dt * dt / dx / dx;

    for (int k = 0; k < Nt; k++)
    {
        if (t < dt)
        {
            for (int j = 1; j < N - 1; j++)
            {
                temp1[j] = u[j];
                temp2[j] = dt * v[j] + (1 - alpha2) * u[j] + 0.5 * alpha2 * (u[j - 1] + u[j + 1]);
                /*strictly the BC function is not yet correctly implemented for the first step
                 -- but doesn't matter for a Gaussian. */
                imposeBoundaryConditionsLeapFrog(BC, N, dx, alpha2, u, v, temp1, temp2); 
            }
            for (int j = 0; j < N; j++)
            {
                v[j] = temp1[j];
                u[j] = temp2[j];
            }
        }
        else
        {
            for (int j = 1; j < N - 1; j++)
            {
                temp1[j] = u[j];
                temp2[j] = -v[j] + 2 * (1 - alpha2) * u[j] + alpha2 * (u[j + 1] + u[j - 1]);
                imposeBoundaryConditionsLeapFrog(BC, N, dx, alpha2, u, v, temp1, temp2);
            }
            for (int j = 0; j < N; j++)
            {
                v[j] = temp1[j];
                u[j] = temp2[j];
            }
        }
        t += dt;
    }
}

void rk4(int BC, int N, int Nt,
         double dx, double &t, double dt, double *u, double *v)
{

    for (int k = 0; k < Nt; ++k)
    {
        // Note:  v = du_dt;
        double k1[N], k2[N], k3[N], k4[N];
        double p1[N], p2[N], p3[N], p4[N];
        double u_temp[N], v_temp[N];
        waveEquation(BC, N, dx, u, v, k1, p1);
        for (int j = 0; j < N; j++)
        {
            u_temp[j] = u[j] + 0.5 * dt * k1[j];
            v_temp[j] = v[j] + 0.5 * dt * p1[j];
        }
        waveEquation(BC, N, dx, u_temp, v_temp, k2, p2);
        for (int j = 0; j < N; j++)
        {
            u_temp[j] = u[j] + 0.5 * dt * k2[j];
            v_temp[j] = v[j] + 0.5 * dt * p2[j];
        }
        waveEquation(BC, N, dx, u_temp, v_temp, k3, p3);
        for (int j = 0; j < N; j++)
        {
            u_temp[j] = u[j] + dt * k3[j];
            v_temp[j] = v[j] + dt * p3[j];
        }
        waveEquation(BC, N, dx, u_temp, v_temp, k4, p4);
        for (int j = 0; j < N; j++)
        {
            u[j] += dt / 6.0 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
            v[j] += dt / 6.0 * (p1[j] + 2 * p2[j] + 2 * p3[j] + p4[j]);
        }
    }
    t += (Nt * dt);
    return;
}

double calculateEnergy(int N, double dx, double *u, double *v, double &KE, double &PE)
{
    // Calculates the energy using trapezoidal rule -- currently this only works for Dirichlet and Euler / RK4 steppers

    double E = 0.0;

    // Interior Points
    for (int j = 1; j < N - 1; j++)
    {
        KE += v[j] * v[j];
        double ux = (u[j + 1] - u[j - 1]) / 2.0 / dx;
        PE += ux * ux;
    }

    // End Values are weighted by 1/2 compared to interior points for trapz ..
    KE += v[0] * v[0] / 2.0;
    KE += v[N - 1] * v[N - 1] / 2.0; // these should be zero anyway for Dirichlet

    double ux0, uxL;
    ux0 = 2 * u[1] / 2.0 / dx; // happens to work out same as doing first-order differences
    uxL = 2 * u[N - 2] / 2.0 / dx;
    PE += ux0 * ux0 / 2.0;
    PE += uxL * uxL / 2.0;

    PE *= 0.5 * dx;
    KE *= 0.5 * dx;
    E = KE + PE;

    return E;
}