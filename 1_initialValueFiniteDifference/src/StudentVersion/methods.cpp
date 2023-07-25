#include "methods.h"

void imposeBoundaryConditions(int BC_type, int N, double dx,
                              double *u, double *du_dt, double *U_t, double *U_tt)
{
    int L = 0, R = N - 1;
    switch (BC_type)
    {
    case 0: // Periodic Boundary Conditions
        /* Question 2: */
        break;
    case 1: // Dirichlet Boundary Conditions
        U_t[L] = 0;
        U_t[R] = 0;

        U_tt[L] = 0;
        U_tt[R] = 0;
        break;
    case 2: // "Neumann Boundary Conditions"
        /* Question 2: */
        break;
    case 3: // Mixed Boundary Conditions
        /* Question 2: */
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

