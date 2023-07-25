#ifndef METHODS_H
#define METHODS_H
    void imposeBoundaryConditions(int BC_type, int N, double dx, double* u, double* du_dt, double* U_t, double* U_tt );
    void waveEquation(int BC_type, int N, double dx, double* u, double* du_dt, double* U_t, double* U_tt );
    void forwardEuler(int BC_type, int N, int Nt,  double dx,  double & t, double dt, double* u, double* du_dt);
    //void rk4(int BC_type, int N, int Nt,  double dx,  double & t, double dt, double* u, double* du_dt);
#endif

