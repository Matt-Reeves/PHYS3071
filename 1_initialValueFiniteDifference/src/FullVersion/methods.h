#ifndef METHODS_H
#define METHODS_H
    void imposeBoundaryConditions(int BC_type, int N, double dx, double* u, double* du_dt, double* U_t, double* U_tt );
    void imposeBoundaryConditionsLeapFrog(int BC_type, int N, double dx, double alpha2, double *u, double *v, double* temp1, double* temp2)
    void waveEquation(int BC_type, int N, double dx, double* u, double* du_dt, double* U_t, double* U_tt );
    void forwardEuler(int BC_type, int N, int Nt,  double dx,  double & t, double dt, double* u, double* du_dt);
    void rk4(int BC_type, int N, int Nt,  double dx,  double & t, double dt, double* u, double* du_dt);
    void leapFrog(int BC, int N, int Nt,double dx, double &t,double dt, double* u, double* v);
    double calculateEnergy(int N, double dx, double* u, double* v, double& KE, double & PE);
#endif

