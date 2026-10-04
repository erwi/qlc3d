#ifndef SOLVER_H
#define SOLVER_H
#include <stdio.h>

const int DEFAULT_N_THREADS = 1;
const int DEFAULT_Q_Newton_Panic_Iter = 10;
const double DEFAULT_Q_Newton_Panic_Coeff = 0.1;
const int DEFAULT_Matrix_Maxiter = 2000;
const int DEFAULT_GMRES_Restart  = 100;
const double DEFAULT_Matrix_Toler = 1e-7;

class SolverSettings {
private:
    int         nThreads;   // Sets number of threads used by openMP
    int         Q_Newton_Panic_Iter;
    double      Q_Newton_Panic_Coeff;
    int         Q_GMRES_Maxiter;
    int         Q_GMRES_Restart;
    double      Q_GMRES_Toler;
    int         V_GMRES_Maxiter;
    int         V_GMRES_Restart;
    double      V_GMRES_Toler;
public:
    SolverSettings();
    void    setnThreads(int num);
    void    setQ_Newton_Panic_Iter(int i);
    void    setQ_Newton_Panic_Coeff(double c);
    [[nodiscard]] unsigned int getnThreads() const;
    int     getQ_Newton_Panic_Iter() const;
    double  getQ_Newton_Panic_Coeff() const;
    void setQ_GMRES_Maxiter(int maxiter);
    void setQ_GMRES_Toler(double toler);
    void setQ_GMRES_Restart(int restart);
    int     getQ_GMRES_Maxiter() const;
    int     getQ_GMRES_Restart() const;
    double  getQ_GMRES_Toler() const;
    void setV_GMRES_Maxiter(int maxiter);
    void setV_GMRES_Toler(double toler);
    void setV_GMRES_Restart(int restart);
    int     getV_GMRES_Maxiter() const;
    int     getV_GMRES_Restart() const;
    double  getV_GMRES_Toler() const;
};
#endif

