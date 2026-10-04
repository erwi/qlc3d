#include <solver-settings.h>
#include <reader.h>
#include <iostream>
SolverSettings::SolverSettings():
    nThreads(DEFAULT_N_THREADS),
    Q_Newton_Panic_Iter(DEFAULT_Q_Newton_Panic_Iter),
    Q_Newton_Panic_Coeff(DEFAULT_Q_Newton_Panic_Coeff),
    Q_GMRES_Maxiter(DEFAULT_Matrix_Maxiter),
    Q_GMRES_Restart(DEFAULT_GMRES_Restart),
    Q_GMRES_Toler(DEFAULT_Matrix_Toler),
    V_GMRES_Maxiter(DEFAULT_Matrix_Maxiter),
    V_GMRES_Restart(DEFAULT_GMRES_Restart),
    V_GMRES_Toler(DEFAULT_Matrix_Toler) {
}

void SolverSettings::setnThreads(int num) {
  if (num < 0) {
    throw std::runtime_error("Number of threads must be 0 or positive");
  }
  nThreads = num;
}

unsigned int SolverSettings::getnThreads() const           {
  if (nThreads < 0) {
    throw std::runtime_error("Number of threads must be 0 or positive");
  }
  return (unsigned int) nThreads;
}

void    SolverSettings::setQ_Newton_Panic_Iter(int i)     {
    Q_Newton_Panic_Iter = i;
}
void    SolverSettings::setQ_Newton_Panic_Coeff(double c) {
    Q_Newton_Panic_Coeff = c;
}
int     SolverSettings::getQ_Newton_Panic_Iter() const    {
    return Q_Newton_Panic_Iter;
}
double  SolverSettings::getQ_Newton_Panic_Coeff() const   {
    return Q_Newton_Panic_Coeff;
}

void SolverSettings::setQ_GMRES_Maxiter(int maxiter)  {
    Q_GMRES_Maxiter = maxiter;
}
void SolverSettings::setQ_GMRES_Toler(double toler)   {
    Q_GMRES_Toler = toler;
}
void SolverSettings::setQ_GMRES_Restart(int rest)     {
    Q_GMRES_Restart = rest;
}

int SolverSettings::getQ_GMRES_Maxiter() const                {
    return Q_GMRES_Maxiter;
}
int SolverSettings::getQ_GMRES_Restart() const                {
    return Q_GMRES_Restart;
}
double SolverSettings::getQ_GMRES_Toler() const               {
    return Q_GMRES_Toler;
}

void SolverSettings::setV_GMRES_Maxiter(int maxiter)  {
    V_GMRES_Maxiter = maxiter;
}
void SolverSettings::setV_GMRES_Toler(double toler)   {
    V_GMRES_Toler = toler;
}
void SolverSettings::setV_GMRES_Restart(int rest)     {
    V_GMRES_Restart = rest;
}

int SolverSettings::getV_GMRES_Maxiter() const {
    return V_GMRES_Maxiter;
}
int SolverSettings::getV_GMRES_Restart() const {
    return V_GMRES_Restart;
}
double SolverSettings::getV_GMRES_Toler() const {
    return V_GMRES_Toler;
}
