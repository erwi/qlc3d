#include <solver-settings.h>
#include <reader.h>
#include <iostream>
SolverSettings::SolverSettings():
    nThreads(DEFAULT_N_THREADS),
    Q_Solver(Q_Solver_GMRES),
    V_Solver(V_Solver_GMRES),
    Q_Newton_Panic_Iter(DEFAULT_Q_Newton_Panic_Iter),
    Q_Newton_Panic_Coeff(DEFAULT_Q_Newton_Panic_Coeff),

    Q_PCG_Preconditioner(Diagonal),
    Q_PCG_Maxiter(DEFAULT_Matrix_Maxiter),
    Q_PCG_Toler(DEFAULT_Matrix_Toler),

    Q_GMRES_Preconditioner(LUinc),
    Q_GMRES_Maxiter(DEFAULT_Matrix_Maxiter),
    Q_GMRES_Restart(DEFAULT_GMRES_Restart),
    Q_GMRES_Toler(DEFAULT_Matrix_Toler),

    V_PCG_Preconditioner(Cholinc),
    V_PCG_Maxiter(DEFAULT_Matrix_Maxiter),
    V_PCG_Toler(DEFAULT_Matrix_Toler),

    V_GMRES_Preconditioner(LUinc),
    V_GMRES_Maxiter(DEFAULT_Matrix_Maxiter),
    V_GMRES_Restart(DEFAULT_GMRES_Restart),
    V_GMRES_Toler(DEFAULT_Matrix_Toler) {
}

void SolverSettings::setnThreads(int num) {
  if (nThreads < 0) {
    throw std::runtime_error("Number of threads must be 0 or positive");
  }
  nThreads = num;
}
void SolverSettings::setQ_Solver(int s)           {
    Q_Solver = s;
}
void SolverSettings::setV_Solver(int s)           {
    V_Solver = s;
}

unsigned int SolverSettings::getnThreads() const           {
  if (nThreads < 0) {
    throw std::runtime_error("Number of threads must be 0 or positive");
  }
  return (unsigned int) nThreads;
}
int SolverSettings::getQ_Solver() const           {
    return Q_Solver;
}
int SolverSettings::getV_Solver() const           {
    return V_Solver;
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

//  Q-PCG SETTINGS
void SolverSettings::setQ_PCG_Maxiter(int maxiter)    {
    Q_PCG_Maxiter = maxiter;
}
void SolverSettings::setQ_PCG_Toler(double toler)     {
    Q_PCG_Toler = toler;
}
void SolverSettings::setQ_PCG_Preconditioner(int pc)  {
    Q_PCG_Preconditioner = pc;
}

int SolverSettings::getQ_PCG_Maxiter() const              {
    return Q_PCG_Maxiter;
}
int SolverSettings::getQ_PCG_Preconditioner() const           {
    return Q_PCG_Preconditioner;
}
double SolverSettings::getQ_PCG_Toler() const             {
    return Q_PCG_Toler;
}

// Q-GMRES SETTINGS
void SolverSettings::setQ_GMRES_Preconditioner(int pc) {
    Q_GMRES_Preconditioner = pc;
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

int SolverSettings::getQ_GMRES_Preconditioner() const     {
    return Q_GMRES_Preconditioner;
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
//
//      POTENTIAL  --- V
//
//  V-PCG SETTINGS
void SolverSettings::setV_PCG_Maxiter(int maxiter)    {
    V_PCG_Maxiter = maxiter;
}
void SolverSettings::setV_PCG_Toler(double toler)     {
    V_PCG_Toler = toler;
}
void SolverSettings::setV_PCG_Preconditioner(int pc)  {
    V_PCG_Preconditioner = pc;
}

int SolverSettings::getV_PCG_Maxiter() const {
    return V_PCG_Maxiter;
}
int SolverSettings::getV_PCG_Preconditioner() const {
    return V_PCG_Preconditioner;
}
double SolverSettings::getV_PCG_Toler() const {
    return V_PCG_Toler;
}

// V-GMRES SETTINGS
void SolverSettings::setV_GMRES_Preconditioner(int pc) {
    V_GMRES_Preconditioner = pc;
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

int SolverSettings::getV_GMRES_Preconditioner() const {
    return V_GMRES_Preconditioner;
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
