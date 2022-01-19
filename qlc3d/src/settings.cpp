#include <settings.h>
#include <reader.h>
#include <iostream>
Settings::Settings():
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

void Settings::setnThreads(int num)         {
    nThreads = num;
}
void Settings::setQ_Solver(int s)           {
    Q_Solver = s;
}
void Settings::setV_Solver(int s)           {
    V_Solver = s;
}

int Settings::getnThreads() const           {
    return nThreads;
}
int Settings::getQ_Solver() const           {
    return Q_Solver;
}
int Settings::getV_Solver() const           {
    return V_Solver;
}

void    Settings::setQ_Newton_Panic_Iter(int i)     {
    Q_Newton_Panic_Iter = i;
}
void    Settings::setQ_Newton_Panic_Coeff(double c) {
    Q_Newton_Panic_Coeff = c;
}
int     Settings::getQ_Newton_Panic_Iter() const    {
    return Q_Newton_Panic_Iter;
}
double  Settings::getQ_Newton_Panic_Coeff() const   {
    return Q_Newton_Panic_Coeff;
}

//  Q-PCG SETTINGS
void Settings::setQ_PCG_Maxiter(int maxiter)    {
    Q_PCG_Maxiter = maxiter;
}
void Settings::setQ_PCG_Toler(double toler)     {
    Q_PCG_Toler = toler;
}
void Settings::setQ_PCG_Preconditioner(int pc)  {
    Q_PCG_Preconditioner = pc;
}

int Settings::getQ_PCG_Maxiter() const              {
    return Q_PCG_Maxiter;
}
int Settings::getQ_PCG_Preconditioner() const           {
    return Q_PCG_Preconditioner;
}
double Settings::getQ_PCG_Toler() const             {
    return Q_PCG_Toler;
}

// Q-GMRES SETTINGS
void Settings::setQ_GMRES_Preconditioner(int pc) {
    Q_GMRES_Preconditioner = pc;
}
void Settings::setQ_GMRES_Maxiter(int maxiter)  {
    Q_GMRES_Maxiter = maxiter;
}
void Settings::setQ_GMRES_Toler(double toler)   {
    Q_GMRES_Toler = toler;
}
void Settings::setQ_GMRES_Restart(int rest)     {
    Q_GMRES_Restart = rest;
}

int Settings::getQ_GMRES_Preconditioner() const     {
    return Q_GMRES_Preconditioner;
}
int Settings::getQ_GMRES_Maxiter() const                {
    return Q_GMRES_Maxiter;
}
int Settings::getQ_GMRES_Restart() const                {
    return Q_GMRES_Restart;
}
double Settings::getQ_GMRES_Toler() const               {
    return Q_GMRES_Toler;
}
//
//      POTENTIAL  --- V
//
//  V-PCG SETTINGS
void Settings::setV_PCG_Maxiter(int maxiter)    {
    V_PCG_Maxiter = maxiter;
}
void Settings::setV_PCG_Toler(double toler)     {
    V_PCG_Toler = toler;
}
void Settings::setV_PCG_Preconditioner(int pc)  {
    V_PCG_Preconditioner = pc;
}

int Settings::getV_PCG_Maxiter() const {
    return V_PCG_Maxiter;
}
int Settings::getV_PCG_Preconditioner() const {
    return V_PCG_Preconditioner;
}
double Settings::getV_PCG_Toler() const {
    return V_PCG_Toler;
}

// V-GMRES SETTINGS
void Settings::setV_GMRES_Preconditioner(int pc) {
    V_GMRES_Preconditioner = pc;
}
void Settings::setV_GMRES_Maxiter(int maxiter)  {
    V_GMRES_Maxiter = maxiter;
}
void Settings::setV_GMRES_Toler(double toler)   {
    V_GMRES_Toler = toler;
}
void Settings::setV_GMRES_Restart(int rest)     {
    V_GMRES_Restart = rest;
}

int Settings::getV_GMRES_Preconditioner() const {
    return V_GMRES_Preconditioner;
}
int Settings::getV_GMRES_Maxiter() const {
    return V_GMRES_Maxiter;
}
int Settings::getV_GMRES_Restart() const {
    return V_GMRES_Restart;
}
double Settings::getV_GMRES_Toler() const {
    return V_GMRES_Toler;
}
