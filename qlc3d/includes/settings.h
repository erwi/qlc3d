
#ifndef SOLVER_H
#define SOLVER_H

#include <stdio.h>

#define Q_SOLVER_PCG		0
#define Q_SOLVER_GMRES		1

#define V_SOLVER_PCG		0
#define V_SOLVER_GMRES		1

#define DIAG_PRECONDITIONER 0
#define IC_PRECONDITIONER   1
#define ILU_PRECONDITIONER  2

class Settings
{
	private:
		
		
		int 		nThreads;	// Sets number of threads used by openMP
		int 		Q_Solver;	// Selects which solver to use for Q-tensor solution
		int 		V_Solver;	// Selects which solver to use for potential solution
//
// QSOLVER
//		
		int 		Q_Newton_Panic_Iter;
		double		Q_Newton_Panic_Coeff;
	
	// Preconditioned Conjugate Gradient Q-tensor solver settings
		int 		Q_PCG_Preconditioner;
		int 		Q_PCG_Maxiter;
		double		Q_PCG_Toler;
		
	// GMRES Q-tensor solver settings	
		int 		Q_GMRES_Preconditioner;
		int 		Q_GMRES_Maxiter;
		int			Q_GMRES_Restart;
		double		Q_GMRES_Toler;
//
// VSOLVER
//		
	// Preconditioned Conjugate Gradient potential solver settings
		int 		V_PCG_Preconditioner;
		int 		V_PCG_Maxiter;
		double		V_PCG_Toler;
		
	// GMRES potential solver settings	
		int 		V_GMRES_Preconditioner;
		int 		V_GMRES_Maxiter;
		int			V_GMRES_Restart;
		double		V_GMRES_Toler;
		
	
	public:
		Settings();
		void PrintSettings(); // Prints the settings
		
		void 	setnThreads(int num);
		void 	setQ_Solver(int s);
		void	setV_Solver(int s);
		void 	setQ_Newton_Panic_Iter(int i);
		void	setQ_Newton_Panic_Coeff(double c);
		
		int 	getnThreads();
		int 	getQ_Solver();
		int 	getV_Solver();
		int 	getQ_Newton_Panic_Iter();
		
		double  getQ_Newton_Panic_Coeff();
		
	// Preconditioned Conjugate Gradient Q-tensor solver settings	
		void setQ_PCG_Maxiter(int maxiter);
		void setQ_PCG_Toler(double toler);
		void setQ_PCG_Preconditioner(int pc);
		
		int  	getQ_PCG_Maxiter();
		int		getQ_PCG_Preconditioner();
		double	getQ_PCG_Toler();
		
	// GMRES Q-tensor solver settings	
		void setQ_GMRES_Maxiter(int maxiter);
		void setQ_GMRES_Toler(double toler);
		void setQ_GMRES_Preconditioner(int pc);
		void setQ_GMRES_Restart(int restart);
		
		int  	getQ_GMRES_Maxiter();
		int     getQ_GMRES_Restart();
		int		getQ_GMRES_Preconditioner();
		double	getQ_GMRES_Toler();

//
//		Potential solvers
//
		
	// Preconditioned Conjugate Gradient potential solver settings	
		void setV_PCG_Maxiter(int maxiter);
		void setV_PCG_Toler(double toler);
		void setV_PCG_Preconditioner(int pc);
		
		int  	getV_PCG_Maxiter();
		int		getV_PCG_Preconditioner();
		double	getV_PCG_Toler();
		
	// GMRES potential solver settings	
		void setV_GMRES_Maxiter(int maxiter);
		void setV_GMRES_Toler(double toler);
		void setV_GMRES_Preconditioner(int pc);
		void setV_GMRES_Restart(int restart);
		
		int  	getV_GMRES_Maxiter();
		int     getV_GMRES_Restart();
		int		getV_GMRES_Preconditioner();
		double	getV_GMRES_Toler();
};


#endif

