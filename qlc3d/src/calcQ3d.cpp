
#include <math.h>
#include <qlc3d.h>
#include <simulation-state.h>
#include <util/logging.h>
#include <util/hash.h>
#include <thread>

// SpaMtrix INCLUDES
#include <spamtrix_ircmatrix.hpp>
#include <spamtrix_vector.hpp>
#include <spamtrix_tickcounter.hpp>
#include <spamtrix_iterativesolvers.hpp>
#include <spamtrix_diagpreconditioner.hpp>

Simu::QMatrixSolvers selectQMatrixSolver(const Simu &simu, const LC &lc) {
    // SELECTS WHICH MATRIX SOLVER TO USE FOR Q-TENSOR
    // EQUATIONS. IF MATRIX IS SYMMETRIC USE PCG, ELSE
    // GMRES

    // IF SOLVER HAS ALREADY BEEN CHOSEN IN SETTINGS
    // FILE, DON'T DO ANYTHING
    if (simu.getQMatrixSolver()!=Simu::Auto) {
        return simu.getQMatrixSolver();
    }
    // SINGLE ELASTIC COEFF. EQUATIONS -> SYMMETRIC
    // CHIRALITY -> NON-SYMMETRIC
    bool isSymmetric = true;
    if (lc.K11() != lc.K22() ) { isSymmetric = false; }
    if (lc.K11() != lc.K33() ) { isSymmetric = false; }
    if (lc.p0() != 0.0) { isSymmetric = false; }

    return isSymmetric ? Simu::PCG : Simu::GMRES;
}

void solve_QTensor(SpaMtrix::IRCMatrix &K,
                   SpaMtrix::Vector &b,
                   SpaMtrix::Vector &x,
                   const Simu &simu,
                   const SolverSettings &settings,
                   const LC &lc) {
    SpaMtrix::DiagPreconditioner M(K);
    // Iterative solvers' settings...
    idx maxiter 	= settings.getQ_GMRES_Maxiter();
    idx restart 	= settings.getQ_GMRES_Restart();
    real toler      = settings.getQ_GMRES_Toler();
    SpaMtrix::IterativeSolvers solver(maxiter, restart, toler);
    Simu::QMatrixSolvers solverType = selectQMatrixSolver(simu, lc);

    if (solverType == Simu::PCG ) {
        if (!solver.pcg(K,x,b,M) ) {
            Log::warn("PCG did not converge in {} iterations. Tolerance achieved is {}.", solver.maxIter, solver.toler);
        }
    }
    else if (solverType == Simu::GMRES) {
        if (!solver.gmres(K,x,b,M)) {
            Log::warn("GMRES did not converge in {} iterations. Tolerance achieved is {}.", solver.maxIter, solver.toler);
        }
    }
}

void setThreadCount(unsigned int nt)
{
#ifdef NDEBUG
    omp_set_num_threads(nt);
#endif
}

/**
 * Updates the current Q-tensor solution. q = q + dq
 * @param q current Q-tensor solution
 * @param dq change in Q-tensor solution
 * @param maxdqOut return value magnitude of largest change in Q-tensor
 * @param simu
 * @return damping value, which may be < 1.0, if steady state solver solution change is too large (as configured in Simu)
 */
double updateSolutionVector(SolutionVector &q,
                          const SpaMtrix::Vector &dq,
                          double &maxdqOut,
                          const Simu &simu ) {
    // CALCULATE DAMPING COEFFICIENT IF GOING FOR STEADY STATE WITH NEWTON METHOD
    maxdqOut = fabs(dq[0]);
    double damping = 1.0;
    if (simu.simulationMode() == SteadyState) {
        // FIND MAXIMUM CHANGE
        for (idx i = 0 ; i < dq.getLength() ; i++) {
            maxdqOut = maxdqOut > fabs(dq[i]) ? maxdqOut : fabs(dq[i]);
        }

        if (maxdqOut > simu.getTargetdQ() ) { // IF DAMPING IS NEEDED
            damping = simu.getTargetdQ() / maxdqOut;
        }
    }

    const idx npLC = q.getnDoF();
    for (unsigned int i = 0 ; i < 5 ; i++) {   // LOOP OVER DIMENSIONS
        for (idx j = 0; j < npLC ; j++) { // LOOP OVER EACH NODE IN DIMENSION i
            const idx n = j + i*npLC;
            const idx effDoF = q.getEquNode(n);

            // EQUIVALENT DOF OF FIXED NODES ARE LABELLED AS "NOT_AN_INDEX"
            if (effDoF < NOT_AN_INDEX ) {
                const double dqj = dq[ effDoF ]; // TODO: should we multiply by "damping" here, it looks like it's never used !!
                //q.Values[n] -= dqj ;
                q[n] -= dqj;
                // KEEP TRACK OF LARGEST CHANGE IN Q-TENSOR
                maxdqOut = fabs(dqj) > maxdqOut ? fabs(dqj) : maxdqOut;
            }
        }
    }
    return damping;
}

double calcQ3d(SolutionVector *q,   // current Q-tensor
               SolutionVector *qn,  // previous Q-tensor
               SolutionVector *v,   // potential
               Geometry& geom,
               LC *mat_par,
               Simu* simu,
               SimulationState &simulationState,
               SpaMtrix::IRCMatrix &K,
               SolverSettings* settings,
               Alignment* alignment) {
    const idx numCols = K.getNumCols();
    double maxdq = 10;
    double maxdq_initial = 0;

    const idx isTimeStepping = simu->simulationMode() == TimeStepping ? 1:0;
    double timeStep = simulationState.dt();
    SpaMtrix::Vector RHS(isTimeStepping * numCols); // not used in steady-state solution
    SpaMtrix::Vector dq(numCols);
    SpaMtrix::Vector L(numCols);

    int newton_iter = 0;    // COUNTER FOR NEWTON ITERATIONS

    // SAVE Q FROM PREVIOUS TIME STEP
    if (simu->simulationMode() == TimeStepping) {
        qn->setValuesTo(*q);
    }

    // CREATE MILLISECOND ACCURACY TIMER
    TickCounter <std::chrono::milliseconds> timer;
    timer.start();

    bool LOOP = true;
    while ( LOOP ){ // Q-TENSOR SOLUTION LOOP
        newton_iter++;
        timer.reset();
        // SET NUMBER OF THREADS USED IN ASSEMBLY
        setThreadCount(simu->getAssemblyThreadCount() );

        // ASSEMBLE RHS CONTRIBUTION FROM PREVIOUS TIME-STEP, IF NEEDED
        if ( (isTimeStepping) && (newton_iter == 1) ){
            RHS = 0;
            assemble_prev_rhs(RHS, *qn, *v, *mat_par, timeStep, geom );
        }
        // CLEAR MATRIX, L AND dq
        K = 0;
        L = 0;
        dq = 0;

        //======================================
        //  MATRIX ASSEMBLY
        //======================================
#ifdef NDEBUG
      int numThreads = simu->getAssemblyThreadCount() == 0 ? (int) std::thread::hardware_concurrency() : (int) simu->getAssemblyThreadCount();
        omp_set_num_threads(numThreads);
#endif
        assembleQ(K, L, q, v, &geom.getTetrahedra(), &geom.getTriangles(), geom.getCoordinates(), mat_par, timeStep, alignment, geom.getNodeNormals());

#ifdef LOG_DEBUG_HASH
        int64_t lHash = hashCode64(&L[0], &L[0] + L.getLength());
        int64_t rhsHash = isTimeStepping ? hashCode64(&RHS[0], &RHS[0] + RHS.getLength()) : 0;
        int64_t vHash = v->hashCode();
        int64_t qHash = q->hashCode();
        int64_t qnHash = qn->hashCode();
        Log::info("lHash={:X}, rhsHash={:X}, vHash={:X}, qHash={:X}, qnHash{:}", lHash, rhsHash, vHash, qHash, qnHash);
#endif

        if (isTimeStepping) { // make Non-linear Crank-Nicholson RHS
            for (idx i = 0; i < numCols; i++) {
                L[i] += RHS[i];
            }
        }
        double matrixAssemblyTimeSeconds = (double) timer.getElapsed() / 1000.0;
        timer.reset();

        //======================================
        //  MATRIX SOLUTION
        //======================================
        setThreadCount(simu->getMatrixSolverThreadCount());
        // SOLVES Ax = b MATRIX PROBLEM
        solve_QTensor(K, L, dq, *simu, *settings, *mat_par);
        double damping = updateSolutionVector(*q, dq, maxdq, *simu); // q += damping*dq , taking into account periodic and fixed nodes
        double matrixSolutionTimeSeconds = (double) timer.getElapsed() / 1000.0;

        Log::info("Newton iteration {} for Q-tensor. Matrix assembly time = {}s. Matrix solution time = {}s. maxdq = {}.",
                  newton_iter, matrixAssemblyTimeSeconds, matrixSolutionTimeSeconds, maxdq);

        if (damping < 1.0) { // if damped, display by how much
            Log::warn("Damping by {}.", damping);
        }
/*
        // PANIC!! if looks like no convergence
        if (newton_iter > settings->getQ_Newton_Panic_Iter() ) {

            double newTimeStep = std::max(settings->getQ_Newton_Panic_Coeff() * simulationState.dt(), simu->getMindt());

            Log::warn("Newton iteration count {} exceeds threshold {}. Reducing time-step, new time-step is {}s.",
                      newton_iter, settings->getQ_Newton_Panic_Iter(), newTimeStep);
            simulationState.dt(newTimeStep);
            newton_iter = 0;
            q->setValuesTo(*qn);
        }

        // maxdq_initial is needed elsewhere to adjust time-step size. The time step is adjusted up or down in an
        // attempt to keep maxdq_initial at some constant value configured by the user.
        if (newton_iter == 1) {
            maxdq_initial = maxdq;
        }
        */
        // DETERMINE WHETHER NEWTON LOOP IS DONE
        if (!isTimeStepping) {
            LOOP = false;
        } else if ( fabs(maxdq) < simu->getMaxError() ) { // EXIT IF ACCURATE ENOUGH
            LOOP = false;
        }


    }
    return maxdq_initial;
}

