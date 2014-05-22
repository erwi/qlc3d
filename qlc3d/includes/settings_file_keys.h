#ifndef SETTINGS_FILE_KEYS_H
#define SETTINGS_FILE_KEYS_H

/*! All the settings file keys are declared and defined in this file
 * as const string variables.
 *
 * Start each variable name with SFK_
 * Keys are not case sensitive.
 *
*/
const char SFK_WILDCARD = '*';
const string SFK_MESH_NAME = "MeshName";
const string SFK_LOAD_Q    = "LoadQ";
const string SFK_SAVE_DIR  = "SaveDir";
const string SFK_Q_MATRIX_SOLVER = "QMatrixSolver";
const string SFK_SAVE_FORMAT = "SaveFormat";
const string SFK_END_VALUE = "EndValue";
const string SFK_DT = "dt";
const string SFK_TARGET_DQ = "TargetdQ";
const string SFK_MAX_DT  = "Maxdt";
const string SFK_MAX_ERROR = "MaxError";
const string SFK_OUTPUT_ENERGY = "OutputEnergy";
const string SFK_OUTPUT_FORMAT = "OutputFormat";
const string SFK_SAVE_ITER ="SaveIter";
const string SFK_NUM_ASSEMBLY_THREADS = "NumAssemblyThreads";
const string SFK_NUM_MATRIX_SOLVER_THREADS = "NumMatrixSolverThreads";
const string SFK_STRETCH_VECTOR = "StretchVector";
const string SFK_DT_LIMITS = "dtLimits";
const string SFK_DT_FUNCTION = "dtFunction";
const string SFK_REGULAR_GRID_SIZE = "RegularGridSize";
// LC Parameter keys
const string SFK_K11 = "K11";
const string SFK_K22 = "K22";
const string SFK_K33 = "K33";
const string SFK_p0  = "p0";
const string SFK_A = "A";
const string SFK_B = "B";
const string SFK_C = "C";
const string SFK_EPS_PAR = "eps_par";
const string SFK_EPS_PER = "eps_per";
const string SFK_E1 = "e1";
const string SFK_E3 = "e3";
const string SFK_GAMMA1 = "gamma1";
// Electrodes and potential calculation keys
const string SFK_E_TIME = "E*.Time";
const string SFK_E_POTS = "E*.Pot";
const string SFK_E_FIELD = "EField";
const string SFK_EPS_DIELECTRIC = "eps_dielectric";
// Initial LC orientation boxes
const string SFK_BOX_TYPE = "Box*.Type";
const string SFK_BOX_PARAMS = "Box*.Params";
const string SFK_BOX_X = "Box*.X";
const string SFK_BOX_Y = "Box*.Y";
const string SFK_BOX_Z = "Box*.Z";
const string SFK_BOX_TILT = "Box*.Tilt";
const string SFK_BOX_TWIST = "Box*.Twist";
// Anchoring surfaces keys
const string SFK_FIXLC_ANCHORING = "FIXLC*.Anchoring";
const string SFK_FIXLC_STRENGTH = "FIXLC*.Strength";
const string SFK_FIXLC_EASY = "FIXLC*.Easy";
const string SFK_FIXLC_K1 = "FIXLC*.K1";
const string SFK_FIXLC_K2 = "FIXLC*.K2";
const string SFK_FIXLC_PARAMS = "FIXLC*.Params";

// Solver settings keys
const string SFK_Q_SOLVER = "Q_Solver";
const string SFK_V_SOLVER = "V_Solver";
const string SFK_Q_NEWTON_PANIC_ITER = "Q_Newton_Panic_Iter";
const string SFK_Q_NEWTON_PANIC_COEFF = "Q_Newton_Panic_Coeff";
const string SFK_Q_PCG_PECONDITIONER = "Q_PCG_Preconditioner";
const string SFK_Q_PCG_MAXITER = "Q_PCG_Maxiter";
const string SFK_Q_PCG_TOLER = "Q_PCG_Toler";
const string SFK_Q_GMRES_PRECONDITIONER = "Q_GMRES_Preconditioner";
const string SFK_Q_GMRES_MAXITER = "Q_GMRES_Maxiter";
const string SFK_Q_GMRES_RESTART = "Q_GMRES_Restart";
const string SFK_Q_GMRES_TOLER = "Q_GMRES_Toler";
const string SFK_V_PCG_PRECONDITIONER = "V_PCG_Preconditioner";
const string SFK_V_PCG_MAXITER = "V_PCG_Maxiter";
const string SFK_V_PCG_TOLER = "V_PCG_Toler";
const string SFK_V_GMRES_PRECONDITIONER = "V_GMRES_Preconditioner";
const string SFK_V_GMRES_MAXITER = "V_GMRES_Maxiter";
const string SFK_V_GMRES_RESTART = "V_GMRES_Restart";
const string SFK_V_GMRES_TOLER = "V_GMRES_Toler";
#endif
