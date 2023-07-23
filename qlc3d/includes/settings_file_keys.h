#ifndef SETTINGS_FILE_KEYS_H
#define SETTINGS_FILE_KEYS_H

/*! All the settings file keys are declared and defined in this file
 * as const std::string variables.
 *
 * Start each variable name with SFK_
 * Keys are not case sensitive.
 *
*/

const char SFK_WILDCARD = '*';
const std::string SFK_MESH_NAME = "MeshName";
const std::string SFK_LOAD_Q    = "loadQ";
const std::string SFK_SAVE_DIR  = "saveDir";
const std::string SFK_Q_MATRIX_SOLVER = "QMatrixSolver";
const std::string SFK_SAVE_FORMAT = "SaveFormat";
const std::string SFK_END_CRITERION = "EndCriterion";
const std::string SFK_END_VALUE = "EndValue";
const std::string SFK_DT = "dt";
const std::string SFK_TARGET_DQ = "TargetdQ";
const std::string SFK_MAX_DT  = "Maxdt";
const std::string SFK_MAX_ERROR = "MaxError";
const std::string SFK_OUTPUT_ENERGY = "outputEnergy";
const std::string SFK_OUTPUT_FORMAT = "outputFormat";
const std::string SFK_SAVE_ITER ="SaveIter";
const std::string SFK_SAVE_TIME ="SaveTime";
const std::string SFK_NUM_ASSEMBLY_THREADS = "NumAssemblyThreads";
const std::string SFK_NUM_MATRIX_SOLVER_THREADS = "NumMatrixSolverThreads";
const std::string SFK_STRETCH_VECTOR = "stretchVector";
const std::string SFK_DT_LIMITS = "dtLimits";
const std::string SFK_DT_FUNCTION = "dtFunction";
const std::string SFK_REGULAR_GRID_SIZE = "RegularGridSize";
// LC Parameter keys
const std::string SFK_K11 = "K11";
const std::string SFK_K22 = "K22";
const std::string SFK_K33 = "K33";
const std::string SFK_K24 = "K24";
const std::string SFK_p0  = "p0";
const std::string SFK_A = "A";
const std::string SFK_B = "B";
const std::string SFK_C = "C";
const std::string SFK_EPS_PAR = "eps_par";
const std::string SFK_EPS_PER = "eps_per";
const std::string SFK_E1 = "e1";
const std::string SFK_E3 = "e3";
const std::string SFK_GAMMA1 = "gamma1";
// Electrodes and potential calculation keys
const std::string SFK_E_TIME = "E*.Time";
const std::string SFK_E_POTS = "E*.Pot";
const std::string SFK_E_FIELD = "EField";
const std::string SFK_EPS_DIELECTRIC = "eps_dielectric";
// Initial LC orientation boxes
const std::string SFK_BOX_TYPE = "Box*.Type";
const std::string SFK_BOX_PARAMS = "Box*.Params";
const std::string SFK_BOX_X = "Box*.X";
const std::string SFK_BOX_Y = "Box*.Y";
const std::string SFK_BOX_Z = "Box*.Z";
const std::string SFK_BOX_TILT = "Box*.Tilt";
const std::string SFK_BOX_TWIST = "Box*.Twist";
// Anchoring surfaces keys
const std::string SFK_FIXLC_ANCHORING = "FIXLC*.Anchoring";
const std::string SFK_FIXLC_STRENGTH = "FIXLC*.Strength";
const std::string SFK_FIXLC_EASY = "FIXLC*.Easy";
const std::string SFK_FIXLC_K1 = "FIXLC*.K1";
const std::string SFK_FIXLC_K2 = "FIXLC*.K2";
const std::string SFK_FIXLC_PARAMS = "FIXLC*.Params";

// Solver settings keys
const std::string SFK_Q_SOLVER = "Q_Solver";
const std::string SFK_V_SOLVER = "V_Solver";
const std::string SFK_Q_NEWTON_PANIC_ITER = "Q_Newton_Panic_Iter";
const std::string SFK_Q_NEWTON_PANIC_COEFF = "Q_Newton_Panic_Coeff";
const std::string SFK_Q_PCG_PECONDITIONER = "Q_PCG_Preconditioner";
const std::string SFK_Q_PCG_MAXITER = "Q_PCG_Maxiter";
const std::string SFK_Q_PCG_TOLER = "Q_PCG_Toler";
const std::string SFK_Q_GMRES_PRECONDITIONER = "Q_GMRES_Preconditioner";
const std::string SFK_Q_GMRES_MAXITER = "Q_GMRES_Maxiter";
const std::string SFK_Q_GMRES_RESTART = "Q_GMRES_Restart";
const std::string SFK_Q_GMRES_TOLER = "Q_GMRES_Toler";
const std::string SFK_V_PCG_PRECONDITIONER = "V_PCG_Preconditioner";
const std::string SFK_V_PCG_MAXITER = "V_PCG_Maxiter";
const std::string SFK_V_PCG_TOLER = "V_PCG_Toler";
const std::string SFK_V_GMRES_PRECONDITIONER = "V_GMRES_Preconditioner";
const std::string SFK_V_GMRES_MAXITER = "V_GMRES_Maxiter";
const std::string SFK_V_GMRES_RESTART = "V_GMRES_Restart";
const std::string SFK_V_GMRES_TOLER = "V_GMRES_Toler";

inline std::string wildcardToNum(const std::string& base, int num) {
    /*!Replaces wildcard character in input std::string with input
number, returning a copy*/
    size_t idx = base.find_first_of(SFK_WILDCARD);
    if (idx == std::string::npos) {
        std::cerr << "error in " << __PRETTY_FUNCTION__ << std::endl;
        std::exit(1);
    }
    std::string newStr = base.substr(0,idx);
    newStr += std::to_string(num);
    newStr += base.substr(idx+1 , std::string::npos);
    return newStr;
}

#endif
