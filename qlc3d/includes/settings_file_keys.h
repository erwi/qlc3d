#ifndef SETTINGS_FILE_KEYS_H
#define SETTINGS_FILE_KEYS_H

/*! All the settings file keys are declared and defined in this file
 * as const char variables.
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
#endif
