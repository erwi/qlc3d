#ifndef SETTINGS_FILE_KEYS_H
#define SETTINGS_FILE_KEYS_H

/*! All the settings file keys are declared and defined in this file
 * as const char variables.
 *
 * Start each variable name with SFK_
 * Keys are not case sensitive.
 *
*/
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



#endif
