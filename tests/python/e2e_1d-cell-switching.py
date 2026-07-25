import os
import tempfile
import shutil
import math
from qlc3d_test_utils import *

# The expected result has a mid-cell tilt angle of 44.8007 degrees, obtained using Richard's 1D MATLAB solver
EXPECTED_MID_CELL_TILT_DEGREES = 44.8007


def check_results(result_dir):
    """
    Check the results match expected.
    :param result_dir: the directory where results are written to
    :return:
    """
    print("resultDir:" + result_dir)

    # Check that first and last result files exist in the results directory
    initial_result_file = 'dirstacksz00000000.csv'
    final_result_file = 'dirstacksz-final.csv'
    files = os.listdir(result_dir)

    assert_true(initial_result_file in files, "{} should contain the file {}".format(str(files), initial_result_file))
    assert_true(final_result_file in files, "{} should contain the file {}".format(str(files), final_result_file))

    # Read the final result file. The first line contains header info and the second line contains the
    # mid-cell director x, y, z values.
    fid = open(result_dir + "/" + final_result_file, "r")
    fid.readline()  # ignore header line
    director_string = fid.readline().split(",")  # read the director line
    fid.close()
    nz = float(director_string[2])  # z-component of the director
    tilt_degrees = 180. * math.asin(nz) / math.pi
    print("mid-cell tilt angle degrees = " + str(tilt_degrees))
    print("expected mid-cell tilt angle degrees = " + str(EXPECTED_MID_CELL_TILT_DEGREES))
    tilt_error = math.fabs(EXPECTED_MID_CELL_TILT_DEGREES - tilt_degrees)
    assert_true(tilt_error < 0.05, "mid-cell tilt error {} >= {}".format(str(tilt_error), str(0.05)))


def run_test(qlc3d_executable):
    """
    In this test, a potential is applied across a 1 micron thick 1D cell and the steady-state mid-cell tilt
    angle is measured and compared to a known (assumed) correct value.

    Material parameters are for 5CB and the cell has strong anchoring with 5 degrees pre-tilt on both surfaces.

    """
    print("current working directory=" + os.getcwd())
    # Setup
    # 1. create temp directory for project
    project_dir = tempfile.mkdtemp()

    # 2. copy mesh and settings file to project dir
    settings_file = project_dir + "/settings.txt"
    shutil.copy("./resources/settings-switch-1d.txt", settings_file)
    shutil.copy("./resources/thin.msh", project_dir + "/thin.msh")

    # 3. run qlc3d executable
    run_executable(qlc3d_executable, settings_file, project_dir)

    # 4. check results
    check_results(project_dir + '/res')

    # 5. delete the project directory, ignoring any errors
    shutil.rmtree(project_dir, True)


if __name__ == "__main__":
    print("args:" + str(sys.argv))
    run_test(sys.argv[1])
