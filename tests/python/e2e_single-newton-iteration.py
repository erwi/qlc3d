import sys
import os
import tempfile as tf
import shutil
from qlc3d_test_utils import *
'''
In this test a single iteration of Newton method is run with a potential turned on and the result file
contents are compared to expected values.  
'''
RESULT_BACKUP_CONFIG_FILE = 'settings.qfg'
RESULT_FILE_0 = 'dirstacksz00000000.csv'
RESULT_FILE_1 = 'dirstacksz00000001.csv'
RESULT_FILE_FINAL = 'dirstacksz-final.csv'
# first file line is header with info about grid. 1 x, 1 y, 3 z points
EXPECTED_HEADER_LINE = '1,1,3,0'
# second line has three director x,y,z components along a single column of 3 points
# values observed from a presumed correctly run. This is not the expected final correct result, just the first step result
# and may depend on many configuration settings
EXPECTED_DIRECTORS = [
    (0.998782, 0.0348782, 0.0348995),
    (0.958874,0.00504492,0.283788),
    (0.998782, -0.0348782, 0.0348995)
]

EPSILON = 1e-5

def check_results(result_dir):
    print("resultDir:" + result_dir)
    files = os.listdir(result_dir)
    print("resultDir contents =", str(files))
    assert_equals(4, len(files)) # expect 4 files in the result directory
    assert_true(RESULT_BACKUP_CONFIG_FILE in files, RESULT_BACKUP_CONFIG_FILE + " not found")
    assert_true(RESULT_FILE_0 in files, RESULT_FILE_0 + " not found")
    assert_true(RESULT_FILE_1 in files, RESULT_FILE_1 + " not found")
    assert_true(RESULT_FILE_FINAL in files, RESULT_FILE_FINAL + " not found")

    # Check the result file contents
    fid = open(result_dir + "/" + RESULT_FILE_1)
    count = 0
    for line in fid:
        if count == 0:
            assert_equals(EXPECTED_HEADER_LINE, line.rstrip())
        elif count == 1:
            # split the line into three tuples, one per director
            nums = [float(x) for x in line.rstrip().split(',')]
            directors = [tuple(nums[i:i+3]) for i in range(0, 9, 3)]

            for i, (d, expected) in enumerate(zip(directors, EXPECTED_DIRECTORS)):
                dot = sum(a * b for a, b in zip(d, expected))
                assert_true(abs(dot - 1.0) < EPSILON,
                            "director %d mismatch, dot product = %f" % (i, dot))
            break;
        count += 1
    assert_true(count == 1, "unexpected line count in result file, may have skipped assertions")
    fid.close()


def run_test(executable):
    print("current directory=" + os.getcwd())
    # setup
    # 1. create temp directory for project
    project_dir = tf.mkdtemp()

    # 2. copy mesh and settings file to project dir
    settings_file = project_dir + "/settings.txt"
    shutil.copy("./resources/settings.txt", settings_file)
    shutil.copy("./resources/thin.msh", project_dir + "/thin.msh")

    # 3. run qlc3d executable
    run_executable(executable, settings_file, project_dir)

    # 4. check results
    check_results(project_dir + '/res')


if __name__ == "__main__":
    print("args:" + str(sys.argv))
    run_test(sys.argv[1])
