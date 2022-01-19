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
EXPECTED_RESULT_FILE_LINES = ['1,1,3,0',
                              '0.998782,0.0348782,0.0348995,0.958484,-6.13272e-08,0.285145,0.998782,-0.0348782,0.0348995']


def check_results(result_dir):
    print("resultDir:" + result_dir)
    files = os.listdir(result_dir)
    print("resultDir contents =", str(files))
    assert_equals(4, len(files)) # expect 4 files in the result directory
    assert_true(RESULT_BACKUP_CONFIG_FILE in files, RESULT_BACKUP_CONFIG_FILE + " not found")
    assert_true(RESULT_FILE_0 in files, RESULT_FILE_0 + " not found")
    assert_true(RESULT_FILE_1 in files, RESULT_FILE_1 + " not found")
    assert_true(RESULT_FILE_FINAL in files, RESULT_FILE_FINAL + " not found")

    # make sure the file contents are as expected, line-by-line
    fid = open(result_dir + "/" + RESULT_FILE_1)
    count = 0
    for line in fid:
        assert_equals(EXPECTED_RESULT_FILE_LINES[count], line.rstrip())
        count += 1
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
    command = executable + " " + settings_file + " " + project_dir
    print("command=" + command)
    sys.stdout.flush()
    exitCode = os.system(command)
    assert_equals(0, exitCode)
    print("end")
    sys.stdout.flush()

    # 4. check results
    check_results(project_dir + '/res')


if __name__ == "__main__":
    print("args:" + str(sys.argv))
    run_test(sys.argv[1])
