import sys
import os
import tempfile as tf
import shutil
'''
In this test a single iteration of Newton method is run with a potential turned on and the result file
contents are compared to expected values.  
'''
RESULT_BACKUP_CONFIG_FILE = 'settings.qfg'
RESULT_FILE = 'dirstackz1.csv'
EXPECTED_RESULT_FILE_LINES = ['1,1,3,0',
                              '0.998782,0.0348782,0.0348995,0.958484,-4.57858e-08,0.285145,0.998782,-0.0348782,0.0348995']


def assert_true(should_be_true, msg=''):
    if not should_be_true:
        print("Failing because " + msg)
        sys.exit(1)


def assert_equals(expected, actual):
    print(expected)
    print(actual)
    assert_true(expected == actual, '"' + str(actual) + '" should be "' + str(expected) + '"')


def check_results(result_dir):
    print("resultDir:" + result_dir)
    files = os.listdir(result_dir)
    print(str(files))
    assert_true(RESULT_BACKUP_CONFIG_FILE in files, RESULT_BACKUP_CONFIG_FILE + " not found")
    assert_true(RESULT_FILE in files, RESULT_FILE + " not found")

    # make sure the file contents are as expected
    fid = open(result_dir + "/" + RESULT_FILE)
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
