import sys
import os


def assert_true(should_be_true: bool, msg: str = '') -> None:
    if not should_be_true:
        print("Failing because " + msg)
        sys.exit(1)


def assert_equals(expected, actual) -> None:
    assert_true(expected == actual, '"' + str(actual) + '" should be "' + str(expected) + '"')

def run_executable(executable: str, settings_file: str, project_dir:str) -> None:
    """Launches the provided executable and waits for it to return successfully"""
    command = executable + " " + settings_file + " " + project_dir
    print("command=" + command)
    sys.stdout.flush()
    exit_code = os.system(command)
    assert_equals(0, exit_code)
    print("end")
    sys.stdout.flush()
