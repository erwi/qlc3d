import sys


def assert_true(should_be_true, msg=''):
    if not should_be_true:
        print("Failing because " + msg)
        sys.exit(1)


def assert_equals(expected, actual):
    assert_true(expected == actual, '"' + str(actual) + '" should be "' + str(expected) + '"')

