import sys, os
import tempfile as tf
import shutil


def runTest(executable):
    print("current directory=" + os.getcwd());
    # setup
    # 1. create temp directory for project
    projectDir = tf.mkdtemp()
    print("projectDir=" + projectDir)
    # 2. copy mesh and settings file to project dir
    print(os.listdir("./resources"))
    dst = projectDir + "/settings.txt"
    print(dst)
    shutil.copy("./resources/settings.txt", dst)
    shutil.copy("./resources/thin.msh", projectDir + "/thin.msh")

    print(os.listdir(projectDir))

    # 3. run qlc3d executable
    command = executable + " " + dst + " " + projectDir;
    print("command=" + command)
    sys.stdout.flush()
    os.system(command)
    print("end")
    sys.stdout.flush()
    #sys.exit(0)


if __name__ == "__main__":

    print("args:" + str(sys.argv))
    runTest(sys.argv[1])
