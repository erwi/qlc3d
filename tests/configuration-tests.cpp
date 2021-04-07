//
// Created by eero on 07/04/2021.
//
#include <catch.h>

#include <settings-reader.h>
#include <reader.h>
#include <test-util.h>

TEST_CASE("Catch library should work") {
    REQUIRE(true);
}

TEST_CASE("MeshName is required in settings file") {
    try {
        auto settingsFile = TestUtil::TemporaryFile::empty();
        SettingsReader reader(settingsFile.name());
    } catch (ReaderError &e) { // expect
        REQUIRE(e.errorMessage.find("Key not found :MeshName") == 0);
        return;
    }
    REQUIRE_FALSE(true); // should not reach this
}

TEST_CASE("Read Simu from settings file") {
// ARRANGE
    std::string contents;
    contents += "MeshName= wowowoo.txt\n";
    contents += "dt=123\n";
    contents += "qmatrixsolver=PCG\n";
    contents += "targetdq=1.2\n";
    contents += "maxError=2.2\n";
    contents += "endcriterion=CHANGE\n";
    contents += "loadQ=some-file.abc\n";
    contents += "savedir=someDir\n";
    contents += "endValue=1.5\n";
    contents += "outputEnergy=1\n";
    contents += "outputformat=123\n";
    contents += "saveiter=13\n";
    contents += "saveFormat=[regularvecmat, LCviewTXT]\n";
    contents += "numAssemblyThreads=99\n";
    contents += "numMatrixSolverThreads=98\n";
    contents += "dtlimits=[3.2, 4.2]\n";
    contents += "dtfunction=[1.2, 2.2, 3.2, 4.2]\n";
    contents += "stretchvector=[1.2, 2.2, 3.2]\n";
    contents += "regularGridSize=[2,3,4]\n";

    auto settingsFile = TestUtil::TemporaryFile::withContents(contents);

// ACT
    SettingsReader reader(settingsFile.name());
    auto simu = reader.simu();

// ASSERT - check that the read simu matches the file contents specified above
    REQUIRE(simu->meshName() == "wowowoo.txt");
    REQUIRE(simu->initialTimeStep() == 123);
    REQUIRE(simu->getQMatrixSolver() == Simu::QMatrixSolvers::PCG);
    REQUIRE(simu->getTargetdQ() == 1.2);
    REQUIRE(simu->getMaxError() == 2.2);
    REQUIRE(simu->getMindt() == 3.2);
    REQUIRE(simu->getMaxdt() == 4.2);
    REQUIRE(simu->getEndCriterion() == Simu::EndCriteria::Change);
    REQUIRE(simu->getLoadQ() == "some-file.abc");
    REQUIRE(simu->getSaveDir() == "somedir");
    REQUIRE(simu->getEndValue() == 1.5);
    REQUIRE(simu->getOutputEnergy() == 1);
    REQUIRE(simu->getOutputFormat() == 123);
    REQUIRE(simu->getSaveIter() == 13);
    REQUIRE(simu->getAssemblyThreadCount() == 99);
    REQUIRE(simu->getMatrixSolverThreadCount() == 98);

// check lists
// dtLimits
    REQUIRE(simu->getMindt() == 3.2);
    REQUIRE(simu->getMaxdt() == 4.2);

// dtFunction
    double dtf[4]{0, 0, 0, 0};
    simu->getdtFunction(dtf);
    REQUIRE((dtf[0] == 1.2 && dtf[1] == 2.2 && dtf[2] == 3.2 && dtf[3] == 4.2));

// strech vector
    REQUIRE(simu->getStretchVectorX() == 1.2);
    REQUIRE(simu->getStretchVectorY() == 2.2);
    REQUIRE(simu->getStretchVectorZ() == 3.2);

// regula grid size
    REQUIRE(simu->getRegularGridXCount() == 2);
    REQUIRE(simu->getRegularGridYCount() == 3);
    REQUIRE(simu->getRegularGridZCount() == 4);

// specified save formats
    REQUIRE(simu->getSaveFormat().count(Simu::RegularVecMat));
    REQUIRE(simu->getSaveFormat().count(Simu::LCviewTXT));
}