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
    contents += "saveTime=6e-6\n";
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
    REQUIRE(simu->getSaveTime() == 6e-6);
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

TEST_CASE("read LC from settings file") {
    // ARRANGE
    std::string contents;
    contents += "MeshName = ""wowowoo.msh\n"; // required in every settings file
    contents += "K11 = 1.23e-12\n";
    contents += "K22 = 2.23e-12\n";
    contents += "K33 = 3.33e-12\n";
    contents += "p0 = 0.01\n";
    contents += "A = -0.12\n";
    contents += "B = -2.1333\n";
    contents += "C = 1.733\n";
    contents += "eps_par = 12\n";
    contents += "eps_per = 13\n";
    contents += "e1 = 1e-11\n";
    contents += "e3 = -1e-11\n";
    contents += "gamma1 = 0.101\n";
    auto settingsFile = TestUtil::TemporaryFile::withContents(contents);

    // ACT:
    SettingsReader reader(settingsFile.name());
    auto lc = reader.lc();

    // ASSERT:
    const double eps = 1e-9;
    // Explicitly defined values
    REQUIRE(lc->K11() == Approx(1.23e-12).margin(eps));
    REQUIRE(lc->K22() == Approx(2.23e-12).margin(eps));
    REQUIRE(lc->K33() == Approx(3.33e-12).margin(eps));
    REQUIRE(lc->p0() == Approx(0.01).margin(eps));
    REQUIRE(lc->A() == Approx(-0.12).margin(eps));
    REQUIRE(lc->B() == Approx(-2.1333).margin(eps));
    REQUIRE(lc->C() == Approx(1.733).margin(eps));
    REQUIRE(lc->eps_par() == Approx(12).margin(eps));
    REQUIRE(lc->eps_per() == Approx(13).margin(eps));
    REQUIRE(lc->e1() == Approx(1e-11).margin(eps));
    REQUIRE(lc->e3() == Approx(-1e-11).margin(eps));
    REQUIRE(lc->gamma1() == Approx(0.101).margin(eps));

    // Implicit values, calculated in LC constructor from the above explicit values
    REQUIRE(lc->S0() == Approx(0.50224218371360507).margin(eps));
    REQUIRE(lc->L1() == Approx(2.5812420611831688e-12).margin(eps));
    REQUIRE(lc->L2() == Approx(-1.7619399735038694e-12).margin(eps));
    REQUIRE(lc->L3() == 0);
    REQUIRE(lc->L4() == Approx(4.9374855277287484e-09).margin(eps));
    REQUIRE(lc->L5() == 0);
    REQUIRE(lc->L6() == Approx(2.4557036852882312e-12).margin(eps));
    REQUIRE(lc->u1() == Approx(0.089376978566352377).margin(eps));
}