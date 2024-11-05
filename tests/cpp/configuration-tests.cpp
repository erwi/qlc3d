#include <catch.h>

#include <settings-reader.h>
#include <reader.h>
#include <test-util.h>
#include "geom/vec3.h"
#include "lc-representation.h"

TEST_CASE("Catch library should work") {
    REQUIRE(true);
}

TEST_CASE("MeshName is required in settings file") {
    try {
        auto settingsFile = TestUtil::TemporaryFile::empty();
        SettingsReader reader(settingsFile.name());
    } catch (ReaderError &e) { // expect
        REQUIRE(e.errorMessage.find("Key not found: MeshName") == 0);
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
    Vec3 sv = simu->getStretchVector();
    REQUIRE(sv.x() == 1.2);
    REQUIRE(sv.y() == 2.2);
    REQUIRE(sv.z() == 3.2);

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
    contents += "K24 = 4.44e-12\n";
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
    REQUIRE(lc->K24() == Approx(4.44e-12).margin(eps));
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
    REQUIRE(lc->L3() == Approx(7.8230134823571808e-12).margin(eps)); // 4.0 * K24 / (9.0 * S0 * S0);
    REQUIRE(lc->L4() == Approx(4.9374855277287484e-09).margin(eps));
    REQUIRE(lc->L5() == 0);
    REQUIRE(lc->L6() == Approx(2.4557036852882312e-12).margin(eps));
    REQUIRE(lc->u1() == Approx(0.08897796866194542).margin(eps)); // 2  * gamma1 / (9*So*S0)
}

TEST_CASE("read REFINEMENT from settings file") {
    std::string contents;
    contents += "MeshName = test.msh\n"; // required in every settings file
    contents += "REFINEMENT1.TYPE = Sphere\n";
    contents += "REFINEMENT1.X=[0.0]\n";
    contents += "REFINEMENT1.Y=[0.2]\n";
    contents += "REFINEMENT1.Z=[0.4]\n";
    contents += "REFINEMENT1.Iterations = [3]\n";
    contents += "REFINEMENT1.Times = [1, 2]\n";
    contents += "REFINEMENT1.Values = [0.12]\n";

    contents += "RepRefIter=10 \n";
    contents += "RepRefTime = 2e-9\n";

    auto settingsFile = TestUtil::TemporaryFile::withContents(contents);

    SettingsReader reader(settingsFile.name());
    auto refinement = reader.refinement();
    const vector<RefinementConfig> &refs = refinement->getRefinementConfig();


    REQUIRE(refinement->getRepRefIter() == 10);
    REQUIRE(refinement->getRepRefTime() == 2e-9);

    REQUIRE(refs.size() == 1);

    auto ref = refs[0];

    REQUIRE(ref.type_ == "sphere");

    REQUIRE((ref.x_.size() == 1 && ref.x_[0] == 0));
    REQUIRE((ref.y_.size() == 1 && ref.y_[0] == 0.2));
    REQUIRE((ref.z_.size() == 1 && ref.z_[0] == 0.4));

    REQUIRE((ref.iterations_.size() == 1 && ref.iterations_[0] == 3));
    REQUIRE((ref.times_.size() == 2 && ref.times_[0] == 1 && ref.times_[1] == 2));
    REQUIRE((ref.iterations_.size() == 1 && ref.iterations_[0] == 3));
    REQUIRE((ref.values_.size() == 1 && ref.values_[0] == 0.12));
}

TEST_CASE("read electrodes from settings file") {
  std::string contents;
  contents += "MeshName= test.msh\n"; // required in every settings file

  contents += "E1.Time = [0]\n";
  contents += "E1.Pot = [0]\n";

  contents += "E3.Time = [1, 2]\n";
  contents += "E3.Pot = [3, 4]\n";

  auto settingsFile = TestUtil::TemporaryFile::withContents(contents);

  SettingsReader reader(settingsFile.name());

  auto electrodes = reader.electrodes();

  REQUIRE(2 == electrodes->getnElectrodes());

  // potentials at time 0
  auto potByElectrode = electrodes->getCurrentPotentials(0);
  REQUIRE(0 == potByElectrode[1]);
  REQUIRE(0 == potByElectrode[3]);

  // potentials at time 1.1
  potByElectrode = electrodes->getCurrentPotentials(1.1);
  REQUIRE(0 == potByElectrode[1]);
  REQUIRE(3 == potByElectrode[3]);

  // potentials at time 1000
  potByElectrode = electrodes->getCurrentPotentials(1000);
  REQUIRE(0 == potByElectrode[1]);
  REQUIRE(4 == potByElectrode[3]);
}

TEST_CASE("Read solver settings from settings file") {
  std::string contents;
  contents += "MeshName= test.msh\n"; // required in every settings file
  contents += "NumAssemblyThreads = 1\n";
  contents += "Q_Solver = 2\n";
  contents += "V_Solver = 3\n";
  contents += "Q_Newton_Panic_Iter = 4\n";
  contents += "Q_Newton_Panic_Coeff = 0.5\n";
  contents += "Q_PCG_Preconditioner = 6\n";
  contents += "Q_PCG_Maxiter = 7\n";
  contents += "Q_PCG_Toler = 8\n";
  contents += "Q_GMRES_Preconditioner = 9\n";
  contents += "Q_GMRES_Maxiter = 10\n";
  contents += "Q_GMRES_Restart = 11\n";
  contents += "Q_GMRES_Toler = 12\n";
  contents += "V_PCG_Preconditioner = 13\n";
  contents += "V_PCG_Maxiter = 14\n";
  contents += "V_PCG_Toler = 15\n";
  contents += "V_GMRES_Preconditioner = 16\n";
  contents += "V_GMRES_Maxiter = 17\n";
  contents += "V_GMRES_Restart = 18\n";
  contents += "V_GMRES_Toler = 19\n";

  auto settingsFile = TestUtil::TemporaryFile::withContents(contents);
  SettingsReader reader(settingsFile.name());
  auto solverSettings = reader.solverSettings();

  REQUIRE(1 == solverSettings->getnThreads());
  REQUIRE(2 == solverSettings->getQ_Solver());
  REQUIRE(3 == solverSettings->getV_Solver());
  REQUIRE(4 == solverSettings->getQ_Newton_Panic_Iter());
  REQUIRE(0.5 == solverSettings->getQ_Newton_Panic_Coeff());
  REQUIRE(6 == solverSettings->getQ_PCG_Preconditioner());
  REQUIRE(7 == solverSettings->getQ_PCG_Maxiter());
  REQUIRE(8 == solverSettings->getQ_PCG_Toler());
  REQUIRE(9 == solverSettings->getQ_GMRES_Preconditioner());
  REQUIRE(10 == solverSettings->getQ_GMRES_Maxiter());
  REQUIRE(11 == solverSettings->getQ_GMRES_Restart());
  REQUIRE(12 == solverSettings->getQ_GMRES_Toler());
  REQUIRE(13 == solverSettings->getV_PCG_Preconditioner());
  REQUIRE(14 == solverSettings->getV_PCG_Maxiter());
  REQUIRE(15 == solverSettings->getV_PCG_Toler());
  REQUIRE(16 == solverSettings->getV_GMRES_Preconditioner());
  REQUIRE(17 == solverSettings->getV_GMRES_Maxiter());
  REQUIRE(18 == solverSettings->getV_GMRES_Restart());
  REQUIRE(19 == solverSettings->getV_GMRES_Toler());
}

TEST_CASE("Read alignment from settings file") {
  std::string contents;
  contents += "MeshName= test.msh\n"; // required in every settings file

  // add case for strong anchoring
  contents += "FIXLC1.Anchoring = Strong\n";
  contents += "FIXLC1.Strength = 1e-4\n";
  contents += "FIXLC1.Easy = [80.0000, 45.0000, 0.0000]\n";
  contents += "FIXLC1.K1 = 1.0000\n";
  contents += "FIXLC1.K2 = 1.0000\n";

  // add case for homeotropic anchoring
  contents += "FIXLC2.Anchoring = HomeOtropic\n";
  contents += "FIXLC2.Strength = 1e-4\n";
  contents += "FIXLC2.Easy = [80.0000, 45.0000, 0.0000]\n";
  contents += "FIXLC2.K1 = 1.0000\n";
  contents += "FIXLC2.K2 = 1.0000\n";

  // add case for weak anchoring
  contents += "FIXLC3.Anchoring = weak\n";
  contents += "FIXLC3.Strength = 1e-4\n";
  contents += "FIXLC3.Easy = [80.0000, 45.0000, 0.0000]\n";
  contents += "FIXLC3.K1 = 1.0000\n";
  contents += "FIXLC3.K2 = 2.0000\n";

  // add case for degenerate anchoring
  contents += "FIXLC4.Anchoring = Degenerate\n";
  contents += "FIXLC4.strength = 1e-4\n";


  auto settingsFile = TestUtil::TemporaryFile::withContents(contents);

  SettingsReader reader(settingsFile.name());

  auto alignment = reader.alignment();

  REQUIRE(alignment->getnSurfaces() == 4);

  // Check that the correct types were added. Other properties are tested in the Surface class tests
  REQUIRE(alignment->getSurface(0).getAnchoringType() == AnchoringType::Strong);
  REQUIRE(alignment->getSurface(1).getAnchoringType() == AnchoringType::Homeotropic);
  REQUIRE(alignment->getSurface(2).getAnchoringType() == AnchoringType::Weak);
  REQUIRE(alignment->getSurface(3).getAnchoringType() == AnchoringType::Degenerate);
}

TEST_CASE("Read alignment with analytic expressions for tilt and twist") {
  std::string contents;
  contents += "MeshName= test.msh\n"; // required in every settings file

  contents += "FIXLC1.Anchoring = Strong\n";
  contents += "FIXLC1.Strength = 1e-4\n";
  contents += "FIXLC1.Easy = [\"X * 90\", \"X * 45 \", \" Z + Y + Z\"]\n";

  auto settingsFile = TestUtil::TemporaryFile::withContents(contents);
  SettingsReader reader(settingsFile.name());

  auto alignment = reader.alignment();

  auto &surface = alignment->getSurface(0);

  REQUIRE(surface.getAnchoringType() == AnchoringType::Strong);

  REQUIRE(surface.getEasyTiltAngleAt(Vec3(0, 0, 0)) == Approx(0.0).margin(1e-15));
  REQUIRE(surface.getEasyTwistAngleAt(Vec3(0, 0, 0)) == Approx(0.0).margin(1e-15));

  REQUIRE(surface.getEasyTiltAngleAt(Vec3(1, 0, 0)) == Approx(90.0).margin(1e-15));
  REQUIRE(surface.getEasyTwistAngleAt(Vec3(1, 0, 0)) == Approx(45.0).margin(1e-15));
}

TEST_CASE("Read boxes from settings file") {
  // ARRANGE
  std::string contents;
  contents += "MeshName= test.msh\n"; // required in every settings file

  contents += "BOX1.Type = Normal\n";
  contents += "BOX1.X = [0.0, 1.0]\n";
  contents += "BOX1.Y = [-1.0, 2.0]\n";
  contents += "BOX1.Z = [0.0, 0.1]\n";
  contents += "BOX1.Params = [1.0, 2., 3.]\n";
  contents += "BOX1.Tilt = [-0.1, 0.5]\n";
  contents += "BOX1.Twist = [-0.2, 0.5]\n";

  contents += "BOX2.Type = Random\n";
  contents += "BOX2.X = [0.0, 1.0]\n";
  contents += "BOX2.Y = [0.0, 1.0]\n";
  contents += "BOX2.Z = [0.0, 1.0]\n";

  contents += "BOX3.Type = Hedgehog\n";
  contents += "BOX3.X = [0.0, 1.0]\n";
  contents += "BOX3.Y = [0.0, 1.0]\n";
  contents += "BOX3.Z = [0.0, 1.0]\n";

  auto settingsFile = TestUtil::TemporaryFile::withContents(contents);

  // ACT
  SettingsReader reader(settingsFile.name());

  // ASSERT
  auto boxes = reader.initialVolumeOrientation();

  REQUIRE(3 == boxes->getBoxCount());

  auto& box1 = boxes->getBox(0);
  REQUIRE(box1.getType() == BoxType::Normal);
  auto bbox = box1.getBoundingBox();

  double defaultParam = -1.;
  REQUIRE(box1.getParam(0, defaultParam) == 1.0);
  REQUIRE(box1.getParam(1, defaultParam) == 2.0);
  REQUIRE(box1.getParam(2, defaultParam) == 3.0);

  auto dirBottom = box1.getDirectorAt(Vec3(0.5, 0, 0.));
  REQUIRE(dirBottom.tiltDegrees() == -.1);
  REQUIRE(dirBottom.twistDegrees() == -.2);

  auto dirTop = box1.getDirectorAt(Vec3(0.5, 0, 0.1));
  REQUIRE(dirTop.tiltDegrees() == 0.4); // bottom + delta
  REQUIRE(dirTop.twistDegrees() == 0.3); // bottom + delta

  REQUIRE(bbox.getXMin() == 0.0);
  REQUIRE(bbox.getXMax() == 1.0);
  REQUIRE(bbox.getYMin() == -1.0);
  REQUIRE(bbox.getYMax() == 2.0);
  REQUIRE(bbox.getZMin() == 0.0);
  REQUIRE(bbox.getZMax() == 0.1);

  // Check the other boxe types.
  auto& box2 = boxes->getBox(1);
  REQUIRE(box2.getType() == BoxType::Random);

  auto& box3 = boxes->getBox(2);
  REQUIRE(box3.getType() == BoxType::Hedgehog);
}

TEST_CASE("Read normal box with expression tilt and twist") {
  std::string contents;
  contents += "MeshName= test.msh\n"; // required in every settings file

  contents += "BOX1.Type = Normal\n";
  contents += "BOX1.X = [0.0, 1.0]\n";
  contents += "BOX1.Y = [-1.0, 2.0]\n";
  contents += "BOX1.Z = [0.0, 0.1]\n";
  contents += "BOX1.Tilt = \"x + y + z\"\n";
  contents += "BOX1.Twist = \"x - y - z\"\n";

  auto settingsFile = TestUtil::TemporaryFile::withContents(contents);
  SettingsReader reader(settingsFile.name());

  auto boxes = reader.initialVolumeOrientation();
  auto &box = boxes->getBox(0);
}