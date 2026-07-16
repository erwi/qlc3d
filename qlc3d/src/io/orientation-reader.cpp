#include <io/orientation-reader.h>
#include <io/director-csv-reader.h>
#include <cstdio>
#include <string>
#include <filesystem>
#include <algorithm>
#include <util/logging.h>
#include <util/exception.h>
#include <lc-representation.h>
#include <geom/vec3.h>
#include <globals.h>

namespace qlc3d {
  namespace {
    const char LCVIEW_TEXT_FORMAT_STRING[] = "%i %f %f %f %f %f %f\n";
  }

  std::vector<OrientationSample> LcViewTextReader::read(const std::string &fileName, double s0) const {
    FILE *fid = fopen(fileName.c_str(), "rt");
    if (!fid) {
      RUNTIME_ERROR("Could not open file " + fileName + ".");
    }
    const int lineLen = 256;
    char *line = new char[lineLen];
    // READ 3 LINES OF HEADER DATA
    line = fgets(line, lineLen, fid);
    line = fgets(line, lineLen, fid);
    line = fgets(line, lineLen, fid);
    delete [] line;

    std::vector<OrientationSample> samples;
    int id;
    float n[3], S[2], v;
    // READ FROM FILE UNTIL EOF OR END OF LC REGION (|n| < 1)
    while (fscanf(fid, LCVIEW_TEXT_FORMAT_STRING,
                  &id, &n[0], &n[1], &n[2], &v, &S[0], &S[1]) != EOF) {

      // All zero director vector indicated end of LC region, and we can stop reading LC data
      Vec3 vector(n[0], n[1], n[2]);
      if (vector.norm() == 0.) {
        break;
      }

      auto dir = qlc3d::Director(vector.normalized(), S[0]);
      samples.push_back(OrientationSample{TTensor::fromDirector(dir), std::nullopt});
    }
    fclose(fid);
    return samples;
  }

  std::vector<OrientationSample> LcViewBinaryReader::read(const std::string &fileName, double s0) const {
    FILE *fid = fopen(fileName.c_str(), "rb");
    if (!fid) {
      RUNTIME_ERROR("Could not open file " + fileName);
    }
    // keep a return value ptr to suppress warnings
    char *str;
    const int tempLineLength = 100;
    str = new char[tempLineLength];
    float S0;
    idx np, nsol;
    // READS 5 LINES DISCARDING DATA
    str = fgets(str, tempLineLength, fid);
    str = fgets(str, tempLineLength, fid);
    str = fgets(str, tempLineLength, fid);
    str = fgets(str, tempLineLength, fid);
    str = fgets(str, tempLineLength, fid);
    delete[] str;
    size_t numRead = (size_t) fscanf(fid, "%f %i %i\n", &S0, &np, &nsol);
    numRead++;

    std::vector<OrientationSample> samples;
    samples.reserve(np);
    float q1, q2, q3, q4, q5, temp;
    for (idx i = 0; i < np; i++) {
      numRead = fread(&q1, sizeof(float), 1, fid);
      numRead = fread(&q2, sizeof(float), 1, fid);
      numRead = fread(&q3, sizeof(float), 1, fid);
      numRead = fread(&q5, sizeof(float), 1, fid);
      numRead = fread(&q4, sizeof(float), 1, fid);
      samples.push_back(OrientationSample{TTensor{q1, q2, q3, q4, q5}, std::nullopt});
      for (idx j = 0; j < nsol - 5; j++) // READ&DISCARD POTENTIAL AND FLOW
        numRead = fread((void *) &temp, sizeof(float), 1, fid);
    }
    fclose(fid);
    return samples;
  }

  std::unique_ptr<OrientationReader> createLcViewReader(const std::string &fileName) {
    // check if file exists
    namespace fs = std::filesystem;
    auto filePath = fs::path(fileName);
    if (!fs::exists(filePath)) {
      throw std::invalid_argument("Result file " + fileName + " does not exist.");
    }

    Log::info("Reading result file {}.", fileName);
    FILE *fid = fopen(fileName.c_str(), "rt");
    if (!fid) {
      RUNTIME_ERROR("Could not open result file " + fileName);
    }
    // READ SOME LINES FROM THE FILE AND TRY TO FIND OUT WHICH TYPE IT IS
    bool isBinary = false;
    const int lineLength = 256;
    char *line = new char[lineLength];
    // IF FILE CONTAINS BELOW MAGIC TEXT, IT IS IN BINARY MODE
    const char binaryMarker[] = "RAW FLOAT TRI";
    for (int i = 0; i < 5; i++) {
      line = fgets(line, lineLength, fid); // returns null pointer if fails to read
      std::string sline = line;
      if (sline.find(binaryMarker) < std::string::npos) {
        isBinary = true;
        break;
      }
    }
    delete [] line;
    fclose(fid);

    if (isBinary) {
      Log::info("Result file format is binary.");
      return std::make_unique<LcViewBinaryReader>();
    } else {
      Log::info("Result file format is text.");
      return std::make_unique<LcViewTextReader>();
    }
  }

  std::unique_ptr<OrientationReader> createOrientationReader(const std::string &fileName) {
    namespace fs = std::filesystem;
    std::string extension = fs::path(fileName).extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(),
                    [](unsigned char c) { return std::tolower(c); });

    if (extension == ".csv") {
      return std::make_unique<DirectorCsvReader>();
    }
    return createLcViewReader(fileName);
  }
}
