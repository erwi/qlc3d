#include <io/director-csv-reader.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <optional>
#include <unordered_map>
#include <lc-representation.h>
#include <geom/vec3.h>
#include <util/exception.h>
#include <fmt/format.h>

namespace qlc3d {
  namespace {
    std::string trim(const std::string &s) {
      size_t start = s.find_first_not_of(" \t\r\n");
      if (start == std::string::npos) {
        return "";
      }
      size_t end = s.find_last_not_of(" \t\r\n");
      return s.substr(start, end - start + 1);
    }

    std::string toLower(const std::string &s) {
      std::string out = s;
      std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) { return std::tolower(c); });
      return out;
    }

    std::vector<std::string> splitCsvLine(const std::string &line) {
      std::vector<std::string> fields;
      std::stringstream ss(line);
      std::string field;
      while (std::getline(ss, field, ',')) {
        fields.push_back(trim(field));
      }
      return fields;
    }

    // recognized column names -> whether they are required
    const std::unordered_map<std::string, bool> RECOGNIZED_COLUMNS = {
      {"x", true}, {"y", true}, {"z", true},
      {"nx", true}, {"ny", true}, {"nz", true},
      {"s", false}
    };
  }

  /**
   * read line ignoring anything after a hash character # which is interpreted as start of a comment
   * @param fin
   * @param line
   * @return
   */
  std::ifstream& readLine(std::ifstream &fin, std::string &line) {
    std::getline(fin, line);
    size_t hashPos = line.find_first_of('#');
    if (hashPos != std::string::npos) {
      line = line.substr(0, hashPos);
    }
    return fin;
  }

  std::vector<OrientationSample> DirectorCsvReader::read(const std::string &fileName, double s0) const {
    std::ifstream fin(fileName);
    if (!fin.is_open()) {
      RUNTIME_ERROR("Could not open file " + fileName + ".");
    }

    // Find first non-empty line as the header.
    std::string headerLine;
    while (readLine(fin, headerLine)) {
      if (!trim(headerLine).empty()) {
        break;
      }
      headerLine.clear();
    }
    if (trim(headerLine).empty()) {
      RUNTIME_ERROR("Director CSV file " + fileName + " is empty.");
    }

    std::vector<std::string> headerFields = splitCsvLine(headerLine);
    // maps column name -> index in each data row
    std::unordered_map<std::string, size_t> columnIndex;
    for (size_t i = 0; i < headerFields.size(); i++) {
      std::string name = toLower(headerFields[i]);
      if (RECOGNIZED_COLUMNS.find(name) == RECOGNIZED_COLUMNS.end()) {
        RUNTIME_ERROR(fmt::format("Unrecognized column \"{}\" in director CSV file {}.", headerFields[i], fileName));
      }
      if (columnIndex.count(name)) {
        RUNTIME_ERROR(fmt::format("Duplicate column \"{}\" in director CSV file {}.", headerFields[i], fileName));
      }
      columnIndex[name] = i;
    }

    std::vector<std::string> requiredColumns = {"x", "y", "z", "nx", "ny", "nz"};
    std::vector<std::string> missingColumns;
    for (const auto &col : requiredColumns) {
      if (!columnIndex.count(col)) {
        missingColumns.push_back(col);
      }
    }
    if (!missingColumns.empty()) {
      std::string joined;
      for (size_t i = 0; i < missingColumns.size(); i++) {
        if (i > 0) joined += ", ";
        joined += missingColumns[i];
      }
      RUNTIME_ERROR(fmt::format("Director CSV file {} is missing required column(s): {}.", fileName, joined));
    }
    bool hasS = columnIndex.count("s") > 0;

    std::vector<OrientationSample> samples;
    std::string line;
    size_t lineNumber = 1; // header was line 1
    while (readLine(fin, line)) {
      lineNumber++;
      if (trim(line).empty()) {
        continue;
      }
      std::vector<std::string> fields = splitCsvLine(line);
      if (fields.size() != headerFields.size()) {
        RUNTIME_ERROR(fmt::format("Director CSV file {} line {} has {} column(s), expected {}.",
                                   fileName, lineNumber, fields.size(), headerFields.size()));
      }

      auto parseDouble = [&](const std::string &colName) {
        const std::string &field = fields[columnIndex.at(colName)];
        try {
          return std::stod(field);
        } catch (const std::exception &) {
          throw std::runtime_error(fmt::format("Director CSV file {} line {}: could not parse \"{}\" as a number for column \"{}\".",
                                                fileName, lineNumber, field, colName));
        }
      };

      double x = parseDouble("x");
      double y = parseDouble("y");
      double z = parseDouble("z");
      double nx = parseDouble("nx");
      double ny = parseDouble("ny");
      double nz = parseDouble("nz");

      Vec3 vector(nx, ny, nz);
      if (vector.norm() == 0.) {
        RUNTIME_ERROR(fmt::format("Director CSV file {} line {}: director vector (nx, ny, nz) has zero length.",
                                   fileName, lineNumber));
      }

      double S = hasS ? parseDouble("s") : s0;

      auto dir = qlc3d::Director(vector.normalized(), S);
      samples.push_back(OrientationSample{TTensor::fromDirector(dir), Vec3(x, y, z)});
    }

    if (samples.empty()) {
      RUNTIME_ERROR("Director CSV file " + fileName + " contains no data rows.");
    }

    return samples;
  }
}
