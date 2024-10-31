#include <cstdio>
#include <string>
#include <iostream>
#include <simu.h>
#include <resultio.h>
#include <util/logging.h>
#include <util/exception.h>
#include <lc-representation.h>
#include <filesystem>

namespace ResultIO {
    void ReadResult(const std::string &fileName, SolutionVector &q) {
        /*!
        * Tries to figure out whether a LCView result file on disk
        * is in text or binary format and then load the data using appropriate
        * loading function
        */

        // check if file exists
        namespace fs = std::filesystem;
        auto filePath = fs::path(fileName);
        if (!fs::exists(filePath)) {
            throw std::invalid_argument("Result file " + fileName + " does not exist.");
        }

        Log::info("Reading result file {}.", fileName);
        FILE *fid = fopen(fileName.c_str() , "rt");
        if (!fid) {
            RUNTIME_ERROR("Could not open result file " + fileName);
        }
        /// RAD SOME LINES FROM THE FILE AND TRY TO FIND OUT WHICH TYPE IT IS
        bool isBinary = false;
        const int lineLength = 256;
        char *line = new char[lineLength];
        // IF FILE CONTAINS BELOW MAGIC TEXT, IT IS IN BINARY MODE
        const char binaryMarker[] = "RAW FLOAT TRI";
        for (int i = 0 ; i < 5 ; i++) {
            line = fgets(line, lineLength, fid); // returns null pointer if fails to read
            string sline = line;
            if (sline.find(binaryMarker) < std::string::npos) {
                isBinary = true;
                break;
            }
        }
        delete [] line;
        fclose(fid);
        if (isBinary) {
            Log::info("Result file format is binary.");
          ResultIO::readBinaryLcViewResultFile(fileName, q);
        } else  {
            Log::info("Result file format is text.");
          ResultIO::readTextLcViewResultFile(fileName, q);
        }
    }

    void readTextLcViewResultFile(const std::string &fileName, SolutionVector &q) {
        FILE *fid = fopen(fileName.c_str(), "rt");
        if (!fid) {
            RUNTIME_ERROR("Could not open file " + fileName + ".");
        }
        const int lineLen = 256;
        char *line = new char[lineLen];
        // READ 3  LINES OF HEADER DATA
        line = fgets(line, lineLen, fid);
        line = fgets(line, lineLen, fid);
        line = fgets(line, lineLen, fid);
        delete [] line;
        int id;
        float n[3], S[2], v;
        const unsigned int npLC = q.getnDoF();
        // READ FROM FILE UNTIL EOF OR END OF LC REGION (|n| < 1)
        unsigned int i = 0;
        while (fscanf(fid, LCVIEW_TEXT_FORMAT_STRING,
                      &id, &n[0], &n[1], &n[2], &v, &S[0], &S[1]) != EOF) {

            // All zero director vector indicated end of LC region, and we can stop reading LC data
            Vec3 vector(n[0], n[1], n[2]);
            if (vector.norm() == 0.) {
              break;
            }

            auto dir = qlc3d::Director(vector.normalized(), S[0]);
            q.setValue(i, dir);
            i++;
        }
        fclose(fid);
        if (i != npLC) {
            RUNTIME_ERROR(fmt::format("The loaded result file size {} does not match the expected size {}", i, npLC));
        }
    }

    void readBinaryLcViewResultFile(const std::string &fileName, SolutionVector &q) {
        FILE *fid = fopen(fileName.c_str(), "rb");
        // keep a return value ptr to suppress warnings
        char *str;
        const int tempLineLength = 100;
        str = new char[tempLineLength];
        if (fid) {
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
            if (np < q.getnDoF()) {
                RUNTIME_ERROR("The loaded result file does not match the mesh size.");
            }
            float q1, q2, q3, q4, q5, temp;
            for (idx i = 0; i < q.getnDoF(); i++) {
                numRead = fread(&q1, sizeof(float), 1, fid);
                numRead = fread(&q2, sizeof(float), 1, fid);
                numRead = fread(&q3, sizeof(float), 1, fid);
                numRead = fread(&q5, sizeof(float), 1, fid);
                numRead = fread(&q4, sizeof(float), 1, fid);
                q.setValue(i, 0, q1);
                q.setValue(i, 1, q2);
                q.setValue(i, 2, q3);
                q.setValue(i, 3, q4);
                q.setValue(i, 4, q5);
                for (idx j = 0; j < nsol - 5; j++) // READ&DISCARD POTENTIAL AND FLOW
                    numRead = fread((void *) &temp, sizeof(float), 1, fid);
            }
            fclose(fid);
        } else {
            RUNTIME_ERROR("Could not open file " + fileName);
        }
    }
} // end namespace // WriteResults
