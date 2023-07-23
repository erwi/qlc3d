#include <stdio.h>
#include <string>
#include <iostream>
#include <simu.h>
#include <resultio.h>
#include <util/logging.h>
#include <util/exception.h>

namespace ResultIO {
    void ReadResult(Simu &simu, SolutionVector &q) {
        /*!
        * Tries to figure out whether a LCView result file on disk
        * is in text or binary format and then load the data using appropriate
        * loading function
        */
        string filename = simu.getLoadQ();
        Log::info("Reading result file {}.", filename);
        FILE *fid = fopen(filename.c_str() , "rt");
        if (!fid) {
            RUNTIME_ERROR("Could not open file " + filename + ".");
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
            ResultIO::ReadLCD_B(&simu, &q);
        } else  {
            Log::info("Result file format is text.");
            ResultIO::ReadLCD_T(simu, q);
        }
    }

    void ReadLCD_T(Simu &simu, SolutionVector &q) {
        /*!
         * loads Q-tensor from a result file, assuming the file is written as a text file
         */
        string filename = simu.getLoadQ();
        FILE *fid = fopen(filename.c_str(), "rt");
        if (!fid) {
            RUNTIME_ERROR("Could not open file " + filename + ".");
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
        //fgets(line, lineLen, fid);
        std::vector<float> q1;
        std::vector<float> q2;
        std::vector<float> q3;
        std::vector<float> q4;
        std::vector<float> q5;
        double rt6 = sqrt(6.0);
        double rt2 = sqrt(2.0);
        // READ FROM FILE UNTIL EOF OR END OF LC REGION (|n| < 1)
        while (fscanf(fid, LCVIEW_TEXT_FORMAT_STRING,
                      &id, &n[0], &n[1], &n[2], &v, &S[0], &S[1]) != EOF) {
            // MAKE SURE DIRECTOR LENGTH ~= 1.
            // ALL ZERO DIRECTOR MEANS DIELECTRIC REGION, WHICH WE DONT WANT TO READ
            float dirLen = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
            if (dirLen < 0.95)
                break;
            // CONVERT DIRECTOR TO Q-TENSOR
            double a1, a2, a3, a4, a5; // "NORMAL" Q-TENSOR COMPONENTS
            a1 = S[0] * (3 * n[0] * n[0] - 1) / 2.0; // Qxx
            a2 = S[0] * (3 * n[1] * n[1] - 1) / 2.0; // Qyy
            a3 = S[0] * (3 * n[0] * n[1]) / 2.0; // Qxy
            a4 = S[0] * (3 * n[1] * n[2]) / 2.0; // Qyz
            a5 = S[0] * (3 * n[0] * n[2]) / 2.0; // Qxz
            // CONVERT TO TRACELESS BASIS
            q1.push_back(0.5 * (a1 + a2) *rt6);
            q2.push_back((a1 + (a1 + a2) / 2.0)*rt2);
            q3.push_back(a3 * rt2);
            q4.push_back(a4 * rt2);
            q5.push_back(a5 * rt2);
        }
        fclose(fid);
        if ((idx) q1.size() != q.getnDoF()) {
            RUNTIME_ERROR("The loaded result file does not match the msh size.");
            exit(1);
        }
        for (idx i = 0 ; i < q.getnDoF() ; ++i) {
            q.setValue(i, 0, q1[i]);
            q.setValue(i, 1, q2[i]);
            q.setValue(i, 2, q3[i]);
            q.setValue(i, 3, q4[i]);
            q.setValue(i, 4, q5[i]);
        }
    }

    void ReadLCD_B(Simu *simu, SolutionVector *q) {
        // READS BINARY FORMATED RESULT FILE
        string filename = simu->getLoadQ();
        FILE *fid = fopen(filename.c_str(), "rb");
        // keep a return value ptr to suppress warnings
        char *str;
        const int tempLineLength = 100;
        str = new char[tempLineLength];
        if (fid) {
            float S0;
            idx np, nsol;
            // READS 5 LINES DISCARIDING DATA
            str = fgets(str, tempLineLength, fid);
            str = fgets(str, tempLineLength, fid);
            str = fgets(str, tempLineLength, fid);
            str = fgets(str, tempLineLength, fid);
            str = fgets(str, tempLineLength, fid);
            delete[] str;
            size_t numRead = (size_t) fscanf(fid, "%f %i %i\n", &S0, &np, &nsol);
            numRead++;
            if (np < q->getnDoF()) {
                RUNTIME_ERROR("The loaded result file does not match the mesh size.");
            }
            float q1, q2, q3, q4, q5, temp;
            for (idx i = 0; i < q->getnDoF(); i++) {
                numRead = fread(&q1, sizeof(float), 1, fid);
                numRead = fread(&q2, sizeof(float), 1, fid);
                numRead = fread(&q3, sizeof(float), 1, fid);
                numRead = fread(&q5, sizeof(float), 1, fid);
                numRead = fread(&q4, sizeof(float), 1, fid);
                q->setValue(i, 0, q1);
                q->setValue(i, 1, q2);
                q->setValue(i, 2, q3);
                q->setValue(i, 3, q4);
                q->setValue(i, 4, q5);
                for (idx j = 0; j < nsol - 5; j++) // READ&DISCARD POTENTIAL AND FLOW
                    numRead = fread((void *) &temp, sizeof(float), 1, fid);
            }
            fclose(fid);
        } else {
            RUNTIME_ERROR("Could not open file " + filename + ".");
        }
    }
} // end namespace // WriteResults
