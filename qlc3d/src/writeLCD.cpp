#include <geometry.h>
#include <stdio.h>
#include <qlc3d.h>
#include <string>
#include <iostream>
#include <simu.h>
#include <filesysfun.h>
#include <globals.h>
#include <resultio.h>
#include <util/logging.h>
#include <util/exception.h>

namespace ResultIO {
    void writeMesh(double *p,
                   Mesh *t,
                   Mesh *e,
                   idx np,
                   const std::string &fileName) {
        idx i;
        FILE *fid = fopen(fileName.c_str(), "w");
        if (fid != NULL) {
            fputs("MESH    dimension 3 ElemType Tetrahedra  Nnode 4\nCoordinates\n", fid);
            for (i = 0; i < np; i++) {
                fprintf(fid, "%i\t%f\t%f\t%f\n", i + 1, p[3 * i], p[3 * i + 1], p[3 * i + 2]);
            }
            fprintf(fid, "end coordinates\n\nElements\n");
            for (i = 0 ; i < t->getnElements() ; i++)
                fprintf(fid, "%i\t%i\t%i\t%i\t%i\t%i\n", i + 1, t->getNode(i, 0) + 1, t->getNode(i, 1) + 1, t->getNode(i, 2) + 1, t->getNode(i, 3) + 1, t->getMaterialNumber(i));
            fprintf(fid, "end elements\nMESH    dimension 3 ElemType Triangle  Nnode 3\nCoordinates\nend coordinates\n\nElements\n");
            for (i = 0 ; i < e->getnElements() ; i++)
                fprintf(fid, "%i\t%i\t%i\t%i\t%i\n", i + 1, e->getNode(i, 0) + 1, e->getNode(i, 1) + 1, e->getNode(i, 2) + 1, e->getMaterialNumber(i));
            fprintf(fid, "end elements\n");
            fclose(fid);
        } else {
            RUNTIME_ERROR("Could not open file for output mesh: " + fileName + ".");
        }
    }

/*!
  Writes LCview result in text format.
*/
    void writeLCD_T(double *p,
                    Mesh *t,
                    Mesh *e,
                    SolutionVector *v,
                    SolutionVector *q,
                    int currentIteration,
                    double currentTime,
                    const std::string &meshFileName) {
        int np = v->getnDoF();
        FILE *fid;
        int i;
        char str[15];
        // check whether final result in simulation, if yes, use special filename
        if (currentIteration != SIMU_END_SIMULATION) { // TODO: should not be responsibility of this function
            sprintf(str, "result_t%05i", currentIteration);
        } else {
            sprintf(str, "result_t_final");
        }
        string resname = str;
        resname.append(".dat");
        fid = fopen(resname.c_str(), "w");
        if (fid != nullptr) {
            int npLC = q->getnDoF();
            double *n = tensortovector(q->Values, npLC); // get vector data
            sprintf(str, "%f", currentTime);
            std::string text = "** Result Time :    ";
            text.append(str);
            text.append("\n** z Compression Ratio :  1.00000\n");
            text.append(meshFileName + "\n");
            fprintf(fid, "%s", text.c_str()); //** Result Time :    0.00000000\n** z Compression Ratio :  1.00000\nmeshout.txt\n");
            for (i = 0; i < np; i++) {
                if (i < npLC)
                    fprintf(fid, LCVIEW_TEXT_FORMAT_STRING, i + 1, n[i], n[i + npLC], n[i + 2 * npLC], v->Values[i], n[i + 3 * npLC], n[i + 4 * npLC]);
                else
                    fprintf(fid, LCVIEW_TEXT_FORMAT_STRING, i + 1, 0., 0., 0., v->Values[i], 0., 0.);
            }
            fclose(fid);
            delete [] n;
        } else {
            RUNTIME_ERROR("Could not open result file: " + resname + ".");
        }
    }
// end writeLCD_T

    void writeLCD_B(double *p,
                    Mesh *t, Mesh *e,
                    SolutionVector *v, SolutionVector *q,
                    int currentIteration,
                    double currentTime,
                    double S0,
                    const std::string &meshFileName) {
        int npLC = q->getnDoF();
        int np = v->getnDoF();
        FILE *fid;
        int i;
        char str[15];
        // check whether final result in simulation, if yes, use special filename
        // TODO: pass filename as input arg
        if (currentIteration != SIMU_END_SIMULATION) {
            sprintf(str, "result%05i", currentIteration);
        } else {
            sprintf(str,"result_final");
        }
        string resname = str;
        resname.append(".dat");
        fid = fopen(resname.c_str(), "wb");
        char time[20];
        sprintf(time, "%1.9f\n", currentTime);
        string text = "** Result Time :   ";
        text.append(time);
        fprintf(fid, "%s\n", text.c_str());
        fprintf(fid, "** z Compression Ratio :  1.00000\n");
        fprintf(fid, "%s\n", meshFileName.c_str());
        fprintf(fid, "RAW FLOAT TRI - S0, np, nsols\n");
        fprintf(fid, "%g %d %d\r\n", S0, np, 6);
        for (i = 0;  i < np; i++) {
            if (i < npLC) {
                // WRITE LC REGIONS
                float value;
                value = (float)q->getValue(i, 0);
                fwrite((void *) &value, sizeof(float), 1, fid);
                value = (float)q->getValue(i, 1);
                fwrite((void *) &value, sizeof(float), 1, fid);
                value = (float)q->getValue(i, 2);
                fwrite((void *) &value, sizeof(float), 1, fid);
                value = (float)q->getValue(i, 4);
                fwrite((void *) &value, sizeof(float), 1, fid);
                value = (float)q->getValue(i, 3);
                fwrite((void *) &value, sizeof(float), 1, fid);
                value = v->getValue(i);
                fwrite((void *) &value , sizeof(float), 1 , fid);
            } else {
                // WRITE DIELECTRIC REGIONS
                float value = 0;
                for (int qc = 0; qc < 5; qc++) // q-component
                    fwrite((void *) &value, sizeof(float), 1, fid);
                value = v->getValue(i);
                fwrite((void *) &value, sizeof(float), 1, fid);
            }
        }
        fclose(fid);
    }


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
    void CreateSaveDir(Simu &simu) {
        // CREATES DIRECTORY FOR RESULTS, IF IT DOES NOT
        // ALREADY EXIST
        // first check if savedit already exists
        if (FilesysFun::dirExists(simu.getSaveDir()) || FilesysFun::createDirectory(simu.getSaveDir())) {
            return;
        } else {
            RUNTIME_ERROR("Could not create directory for results, saveDir=" + simu.getSaveDir() + ".");
        }
    }
} // end namespace // WriteResults
