#include <stdio.h>
#include <stdlib.h>
#include <qlc3d.h>
#include <string>
#include <vector>
#include <iostream>
#include <meshrefinement.h>
#include <refinfo.h>
#include <filesysfun.h>
#include <globals.h>
#include <reader.h>
#include <settings_file_keys.h>
#include <simu.h>
using std::cerr;
using std::endl;




void readLC(LC& lc,Reader& reader) {
    try {
        lc.K11 = reader.get<double>(SFK_K11, LC::DEFAULT_K11);
        lc.K22 = reader.get<double>(SFK_K22, LC::DEFAULT_K22);
        lc.K33 = reader.get<double>(SFK_K33, LC::DEFAULT_K33);
        lc.p0  = reader.get<double>(SFK_p0, LC::DEFAULT_P0);
        lc.A = reader.get<double>(SFK_A, LC::DEFAULT_A);
        lc.B = reader.get<double>(SFK_B, LC::DEFAULT_B);
        lc.C = reader.get<double>(SFK_C, LC::DEFAULT_C);
        lc.eps_par = reader.get<double>( SFK_EPS_PAR, LC::DEFAULT_EPS_PAR);
        lc.eps_per = reader.get<double>(SFK_EPS_PER, LC::DEFAULT_EPS_PER);
        lc.e11 = reader.get<double>(SFK_E1, LC::DEFAULT_E1);
        lc.e33 = reader.get<double>(SFK_E3, LC::DEFAULT_E3);
        lc.gamma1 = reader.get<double>(SFK_GAMMA1, LC::DEFAULT_GAMMA1);
        lc.convert_params_n2Q();
    } catch (ReaderError e) {
        e.printError();
        std::exit(qlc3d_GLOBALS::ERROR_CODE_BAD_SETTINGS_FILE);
    }
}//end void readLC

void readSimu(Simu &simu, EventList &eventList, Reader &reader) {
    /*! loads values from settings file reader */
    try {
        simu.MeshName = reader.getValueByKey<string>(SFK_MESH_NAME);
        std::string key;
        // Read optional string values
        simu.setLoadQ(reader.get<string>(SFK_LOAD_Q, Simu::DEFAULT_LOAD_Q));
        simu.setSaveDir(reader.get<string>(SFK_SAVE_DIR, Simu::DEFAULT_SAVE_DIR));
        simu.setQMatrixSolver(reader.get<string>(SFK_Q_MATRIX_SOLVER, Simu::DEFAULT_Q_MATRIX_SOLVER));
        simu.setEndCriterion(reader.get<string>(SFK_END_CRITERION, Simu::DEFAULT_END_CRITERION));
        // string arrays
        simu.setSaveFormats(reader.get<vector<string> >(SFK_SAVE_FORMAT, Simu::DEFAULT_SAVE_FORMATS));
        // Read optional scalar values
        simu.setEndValue(reader.get<double>(SFK_END_VALUE, Simu::DEFAULT_END_VALUE));
        simu.setdt(reader.get<double>(SFK_DT, Simu::DEFAULT_DT));
        simu.setTargetdQ(reader.get<double>(SFK_TARGET_DQ, Simu::DEFAULT_TARGET_DQ));
        simu.setMaxdt(reader.get<double>(SFK_MAX_DT, Simu::DEFAULT_MAX_DT));
        simu.setMaxError(reader.get<double>(SFK_MAX_ERROR, Simu::DEFAULT_MAX_ERROR));
        // int values
        simu.setOutputEnergy(reader.get<int>(SFK_OUTPUT_ENERGY, Simu::DEFAULT_OUTPUT_ENERGY));
        simu.setOutputFormat(reader.get<int>(SFK_OUTPUT_FORMAT, Simu::DEFAULT_OUTPUT_FORMAT));
        simu.setSaveIter(reader.get<int>(SFK_SAVE_ITER, Simu::DEFAULT_SAVE_ITER));
        simu.setAsseblyThreadCount(reader.get<unsigned int>(SFK_NUM_ASSEMBLY_THREADS, Simu::DEFAULT_NUM_ASSEMBLY_THREADS));
        simu.setMatrixSolverThreadCount(reader.get<unsigned int>(SFK_NUM_MATRIX_SOLVER_THREADS, Simu::DEFAULT_NUM_MATRIX_SOLVER_THREADS));
        // Number arrays
        simu.setStretchVector(reader.get<vector<double> >(SFK_STRETCH_VECTOR, Simu::DEFAULT_STRETCH_VECTOR));
        simu.setdtLimits(reader.get<vector<double> > (SFK_DT_LIMITS, Simu::DEFAULT_DT_LIMITS));
        simu.setdtFunction(reader.get<vector<double> >(SFK_DT_FUNCTION, Simu::DEFAULT_DT_FUNCTION));
        simu.setRegularGridSize(reader.get<vector<idx> >(SFK_REGULAR_GRID_SIZE, Simu::DEFAULT_REGULAR_GRID_SIZE));

        eventList.setSaveIter(simu.getSaveIter());
    } catch (ReaderError e) {
        e.printError();
    }
}//end void readSimu

void readBoxes(Boxes &boxes, Reader& reader) {
    using std::string;
    const int MAX_NUM_BOXES = 100;
    // Loop over all possible boxes ad try to read
    for (int boxNum = 1; boxNum < MAX_NUM_BOXES; boxNum++) {
        string typeKey = wildcardToNum(SFK_BOX_TYPE, boxNum);
        if (reader.containsKey(typeKey)) {
            // found box with current number
            // Make keys for current box number
            string paramsKey = wildcardToNum(SFK_BOX_PARAMS, boxNum);
            string xKey = wildcardToNum(SFK_BOX_X, boxNum);
            string yKey = wildcardToNum(SFK_BOX_Y, boxNum);
            string zKey = wildcardToNum(SFK_BOX_Z, boxNum);
            string tiltKey = wildcardToNum(SFK_BOX_TILT, boxNum);
            string twistKey = wildcardToNum(SFK_BOX_TWIST, boxNum);
            // add box
            boxes.addBox(boxNum,
                         reader.get<string>(typeKey, Box::DEFAULT_TYPE),
                         reader.get<vector<double>>(paramsKey, Box::DEFAULT_PARAMS),
                         reader.get<vector<double>>(xKey, Box::DEFAULT_X_Y_Z),
                         reader.get<vector<double>>(yKey, Box::DEFAULT_X_Y_Z),
                         reader.get<vector<double>>(zKey, Box::DEFAULT_X_Y_Z),
                         reader.get<vector<double>>(tiltKey, Box::DEFAULT_TILT_TWIST),
                         reader.get<vector<double>>(twistKey, Box::DEFAULT_TILT_TWIST));
        }
    } // end for boxNum
}// end void readBoxes



void readAlignment(Alignment& alignment, Reader& reader) {
    /*!Reads surface anchoring settings from file*/
    const int MAX_NUM_SURFACES = 99;
    // Loop over all possible FixLC numbers
    for (int i = 0; i < MAX_NUM_SURFACES; i++) {
        //
        // Use anchoring type key to determine whether FixLC with current
        // number has been defined
        string anchoringKey = wildcardToNum(SFK_FIXLC_ANCHORING, i);
        if (reader.containsKey(anchoringKey)) {
            //
            // create rest of key for this FixLC surface number
            string strengthKey = wildcardToNum(SFK_FIXLC_STRENGTH, i);
            string easyKey = wildcardToNum(SFK_FIXLC_EASY, i);
            string k1Key = wildcardToNum(SFK_FIXLC_K1, i);
            string k2Key = wildcardToNum(SFK_FIXLC_K2, i);
            string paramsKey = wildcardToNum(SFK_FIXLC_PARAMS, i);
            // add surface to collection of surfaces
            alignment.addSurface(i,
                                 reader.get<string>(anchoringKey, Surface::DEFAULT_ANCHORING_TYPE),
                                 reader.get<double>(strengthKey, Surface::DEFAULT_ANCHORING_STRENGTH),
                                 reader.get<vector<double>>(easyKey, Surface::DEFAULT_ANCHORING_EASY),
                                 reader.get<double>(k1Key, Surface::DEFAULT_ANCHORING_K1),
                                 reader.get<double>(k2Key, Surface::DEFAULT_ANCHORING_K2),
                                 reader.get<vector<double>>(paramsKey, Surface::DEFAULT_ANCHORING_PARAMS));
        }
    }
} //end void readAlignment


void readElectrodes(Electrodes &electrodes,
                    EventList &evli,
                    Reader &reader) {
    //
    // Potential calculation related settings are read here.
    // This includes electrode switching, uniform E-field and
    // relative dielectric permittivity values
    //
    // Loop over range of possible (1-99) Electrodes
    const int MAX_NUM_ELECTRODES = 100;
    for (int electrodeNumber = 1; electrodeNumber < MAX_NUM_ELECTRODES; electrodeNumber++) {
        std::string keyTime = wildcardToNum(SFK_E_TIME, electrodeNumber);
        std::string keyPots = wildcardToNum(SFK_E_POTS, electrodeNumber);
        vector<double> times;
        vector<double> pots;
        if (reader.containsKey(keyTime))
            times = reader.get<vector<double>>();
        if (reader.containsKey(keyPots))
            pots = reader.get<vector<double>>();
        //
        // Make sure number of switching times matches number
        // of switching potentials
        if (times.size() != pots.size()) {
            cerr << "error reading Electrode " << electrodeNumber <<
                    ". Swithing times don't match switching potentials - bye!" << endl;
            std::exit(qlc3d_GLOBALS::ERROR_CODE_BAD_SETTINGS_FILE);
        }
        //
        // Decompose switching times and potentials to separate switching events
        // TODO: move this loop inside Event
        for (size_t i = 0; i < times.size(); i++) {
            SwitchingInstance *switchingInstance = new SwitchingInstance(times.at(i),
                                                                         pots.at(i),
                                                                         electrodeNumber - 1);
            Event *event = new Event(EVENT_SWITCHING,
                                     times.at(i));
            event->setEventDataPtr(static_cast<void *>(switchingInstance));
            evli.insertTimeEvent(event);
        }
    }// end for electrodeNumber
    //
    // Read uniform bulk E-fields here
    if (reader.containsKey(SFK_E_FIELD)) {
        electrodes.setEField(reader.get<vector<double>>());
        // A dummy switching event at time 0 is needed for uniform E-fields.
        // TODO: make this a part of the Event class
        Event *swEvent = new Event(EVENT_SWITCHING , 0.0);
        SwitchingInstance *si = new SwitchingInstance(0 , 0,
                                                      SwitchingInstance::UNIFORM_E_FIELD);
        swEvent->setEventDataPtr(static_cast <void *>(si));
        evli.insertTimeEvent(swEvent);
    }
    //
    // Read dielectric permittivities
    if (reader.containsKey(SFK_EPS_DIELECTRIC))
        electrodes.eps_dielectric = reader.get<vector<double>>();
    // Finally set internal flags when all other settings are done
    // electrodes.setImplicitVariables();
}


void readRefinement(Reader &reader,
                    EventList &evli) {
    // Reads mesh refinement settings and adds any refinement
    // event to event list
    // First read periodically occurring refinement settings
    if (reader.containsKey("RefRefIter"))
        evli.setRepRefIter(reader.get<int>());
    if (reader.containsKey("RepRefTime"))
        evli.setRepRefTime(reader.get<double>());
    // Then read any numbered refinement "objects"
    // Loop over each possible one
    const int MAX_REF_COUNT = 100;
    for (int refNum = 1; refNum < MAX_REF_COUNT; refNum++) {
        string keyBase = "REFINEMENT" + std::to_string(refNum) + ".";
        string key = keyBase + "Type";
        //
        // If found refinement event with current number
        if (reader.containsKey(key)) {
            // Read all fields for this event
            string typeValue;
            vector<long int> iterationsValue;
            vector<double>   timesValue;
            vector<double>   valuesValue;
            vector<double>   xValue;
            vector<double>   yValue;
            vector<double>   zValue;
            typeValue = reader.getValueByKey<string>(keyBase + "Type");
            valuesValue = reader.getValueByKey<vector<double>>(keyBase + "Values");
            xValue = reader.getValueByKey<vector<double>>(keyBase + "X");
            yValue = reader.getValueByKey<vector<double>>(keyBase + "Y");
            zValue = reader.getValueByKey<vector<double>>(keyBase + "Z");
            if (reader.containsKey(keyBase + "Iterations"))
                iterationsValue = reader.get<vector<long int>>();
            if (reader.containsKey(keyBase + "Times"))
                timesValue = reader.get<vector<double>>();
            // Can't have both times and iterations defined for same
            // refinement onject
            if (iterationsValue.size() > 0 && timesValue.size() > 0) {
                std::cerr << "error reading REFINEMENT" + std::to_string(refNum) << std::endl;
                std::cerr << "Can not define both Iterations and Times for same event - bye!" << std::endl;
                std::exit(qlc3d_GLOBALS::ERROR_CODE_BAD_SETTINGS_FILE);
            }
            //
            // Now have all info for current refinement event.
            // Create it and Add to events list
            if (!iterationsValue.empty()) {
                for (auto iter : iterationsValue) {
                    Event *e = new Event(EVENT_REFINEMENT, (unsigned int) iter);
                    e->setEventDataPtr(static_cast<void *>
                                       (RefInfo::make(typeValue,
                                                      iter, 0,
                                                      valuesValue,
                                                      xValue, yValue, zValue)));
                }
            } else if (!timesValue.empty()) {
                for (auto time : timesValue) {
                    Event *e = new Event(EVENT_REFINEMENT, time);
                    e->setEventDataPtr(static_cast<void *>
                                       (RefInfo::make(typeValue,
                                                      0, time,
                                                      valuesValue,
                                                      xValue, yValue, zValue)));
                    evli.insertTimeEvent(e);
                }
            } else {
                RefInfo *refInfo = RefInfo::make(typeValue,
                                                 0, 0, valuesValue,
                                                 xValue, yValue, zValue);
                // REPETITION BY ITERATION ONLY CURRENTLY SUPPORTED
                // TODO Also support repetition by fixed time period
                Event *e = new Event(EVENT_REFINEMENT, refInfo->getRefIter());
                e->setEventDataPtr(static_cast<void *>
                                   (refInfo));
                evli.addRepRefInfo(e);
            }
        } // end if found event by number
    }// end for event numbers
} // end readRefinement

void readSolverSettings(Settings &settings, Reader &reader) {
    if (reader.containsKey(SFK_NUM_ASSEMBLY_THREADS))
        settings.setnThreads(reader.get<int>());
    if (reader.containsKey(SFK_Q_SOLVER))
        settings.setQ_Solver(reader.get<int>());
    if (reader.containsKey(SFK_V_SOLVER))
        settings.setV_Solver(reader.get<int>());
    if (reader.containsKey(SFK_Q_NEWTON_PANIC_ITER))
        settings.setQ_Newton_Panic_Iter(reader.get<int>());
    if (reader.containsKey(SFK_Q_NEWTON_PANIC_COEFF))
        settings.setQ_Newton_Panic_Coeff(reader.get<double>());
    if (reader.containsKey(SFK_Q_PCG_PECONDITIONER))
        settings.setQ_PCG_Preconditioner(reader.get<int>());
    if (reader.containsKey(SFK_Q_PCG_MAXITER))
        settings.setQ_PCG_Maxiter(reader.get<int>());
    if (reader.containsKey(SFK_Q_PCG_TOLER))
        settings.setQ_PCG_Toler(reader.get<double>());
    if (reader.containsKey(SFK_Q_GMRES_PRECONDITIONER))
        settings.setQ_GMRES_Preconditioner(reader.get<int>());
    if (reader.containsKey(SFK_Q_GMRES_MAXITER))
        settings.setQ_GMRES_Maxiter(reader.get<int>());
    if (reader.containsKey(SFK_Q_GMRES_RESTART))
        settings.setQ_GMRES_Restart(reader.get<int>());
    if (reader.containsKey(SFK_Q_GMRES_TOLER))
        settings.setQ_GMRES_Toler(reader.get<double>());
    if (reader.containsKey(SFK_V_PCG_PRECONDITIONER))
        settings.setV_PCG_Preconditioner(reader.get<int>());
    if (reader.containsKey(SFK_V_PCG_MAXITER))
        settings.setV_PCG_Maxiter(reader.get<int>());
    if (reader.containsKey(SFK_V_PCG_TOLER))
        settings.setV_PCG_Toler(reader.get<double>());
    if (reader.containsKey(SFK_V_GMRES_PRECONDITIONER))
        settings.setV_GMRES_Preconditioner(reader.get<int>());
    if (reader.containsKey(SFK_V_GMRES_MAXITER))
        settings.setV_GMRES_Maxiter(reader.get<int>());
    if (reader.containsKey(SFK_V_GMRES_RESTART))
        settings.setV_GMRES_Restart(reader.get<int>());
    if (reader.containsKey(SFK_V_GMRES_TOLER))
        settings.setV_GMRES_Toler(reader.get<double>());
} // end readSolverSettings


void ReadSettings(string settingsFileName,
                  Simu &simu,
                  LC &lc,
                  Boxes &boxes,
                  Alignment &alignment,
                  Electrodes &electrodes,
                  MeshRefinement &meshrefinement, // <--- unused param.
                  EventList &eventList, Settings &settings) {
    Reader reader;
    reader.setCaseSensitivity(false);
    reader.readSettingsFile(settingsFileName);
    try {
        readSimu(simu, eventList,reader);
        readLC(lc, reader);
        readBoxes(boxes, reader);
        readSolverSettings(settings, reader);
        readAlignment(alignment, reader);
        readElectrodes( electrodes, eventList, reader);
        readRefinement(reader, eventList);
    } catch (ReaderError e) {
        e.printError();
    }
    //std::cout << "EXIT OK " << std::endl;
    //exit(1);
}
// end ReadSettings




