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
using std::cerr;
using std::endl;

std::string wildcardToNum(const std::string& base, int num) {
    /*!Replaces wildcard character in input string with input
number, returning a copy*/
    size_t idx = base.find_first_of(SFK_WILDCARD);
    if (idx == std::string::npos) {
        cerr << "error in " << __PRETTY_FUNCTION__ << endl;
        std::exit(1);
    }
    std::string newStr = base.substr(0,idx);
    newStr+= std::to_string(num);
    newStr+= base.substr(idx+1 , std::string::npos);
    return newStr;
}


void readLC(LC& lc,Reader& reader) {
    try {
        if (reader.containsKey(SFK_K11))
            lc.K11 = reader.get<double>();
        if (reader.containsKey(SFK_K22))
            lc.K22 = reader.get<double>();
        if (reader.containsKey(SFK_K33))
            lc.K33 = reader.get<double>();
        if (reader.containsKey(SFK_p0))
            lc.p0 = reader.get<double>();
        if (reader.containsKey(SFK_A))
            lc.A = reader.get<double>();
        if (reader.containsKey(SFK_B))
            lc.B = reader.get<double>();
        if (reader.containsKey(SFK_C))
            lc.C = reader.get<double>();
        if (reader.containsKey(SFK_EPS_PAR))
            lc.eps_par = reader.get<double>();
        if (reader.containsKey(SFK_EPS_PER))
            lc.eps_per = reader.get<double>();
        if (reader.containsKey(SFK_E1))
            lc.e11 = reader.get<double>();
        if (reader.containsKey(SFK_E3))
            lc.e33 = reader.get<double>();
        if (reader.containsKey(SFK_GAMMA1))
            lc.gamma1 = reader.get<double>();

        lc.convert_params_n2Q();
    } catch (ReaderError e) {
        e.printError();
        std::exit(qlc3d_GLOBALS::ERROR_CODE_BAD_SETTINGS_FILE);
    }

}//end void readLC

void readSimu(Simu &simu, Reader &reader) {
    /*! loads values from settings file reader */
    try {
        simu.MeshName = reader.getValueByKey<string>(SFK_MESH_NAME);
        std::string key;
        // Read optional string values
        if (reader.containsKey(SFK_LOAD_Q))
            simu.LoadQ = reader.get<string>();
        if (reader.containsKey(SFK_SAVE_DIR))
            simu.SaveDir = reader.get<string>();
        if (reader.containsKey(SFK_Q_MATRIX_SOLVER))
            simu.setQMatrixSolver(reader.get<string>());
        // string arrays
        if (reader.containsKey(SFK_SAVE_FORMAT))
            simu.setSaveFormats(reader.get<vector<string> >());

        // Read optional scalar values
        if (reader.containsKey(SFK_END_VALUE))
            simu.setEndValue(reader.get<double>());
        if (reader.containsKey(SFK_DT))
            simu.setdt(reader.get<double>());
        if (reader.containsKey(SFK_TARGET_DQ))
            simu.setTargetdQ(reader.get<double>());
        if (reader.containsKey(SFK_MAX_DT))
            simu.setMaxdt(reader.get<double>());
        if (reader.containsKey(SFK_MAX_ERROR))
            simu.setMaxError(reader.get<double>());

        if (reader.containsKey(SFK_OUTPUT_ENERGY))
            simu.setOutputEnergy(reader.get<int>());
        if (reader.containsKey(SFK_OUTPUT_FORMAT))
            simu.setOutputFormat(reader.get<int>());
        if (reader.containsKey(SFK_SAVE_ITER))
            simu.setSaveIter(reader.get<int>());
        if (reader.containsKey(SFK_NUM_ASSEMBLY_THREADS))
            simu.setAsseblyThreadCount(reader.get<unsigned int>());
        if (reader.containsKey(SFK_NUM_MATRIX_SOLVER_THREADS))
            simu.setMatrixSolverThreadCount(reader.get<unsigned int>());

        // Number arrays
        if (reader.containsKey(SFK_STRETCH_VECTOR))
            simu.setStretchVector(reader.get<vector<double> >());
        if (reader.containsKey(SFK_DT_LIMITS))
            simu.setdtLimits(reader.get<vector<double> > ());
        if (reader.containsKey(SFK_DT_FUNCTION))
            simu.setdtFunction(reader.get<vector<double> >());
        if (reader.containsKey(SFK_REGULAR_GRID_SIZE))
            simu.setRegularGridSize(reader.get<vector<idx> >());
    } catch (ReaderError e) {
        e.printError();
    }
}//end void readSimu



void readBoxes(Boxes &boxes, Reader& reader) {
    using std::string;
    const int MAX_NUM_BOXES = 100;
    // Loop over all possible boxes ad try to read
    for (int boxNum = 1; boxNum < MAX_NUM_BOXES; boxNum++) {
        string key = "BOX"+std::to_string(boxNum) + ".Params";
        // if box with current number found, read it fully
        if (reader.containsKey(key)) {
            Box *box = new Box(boxNum);
            string keyBase = "BOX"+std::to_string(boxNum) + ".";
            string type = reader.getValueByKey<string>(keyBase+"Type");
            box->setBoxType(type);
            box->Params = reader.getValueByKey<vector<double>>(keyBase+"Params");
            box->X = reader.getValueByKey<vector<double>>(keyBase+"X");
            box->Y = reader.getValueByKey<vector<double>>(keyBase+"Y");
            box->Z = reader.getValueByKey<vector<double>>(keyBase+"Z");
            box->Tilt = reader.getValueByKey<vector<double>>(keyBase+"Tilt");
            box->Twist = reader.getValueByKey<vector<double>>(keyBase+"Twist");
            boxes.addBox(box);
        }
    }
}// end void readBoxes



void readAlignment(Alignment& alignment, Reader& reader) {
    /*!Reads surface anchoring settings from file*/
    const int MAX_NUM_SURFACES = 99;
    // Loop over all possible FixLC numbers
    for (int i = 0; i < MAX_NUM_SURFACES; i++) {
        string keyBase = "FIXLC"+std::to_string(i) + ".";
        string key = keyBase + "Anchoring";
        if (reader.containsKey(key)) {
            string name = reader.getValueByKey<string>(key);
            // Create reasonable default values for surface parameters
            double strength  = Surface::DEFAULT_ANCHORING_STRENGTH;
            double K1 = 1.0;
            double K2 = 1.0;
            vector <double> easyAngle = {0,0,0};
            vector<double> params;
            // Read optional values to overwrite defaults
            if (reader.containsKey(keyBase + "Strength"))
                strength = reader.get<double>();
            if (reader.containsKey(keyBase + "K1"))
                K1 = reader.get<double>();
            if (reader.containsKey(keyBase + "K2"))
                K2 = reader.get<double>();
            if (reader.containsKey(keyBase + "Easy"))
                easyAngle = reader.get<vector<double>>();
            if (reader.containsKey(keyBase+"Params"))
                params = reader.get<vector<double>>();
            // Create new Surface object and set all values
            Surface *s = new Surface(i);
            s->setAnchoringType(name);
            s->setEasyAngles(easyAngle);
            s->setStrength(strength);
            s->setK1(K1);
            s->setK2(K2);
            s->Params = params;
            alignment.addSurface(s);
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
    electrodes.setImplicitVariables();
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
    if (reader.containsKey("nThreads"))
        settings.setnThreads(reader.get<int>());
    if (reader.containsKey("Q_Solver"))
        settings.setQ_Solver(reader.get<int>());
    if (reader.containsKey("V_Solver"))
        settings.setV_Solver(reader.get<int>());
    if (reader.containsKey("Q_Newton_Panic_Iter"))
        settings.setQ_Newton_Panic_Iter(reader.get<int>());
    if (reader.containsKey("Q_Newton_Panic_Coeff"))
        settings.setQ_Newton_Panic_Coeff(reader.get<double>());
    if (reader.containsKey("Q_PCG_Preconditioner"))
        settings.setQ_PCG_Preconditioner(reader.get<int>());
    if (reader.containsKey("Q_PCG_Maxiter"))
        settings.setQ_PCG_Maxiter(reader.get<int>());
    if (reader.containsKey("Q_PCG_Toler"))
        settings.setQ_PCG_Toler(reader.get<double>());
    if (reader.containsKey("Q_GMRES_Preconditioner"))
        settings.setQ_GMRES_Preconditioner(reader.get<int>());
    if (reader.containsKey("Q_GMRES_Maxiter"))
        settings.setQ_GMRES_Maxiter(reader.get<int>());
    if (reader.containsKey("Q_GMRES_Restart"))
        settings.setQ_GMRES_Restart(reader.get<int>());
    if (reader.containsKey("Q_GMRES_Toler"))
        settings.setQ_GMRES_Toler(reader.get<double>());
    if (reader.containsKey("V_PCG_Preconditioner"))
        settings.setV_PCG_Preconditioner(reader.get<int>());
    if (reader.containsKey("V_PCG_Maxiter"))
        settings.setV_PCG_Maxiter(reader.get<int>());
    if (reader.containsKey("V_PCG_Toler"))
        settings.setV_PCG_Toler(reader.get<double>());
    if (reader.containsKey("V_GMRES_Preconditioner"))
        settings.setV_GMRES_Preconditioner(reader.get<int>());
    if (reader.containsKey("V_GMRES_Maxiter"))
        settings.setV_GMRES_Maxiter(reader.get<int>());
    if (reader.containsKey("V_GMRES_Restart"))
        settings.setV_GMRES_Restart(reader.get<int>());
    if (reader.containsKey("V_GMRES_Toler"))
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
        readSimu(simu, reader);
        readLC(lc, reader);
        readBoxes(boxes, reader);
        readSolverSettings(settings, reader);
        readAlignment(alignment, reader);
        readElectrodes( electrodes, eventList, reader);
        readRefinement(reader, eventList);
        // readAlignment(*alignment, reader);
    } catch (ReaderError e) {
        e.printError();
    }
    std::cout << "EXIT OK " << std::endl;
    exit(1);
}
// end ReadSettings




