//
// Created by eero on 04/04/2021.
//

#ifndef PROJECT_QLC3D_SETTINGS_READER_H
#define PROJECT_QLC3D_SETTINGS_READER_H

#include <filesystem>
#include <string>
#include <memory>
#include "simu.h"
#include "lc.h"
#include "meshrefinement.h"
#include "electrodes.h"
//class MeshRefinement;

class SettingsReader {
    std::filesystem::path fileName_;
    std::unique_ptr<Simu> simu_;
    std::unique_ptr<LC> lc_;

    std::unique_ptr<MeshRefinement> meshRefinement_;

    /*!
     * Reads the contents of the settings file. Be sure to call this before accessing any of the read settings
     * data object.
     * @throws ReaderError (defined in reader.h) in case anything unexpected, like bad file format etc.
     */
    void read();
    void readSimu(Reader &reader);
    void readLC(Reader &reader);
    void readElectrodes(Reader &reader);
    /** Reads optional mesh refinement configuration */
    void readRefinement(Reader &reader);
    //! utility assertion for checking some input file format related stuff. Throws ReaderError
    void assertTrue(bool condition, const std::string &msg);
public:
    SettingsReader(const std::filesystem::path &fileName);

    std::unique_ptr<Simu> simu();
    std::unique_ptr<LC> lc();
    std::unique_ptr<MeshRefinement> refinement();
    std::unique_ptr<Electrodes> electrodes();
};
#endif //PROJECT_QLC3D_SETTINGS_READER_H
