//
// Created by eero on 07/04/2021.
//
#include "test-util.h"

#include <filesystem>
#include <vector>
#include <fstream>

using namespace TestUtil;

//<editor-fold desc=TemporaryFile>
TemporaryFile::TemporaryFile() {
    name_ = std::tmpnam(nullptr);
}

TemporaryFile::~TemporaryFile() {
    using namespace std::filesystem;
    if (exists(name_)) {
        remove(name_);
    }
}

TestUtil::TemporaryFile TemporaryFile::empty() {
    return TemporaryFile::withContents("");
}

TestUtil::TemporaryFile TemporaryFile::withContents(const std::string &fileContents) {
    TemporaryFile f;
    FILE *fid = fopen(f.name().c_str(), "wt");
    fprintf(fid, "%s", fileContents.c_str());
    fclose(fid);
    return f;
}

std::vector<std::string> TemporaryFile::readContentsAsText() const {
    using namespace std;
    ifstream fin(name());

    vector<string> lines;
    string line;
    while (getline(fin, line)) {
        lines.emplace_back(string(line));
    }
    fin.close();
    return lines;
}
//</editor-fold>