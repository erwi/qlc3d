//
// Created by eero on 31/03/2021.
//
#include <iostream>
#include <qlc3d.h>
#include <filesystem>
int main(int narg, char** args) {
    using namespace std;


    //cout <<"temp dir" << std::filesystem

    for (int i = 0; i < narg; i++) {
        cout << args[i] << endl;
    }
    try {
        runQlc3d(narg, args);
    } catch (...) {
        cout << "booo" << endl;
    }

    cout << "EXIT OK " << endl;
}