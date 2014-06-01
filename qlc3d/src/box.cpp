#include <box.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <reader.h>
#include <vector>

const std::string Box::DEFAULT_TYPE = "Normal";
const std::vector<double> Box::DEFAULT_PARAMS = {};
const std::vector<double> Box::DEFAULT_X_Y_Z = {0,0};
const std::vector<double> Box::DEFAULT_TILT_TWIST = {0,0};

using std::vector;
using std::cerr;
Box::Box(int boxnum) {
    //Type = Normal;
    setBoxType(Box::DEFAULT_TYPE);
    BoxNumber = boxnum;
    Params = Box::DEFAULT_PARAMS;
    X = Box::DEFAULT_X_Y_Z;
    Y = Box::DEFAULT_X_Y_Z;
    Z = Box::DEFAULT_X_Y_Z;;
    Tilt = Box::DEFAULT_TILT_TWIST;
    Twist = Box::DEFAULT_TILT_TWIST;
}

void Box::setX(std:: vector<double> x) {
    if (x.size() == 2) {
        X[0] = x[0];
        X[1] = x[1];
    } else {
        std::cerr << "error, Box::setX, invalid X length - bye!" << std::endl;
        exit(1);
    }
}

void Box::setY(std:: vector<double> y) {
    if (y.size() == 2) {
        Y[0] = y[0];
        Y[1] = y[1];
    } else {
        std::cerr << "error, Box::setY, invalid Y length - bye!" << std::endl;
        exit(1);
    }
}

void Box::setZ(std:: vector<double> z) {
    if (z.size() == 2) {
        Z[0] = z[0];
        Z[1] = z[1];
    } else {
        std::cerr << "error, Box::setZ, invalid Z length - bye!" << std::endl;
        exit(1);
    }
}

void Box::setTilt(std:: vector<double> tlt) {
    if (tlt.size() == 2) {
        Tilt[0] = tlt[0];
        Tilt[1] = tlt[1];
    } else {
        std::cerr << "error, Box::setTilt, invalid Tilt length - bye!" << std::endl;
        exit(1);
    }
}
void Box::setTwist(std:: vector<double> twt) {
    if (twt.size() == 2) {
        Twist[0] = twt[0];
        Twist[1] = twt[1];
    } else {
        std::cerr << "error, Box::setTwist, invalid Twist length - bye!" << std::endl;
        exit(1);
    }
}
bool Box::isInBox(double *coords) {
    double x = coords[0];
    double y = coords[1];
    double z = coords[2];
    //  test if coordinate is outside of box and return false
    if ((x < this->X[0]) || (x > this->X[1])) return false;      // if smaller than minimum or larger than maximum...
    if ((y < this->Y[0]) || (y > this->Y[1])) return false;
    if ((z < this->Z[0]) || (z > this->Z[1])) return false;
    // otherwise return true
    return true;
}

void Box::setBoxType(const std::string &bt) {
    this->TypeString = bt;
    //std::transform(bt.begin(), bt.end(), bt.begin(), ::tolower);
    if (bt.compare("normal") == 0)
        Type = Normal;
    else if (bt.compare("random") == 0)
        Type = Random;
    else if (bt.compare("hedgehog") == 0)
        Type = Hedgehog;
    else {
        using std::cout;
        using std::endl;
        cerr << "error specyfying Box" << BoxNumber << ".Type as :\"" << TypeString << "\"" << endl;
        cerr << "possible types are: \n" << "\tnormal\n" << "\trandom\n"
             << "\thedgehog" << endl;
        exit(1);
    }
}
//===================================================================
Boxes::Boxes() {
    n_Boxes = 0;
}

Boxes::~Boxes() {
    std::vector<Box *>::iterator itr;
    for (itr = box.begin() ; itr != box.end() ; ++itr)
        delete(*itr);
}

void Boxes::addBox(Box *b) {
    box.push_back(b);
    n_Boxes ++;
}

void Boxes::addBox(const int &boxNum,
                   const std::string &boxType,
                   const std::vector<double> &params,
                   const std::vector<double> &x,
                   const std::vector<double> &y,
                   const std::vector<double> &z,
                   const std::vector<double> &tilt,
                   const std::vector<double> &twist) {
    Box *b = new Box(boxNum);
    b->Params = params;
    b->setX(x);
    b->setY(y);
    b->setZ(z);
    b->setTilt(tilt);
    b->setTwist(twist);
    this->addBox(b);
}


