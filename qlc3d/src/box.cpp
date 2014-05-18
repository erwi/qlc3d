#include <box.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <reader.h>
#include <vector>

using std::vector;

Box::Box(int boxnum) {

    Type = Normal;
    BoxNumber = boxnum;
    Params = {0,0};
    X = {0,0};
    Y = {0,0};
    Z = {0,0};
    Tilt = {0,0};
    Twist = {0,0};

}

void Box::printBox() {
    printf("\tType: ");
    printf("\n\tparams : %f, %f\n", Params[0], Params[1]);
    printf("\tX\t= [%1.1f, %1.1f]\n\tY\t= [%1.1f, %1.1f]\n\tZ\t= [%1.1f, %1.1f]\n", X[0], X[1], Y[0], Y[1], Z[0], Z[1]);
    printf("\tTilt\t= [%1.1f, %1.1f]\n", Tilt[0], Tilt[1]);
    printf("\tTwist\t= [%1.1f, %1.1f]\n", Twist[0], Twist[1]);
}

/*
void Box::setParams(std::vector<double> p) {
    if (p.size() == 2) {
        Params[0] = p[0];
        Params[1] = p[1];
    } else {
        std::cout << "error, Box::setParams, invalid Params length - bye!" << std::endl;
        exit(1);
    }
}
*/
void Box::setX(std:: vector<double> x) {
    if (x.size() == 2) {
        X[0] = x[0];
        X[1] = x[1];
    } else {
        std::cout << "error, Box::setX, invalid X length - bye!" << std::endl;
        exit(1);
    }
}

void Box::setY(std:: vector<double> y) {
    if (y.size() == 2) {
        Y[0] = y[0];
        Y[1] = y[1];
    } else {
        std::cout << "error, Box::setY, invalid Y length - bye!" << std::endl;
        exit(1);
    }
}

void Box::setZ(std:: vector<double> z) {
    if (z.size() == 2) {
        Z[0] = z[0];
        Z[1] = z[1];
    } else {
        std::cout << "error, Box::setZ, invalid Z length - bye!" << std::endl;
        exit(1);
    }
}

void Box::setTilt(std:: vector<double> tlt) {
    if (tlt.size() == 2) {
        Tilt[0] = tlt[0];
        Tilt[1] = tlt[1];
    } else {
        std::cout << "error, Box::setTilt, invalid Tilt length - bye!" << std::endl;
        exit(1);
    }
}
void Box::setTwist(std:: vector<double> twt) {
    if (twt.size() == 2) {
        Twist[0] = twt[0];
        Twist[1] = twt[1];
    } else {
        std::cout << "error, Box::setX, invalid Twist length - bye!" << std::endl;
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

void Box::setBoxType(std::string &bt) {
    this->TypeString = bt;
    std::transform(bt.begin(), bt.end(), bt.begin(), ::tolower);
    if (bt.compare("normal") == 0)
        Type = Normal;
    else if (bt.compare("random") == 0)
        Type = Random;
    else if (bt.compare("hedgehog") == 0)
        Type = Hedgehog;
    else {
        using std::cout;
        using std::endl;
        cout << "error specyfying Box" << BoxNumber << ".Type as :\"" << TypeString << "\"" << endl;
        cout << "possible types are: \n" << "\tnormal\n" << "\trandom\n"
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

void Boxes::printBoxes() {
    std::vector<Box *>::iterator i;
    printf("%i Boxes:\n", n_Boxes);
    for (int i = 0 ; i < n_Boxes; i++) {
        printf("Box%i :\n", i + 1);
        box[i]->printBox();
    }
}



