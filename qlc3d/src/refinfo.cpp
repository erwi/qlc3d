#include <refinfo.h>
#include <assert.h>
#include <globals.h>
void RefInfo::setType(std::string s) {
// TESTS WHETHER "REFINEMENT.Type" STRING IN SETTNGS FILE IS VALID.
// IF YES, SETS type_ ENUMERATOR FOR THIS OBJECT
// IF NOT, TERMINATES PROGRAM
    // MAKE LOWERCASE
    std::transform(s.begin(), s.end(), s.begin(), std::ptr_fun<int, int>(std::tolower));
    if (! s.compare("change")) {
        _type = Change;
        return;
    } else if (! s.compare("sphere")) {
        _type = Sphere;
        return;
    } else if (! s.compare("box")) {
        _type = Box;
        return;
    } else {
        std::cout << "error - bad refinement type: " << s << " bye!" << std::endl;
        exit(1);
    }
}
void RefInfo::setRefIter() {
// DETERMINES NUMBER OF REFINEMENT ITERATIONS THIS OBJECT DESCRIBES
// AND SETS THE refIter_ VARIABLE
// DEPENDING ON THE type_, NUMBER OF REFINEMENT ITERATIONS IS DETERMINED BY
// DIFFERENT CONDITIONS
    switch (_type) {
    case (Change): {
        _refIter = (int) _values.size();
        break;
    }
    case (Sphere): {
        _refIter = (int) _values.size();
        // MAKE SURE NOT EMPTY AND EQUAL NUMBER OF COORDS HAVE BEEN SPECIFIED
        break;
    }
    case (Box): {
        // FOR A BOX, refIter DEPENDS ON NUMBER OF COORDS SPECIFIED
        size_t numx = _x.size();
        _refIter = (int) numx / 2;
        break;
    }
    default:
        printf("error in %s, unhandled refinement type - bye\n", __func__);
        exit(1);
    }
}


RefInfo::RefInfo(const std::string &Type):
    _type(None),
    _iter(0),
    _time(0),
    _refIter(0),
    _x(1, 0),
    _y(1, 0),
    _z(1, 0)

{
    setType(Type);
}

RefInfo::RefInfo(const RefInfo &other):
    _type(other._type),
    _iter(other._iter),
    _time(other._time),
    _refIter(other._refIter),
    _values(other._values),
    _x(other._x),
    _y(other._y),
    _z(other._z) {
}



void RefInfo::setValues(std::vector<double> &values) {
    _values.clear();
    _values.insert(_values.begin(), values.begin(), values.end());
    this->setRefIter();
}



void RefInfo::setCoords(const std::vector<double> &x,
                        const std::vector<double> &y,
                        const std::vector<double> &z) {
    _x.clear();
    _y.clear();
    _z.clear();
    _x.insert(_x.begin(), x.begin(), x.end());
    _y.insert(_y.begin(), y.begin(), y.end());
    _z.insert(_z.begin(), z.begin(), z.end());
    this->setRefIter();
}


void RefInfo::getCoords(std::vector<double> &x,
                        std::vector<double> &y,
                        std::vector<double> &z) const {
    x = _x;
    y = _y;
    z = _z;
}

void RefInfo::getCoord(double &x, double &y, double &z) const {
    assert(_x.size());
    assert(_y.size());
    assert(_z.size());
    x = _x[0];
    y = _y[0];
    z = _z[0];
}

double RefInfo::getValue(const size_t i) const {
#ifdef DEBUG
    assert(i <= _values.size());
#endif
    return _values[i];
}

void RefInfo::printRefInfo(FILE *fid) const {
    std::cout << "Iteration = " << this->_iter << std::endl;
    std::cout << "Time = " << this->_time << std::endl;
    std::cout << "Type = " << this->_type << std::endl;
    for (size_t i = 0 ; i < _values.size() ; i++)
        std::cout << "values[" << i << "]=" << _values[i] << std::endl;
}

void RefInfo::validate(const RefInfo &refinfo) {
#define ERRORMSG printf("\n\nerror, bad or missing REFINEMENT data - bye!\n")
    switch (refinfo.getType()) {
    case (RefInfo::Change):
        // MAKE SURE VALUES EXIST
        if (!refinfo._values.size()) {
            ERRORMSG;
            exit(1);
        }
        break;
    case (RefInfo::Sphere): {
        // MAKE SURE VALUES AND COORDINATES EXIST
        if ((!refinfo._values.size()) ||
                (!refinfo._x.size()) ||
                (!refinfo._y.size()) ||
                (!refinfo._z.size())) {
            ERRORMSG;
            exit(1);
        }
        break;
    }
    case (RefInfo::Box): {
        // MAKE SURE CORRECT NUMBER OF X,Y AND Z COORDINATES HAVE BEEN SPECIFIED
        size_t numx = refinfo._x.size();
        size_t numy = refinfo._y.size();
        size_t numz = refinfo._z.size();
        if ((!numx) ||  // COORDINATES MUST BE DEFINED
                (!numy) ||
                (!numz)) {
            ERRORMSG;
            exit(1);
        }
        if ((numx != numy) ||   // EQUAL NUMBER OF COORDS. DEFINED
                (numx != numz)) {
            ERRORMSG;
            exit(1);
        }
        if ((numx % 2 != 0) || // TWO COORDS PER BOX
                (numy % 2 != 0) ||
                (numz % 2 != 0)) {
            ERRORMSG;
            exit(1);
        }
        break;
    }
    case (RefInfo::None):
        printf("RefInfo type is None - bye!");
        exit(1);
        break;
    default:
        printf("error in %s, unknonwn refinement type - bye!\n", __func__);
        exit(1);
    }
}

RefInfo *RefInfo::make(const std::string &type,
                       long iteration,
                       double time,
                       std::vector<double> &values,
                       std::vector<double> &x,
                       std::vector<double> &y,
                       std::vector<double> &z) {
    // convenience factory method for making and validating
    // RefInfo object pointers.
    // TODO: should return a smart pointer instead
    //
    // A refinement event cannot happen both on a given iteration
    // AND given time instance
    if (time != 0 && iteration != 0) {
        std::cerr << "error, cannot have both iteration and time as non-zero in "
                  << __PRETTY_FUNCTION__ << std::endl;
        std::exit(qlc3d_GLOBALS::ERROR_CODE_BAD_SETTINGS_FILE);
    }

    RefInfo *ret = new RefInfo(type);
    ret->setIteration(iteration);
    ret->setTime(time);
    ret->setValues(values);
    ret->setCoords(x,y,z);
    validate(*ret);
    return ret;
}

