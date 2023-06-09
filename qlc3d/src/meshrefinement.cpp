#include <meshrefinement.h>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <util/logging.h>

//<editor-fold RefInfo>
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
        throw std::invalid_argument(fmt::format("Unknown refinement type {}", s));
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
            throw std::runtime_error(fmt::format("Unhandled refinement type in {}, {}", __FILE__, __func__));
    }
}


RefInfo::RefInfo(const std::string &Type):
        _type(None),
        _iter(0),
        _time(0),
        _refIter(0),
        _x(1, 0),
        _y(1, 0),
        _z(1, 0) {
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

void RefInfo::setValues(const std::vector<double> &values) {
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

double RefInfo::getValue(const size_t i) const {
    assert(i <= _values.size());
    return _values[i];
}

void RefInfo::validate(const RefInfo &refinfo) {
    switch (refinfo.getType()) {
        case (RefInfo::Change):
            // MAKE SURE VALUES EXIST
            if (refinfo._values.empty()) {
                throw std::invalid_argument("Refinement values field not set.");
            }
            break;
        case (RefInfo::Sphere): {
            // MAKE SURE VALUES AND COORDINATES EXIST
            if ((!refinfo._values.size()) ||
                (!refinfo._x.size()) ||
                (!refinfo._y.size()) ||
                (!refinfo._z.size())) {
                throw std::invalid_argument("Refinement values or x, y, z bounds fields not set.");
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
                throw std::invalid_argument("Refinement x, y, z bounds fields not set.");
            }
            if ((numx != numy) ||   // EQUAL NUMBER OF COORDS. DEFINED
                (numx != numz)) {
                throw std::invalid_argument("Refinement x, y, z bounds size mismatch.");
            }
            if ((numx % 2 != 0) || // TWO COORDS PER BOX
                (numy % 2 != 0) ||
                (numz % 2 != 0)) {
                throw std::invalid_argument("Unexpected refinement x, y, z bounds sizes.");
            }
            break;
        }
        case (RefInfo::None):
            throw std::runtime_error(fmt::format("RefInfo typw is None in {}, {}.", __FILE__, __func__ ));
        default:
            throw std::runtime_error(fmt::format("Unhandled refinement type in {}, {}.", __FILE__, __func__));
    }
}

RefInfo *RefInfo::make(const std::string &type,
                       long iteration,
                       double time,
                       const std::vector<double> &values,
                       const std::vector<double> &x,
                       const std::vector<double> &y,
                       const std::vector<double> &z) {
    // convenience factory method for making and validating
    // RefInfo object pointers.
    // TODO: should return a smart pointer instead
    //
    // A refinement event cannot happen both on a given iteration
    // AND given time instance
    if (time > 0.0 && iteration > 0) {
        throw std::runtime_error("can not have both iteration and time set for RefInfo");
    }

    RefInfo *ret = new RefInfo(type);
    ret->setIteration(iteration);
    ret->setTime(time);
    ret->setValues(values);
    ret->setCoords(x,y,z);
    validate(*ret);
    return ret;
}

RefInfo* RefInfo::ofPeriodicMeshRefinement(const std::string &type, const std::vector<double> &values,
                                           const std::vector<double> &x, const std::vector<double> &y,
                                           const std::vector<double> &z) {
    return RefInfo::make(type, -1, -1, values, x, y, z);
}
//</editor-fold>

//<editor-fold RefinementConfig>
void RefinementConfig::validate() {
    // TODO:
    // check that Type is known
    // check that correct number of values are given
    // etc.
}

void MeshRefinement::setRefinementConfig(std::vector<RefinementConfig> &&ref) {
    refinementConfigs_ = std::move(ref);
}
//</editor-fold RefinementConfig>
