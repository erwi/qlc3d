#include <refinement/refinement-spec.h>
#include <stdexcept>
#include <algorithm>
#include <cassert>
#include <fmt/format.h>

RefinementSpec::Type RefinementSpec::parseType(const std::string &type) {
    std::string lower = type;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower == "change") return Type::Change;
    if (lower == "sphere") return Type::Sphere;
    if (lower == "box") return Type::Box;
    throw std::invalid_argument(fmt::format("Unknown refinement type '{}'.", type));
}

unsigned int RefinementSpec::calcRefIter(Type type, const std::vector<double> &values, const std::vector<double> &x) {
    switch (type) {
        case Type::Change: return static_cast<unsigned int>(values.size());
        case Type::Sphere: return static_cast<unsigned int>(values.size());
        case Type::Box:    return static_cast<unsigned int>(x.size() / 2);
        default:
            throw std::runtime_error("Unhandled refinement type in RefinementSpec::calcRefIter");
    }
}

void RefinementSpec::validateFields(Type type,
                                    const std::vector<double> &values,
                                    const std::vector<double> &x,
                                    const std::vector<double> &y,
                                    const std::vector<double> &z) {
    switch (type) {
        case Type::Change:
            if (values.empty()) {
                throw std::invalid_argument("Refinement type 'Change' requires non-empty values.");
            }
            break;
        case Type::Sphere:
            if (values.empty() || x.empty() || y.empty() || z.empty()) {
                throw std::invalid_argument("Refinement type 'Sphere' requires non-empty values and x, y, z coordinates.");
            }
            break;
        case Type::Box: {
            size_t nx = x.size(), ny = y.size(), nz = z.size();
            if (nx == 0 || ny == 0 || nz == 0) {
                throw std::invalid_argument("Refinement type 'Box' requires non-empty x, y, z coordinates.");
            }
            if (nx != ny || nx != nz) {
                throw std::invalid_argument("Refinement type 'Box' requires equal numbers of x, y, z coordinates.");
            }
            if (nx % 2 != 0) {
                throw std::invalid_argument("Refinement type 'Box' requires an even number of x, y, z coordinate pairs.");
            }
            break;
        }
        default:
            throw std::runtime_error("Unhandled refinement type in RefinementSpec::validateFields");
    }
}

std::unique_ptr<RefinementSpec> RefinementSpec::makeExplicit(
    const std::string &type,
    long iteration,
    double time,
    std::vector<double> values,
    std::vector<double> x,
    std::vector<double> y,
    std::vector<double> z) {

    // Exactly one of iteration or time must be positive
    if (iteration > 0 && time > 0.0) {
        throw std::invalid_argument("RefinementSpec: cannot specify both iteration and time.");
    }
    if (iteration <= 0 && time <= 0.0) {
        throw std::invalid_argument("RefinementSpec: at least one of iteration or time must be positive.");
    }

    Type t = parseType(type);
    validateFields(t, values, x, y, z);

    auto spec = std::unique_ptr<RefinementSpec>(new RefinementSpec());
    spec->type_ = t;
    spec->periodic_ = false;
    spec->iteration_ = iteration;
    spec->time_ = time;
    spec->values_ = std::move(values);
    spec->x_ = std::move(x);
    spec->y_ = std::move(y);
    spec->z_ = std::move(z);
    spec->refIter_ = calcRefIter(spec->type_, spec->values_, spec->x_);
    return spec;
}

std::unique_ptr<RefinementSpec> RefinementSpec::makePeriodic(
    const std::string &type,
    std::vector<double> values,
    std::vector<double> x,
    std::vector<double> y,
    std::vector<double> z) {

    Type t = parseType(type);
    validateFields(t, values, x, y, z);

    auto spec = std::unique_ptr<RefinementSpec>(new RefinementSpec());
    spec->type_ = t;
    spec->periodic_ = true;
    spec->iteration_ = 0;
    spec->time_ = 0.0;
    spec->values_ = std::move(values);
    spec->x_ = std::move(x);
    spec->y_ = std::move(y);
    spec->z_ = std::move(z);
    spec->refIter_ = calcRefIter(spec->type_, spec->values_, spec->x_);
    return spec;
}

double RefinementSpec::getValue(size_t i) const {
    assert(i < values_.size());
    return values_[i];
}

std::unique_ptr<RefinementSpec> RefinementSpec::clone() const {
    auto copy = std::unique_ptr<RefinementSpec>(new RefinementSpec());
    copy->type_ = type_;
    copy->periodic_ = periodic_;
    copy->iteration_ = iteration_;
    copy->time_ = time_;
    copy->refIter_ = refIter_;
    copy->values_ = values_;
    copy->x_ = x_;
    copy->y_ = y_;
    copy->z_ = z_;
    return copy;
}


