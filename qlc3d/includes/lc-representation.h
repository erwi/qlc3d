//
// Created by eero on 10/04/2021.
//

#ifndef PROJECT_QLC3D_LC_REPRESENTATION_H
#define PROJECT_QLC3D_LC_REPRESENTATION_H

#include <cmath>

namespace qlc3d {
    const static double rt2 { sqrt(2.) };
    const static double rt6 { sqrt(6.) };

    class Director {
        double nx_, ny_, nz_, S_;
    public:
        Director() = delete;
        Director(const double &nx, const double &ny, const double &nz, const double S);
        const double &nx() const { return nx_; }
        const double &ny() const { return ny_; }
        const double &nz() const { return nz_; }
        const double &S() const { return S_; }

        Director static fromRadianAngles(const double &tiltRadians, const double &twistRadians, const double &S);
        Director static fromDegreeAngles(const double &tiltDegrees, const double &twistDegrees, const double &S);
    };

    /*!
     *Traceless symmetric 3x3 tensoratic const
     *
     * <p> q1   q3    q5 </p>
     * <p> q3   q2    q4 </p>
     * <p> q5   q4    (-q1 - q2) </p>
     */

    class QTensor {
        double q1_, q2_, q3_, q4_, q5_;
    public:
        QTensor() = delete;
        QTensor(const double &q1, const double &q2, const double &q3, const double &q4, const double &q5);
        // TODO: needs copy constructor?

        QTensor static fromDirector(const Director &director);
        [[nodiscard]] Director toDirector() const;
        [[nodiscard]] const double &q1() const { return q1_; }
        [[nodiscard]] const double &q2() const { return q2_; }
        [[nodiscard]] const double &q3() const { return q3_; }
        [[nodiscard]] const double &q4() const { return q4_; }
        [[nodiscard]] const double &q5() const { return q5_; }
    };

    /*!
     * TTensor representation, used in simulations instead of the Q-Tensor. TODO: add link to reference.
     */
    class TTensor {
        double t1_, t2_, t3_, t4_, t5_;
    public:
        TTensor(const double &t1, const double &t2, const double &t3, const double &t4, const double &t5);
        [[nodiscard]] TTensor static fromQTensor(const QTensor &q);
        [[nodiscard]] TTensor static fromDirector(const Director &d);
        [[nodiscard]] QTensor toQTensor() const;
        [[nodiscard]] Director toDirector() const;
        [[nodiscard]] const double &t1() const { return t1_; }
        [[nodiscard]] const double &t2() const { return t2_; }
        [[nodiscard]] const double &t3() const { return t3_; }
        [[nodiscard]] const double &t4() const { return t4_; }
        [[nodiscard]] const double &t5() const { return t5_; }
    };
}

#endif //PROJECT_QLC3D_LC_REPRESENTATION_H
