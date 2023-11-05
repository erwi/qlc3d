#ifndef PROJECT_QLC3D_LC_REPRESENTATION_H
#define PROJECT_QLC3D_LC_REPRESENTATION_H

#include <cmath>
#include <string>
class Vec3;

namespace qlc3d {
    const static double rt2 { sqrt(2.) };
    const static double rt6 { sqrt(6.) };

    class Director {
        double nx_, ny_, nz_, S_;
    public:
        Director() = delete;
        Director(const double &nx, const double &ny, const double &nz, const double S);
        Director(const Vec3 &n, double S);
        const double &nx() const { return nx_; }
        const double &ny() const { return ny_; }
        const double &nz() const { return nz_; }
        const double &S() const { return S_; }
        [[nodiscard]] Vec3 vector() const;

        Director static fromRadianAngles(const double &tiltRadians, const double &twistRadians, const double &S);
        Director static fromDegreeAngles(const double &tiltDegrees, const double &twistDegrees, const double &S);

        double tiltRadians() const;
        double tiltDegrees() const;

        double twistRadians() const;
        double twistDegrees() const;
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
     * TTensor: traceless representation of Q-tensor, used in simulations instead of the Q-Tensor.
     */
    class TTensor {
        double t1_, t2_, t3_, t4_, t5_; // todo array of length 5

    public:
      /** Un-initialised TTensor */
      TTensor() = default;
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

      void set(double t1, double t2, double t3, double t4, double t5);

      const double& operator[](int i) const {
        switch (i) {
          case 0: return t1_;
          case 1: return t2_;
          case 2: return t3_;
          case 3: return t4_;
          case 4: return t5_;
          default: throw "DielectricPermittivity index out of range " + std::to_string(i);
        }
      }
    };

    class DielectricPermittivity {
      //double e11_, e22_, e33_, e12_, e13_, e23_;
      double e[6];
    public:
      DielectricPermittivity() = default;
      DielectricPermittivity(double e11, double e22, double e33, double e12, double e13, double e23);
      [[nodiscard]] DielectricPermittivity static fromTTensor(const TTensor &t, double S0, double deleps, double eper);
      [[nodiscard]] const double &e11() const { return e[0]; }
      [[nodiscard]] const double &e22() const { return e[1]; }
      [[nodiscard]] const double &e33() const { return e[2]; }
      [[nodiscard]] const double &e12() const { return e[3]; }
      [[nodiscard]] const double &e13() const { return e[4]; }
      [[nodiscard]] const double &e23() const { return e[5]; }
      void set(double e11, double e22, double e33, double e12, double e13, double e23);

      const double& operator[](int i) const { return e[i]; }
    };
}

#endif //PROJECT_QLC3D_LC_REPRESENTATION_H
