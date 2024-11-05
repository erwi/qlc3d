#include <lc-representation.h>
#include <qlc3d.h> // TODO: deleteme
#include <geom/vec3.h>
using namespace qlc3d;

//<editor-fold desk=QTensor>

// TODO: this should probably be **the** definition, only accessible to QTensor
void eig2(int mv, int n, double *a, double *r) {
    /*!eig calculates eigenvalues and eigenvectors */
    int i, il, ilq, ilr, im, imq, imr, ind, iq, j, l, ll, lm, lq, m, mm, mq;
    double anorm, anrmx, cosx, cosx2, sincs, sinx, sinx2, thr, x, y;
    double theshold = LDBL_EPSILON;
    if (mv) {
        for (i = 1; i < n * n; i++) r[i] = 0.;
        for (i = 0; i < n * n; i += n + 1) r[i] = 1.;
    }
    /* Initial and final norms (anorm & anrmx). */
    anorm = 0.;
    iq = 0;
    for (i = 0; i < n; i++)
        for (j = 0; j <= i; j++) {
            if (j != i) anorm += a[iq] * a[iq];
            ++iq;
        }
    if (anorm > 0.) {
        anorm = sqrt(2. * anorm);
        anrmx = theshold * anorm / n;
        /* Compute threshold and initialise flag. */
        thr = anorm;
        do {
            thr /= n;
            do {
                ind = 0;
                l = 0;
                do {
                    lq = l * (l + 1) / 2;
                    ll = l + lq;
                    m = l + 1;
                    ilq = n * l;
                    do
                        /* Compute sin & cos. */
                    {
                        mq = m * (m + 1) / 2;
                        lm = l + mq;
                        if (fabs(a[lm]) >= thr) {
                            ind = 1;
                            mm = m + mq;
                            x = .5 * (a[ll] - a[mm]);
                            y = -a[lm] / sqrt(a[lm] * a[lm] + x * x);
                            if (x < 0.) y = -y;
                            sinx = y / sqrt(2. * (1. + (sqrt(1. - y * y))));
                            sinx2 = sinx * sinx;
                            cosx = sqrt(1. - sinx2);
                            cosx2 = cosx * cosx;
                            sincs = sinx * cosx;
                            /* Rotate l & m columns. */
                            imq = n * m;
                            for (i = 0; i < n; i++) {
                                iq = i * (i + 1) / 2;
                                if (i != l && i != m) {
                                    if (i < m) im = i + mq;
                                    else im = m + iq;
                                    if (i < l) il = i + lq;
                                    else il = l + iq;
                                    x = a[il] * cosx - a[im] * sinx;
                                    a[im] = a[il] * sinx + a[im] * cosx;
                                    a[il] = x;
                                }
                                if (mv) {
                                    ilr = ilq + i;
                                    imr = imq + i;
                                    x = r[ilr] * cosx - r[imr] * sinx;
                                    r[imr] = r[ilr] * sinx + r[imr] * cosx;
                                    r[ilr] = x;
                                }
                            }
                            x = 2. * a[lm] * sincs;
                            y = a[ll] * cosx2 + a[mm] * sinx2 - x;
                            x = a[ll] * sinx2 + a[mm] * cosx2 + x;
                            a[lm] = (a[ll] - a[mm]) * sincs + a[lm] * (cosx2 - sinx2);
                            a[ll] = y;
                            a[mm] = x;
                        }
                        /* Tests for completion.
                        Test for m = last column. */
                    } while (++m != n);
                    /* Test for l =penultimate column. */
                } while (++l != n - 1);
            } while (ind);
            /* Compare threshold with final norm. */
        } while (thr > anrmx);
    }
}

Director tensortovector(const QTensor &q) {
    const int n = 3;
    double A[6]; // upper diagonal of tensor
    double R[9]; // matrix that will contain 3 x eigenvectors
    A[0] = q.q1();
    A[1] = q.q3();
    A[2] = q.q2();
    A[3] = q.q5();
    A[4] = q.q4();
    A[5] = -q.q1() - q.q2();
    eig2(1, n, A, R);    // find eigenvalue & eigenvectors
    A[1] = A[2];
    A[2] = A[5];
    int iMin, iMax;
    iMin = 0;
    iMax = 0;
    if ((A[0] <= A[1]) && (A[0] <= A[2])) iMin = 0;
    if ((A[1] <= A[0]) && (A[1] <= A[2])) iMin = 1;
    if ((A[2] <= A[0]) && (A[2] <= A[1])) iMin = 2;
    if ((A[0] >= A[1]) && (A[0] >= A[2])) iMax = 0;
    if ((A[1] >= A[0]) && (A[1] >= A[2])) iMax = 1;
    if ((A[2] >= A[0]) && (A[2] >= A[1])) iMax = 2;
    return Director{
            R[3 * iMax + 0],
            R[3 * iMax + 1],
            R[3 * iMax + 2],
            2.0 / 3.0 * (A[iMax] - A[iMin])
    };
}

// <editor-fold desc="Director">
Director::Director(const double &nx, const double &ny, const double &nz, const double S) :
        nx_{nx}, ny_{ny}, nz_{nz}, S_{S} {
    double length = sqrt(nx * nx + ny * ny + nz * nz);
    double lengthError = abs(length - 1);
    if (lengthError > 1e-15) {
        std::string msg = "Non-unit length director. ";
        msg += "Length = " + to_string(length);
        msg += " with (nx, ny, nz) = (" + to_string(nx) + ", " + to_string(ny) + ", " + to_string(nz) + ")";
        throw std::invalid_argument(msg);
    }
}

Director::Director(const Vec3 &n, double S): Director(n.x(), n.y(), n.z(), S) {}

Director Director::fromRadianAngles(const double &tiltRadians, const double &twistRadians, const double &S) {
    // see also definition of vectors v1 and v2 for anchoring in alignment.cpp
    // we should require n = v1 x v2
    return Director{Vec3::fromRadianAngles(tiltRadians, twistRadians), S};
}

Director Director::fromDegreeAngles(const double &tiltDegrees, const double &twistDegrees, const double &S) {
    return fromRadianAngles(M_PI * tiltDegrees / 180., M_PI * twistDegrees / 180., S);
}

double Director::tiltRadians() const {
    return asin(nz_);
}

double Director::tiltDegrees() const {
    return 180. * tiltRadians() / M_PI;
}

double Director::twistRadians() const {
    return atan2(ny_, nx_);
}

double Director::twistDegrees() const {
    return 180. * twistRadians() / M_PI;
}

Vec3 Director::vector() const {
  return Vec3{nx_, ny_, nz_};
}

// </editor-fold>

// <editor-fold desc="QTensor">
QTensor::QTensor(const double &q1, const double &q2, const double &q3, const double &q4, const double &q5) :
        q1_{q1}, q2_{q2}, q3_{q3}, q4_{q4}, q5_{q5} {}

QTensor QTensor::fromDirector(const Director &director) {
    return QTensor{
            (director.S() / 2.) * (3 * director.nx() * director.nx() - 1),
            (director.S() / 2.) * (3 * director.ny() * director.ny() - 1),
            (director.S() / 2.) * (3 * director.nx() * director.ny()),
            (director.S() / 2.) * (3 * director.ny() * director.nz()),
            (director.S() / 2.) * (3 * director.nx() * director.nz())};
}

Director QTensor::toDirector() const {
    return tensortovector(*this);
}
//</editor-fold>

// <editor-fold desc="TTensor">
TTensor::TTensor(const double &t1, const double &t2, const double &t3, const double &t4, const double &t5) :
        t1_{t1}, t2_{t2}, t3_{t3}, t4_{t4}, t5_{t5} {}

TTensor TTensor::fromQTensor(const QTensor &q) {
    return TTensor{
            (q.q1() + q.q2()) * (rt6 / -2.),
            (q.q1() + (q.q1() + q.q2()) / -2.) * rt2,
            q.q3() * rt2,
            q.q4() * rt2,
            q.q5() * rt2
    };
}

TTensor TTensor::fromDirector(const Director &d) {
    return fromQTensor(QTensor::fromDirector(d));
}

QTensor TTensor::toQTensor() const {
    return QTensor{
            -t1_ / rt6 + t2_ / rt2,
            -t1_ / rt6 - t2_ / rt2,
            t3_ / rt2,
            t4_ / rt2,
            t5_ / rt2
    };
    //q[i] =      -a[i] / rt6 + a[i + np] / rt2;
    //q[i + np] =   -a[i] / rt6 - a[i + np] / rt2;
    //q[i + 2 * np] =  a[i + 2 * np] / rt2;
    //q[i + 3 * np] =  a[i + 3 * np] / rt2;
    //q[i + 4 * np] =  a[i + 4 * np] / rt2;
}

Director TTensor::toDirector() const {
    return toQTensor().toDirector();
}

void TTensor::set(double t1, double t2, double t3, double t4, double t5) {
  t1_ = t1;
  t2_ = t2;
  t3_ = t3;
  t4_ = t4;
  t5_ = t5;
}
// </editor-fold>

// <editor-fold desc="DielectricPermittivity">
DielectricPermittivity::DielectricPermittivity(double e11, double e22, double e33, double e12, double e13, double e23) {
  e[0] = e11;
  e[1] = e22;
  e[2] = e33;
  e[3] = e12;
  e[4] = e13;
  e[5] = e23;
}

DielectricPermittivity DielectricPermittivity::fromTTensor(const qlc3d::TTensor &t, double S0, double deleps, double eper) {
  double e11 = (((2.0 / 3.0 / S0) * (-t.t1() / rt6 + t.t2() / rt2) + (1.0 / 3.0)) * deleps + eper); //~nx*nx
  double e22 = (((2.0 / 3.0 / S0) * (-t.t1() / rt6 - t.t2() / rt2) + (1.0 / 3.0)) * deleps + eper); //~ny*ny
  double e33 = (((2.0 / 3.0 / S0) * (2.0 * t.t1() / rt6)        + (1.0 / 3.0)) * deleps + eper); //~nz*nz
  double e12 = (2.0 / 3.0 / S0) * (t.t3() / rt2) * deleps;           //~nx*ny
  double e13 = (2.0 / 3.0 / S0) * (t.t5() / rt2) * deleps;           //~nx*nz
  double e23 = (2.0 / 3.0 / S0) * (t.t4() / rt2) * deleps;       //~ny*nz
  return {e11, e22, e33, e12, e13, e23};
}

void DielectricPermittivity::set(double e11, double e22, double e33, double e12, double e13, double e23) {
  e[0] = e11;
  e[1] = e22;
  e[2] = e33;
  e[3] = e12;
  e[4] = e13;
  e[5] = e23;
}

// </editor-fold>