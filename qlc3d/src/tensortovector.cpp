#include <math.h>
#include <float.h>
#include <algorithm>

/**
 * atoq changes the basis of the Q-tensor
 * from traceless to S/2(3*n*n-I)
 */
void atoq(double *a, double *q, int np) {
    int i;
    double rt2 = sqrtf(2.0), rt6 = sqrtf(6.0);
    for (i = 0; i < np; i++) {
        q[i] =      -a[i] / rt6 + a[i + np] / rt2;
        q[i + np] =   -a[i] / rt6 - a[i + np] / rt2;
        q[i + 2 * np] =  a[i + 2 * np] / rt2;
        q[i + 3 * np] =  a[i + 3 * np] / rt2;
        q[i + 4 * np] =  a[i + 4 * np] / rt2;
    }
} // end void *atoq

/*!
 * eig calculates eigenvalues and eigenvectors
 */
void eig(int mv, int n, double *a, double *r) {

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
    for (i = 0; i < n; i++) for (j = 0; j <= i; j++) {
            if (j != i) anorm += a[iq] * a[iq];
            ++iq;
        }
    if (anorm > 0.) {
        anorm = sqrt(2.*anorm);
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
                            sinx = y / sqrt(2.*(1. + (sqrt(1. - y * y))));
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
                            x = 2.*a[lm] * sincs;
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

double *tensortovector(double *a, int npLC) { // a = Q-tensor in traceless base
    int i;
    // Allocate for return vector and working array
    double *vector =  new double[5 * npLC];
    double *q = new double[5 * npLC];
    atoq(a, q, npLC); // change basis
    const int n = 3;
    double A[6]; // upper diagonal of tensor
    double R[9]; // matrix that will contain 3 x eigenvectors
    // loop over each node and convert Q-tensor to director and order parameter
    for (i = 0; i < npLC; i++) {
        double q1, q2, q3, q4, q5;
        q1 = q[i];
        q2 = q[i + npLC];
        q3 = q[i + 2 * npLC];
        q4 = q[i + 3 * npLC];
        q5 = q[i + 4 * npLC];
        A[0] = q1;
        A[1] = q3;
        A[2] = q2;
        A[3] = q5;
        A[4] = q4;
        A[5] = -q1 - q2;
        eig(1, n, A, R);    // find eigenvalue & eigenvectors
        A[1] = A[2];
        A[2] = A[5];
        int iMin, iMax;
        iMin = 0; iMax = 0;
        if ((A[0] <= A[1]) && (A[0] <= A[2])) iMin = 0;
        if ((A[1] <= A[0]) && (A[1] <= A[2])) iMin = 1;
        if ((A[2] <= A[0]) && (A[2] <= A[1])) iMin = 2;
        if ((A[0] >= A[1]) && (A[0] >= A[2])) iMax = 0;
        if ((A[1] >= A[0]) && (A[1] >= A[2])) iMax = 1;
        if ((A[2] >= A[0]) && (A[2] >= A[1])) iMax = 2;
        vector[i]   = R[3 * iMax + 0];
        vector[i + npLC]  = R[3 * iMax + 1];
        vector[i + 2 * npLC] = R[3 * iMax + 2]; // save director in o/p vector
        // how to get eigenvalues?
        vector[i + 3 * npLC] =  2.0 / 3.0 * (A[iMax] - A[iMin]); // save S1
        vector[i + 4 * npLC] = -2.0 / 3.0 * (A[iMax] + 2 * A[iMin]); // save S2
    }// end for i
    delete [] q;
    return (vector);
} // end tensortovector
