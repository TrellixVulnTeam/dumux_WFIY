// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief This is a program to test the polynomial spline interpolation.
 *
 * It just prints some function to stdout. You can look at the result
 * using the following commands:
 *
----------- snip -----------
./test_spline > spline.csv
gnuplot

gnuplot> plot "spline.csv" using 1:2 w l ti "Curve", \
              "spline.csv" using 1:3 w l ti "Derivative", \
              "spline.csv" using 1:4 w p ti "Monotonical"
----------- snap -----------


*/
#include "config.h"

#include <dumux/common/spline.hh>

#include <dune/common/fvector.hh>

template <class Spline>
void testCommon(const Spline &sp,
                const double *x,
                const double *y)
{
    static double eps = 1e-10;

    int n = sp.numSamples();
    for (int i = 0; i < n; ++i) {
        // sure that we hit all sampling points
        double y0 = (i>0)?sp.eval(x[i]-eps):y[0];
        double y1 = sp.eval(x[i]);
        double y2 = (i<n-1)?sp.eval(x[i]+eps):y[n-1];
        assert(std::abs(y0 - y[i]) < 100*eps);
        assert(std::abs(y1 - y[i]) < eps);
        assert(std::abs(y2 - y[i]) < 100*eps);

        // make sure the derivative is continuous (assuming that the
        // second derivative is smaller than 1000)
        double d1 = sp.evalDerivative(x[i]);
        double d0 = (i>0)?sp.evalDerivative(x[i]-eps):d1;
        double d2 = (i<n-1)?sp.evalDerivative(x[i]+eps):d1;
        assert(std::abs(d1 - d0) < 1000*eps);
        assert(std::abs(d2 - d0) < 1000*eps);

        // make sure the derivative is consistent with the y values
        y0 = sp.eval(x[i] - ((i>0)?eps*1e2:0));
        y2 = sp.eval(x[i] + ((i<n-1)?eps*1e2:0));
        double dC = (y2 - y0)/(2*1e2*eps);
        if (i == 0 || i == n-1)
            dC *= 2;
        assert(std::abs(dC - d0) < 1e-5);
    }
}

template <class Spline>
void testFull(const Spline &sp,
              const double *x,
              const double *y,
              double m0,
              double m1)
{
    // test the common properties of splines
    testCommon(sp, x, y);

    static double eps = 1e-10;
    int n = sp.numSamples();

    // make sure the derivative at both end points is correct
    double d0 = sp.evalDerivative(x[0]);
    double d1 = sp.evalDerivative(x[n-1]);
    assert(std::abs(d0 - m0) < eps);
    assert(std::abs(d1 - m1) < eps);
};

template <class Spline>
void testNatural(const Spline &sp,
                 const double *x,
                 const double *y)
{
    // test the common properties of splines
    testCommon(sp, x, y);

    static double eps = 1e-10;
    int n = sp.numSamples();

    // make sure the second derivatives at both end points are 0
    double d0 = sp.evalDerivative(x[0]);
    double d1 = sp.evalDerivative(x[0] + eps);

    double d2 = sp.evalDerivative(x[n-1] - eps);
    double d3 = sp.evalDerivative(x[n-1]);
    assert(std::abs(d1 - d0)/eps < 1000*eps);
    assert(std::abs(d3 - d2)/eps < 1000*eps);
};

void testAll()
{
    double x[] = { 0, 5, 7.5, 8.75, 9.375 };
    double y[] = { 10, 0, 10, 0, 10 };
    double m0 = 10;
    double m1 = -10;
    double points[][2] =
        {
            {x[0], y[0]},
            {x[1], y[1]},
            {x[2], y[2]},
            {x[3], y[3]},
            {x[4], y[4]},
        };

    std::vector<double> xVec;
    std::vector<double> yVec;
    std::vector<double*> pointVec;
    for (int i = 0; i < 5; ++i) {
        xVec.push_back(x[i]);
        yVec.push_back(y[i]);
        pointVec.push_back(points[i]);
    }

    /////////
    // test fixed length spline, n = 2
    /////////

    // full spline
    { Dumux::Spline<double, 2> sp(x[0], x[1], y[0], y[1], m0, m1); sp.set(x[0],x[1],y[0],y[1],m0, m1); testFull(sp, x, y, m0, m1); };
    { Dumux::Spline<double, 2> sp(x, y, m0, m1); sp.set(x,y,m0, m1); testFull(sp, x, y, m0, m1);  };
    { Dumux::Spline<double, 2> sp(points, m0, m1); sp.set(points,m0, m1); testFull(sp, x, y, m0, m1); };


    /////////
    // test fixed length spline, n > 2
    /////////

    // full spline
    { Dumux::Spline<double, 5> sp(x, y, m0, m1); sp.set(x,y,m0, m1); testFull(sp, x, y, m0, m1);  };
    { Dumux::Spline<double, 5> sp(points, m0, m1); sp.set(points,m0, m1); testFull(sp, x, y, m0, m1); };

    // natural spline
    { Dumux::Spline<double, 5> sp(x, y); sp.set(x, y); testNatural(sp, x, y); };
    { Dumux::Spline<double, 5> sp(points); sp.set(points); testNatural(sp, x, y); };

    /////////
    // test variable length splines
    /////////

    // full spline
    { Dumux::Spline<double, -1> sp(5, x, y, m0, m1); sp.set(5,x,y,m0, m1); testFull(sp, x, y, m0, m1);  };
    { Dumux::Spline<double, -1> sp(5, points, m0, m1); sp.set(5,points,m0, m1); testFull(sp, x, y, m0, m1); };
    { Dumux::Spline<double, -1> sp(xVec, yVec, m0, m1); sp.set(xVec,yVec,m0, m1); testFull(sp, x, y, m0, m1);  };
    { Dumux::Spline<double, -1> sp(pointVec, m0, m1); sp.set(pointVec,m0, m1); testFull(sp, x, y, m0, m1); };

    // natural spline
    { Dumux::Spline<double, -1> sp(5, x, y); sp.set(5,x,y); testNatural(sp, x, y);  };
    { Dumux::Spline<double, -1> sp(5, points); sp.set(5,points); testNatural(sp, x, y); };
    { Dumux::Spline<double, -1> sp(xVec, yVec); sp.set(xVec,yVec); testNatural(sp, x, y); };
    { Dumux::Spline<double, -1> sp(pointVec); sp.set(pointVec); testNatural(sp, x, y); };
}

void plot()
{
    const int numSamples = 5;
    const int n = numSamples - 1;
    typedef Dune::FieldVector<double, numSamples> FV;

    double x_[] = { 0, 5, 7.5, 8.75, 9.375 };
    double y_[] = { 10, 0, 10, 0, 10 };
    double m1 = 10;
    double m2 = -10;
    FV &xs = *reinterpret_cast<FV*>(x_);
    FV &ys = *reinterpret_cast<FV*>(y_);

    Dumux::Spline<double, numSamples> spFull(xs, ys, m1, m2);
    Dumux::Spline<double, numSamples> spNatural(xs, ys);

    spFull.printCSV(-0.01*(x_[n] - x_[0]) + x_[0],
                    0.01*(x_[n] - x_[0]) + x_[n],
                    1000);
    std::cout << "\n";
    spNatural.printCSV(-0.01*(x_[n] - x_[0]) + x_[0],
                       0.01*(x_[n] - x_[0]) + x_[n],
                       1000);
    std::cerr << "Spline is monotonic: " << spFull.monotonic(x_[0], x_[n]) << "\n";
};

int main(int argc, char** argv)
{
    testAll();

    plot();
    return 0;
}
