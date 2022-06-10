#include "SplineAkima.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

void SplineAkima::extrapolate(vector<double*>& x, vector<double*>& y) {
  double dx = *x[2] - *x[0];
  *x[3] = *x[1] + dx;
  *x[4] = *x[2] + dx;

  double s0 = slope(x, y, 0),
         ds = slope(x, y, 1) - s0;

  for (int i = 2; i < 4; i++) {
    double si = s0 + ds * i;

    *y[i + 1] = *y[i] + si * (*x[i + 1] - *x[i]);
  }
}

void SplineAkima::extrapolate() {
  vector<double*> xlow, ylow, xup, yup;
  for (int i = 0; i < 5; i++) {
    xlow.push_back(&(xvec[4 - i]));
    ylow.push_back(&(yvec[4 - i]));

    xup.push_back(&(xvec[i + xvec.size() - 5]));
    yup.push_back(&(yvec[i + xvec.size() - 5]));
  }

  extrapolate(xlow, ylow);
  extrapolate(xup, yup);
}

void SplineAkima::count_slope() {
  vector<double> s(xvec.size() - 1, 0);
  for (int i = 0; i < s.size(); i++) s[i] = slope(i);

  dydx.resize(xvec.size(), 0);
  for (int i = 0; i < xvec.size() - 4; i++) {
    double d32 = fabs(s[i + 3] - s[i + 2]);
    double d10 = fabs(s[i + 1] - s[i]);

    if (d10 == d32 && d32 == 0) dydx[i + 2] = (s[i + 1] + s[i + 2]) / 2;
    else dydx[i + 2] = (d32 * s[i + 1] + d10 * s[i + 2]) / (d32 + d10);
  }
}

SplineAkima::SplineAkima(const vector<double>& xvec_, const vector<double>& yvec_) {
  assert(xvec_.size() >= 5 && yvec_.size() >= 5 && "SplineAkima::There should be at least 5 points for the interpolation");
  assert(xvec_.size() == yvec_.size() && "SplineAkima::The interpolated x, y should have the same size");

  int n = xvec_.size();
  xvec.resize(n + 4, 0);
  yvec.resize(n + 4, 0);

  for (int i = 0; i < n; i++) {
    xvec[i + 2] = xvec_[i];
    yvec[i + 2] = yvec_[i];
  }

  extrapolate();
  count_slope();
}

void SplineAkima::show() const {
  cout << "# SplineAkima status:" << endl
       << "#     x          y         dydx" << endl;
  cout << setprecision(4);
  const int nw = 10;
  for (int i = 0; i < xvec.size(); i++) {
    cout << setw(nw) << xvec[i] << " " << setw(nw) << yvec[i] << " " << setw(nw) << dydx[i] << endl;
  }
}

int get_index(const double* xtab, int n, const double x)
{
  return (upper_bound(xtab, xtab + n, x) - xtab);
}

inline double h00(double t) { return (1 + 2 * t) * (1 - t) * (1 - t); }
inline double h10(double t) { return t * (1 - t) * (1 - t); }
inline double h01(double t) { return t * t * (3 - 2 * t); }
inline double h11(double t) { return t * t * (t - 1); }

double SplineAkima::operator()(double x) const {
  if (x < xvec[2] || xvec[xvec.size() - 3] < x) return 0;

  int upind = get_index(&(xvec[0]), xvec.size(), x);
  if (upind < 3) upind = 3;
  if (upind > xvec.size() - 3) upind = xvec.size() - 3;
  int ind = upind - 1;

  double delta_x = xvec[upind] - xvec[ind];

  double t = (x - xvec[ind]) / (xvec[upind] - xvec[ind]);
  return h00(t) * yvec[ind] + h10(t) * delta_x * dydx[ind]
   + h01(t) * yvec[upind] + h11(t) * delta_x * dydx[upind];
}