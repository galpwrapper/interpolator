#ifndef SPLINEAKIMA_H
#define SPLINEAKIMA_H

#include <vector>

class SplineAkima {
  private:

  std::vector<double> xvec, yvec, dydx;

  inline double slope(const std::vector<double*>& x, const std::vector<double*>& y, int i) const {
    return (*y[i + 1] - *y[i]) / (*x[i + 1] - *x[i]);
  }
  void extrapolate();
  void extrapolate(std::vector<double*>& x, std::vector<double*>& y);

  inline double slope(int i) const { return (yvec[i + 1] - yvec[i]) / (xvec[i + 1] - xvec[i]); }
  void count_slope();

  public:

  SplineAkima(const std::vector<double>& xvec_, const std::vector<double>& yvec_);

  void show() const;

  double operator()(double x) const;
};

#endif /* SPLINEAKIMA_H */
