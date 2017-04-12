#include<cmath>
#include<iostream>
#include<algorithm>
#include"interp.h"

using std::lower_bound;
using std::upper_bound;
using std::vector;
using std::cout;
using std::endl;

#define RANGE(val) val > 700 ? 700 : (val < -700 ? -700 : val)

interp::interp() {}
interp::interp(const double* xtab_, const double* ytab_, unsigned n):
  x_tab(xtab_), y_tab(ytab_), dim(n) {}
interp::interp(const vector <double>& xtab_, const vector <double>& ytab_):
  x_tab(&(xtab_[0])), y_tab(&(ytab_[0])), dim(xtab_.size())
{
  if (xtab_.size() != ytab_.size())
    cout << "interp::interp::Warning::using vectors with different size" << endl;
}
interp::interp(const spectrum& spec_):
  x_tab(&(spec_.E[0])), y_tab(&(spec_.F[0])), dim(spec_.E.size()) {}

int interp::create_lntab()
{
  lnx_tab.resize(dim);
  lny_tab.resize(dim);

  for (unsigned i = 0; i < dim; i++) {
    lnx_tab[i] = log(fmax(x_tab[i], 1e-300));
    lny_tab[i] = log(fmax(y_tab[i], 1e-300));
  }

  return 0;
}

int interp::get_index(const double* xtab, const double x) const
{
  return (upper_bound(xtab, xtab + dim, x) - xtab);
}

/*********************************************************************
  interpolating with Lagrangian interpolating method:
  y(x)=\sum_{j=0}^ny_jl_j(x)  and
  l_j(x)=\Prod_{i=0,i\neq j}^n (x-x_i)/(x_j-x_i)
*********************************************************************/
double interp:: laask(const double x, const int n) const
{
  double l[dim];
  double result = 0;
  int ind = get_index(x_tab, x);

  for (int j = fmax(0, ind - n); j < fmin(dim, ind + n); j++) {
    l[j] = 1;

    for (int i = fmax(0, ind - n); i < fmin(dim, ind + n); i++)
      if (i != j)
        l[j] *= (x - x_tab[i]) / (x_tab[j] - x_tab[i]);

    result += l[j] * y_tab[j];
  }

  return result;
}

double interp::linask(const double* xtab, const double* ytab, const double x) const
{
  int upind = get_index(xtab, x);
  int ind = fmin(fmax(upind - 1, 0), dim - 2);

  return (x - xtab[ind + 1]) * (ytab[ind] - ytab[ind + 1]) / (xtab[ind] - xtab[ind + 1])
         + ytab[ind + 1];
}

double interp::linask(const double x) const
{
  return linask(x_tab, y_tab, x);
}
double interp::lnask(const double x) const
{
  if (x <= 0) {
    cout << "Error::interp::lnask: x out of normal range" << endl;
    exit(0);
  }
  if (lnx_tab.size() == 0) {
    cout << "Error::interp::lnask:Sorry, you should create_lntab before lnask, or use lnask_check instead" << endl;
    exit(0);
  }

  double ln_result = linask(&(lnx_tab[0]), &(lny_tab[0]), log(fmax(x, 1e-300)));
  return exp(RANGE(ln_result));
}

double interp::lnask_check(const double x)
{
  if (lnx_tab.size() == 0) create_lntab();
  return lnask(x);
}
