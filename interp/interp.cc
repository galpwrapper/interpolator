#include<cmath>
#include<iostream>
#include"interp.h"

using std::vector;
using std::cout;
using std::endl;

#define RANGE(val) val > 700 ? 700 : (val < -700 ? -700 : val)

interp::interp(): x_tabv(lnx_tab), y_tabv(lny_tab) {}
interp::interp(const double *xtab_, const double *ytab_, unsigned n):
  x_tab(xtab_), y_tab(ytab_), dim(n),
  x_tabv(lnx_tab), y_tabv(lny_tab) {
  create_lntab();
}
interp::interp(const vector <double> &xtab_, const vector <double> &ytab_):
  x_tab(0), y_tab(0), dim(xtab_.size()),
  x_tabv(xtab_), y_tabv(ytab_) {
  create_lntab();
}
interp::interp(const spectrum &spec_):
  x_tab(0), y_tab(0), dim(spec_.E.size()),
  x_tabv(spec_.E), y_tabv(spec_.F) {
  create_lntab();
}
interp::interp(const double *xtab_, const double *ytab_, unsigned n, bool ln_or_not):
  x_tab(xtab_), y_tab(ytab_), dim(n),
  x_tabv(lnx_tab), y_tabv(lny_tab) {
  if (ln_or_not) create_lntab();
}
interp::interp(const vector <double> &xtab_, const vector <double> &ytab_, bool ln_or_not):
  x_tab(0), y_tab(0), dim(xtab_.size()),
  x_tabv(xtab_), y_tabv(ytab_) {
  if (ln_or_not) create_lntab();
}
interp::interp(const spectrum &spec_, bool ln_or_not):
  x_tab(0), y_tab(0), dim(spec_.E.size()),
  x_tabv(spec_.E), y_tabv(spec_.F) {
  if (ln_or_not) create_lntab();
}

int interp::create_lntab() {
  if (x_tab) return create_lntaba();

  else return create_lntabv();
}

int interp::create_lntaba() {
  for (unsigned i = 0; i < dim; i++) {
    lnx_tab.push_back(log(fmax(x_tab[i], 1e-300)));
    lny_tab.push_back(log(fmax(y_tab[i], 1e-300)));
  }

  return 0;
}

int interp::create_lntabv() {
  for (unsigned i = 0; i < x_tabv.size(); i++) {
    lnx_tab.push_back(log(fmax(x_tabv[i], 1e-300)));
    lny_tab.push_back(log(fmax(y_tabv[i], 1e-300)));
  }

  return 0;
}

int interp::get_index(const double x) const {
  if (x_tab)
    return get_indexa(x);

  else
    return get_indexv(x);
}

int interp::get_indexa(const double x) const {
  extern int get_index(const double *, unsigned, double);
  return get_index(x_tab, dim, x);
}

int interp::get_indexv(const double x) const {
  extern int get_index(const vector <double> &vec, double x);
  return get_index(x_tabv, x);
}

/*********************************************************************
  interpolating with Lagrangian interpolating method:
  y(x)=\sum_{j=0}^ny_jl_j(x)  and
  l_j(x)=\Prod_{i=0,i\neq j}^n (x-x_i)/(x_j-x_i)
*********************************************************************/
double interp:: laask(const double x, const int n) const {
  if (x_tab)
    return laaska(x, n);

  else
    return laaskv(x, n);
}

double interp:: laaska(const double x, const int n) const {
  double l[dim];
  double result = 0;
  int ind = get_indexa(x);

  for (int j = fmax(0, ind - n); j < fmin(dim, ind + n); j++) {
    l[j] = 1;

    for (int i = fmax(0, ind - n); i < fmin(dim, ind + n); i++)
      if (i != j)
        l[j] *= (x - x_tab[i]) / (x_tab[j] - x_tab[i]);

    result += l[j] * y_tab[j];
  }

  return result;
}

double interp:: laaskv(const double x, const int n) const {
  double l[dim];
  double result = 0;
  int ind = get_indexv(x);

  for (int j = fmax(0, ind - n); j < fmin(dim, ind + n); j++) {
    l[j] = 1;

    for (int i = fmax(0, ind - n); i < fmin(dim, ind + n); i++)
      if (i != j) l[j] *= (x - x_tabv[i]) / (x_tabv[j] - x_tabv[i]);

    result += l[j] * y_tabv[j];
  }

  return result;
}

double interp::linask(const double x) const {
  if (x_tab)
    return linaska(x);

  else
    return linaskv(x);
}

double interp::linaska(const double x) const {
  int ind = get_indexa(x);

  if (ind < 0) ind = 0;

  else if (unsigned(ind) > dim - 2) ind = dim - 2;

  double result = (x - x_tab[ind + 1]) * (y_tab[ind] - y_tab[ind + 1]) / (x_tab[ind] - x_tab[ind + 1]) + y_tab[ind + 1];
  return result;
}

double interp::linaskv(const double x) const {
  int ind = get_indexv(x);

  if (ind < 0) ind = 0;

  else if (unsigned(ind) > dim - 2) ind = dim - 2;

  double result = (x - x_tabv[ind + 1]) * (y_tabv[ind] - y_tabv[ind + 1]) / (x_tabv[ind] - x_tabv[ind + 1]) + y_tabv[ind + 1];
  return result;
}

double interp::lnask(double x) const {
  if (x <= 0) cout << "Error::interp::lnask: x out of normal range" << endl;

  int ind = get_index(x);

  if (ind < 0) ind = 0;

  else if (unsigned(ind) > dim - 2) ind = dim - 2;

  x = log(fmax(x, 1e-300));
  const double lnask = (x - lnx_tab[ind + 1]) * (lny_tab[ind] - lny_tab[ind + 1]) / (lnx_tab[ind] - lnx_tab[ind + 1]) + lny_tab[ind + 1];

  return exp(RANGE(lnask));
}
