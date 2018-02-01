#include<cmath>
#include<iostream>
#include<algorithm>
#include<cstring>
#include"interp.h"

using std::string;
using std::lower_bound;
using std::upper_bound;
using std::vector;
using std::cout;
using std::endl;

#define RANGE(val) val > 700 ? 700 : (val < -700 ? -700 : val)

interp::interp() {}
interp::interp(const double* xtab_, const double* ytab_, unsigned n, bool with_spline_):
  x_tab(xtab_), y_tab(ytab_), dim(n), with_spline(with_spline_), ln_created(false)
{
  if (with_spline) set_spline();
}
interp::interp(const vector <double>& xtab_, const vector <double>& ytab_, bool with_spline_):
  x_tab(&(xtab_[0])), y_tab(&(ytab_[0])), dim(xtab_.size()), with_spline(with_spline_), ln_created(false)
{
  if (with_spline) set_spline();

  if (xtab_.size() != ytab_.size())
    cout << "interp::interp::Warning::using vectors with different size" << endl;
}
interp::interp(const spectrum& spec_, bool with_spline_):
  x_tab(&(spec_.E[0])), y_tab(&(spec_.F[0])), dim(spec_.E.size()), with_spline(with_spline_), ln_created(false)
{
  if (with_spline) set_spline();
}


void interp::ini(const double* xtab_, const double* ytab_, unsigned n, bool with_spline_)
{
  x_tab = xtab_;
  y_tab = ytab_;
  dim = n;
  with_spline = with_spline_;
  ln_created = false;

  if (with_spline) set_spline();
}
void interp::ini(const vector<double>& xtab_, const vector<double>& ytab_, bool with_spline_)
{
  return ini(&(xtab_[0]), &(ytab_[0]), xtab_.size(), with_spline_);
}
void interp::ini(const spectrum& spec_, bool with_spline_)
{
  return ini(&(spec_.E[0]), &(spec_.F[0]), spec_.E.size(), with_spline_);
}

void interp::show_vector(const string& k, const double* vec) const
{
  cout << k << " is:" << endl;
  for (unsigned i = 0; i < dim; i++) cout << vec[i] << " "; cout << endl;
}

void interp::show() const
{
  show_vector("x", x_tab);
  show_vector("y", y_tab);

  if (with_spline) show_vector("m", &(m[0]));

  if (ln_created) {
    show_vector("ln(x)", &(lnx_tab[0]));
    show_vector("ln(y)", &(lny_tab[0]));
    if (with_spline) show_vector("m for ln", &(m_log[0]));
  }
}

void interp::set_m(const double* xtab, const double* ytab, vector<double>& m_)
{
  if (dim < 2) return;

  m_.resize(dim);
  m_[0] = (ytab[1] - ytab[0]) / (xtab[1] - xtab[0]);
  m_[dim - 1] = (ytab[dim - 1] - ytab[dim - 2]) / (xtab[1] - xtab[0]);

  for (unsigned i = 1; i < dim - 1; i++)
    m_[i] = 0.5 * ((ytab[i] - ytab[i - 1]) / (xtab[i] - xtab[i - 1]) + (ytab[i + 1] - ytab[i]) / (xtab[i + 1] - xtab[i]));
}

inline double h00(double t) { return (1 + 2 * t) * (1 - t) * (1 - t); }
inline double h10(double t) { return t * (1 - t) * (1 - t); }
inline double h01(double t) { return t * t * (3 - 2 * t); }
inline double h11(double t) { return t * t * (t - 1); }

double interp::spline_ask(const double* xtab, const double* ytab, const double* m_, const double x) const
{
  int upind = get_index(xtab, x);
  if (upind <= 0) upind = 1;
  if (upind >= dim) upind = dim - 1;
  int ind = upind - 1;
  double delta_x = xtab[upind] - xtab[ind];

  const double dh00[2] = { 0, 0 },
        dh10[2] = { 1, 0 },
        dh01[2] = { 0, 0 },
        dh11[2] = { 0, 1 };  // The value of dh/dt at point t = 0 and t = 1;
  if (x < xtab[0] || x > xtab[dim - 1]) {
    unsigned point_ind = x < xtab[0] ? 0 : dim - 1,
             ib = x < xtab[0] ? 0 : 1;

    double k = dh00[ib] * ytab[ind] / delta_x + dh10[ib] * m_[ind] + dh01[ib] * ytab[upind] / delta_x + dh11[ib] * m_[upind];
    return (x - xtab[point_ind]) * k + ytab[point_ind];
  }

  double t = (x - xtab[ind]) / (xtab[upind] - xtab[ind]);
  return h00(t) * ytab[ind] + h10(t) * delta_x * m_[ind] + h01(t) * ytab[upind] + h11(t) * delta_x * m_[upind];
}

void interp::set_spline()
{
  set_m(x_tab, y_tab, m);

  if (ln_created) set_m(&(lnx_tab[0]), &(lny_tab[0]), m_log);
  with_spline = true;
}
void interp::unset_spline() { with_spline = false; }

int interp::create_lntab()
{
  lnx_tab.resize(dim);
  lny_tab.resize(dim);

  for (unsigned i = 0; i < dim; i++) {
    lnx_tab[i] = log(fmax(x_tab[i], 1e-300));
    lny_tab[i] = log(fmax(y_tab[i], 1e-300));
  }

  if (with_spline) set_m(&(lnx_tab[0]), &(lny_tab[0]), m_log);
  ln_created = true;
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
  return with_spline ? spline_ask(x_tab, y_tab, &(m[0]), x) : linask(x_tab, y_tab, x);
}
double interp::lnask(const double x) const
{
  if (x <= 0) {
    cout << "Error::interp::lnask: x out of normal range" << endl;
    exit(0);
  }
  if (!ln_created) {
    cout << "Error::interp::lnask:Sorry, you should create_lntab before lnask, or use lnask_check instead" << endl;
    exit(0);
  }

  double ln_result;
  if (with_spline) ln_result = spline_ask(&(lnx_tab[0]), &(lny_tab[0]), &(m_log[0]), log(fmax(x, 1e-300)));
  else ln_result = linask(&(lnx_tab[0]), &(lny_tab[0]), log(fmax(x, 1e-300)));

  return exp(RANGE(ln_result));
}

double interp::lnask_check(const double x)
{
  if (!ln_created) create_lntab();

  return lnask(x);
}
