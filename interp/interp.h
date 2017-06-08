#ifndef _INTERP_H
#define _INTERP_H
/*********************************************************************
interpolating with linear method or lagrangian method, or interpolating
in the logscale
*********************************************************************/
#include<cstring>

#include"spectrum.h"
#include"vec_utils.h"

class interp
{
private:
  const double* x_tab,
        *y_tab;
  unsigned dim;
  std::vector <double> lnx_tab, lny_tab, m, m_log;
  double linask(const double* xtab, const double* ytab, const double x) const;

  void set_m(const double* xtab, const double* ytab, std::vector<double>& m_);
  double spline_ask(const double* xtab, const double* ytab, const double* m_, const double x) const; // Support the cubic Hermite spline interpolating.
  bool with_spline, ln_created;

  void show_vector(const std::string& key, const double* vec) const;

public:
  interp();
  interp(const double* xtab_, const double* ytab_, unsigned n, bool with_spline_ = false);
  interp(const std::vector <double>& xtab_, const std::vector <double>& ytab_, bool with_spline_ = false);
  interp(const spectrum& spec_, bool with_spline_ = false);

  void set_spline();
  void unset_spline();

  int get_index(const double* xtab, const double x) const;
  double laask(const double x, const int n) const;
  double linask(const double x) const;
  double lnask(const double x) const;
  double lnask_check(const double x);

  int create_lntab();
  void show() const;
};
#endif // for #ifndef _INTERP_H
