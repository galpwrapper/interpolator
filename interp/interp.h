#ifndef _INTERP_H
#define _INTERP_H
/*********************************************************************
interpolating with linear method or lagrangian method, or interpolating
in the logscale
*********************************************************************/
#include"spectrum.h"
#include"vec_utils.h"

class interp
{
private:
  const double* x_tab,
        *y_tab;
  unsigned dim;
  std::vector <double> lnx_tab, lny_tab;
  double linask(const double* xtab, const double* ytab, const double x) const;

public:
  interp();
  interp(const double* xtab_, const double* ytab_, unsigned n);
  interp(const std::vector <double>& xtab_, const std::vector <double>& ytab_);
  interp(const spectrum& spec_);

  int get_index(const double* xtab, const double x) const;
  double laask(const double x, const int n) const;
  double linask(const double x) const;
  double lnask(const double x) const;
  double lnask_check(const double x);

  int create_lntab();
};
#endif // for #ifndef _INTERP_H
