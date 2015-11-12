#ifndef _INTERP_H
#define _INTERP_H
/*********************************************************************
interpolating with linear method or lagrangian method, or interpolating
in the logscale
*********************************************************************/
#include"spectrum.h"
#include"vec_utils.h"

class interp {
private:
  const double *x_tab,
        *y_tab;
  unsigned dim;
  std::vector <double> lnx_tab, lny_tab;
  const std::vector <double> &x_tabv, &y_tabv;
  int create_lntab();
  int create_lntaba();
  int create_lntabv();
  int get_indexa(const double x) const;
  int get_indexv(const double x) const;
  double laaska(const double x, const int n) const;
  double laaskv(const double x, const int n) const;
  double linaska(const double x) const;
  double linaskv(const double x) const;

public:
  interp();
  interp(const double *xtab_, const double *ytab_, unsigned n);
  interp(const std::vector <double> &xtab_, const std::vector <double> &ytab_);
  interp(const spectrum &spec_);
  interp(const double *xtab_, const double *ytab_, unsigned n, bool ln_or_not);
  interp(const std::vector <double> &xtab_, const std::vector <double> &ytab_, bool ln_or_not);
  interp(const spectrum &spec_, bool ln_or_not);

  int get_index(const double x) const;
  double laask(const double x, const int n) const;
  double linask(const double x) const;
  double lnask(double x) const;
};
#endif // for #ifndef _INTERP_H
