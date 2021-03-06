#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<cmath>
#include"Table2D.h"

#include"mydebug.h"

using std::function;
using std::abs;
using std::cout;
using std::endl;
using std::fstream;
using std::ios;
using std::string;
using std::vector;

Table2D::Table2D(const vector <double> &x_, const vector <double> &y_, const vector <vector <double> > &tab_):
  bound({ 0, 0, 0, 0}), empty_x(true), empty_y(true)
{
  for (unsigned i_x = 0; i_x < x_.size(); i_x++)
    for (unsigned i_y = 0; i_y < y_.size(); i_y++)
      insval(x_[i_x], y_[i_y], tab_[i_x][i_y]);
}

Table2D::Table2D(const function<double(double,double)>& func_): bound({ 0, 0, 0, 0}), empty_x(true), empty_y(true),  func(func_) {}

int Table2D::setfunc(const std::function<double(double,double)>& func_) {
  func = func_;
  return 0;
}

inline void Table2D::check_bound(bool &empty, bound_order down, bound_order up, double val) {
  printDebugMsg("Routine", ">>check_bound: empty = %s, down = %d, up = %d, bound size = %d",
                empty ? "yes" : "no", down, up, int(bound.size()));
  if (empty) {
    bound[down] = val, bound[up] = val;
    empty = false;
  } else {
    if (val < bound[down]) bound[down] = val;
    if (val > bound[up]) bound[up] = val;
  }
  printDebugMsg("Routine", "<<check_bound");
}

int Table2D::insline(double x_) {
  printDebugMsg("Routine", ">>insline: x = %f", x_);

  if (xaxis.find(x_) != xaxis.end()) return 0;
  check_bound(empty_x, x_down, x_up, x_);

  value[x_];
  xaxis[x_] = 0;

  TabIter rowmiter = value.find(x_),
          rowliter = rowmiter,
          rowhiter = rowmiter;
  -- rowliter, ++ rowhiter;
  LineIter xmiter = xaxis.find(x_),
           xliter = xmiter,
           xhiter = xmiter;
  -- xliter, ++ xhiter;

  for (LineConsIter yiter = yaxis.begin(); yiter != yaxis.end(); ++yiter)
    rowmiter -> second[yiter -> first] = func(x_, yiter -> first);

  xmiter -> second = dx2(rowmiter);

  if (xmiter != xaxis.begin())  xliter -> second = dx2(rowliter);

  if (xhiter != xaxis.end()) xhiter -> second = dx2(rowhiter);

  LineIter yiter = yaxis.begin(),
           viter = rowmiter -> second.begin(),
           yend = yaxis.end();

  for (; yiter != yend; ++yiter, ++viter)
    yiter -> second = fmax(yiter -> second, abs(dy2_core(rowmiter, viter)));

  printDebugMsg("Routine", "<<insline");
  return 0;
}

int Table2D::inscolm(double y_) {
  printDebugMsg("Routine", ">>inscolm: y = %f", y_);

  if (yaxis.find(y_) != yaxis.end()) return 0;
  check_bound(empty_y, y_down, y_up, y_);

  for (TabIter rowiter = value.begin(); rowiter != value.end(); ++ rowiter)
    rowiter -> second[y_] = func(rowiter -> first, y_);

  yaxis[y_] = 0;

  LineIter ymiter = yaxis.find(y_),
           yliter = ymiter,
           yhiter = ymiter;
  -- yliter, ++ yhiter;

  ymiter -> second = dy2(ymiter);

  if (ymiter != yaxis.begin()) yliter -> second = dy2(yliter);

  if (yhiter != yaxis.end()) yhiter -> second = dy2(yhiter);

  LineIter xiter = xaxis.begin(),
           xend = xaxis.end();
  TabIter riter = value.begin();

  for (; xiter != xend; ++ xiter, ++ riter) {
    xiter -> second = fmax(xiter -> second, abs(dx2(riter, ymiter)));
  }

  printDebugMsg("Routine", "<<inscolm")
  return 0;
}

int Table2D::insval(double x_, double y_, double val_) {
  printDebugMsg("Routine", ">>insval: x, y, val = %f, %f, %f", x_, y_, val_);
  printDebugMsg("State", "empty_x: %s, empty_y: %s", empty_x ? "true" : "false",
                empty_y ? "true" : "false");

  check_bound(empty_x, x_down, x_up, x_);
  check_bound(empty_y, y_down, y_up, y_);

  value[x_][y_] = val_;

  if (xaxis.find(x_) == xaxis.end()) xaxis[x_] = 0;

  if (yaxis.find(y_) == yaxis.end()) yaxis[y_] = 0;

  printDebugMsg("Routine", "<<insval: ");
  return 0;
}

int Table2D::trans() {
  Table valuetmp;
  TabConsIter irow;
  LineConsIter isub;

  for (irow = value.begin(); irow != value.end(); ++irow) {
    double xval = irow -> first;

    for (isub = irow -> second.begin(); isub != irow -> second.end(); ++isub)
      valuetmp[isub -> first][xval] = isub -> second;
  }

  value = valuetmp;
  xaxis.swap(yaxis);

  double tmp;
  tmp = bound[x_down];
  bound[x_down] = bound[y_down];
  bound[y_down] = tmp;

  tmp = bound[x_up];
  bound[x_up] = bound[y_up];
  bound[y_up] = tmp;
  return 0;
}

bool Table2D::inside(double x_, double y_) const {
  if ( x_ < bound[x_down] || x_ > bound[x_up] ||
       y_ < bound[y_down] || y_ > bound[y_up])
    return false;

  return true;
}

int Table2D::list() const {
  printf("                    ");

  for (LineConsIter yiter = yaxis.begin(); yiter != yaxis.end(); ++yiter)
    printf("%4.3e ", yiter -> first);

  printf("\n");
  printf("                    ");

  for (LineConsIter yiter = yaxis.begin(); yiter != yaxis.end(); ++yiter)
    printf("%4.3e ", yiter -> second);

  printf("\n");

  for (LineConsIter xiter = xaxis.begin(); xiter != xaxis.end(); ++xiter) {
    printf("%4.3e %4.3e ", xiter -> first, xiter -> second);

    TabConsIter irow = value.find(xiter -> first);

    for (LineConsIter riter = irow -> second.begin(); riter != irow -> second.end(); ++riter)
      printf("%4.3e ", riter -> second);

    printf("\n");
  }

  return 0;
}

int Table2D::clear() {
  xaxis.clear();
  yaxis.clear();
  value.clear();
  return 0;
}

Table2D::LineConsIter Table2D::getlineiter(double axival, const Line &line) const {
  LineConsIter iter = line.find(axival);

  if (iter == line.end()) notfoundwarn();

  return iter;
}

Table2D::TabConsIter Table2D::gettabiter(double xval) const {
  TabConsIter xiter = value.find(xval);

  if (xiter == value.end()) notfoundwarn();

  return xiter;
}

void Table2D::notfoundwarn() const {
  cout << "Error::Table2D::getlineiter/gettabiter/getiter::cannot find the corresponding pair" << endl;
  exit(1);
}

double Table2D::dx2(TabConsIter rowmiter, LineConsIter yiter) const {
  printDebugMsg("Routine", ">>dx2(sub): x = %f, y = %f", rowmiter -> first, yiter -> first);
  TabConsIter rowliter = rowmiter,
              rowhiter = rowmiter;
  -- rowliter, ++ rowhiter;

  if (rowmiter == value.begin() || rowhiter == value.end()) return 0;

  double x1 = rowliter -> first,
         x2 = rowmiter -> first,
         x3 = rowhiter -> first,
         y = yiter -> first;

  double val1 = getlineiter(y, rowliter -> second) -> second;
  double val2 = getlineiter(y, rowmiter -> second) -> second;
  double val3 = getlineiter(y, rowhiter -> second) -> second;

  //cout << "dx2 is " << val1 << " " << val2 << " " << val3 << " " << (val3 * (x2 - x1) + val2 * (x1 - x3) + val1 * (x3 - x2))*2
  //  /((x3 - x1) * (x2 - x1) * (x3 - x2)) << endl;

  double rslt = (val3 * (x2 - x1) + val2 * (x1 - x3) + val1 * (x3 - x2)) * 2
                / ((x3 - x1) * (x2 - x1) * (x3 - x2));
  printDebugMsg("Result", "%f", rslt);
  printDebugMsg("Routine", "<<dx2(sub)");
  return rslt;
}

double Table2D::dx2(TabConsIter xiter) const {
  printDebugMsg("Routine", ">>dx2(tot): x = %f", xiter -> first);
  double dx2_max = 0,
         dx2tmp;

  for (LineConsIter yiter =  yaxis.begin(); yiter != yaxis.end(); ++yiter) {
    dx2tmp = abs(dx2(xiter, yiter));
    dx2_max = dx2tmp > dx2_max ? dx2tmp : dx2_max;
  }

  printDebugMsg("Result", "%f", dx2_max);
  printDebugMsg("Routine", "<<dx2(tot)");
  return dx2_max;
}

//When you have the ymiter pointing to the exact point in value but not yaxis, you can runing the dy2_core directly
double Table2D::dy2_core(TabConsIter rowiter, LineConsIter ymiter) const {
  LineConsIter yliter = ymiter,
               yhiter = ymiter;
  -- yliter, ++ yhiter;

  if (ymiter == rowiter -> second.begin() || yhiter == rowiter -> second.end()) return 0;

  double y1 = yliter -> first,
         y2 = ymiter -> first,
         y3 = yhiter -> first,
         val1 = yliter -> second,
         val2 = ymiter -> second,
         val3 = yhiter -> second;

  return (val3 * (y2 - y1) + val2 * (y1 - y3) + val1 * (y3 - y2)) * 2
         / ((y3 - y1) * (y2 - y1) * (y3 - y2));
}

double Table2D::dy2(TabConsIter rowiter, LineConsIter ymiter) const {
  ymiter = getlineiter(ymiter -> first, rowiter -> second);
  return dy2_core(rowiter, ymiter);
}

double Table2D::dy2(LineConsIter yiter) const {
  TabConsIter rowiter =  value.begin();
  double dy2_max = 0,
         dy2tmp;

  for (; rowiter != value.end(); ++ rowiter) {
    dy2tmp = abs(dy2(rowiter, yiter));
    dy2_max = dy2tmp > dy2_max ? dy2tmp : dy2_max;
  }

  return dy2_max;
}

int Table2D::Delete() {
  value.clear();
  xaxis.clear();
  yaxis.clear();
  return 0;
}
