#include<iostream>
#include<cstring>
#include<cstdlib>
#include<cmath>
#include<set>
#include"Interp2D.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::log;

#define logwb(expr) log(fmax(expr, 1e-300))
double IntpFunction::operator()(double x, double y) const {
  return ln_or_not ? logwb((*func)(exp(x), exp(y))) : (*func)(x, y);
}
IntpFunction::IntpFunction(gfunction *func_, bool ln_or_not_): func(func_), ln_or_not(ln_or_not_) {}
IntpFunction::IntpFunction(): func(0), ln_or_not(false) {}
int IntpFunction::initial(gfunction *func_, bool ln_or_not_) {
  func = func_;
  ln_or_not = ln_or_not_;
  return 0;
}

int Interp2D::mapping() {
  Table2D::LineConsIter enditer = tab.xaxis.end();

  for (Table2D::LineConsIter iter = tab.xaxis.begin(); iter != enditer; ++ iter)
    lntab.xaxis[logwb(iter -> first)];

  enditer = tab.yaxis.end();

  for (Table2D::LineConsIter iter = tab.yaxis.begin(); iter != enditer; ++ iter)
    lntab.yaxis[logwb(iter -> first)];

  Table2D::TabConsIter riter = tab.value.begin(),
                       renditer = tab.value.end();

  for (; riter != renditer; ++ riter) {
    lntab.value[logwb(riter -> first)];
    Table2D::TabIter cursor = lntab.value.find(logwb(riter -> first));
    enditer = riter -> second.end();
    Table2D::LineConsIter viter = riter -> second.begin(),
                          yiter = lntab.yaxis.begin();

    for (; viter != enditer; ++ viter, ++ yiter)
      cursor -> second[yiter -> first] = logwb(viter -> second);
  }

  return 0;
}

int Interp2D::lnmapping() {
  Table2D::LineConsIter enditer = lntab.xaxis.end();

  for (Table2D::LineConsIter iter = lntab.xaxis.begin(); iter != enditer; ++ iter)
    tab.xaxis[exp(iter -> first)];

  enditer = lntab.yaxis.end();

  for (Table2D::LineConsIter iter = lntab.yaxis.begin(); iter != enditer; ++ iter)
    tab.yaxis[exp(iter -> first)];

  Table2D::TabConsIter riter = lntab.value.begin(),
                       renditer = lntab.value.end();

  for (; riter != renditer; ++ riter) {
    tab.value[exp(riter -> first)];
    Table2D::TabIter cursor = tab.value.find(exp(riter -> first));

    enditer = riter -> second.end();
    Table2D::LineConsIter viter = riter -> second.begin(),
                          yiter = tab.yaxis.begin();

    for (; viter != enditer; ++ viter, ++ yiter)
      cursor -> second[yiter -> first] = exp(viter -> second);
  }

  return 0;
}

#define tri_res(xnum, ynum, x, y)\
  ((*spt[0])[2]+(*spt[1])[2])/2+\
  (y-(*spt[0])[ynum])*( ((*spt[0])[2]+(*spt[1])[2])/2-(*spt[2])[2])/((*spt[0])[ynum]-(*spt[2])[ynum])+\
  (x-(*spt[2])[xnum])*((*spt[0])[2]-(*spt[1])[2])/((*spt[0])[xnum]-(*spt[1])[xnum])

/*********************************************************************
Interpolating in a triangle area from the rectangle area
*********************************************************************/
double Interp2D::tri_intp(double(*spt[3])[3], double x, double y) const {
  if ((*spt[0])[0] == (*spt[1])[0])
    return tri_res(1, 0, y, x);

  else
    return tri_res(0, 1, x, y);
}
/*********************************************************************
To Interpolate in a rectangle area, we first divide it into  four triangle and Interpolating in the corresponding triangle
*********************************************************************/
#define obtiter(itertype, itersm, iterlg, object, num)\
  itertype itersm = object.lower_bound(num);\
  if(itersm == object.begin()) ++ itersm;\
  if(itersm == object.end()) -- itersm;\
  itertype iterlg = itersm --

double Interp2D::linask(const Table2D &table, double x, double y) const {
  obtiter(Table2D::TabConsIter, rowiterup, rowiterdown, table.value, x);
  obtiter(Table2D::LineConsIter, tliter, triter, rowiterup -> second, y);
  obtiter(Table2D::LineConsIter, bliter, briter, rowiterdown -> second, y);

  double points[4][3] = {{rowiterdown -> first, bliter -> first, bliter -> second},
    {rowiterup -> first, tliter -> first, tliter -> second},
    {rowiterup -> first, triter -> first, triter -> second},
    {rowiterdown -> first, briter -> first, briter -> second}
  };

  double(*select[3])[3];
  double k = (points[br][1] - points[tl][1]) / (points[br][0] - points[tl][0]);
  select[0] = (y - points[tl][1] - (x - points[tl][0]) * k > 0) ?
              points + tr : points + bl;

  select[1] = (y - points[bl][1] + (x - points[bl][0]) * k > 0) ?
              points + br : points + tl;

  double midx = (rowiterup -> first + rowiterdown -> first) / 2.0,
         midy = (tliter -> first + triter -> first) / 2.0,
         midval = (points[0][2] + points[1][2] + points[2][2] + points[3][2]) / 4.0;
  double midpoint[1][3] = {{midx, midy, midval}};
  select[2] = midpoint;
  return tri_intp(select, x, y);
}

double Interp2D::linask(double x, double y) const {
  return linask(tab, x, y);
}

#define RANGE(val) val > 700 ? 700 : (val < -700 ? -700 : val)

double Interp2D::lnask(double x, double y) const {
  if (x <= 0 || y <= 0) cout << "Error::Interp2D::lnask: x or y out of range" << endl;

  x = log(x);
  y = log(y);
  return exp(RANGE(linask(lntab, x, y)));
}

#define INTERVAL_MIN 3e-16 //The min interval
#define MAXGNUM 3000
//when checking axis, the axis's size is at least 3
int Interp2D::checkaxis(const Table2D::Line &axis, double err) {
  chlist.clear();
  Table2D::LineConsIter hiter = axis.begin(),
                        liter = hiter ++,
                        miter = hiter ++;
  double midpoint1, midpoint2, interval;

  for (; hiter != axis.end(); liter = miter, miter = hiter, ++hiter) {
    midpoint1 = (liter -> first + miter -> first) / 2;
    midpoint2 = (miter -> first + hiter -> first) / 2;
    interval = (hiter -> first - liter -> first) / 2;

    if (interval < INTERVAL_MIN || axis.size() >= MAXGNUM) return 1;

    if ((miter -> second) * interval * interval > err) {
      if (chlist.empty() || midpoint1 != chlist[chlist.size() - 1])
        chlist.push_back(midpoint1);

      chlist.push_back(midpoint2);
    }
  }

  return 0;
}

//addline to the table until we should not add any more. return false if nothing was done. The size of xaxis should be at least 3
bool Interp2D::tranversx(Table2D &table, double err) {
  double turns;

  for (turns = 0; ; turns++) {
    checkaxis(table.xaxis, err);

    if (chlist.empty()) break;

    for (vector <double>::iterator iter =  chlist.begin(); iter != chlist.end(); ++iter)
      table.insline(*iter);
  }

  //turns == 0 means noline was added.
  return (turns == 0) ? false : true;
}

bool Interp2D::tranversy(Table2D &table, double err) {
  double turns;

  for (turns = 0; ; turns++) {
    checkaxis(table.yaxis, err);

    if (chlist.empty()) break;

    for (vector <double>::iterator iter =  chlist.begin(); iter != chlist.end(); ++iter)
      table.inscolm(*iter);
  }

  return (turns == 0) ? false : true;
}

int Interp2D::create_table(Table2D &table, double range[4], double err) {
  tab.clear();
  lntab.clear();
  table.setfunc(&function);
  double midx = (range[0] + range[2]) / 2,
         midy = (range[1] + range[3]) / 2;
  table.insline(range[0]);
  table.insline(range[2]);
  table.insline(midx);
  table.inscolm(range[1]);
  table.inscolm(range[3]);
  table.inscolm(midy);

  bool xchange = true,
       ychange = true;

  while (xchange || ychange) {
    xchange = tranversx(table, err);
    ychange = tranversy(table, err);
  }

  cout << "x size is " << table.xaxis.size() << " y size is " << table.yaxis.size() << endl;
  return 0;
}

int Interp2D::creating(gfunction *func_, double range[4], double err) {
  function.initial(func_, false);
  create_table(tab, range, err);
  mapping();
  return 0;
}

int Interp2D::lncreating(gfunction *func_, double range[4], double err) {
  function.initial(func_, true);
  double lnrange[4];

  for (int i = 0; i < 4; i++) {
    if (range[i] <= 0) cout << "Error::Interp2D::lncreate_tab:range wrong" << endl;

    lnrange[i] = log(range[i]);
  }

  create_table(lntab, lnrange, err);
  lnmapping();
  return 0;
}

Interp2D::Interp2D(gfunction *func_, double range[4], double err) {
  creating(func_, range, err);
}

Interp2D::Interp2D() {}

Interp2D::Interp2D(const Table2D &tab_): tab(tab_) {
  mapping();
}
int Interp2D::tabling(const Table2D &tab_) {
  tab = tab_;
  mapping();
  return 0;
}
int Interp2D::lntabling(const Table2D &lntab_) {
  lntab = lntab_;
  lnmapping();
  return 0;
}
