#ifndef _TABLE2D_H
#define _TABLE2D_H
/*********************************************************************
2 dimension table which consist of x axis, y axis and the values in
each corresponding point (x, y). The value in each point is determined
by a function you defined before.
 *********************************************************************/
#include<boost/serialization/map.hpp>
#include<vector>
#include<cstring>
#include<fstream>
#include<boost/archive/binary_oarchive.hpp>
#include<boost/archive/binary_iarchive.hpp>
#include"gfunction.h"

class Table2D {
public:
  typedef std::pair <double, double> Pair;
  typedef std::map <double, double> Line;
  typedef std::map <double, double>::iterator LineIter;
  typedef std::map <double, double>::const_iterator LineConsIter;
  typedef std::map <double, Line> Table;
  typedef std::map <double, Line>::iterator TabIter;
  typedef std::map <double, Line>::const_iterator TabConsIter;

  Line xaxis, yaxis;
  Table value;

  Table2D() {};
  Table2D(const std::vector <double> &x_, const std::vector <double> &y_, const std::vector <std::vector <double> > &tab_);
  Table2D(gfunction *func_);

  //int save(const std::string &filename) const;
  //int read(const std::string &filename) const;

  int setfunc(gfunction *func_);
  int insline(double x_);
  int inscolm(double y_);
  int insval(double x_, double y_, double val_);
  int trans();
  int list() const;
  int clear();
private:
  friend class boost::serialization::access;
  template <class Archive>
    void serialize(Archive &ar, const unsigned int version = 0) {
      ar & xaxis;
      ar & yaxis;
      ar & value;
    }

  LineConsIter getlineiter(double axival, const Line &line) const;
  TabConsIter gettabiter(double xval) const;
  void notfoundwarn() const;
  const gfunction *func;

  double dx2(TabConsIter rowmiter, LineConsIter yiter) const;
  double dy2_core(TabConsIter rowmiter, LineConsIter yiter) const;
  double dy2(TabConsIter rowmiter, LineConsIter yiter) const;
  double dx2(TabConsIter rowiter) const;
  double dy2(LineConsIter yiter) const;
};
#endif // for #ifndef _TABLE2D_H
