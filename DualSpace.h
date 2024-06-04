#ifndef DualSpace
#define DualSpace

#include "DualPixel.h"

typedef struct {
  double xMin;//int sizeDomainX;
  double xMax;//int sizeDomainY;
  vector<dualPixel> vecDualPixel; //list of dual pixel (and its preimage)
  vector<myRealPoint> preImage; //region of param space
  bool isClosed;
  bool isEmpty;
} dualSpace;

dualSpace initDualSpace(myPoint pmin, myPoint pmax);
dualSpace initDualSpace(double xmin = -5, double xmax = 5);
dualSpace initDualSpace(myPoint pix, double noise, double xmin = -5, double xmax = 5);
dualSpace getDualSpace(const vector<myPoint>& vecPts, const vector<double>& noise, double xmin = -5, double xmax = 5);

dualSpace addPixel(const dualSpace& sp, myPoint pix, double noise);
dualSpace addDualPixel(const dualSpace& sp, dualPixel dp);

myRealPoint getBarryCenter(const vector<myRealPoint>& polygon);
myRealPoint getBarryCenter(const dualSpace& sp);

pair<myPoint, myPoint> getDomaine(const vector<myPoint>& vecPts);

#endif //DualSpace
