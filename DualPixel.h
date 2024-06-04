#ifndef DualPixel
#define DualPixel

#include <vector>
#include <deque>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/foreach.hpp>

#include "MyPoint.h"

using namespace std;

typedef struct {
  int px;
  int py;
  double taille;
  //P(x,y)-Line(a,b: y=ax+b)-HP(1(>=0), -1(<=0)
  vector<double> line1;//Line(a,b: y=ax+b)
  vector<double> line2;
  vector<double> line3;
  vector<double> line4;
  pair<pair<myRealPoint,myRealPoint>, int> seg1;//Seg(p1->p2) + half-space(HP(1(>=0), -1(<=0))
  pair<pair<myRealPoint,myRealPoint>, int> seg2;
  pair<pair<myRealPoint,myRealPoint>, int> seg3;
  pair<pair<myRealPoint,myRealPoint>, int> seg4;
  vector<myRealPoint> polygon;
} dualPixel;

dualPixel getLinePixel(int x, int y, double noise=1.0);
dualPixel getSegmentPixel(int x, int y, double noise=1.0, double xmin=-5, double xmax=5);
dualPixel getDualPixel(int x, int y, double noise=1.0, double xmin=-5, double xmax=5);
//dualPixel getDualPixel(int x, int y, double noise=0);
//vector<vector<double> > getSegmentPixel(int x, int y, double noise=0);

myRealPoint intersectionTwoLines(vector<double> line1, vector<double> line2);
vector<myRealPoint> polygonIntersection(vector<myRealPoint> poly1, vector<myRealPoint> poly2);
bool isClosedPolygon(vector<myRealPoint> poly, double xmin=-10, double xmax=10);

vector<myRealPoint> drawDualPixel(dualPixel dp);

#endif //DualPixel
