#ifndef MyPoint
#define MyPoint

#include <iostream>
using namespace std;

class myPoint {
public:
  int px;
  int py;
public:
  myPoint(int x = 0, int y = 0) : px(x), py(y) {}
  friend ostream& operator<<(ostream& os, const myPoint& pt);
};


class myRealPoint {
public:
  double px;
  double py;
public:
  myRealPoint(double x = 0, double y = 0) : px(x), py(y) {}
  friend ostream& operator<<(ostream& os, const myRealPoint& pt);
};


#endif //MyPoint
