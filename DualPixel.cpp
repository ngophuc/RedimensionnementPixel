#include <iostream>

#include "DualPixel.h"

//pixel (x,y) -> p1(x-1/2, y+1/2), p2(x+1/2, y+1/2), p3(x-1/2, y-1/2), p4(x+1/2, y-1/2)
//pixel (x,y,n) -> p1(x-1/2-n, y+1/2+n), p2(x+1/2+n, y+1/2+n), p3(x-1/2-n, y-1/2-n), p4(x+1/2+n, y-1/2-n)
//pixel (x,y,n) -> p1(x-n/2, y+n/2), p2(x+n/2, y+n/2), p3(x-n/2, y-n/2), p4(x+n/2, y-n/2)
dualPixel getLinePixel(int x, int y, double noise) {
  dualPixel dPix;
  dPix.px=x;
  dPix.py=y;
  dPix.taille=noise;
  dPix.line1 = {-(x-noise/2.0), y+noise/2.0};//y=ax+b
  dPix.line2 = {-(x+noise/2.0), y+noise/2.0};//y=ax+b
  dPix.line3 = {-(x-noise/2.0), y-noise/2.0};//y=ax+b
  dPix.line4 = {-(x+noise/2.0), y-noise/2.0};//y=ax+b
  return dPix;
}

pair<myRealPoint, myRealPoint> getSegment(vector<double> line, double xmin, double xmax) {
  //line : y=ax+b
  double ymin = line.at(0)*xmin+line.at(1);
  double ymax = line.at(0)*xmax+line.at(1);
  return make_pair(myRealPoint(xmin, ymin), myRealPoint(xmax, ymax));
}

dualPixel getSegmentPixel(int x, int y, double noise, double xmin, double xmax){
  dualPixel dPix = getLinePixel(x, y, noise);
  vector<vector<double> > vecLine = {dPix.line1, dPix.line2, dPix.line3, dPix.line4};
  myRealPoint ip1 = intersectionTwoLines(vecLine.at(0), vecLine.at(1));//intersection between S1 & S2
  myRealPoint ip2 = intersectionTwoLines(vecLine.at(2), vecLine.at(3));//intersection between S3 & S4
  pair<myRealPoint, myRealPoint> s1 = getSegment(vecLine.at(0), xmax, ip1.px);
  pair<myRealPoint, myRealPoint> s2 = getSegment(vecLine.at(1), ip1.px, xmin);
  pair<myRealPoint, myRealPoint> s3 = getSegment(vecLine.at(2), xmin, ip2.px);
  pair<myRealPoint, myRealPoint> s4 = getSegment(vecLine.at(3), ip2.px, xmax);
  dPix.seg1 = make_pair(s1, -1);//h-
  dPix.seg2 = make_pair(s2, -1);//h-
  dPix.seg3 = make_pair(s3, 1);//h+
  dPix.seg4 = make_pair(s4, 1);//h+
  return dPix;
}

dualPixel getDualPixel(int x, int y, double noise, double xmin, double xmax){
  dualPixel dPix = getSegmentPixel(x, y, noise, xmin, xmax);
  //vector<myRealPoint> poly = {dPix.seg1.first.first, dPix.seg2.first.first, dPix.seg2.first.second,
  //  dPix.seg3.first.first, dPix.seg4.first.first, dPix.seg4.first.second, dPix.seg1.first.first};//ACW
  vector<myRealPoint> poly = {dPix.seg1.first.first, dPix.seg4.first.second, dPix.seg4.first.first,
    dPix.seg3.first.first, dPix.seg2.first.second,dPix.seg2.first.first, dPix.seg1.first.first};//CW
  dPix.polygon = poly;
  return dPix;
}
/*
DualPixel getDualPixel(int x, int y, double noise){
  vector<vector<double> > getSegmentPixel;
  
  vector<vector<double> > vecSeg = getSegmentPixel(x, y, noise);
  myRealPoint p1 =intersectionTwoLines(vecSeg[0], vecSeg[1]);//intersection between S1 & S2
  myRealPoint p2 =intersectionTwoLines(vecSeg[2], vecSeg[3]);//intersection between S3 & S4
  
  return getSegmentPixel;
}
*/
/*
vector<vector<double> > getSegmentPixel(int x, int y, double noise){
  vector<double> seg1 = {-(x-0.5-noise), y+0.5+noise};//y=ax+b
  vector<double> seg2 = {-(x+0.5+noise), y+0.5+noise};//y=ax+b
  vector<double> seg3 = {-(x-0.5-noise), y-0.5-noise};//y=ax+b
  vector<double> seg4 = {-(x+0.5+noise), y-0.5-noise};//y=ax+b
  vector<vector<double> > res = {seg1, seg2, seg3, seg4};
  return res;
}
*/

//y1=ax1+c
//y2=bx2+d
myRealPoint intersectionTwoLines(vector<double> line1, vector<double> line2) {
  double a=line1[0];
  double c=line1[1];
  double b=line2[0];
  double d=line2[1];
  double x=(d-c)/(a-b);//if a-b!=0
  double y=a*x+c;//if a-b!=0
  return myRealPoint(x,y);
}

vector<myRealPoint> polygonIntersection(vector<myRealPoint> poly1, vector<myRealPoint> poly2) {
  vector<myRealPoint> res;
  //typedef boost::geometry::model::polygon<boost::geometry::model::d2::point_xy<double> > polygon;
  typedef double coordinate_type;
  typedef boost::geometry::model::d2::point_xy<coordinate_type> point;
  typedef boost::geometry::model::polygon<point> polygon;
  
  std::vector<point> poly1_points;
  for(auto p : poly1)
    poly1_points.push_back(point(p.px,p.py));
  polygon poly1_boost;
  boost::geometry::assign_points(poly1_boost, poly1_points);
  
  std::vector<point> poly2_points;
  for(auto p : poly2)
    poly2_points.push_back(point(p.px,p.py));
  polygon poly2_boost;
  boost::geometry::assign_points(poly2_boost, poly2_points);
  
  std::deque<polygon> output_poly;
  boost::geometry::intersection(poly1_boost, poly2_boost, output_poly);
  
  std::vector<point> output_point;
  boost::geometry::intersection(poly1_boost, poly2_boost,output_point);
  cout<<"output_point.size()="<<output_point.size()<<endl;
  
  //assert(output_boost.size()==1);//region is empty
  //std::cout << "poly1_boost: " << boost::geometry::dsv(poly1_boost)<<endl;
  //std::cout << "poly2_boost: " << boost::geometry::dsv(poly2_boost)<<endl;
  //cout<<"output_boost: "<< boost::geometry::wkt(output_boost.front())<<endl;
  //BOOST_FOREACH(polygon const& p, output_boost)
  //    std::cout << boost::geometry::wkt(p) << std::endl;
  if(output_poly.size()!=0) {
    polygon output = output_poly.front();
    for( const auto& point : output.outer()) {
      double x = point.x();
      double y = point.y();
      //cout<<"x="<<x<<", y="<<y<<endl;
      res.push_back(myRealPoint(x,y));
    }
  }
  return res;
}

bool isClosedPolygon(vector<myRealPoint> poly, double xmin, double xmax) {
  if(poly.size()==0)
    return true;
  //As the polyline is closed, do not consider the last vertex
  int countXmin=0, countXmax=0;
  for(size_t it=0; it<poly.size()-1; it++) {
    if(fabs(poly.at(it).px-xmin)<1e-6)
      countXmin++;
    if(fabs(poly.at(it).px-xmax)<1e-6)
      countXmax++;
  }
  //cout<<"countXmin="<<countXmin<<" vs countXmax="<<countXmax<<endl;
  return countXmin<2 && countXmax<2;
}

vector<myRealPoint> drawDualPixel(dualPixel dp) {
  int px=dp.px;
  int py=dp.py;
  double size = dp.taille/2.0;
  vector<myRealPoint> poly;
  poly.push_back(myRealPoint(px-size, py+size));
  poly.push_back(myRealPoint(px+size, py+size));
  poly.push_back(myRealPoint(px+size, py-size));
  poly.push_back(myRealPoint(px-size, py-size));
  poly.push_back(myRealPoint(px-size, py+size));
  return poly;
}
