#include <iostream>

#include "DualSpace.h"

dualSpace initDualSpace(myPoint pmin, myPoint pmax) {
  dualSpace sp;
  assert(pmin.px<pmax.px);
  sp.xMin = pmin.px;
  sp.xMax = pmax.px;
  sp.isClosed = false;
  sp.isEmpty = false;
  return sp;
}

dualSpace initDualSpace(double xmin, double xmax) {
  dualSpace sp;
  assert(xmin<xmax);
  sp.xMin = xmin;
  sp.xMax = xmax;
  sp.isClosed = false;
  sp.isEmpty = false;
  return sp;
}

dualSpace initDualSpace(myPoint pix, double noise, double xmin, double xmax) {
  dualSpace sp = initDualSpace(xmin, xmax);
  dualPixel dp = getDualPixel(pix.px, pix.py, noise, sp.xMin, sp.xMax);
  sp.vecDualPixel.push_back(dp);
  vector<myRealPoint> poly = dp.polygon;
  sp.preImage = poly;
  sp.isEmpty=false;
  sp.isClosed=false;
  return sp;
}

dualSpace getDualSpace(const vector<myPoint>& vecPts, const vector<double>& noise, double xmin, double xmax) {
  dualSpace sp = initDualSpace(xmin, xmax);
  vector<dualPixel> vecPix;
  assert(vecPts.size()==noise.size());
  
  //list of pixel and its preimage
  for(size_t it=0; it<vecPts.size(); it++) {
    dualPixel p = getDualPixel(vecPts.at(it).px, vecPts.at(it).py, noise.at(it), xmin, xmax);
    vecPix.push_back(p);
  }
  sp.vecDualPixel = vecPix;
  
  //region of param space
  vector<myRealPoint> PreImage;
  myPoint p1 = vecPts.front();
  dualPixel dP1 = getDualPixel(p1.px, p1.py, noise.front(), xmin, xmax);
  vector<myRealPoint> poly1 = dP1.polygon;
  vector<myRealPoint> poly_res = poly1;
  //cout<<"isClosed? "<<isClosedPolygon(poly_res, xmin, xmax)<<endl;
  for(size_t it=1; it<vecPts.size(); it++) {
    myPoint p2 = vecPts.at(it);
    dualPixel dP2 = getDualPixel(p2.px, p2.py, noise.at(it), xmin, xmax);
    vector<myRealPoint> poly2 = dP2.polygon;
    poly_res = polygonIntersection(poly_res, poly2);
    //cout<<it<<" isClosed? "<<isClosedPolygon(poly_res, xmin, xmax)<<endl;
  }
  sp.preImage = poly_res;
  if(poly_res.size()==0)
    sp.isEmpty=true;
  if(!sp.isClosed)
    sp.isClosed = isClosedPolygon(poly_res, xmin, xmax);
  return sp;
}

dualSpace addDualPixel(const dualSpace& sp, dualPixel dp) {
  dualSpace new_sp = sp;
  new_sp.vecDualPixel.push_back(dp);
  vector<myRealPoint> poly = dp.polygon;
  vector<myRealPoint> poly_res = polygonIntersection(new_sp.preImage, poly);
  new_sp.preImage = poly_res;
  if(poly_res.size()==0)
    new_sp.isEmpty=true;
  if(!new_sp.isClosed)
    new_sp.isClosed = isClosedPolygon(poly_res, new_sp.xMin, new_sp.xMax);
  return new_sp;
}

dualSpace addPixel(const dualSpace& sp, myPoint pix, double noise){
  dualPixel dp = getDualPixel(pix.px, pix.py, noise, sp.xMin, sp.xMax);
  dualSpace new_sp = addDualPixel(sp, dp);
  return new_sp;
}

myRealPoint getBarryCenter(const vector<myRealPoint>& polygon){
  double cx=0, cy=0;
  for(size_t it=0; it<polygon.size(); it++) {
    cx += polygon.at(it).px;
    cy += polygon.at(it).py;
  }
  return myRealPoint(cx/polygon.size(), cy/polygon.size());
}

myRealPoint getBarryCenter(const dualSpace& sp){
  return getBarryCenter(sp.preImage);
}

pair<myPoint, myPoint> getDomaine(const vector<myPoint>& vecPts){
  int xmin=vecPts.front().px;
  int xmax=vecPts.front().px;
  int ymin=vecPts.front().py;
  int ymax=vecPts.front().py;
  
  for(size_t it=1; it<vecPts.size(); it++) {
    if(xmin>vecPts.at(it).px)
      xmin=vecPts.at(it).px;
    if(xmax<vecPts.at(it).px)
      xmax=vecPts.at(it).px;
    if(ymin>vecPts.at(it).py)
      ymin=vecPts.at(it).py;
    if(ymax<vecPts.at(it).py)
      ymax=vecPts.at(it).py;
  }
  return make_pair(myPoint(xmin,ymin), myPoint(xmax, ymax));
}

