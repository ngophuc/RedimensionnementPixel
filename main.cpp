#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

#include "DualPixel.h"
#include "DualSpace.h"
#include "UtilityFunctions.h"

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/boards/Board2D.h"

using namespace DGtal;
using namespace Z2i;
using namespace std;

void testDrawDualPixel(vector<Point> vecPts, double noise) {
  Board2D board;
  Point d1(-5,-5);
  Point d2(5,5);
  Domain domain( d1, d2 );
  board << DGtal::SetMode("PointVector", "Both");
  board << domain;// << d1 << d2;
  board.setPenColor(DGtal::Color::Purple);
  board.fillCircle(0,0,0.2);
  for(auto p : vecPts) {
    board << p;
    /*
    //Compute for lines of pixel
    dualPixel dP = getLinePixel(p[0], p[1], noise);
    vector<vector<double> > vecLine = {dP.line1, dP.line2, dP.line3, dP.line4};
    myRealPoint ip1 = intersectionTwoLines(vecSeg.at(0), vecSeg.at(1));//intersection between S1 & S2
    myRealPoint ip2 = intersectionTwoLines(vecSeg.at(2), vecSeg.at(3));//intersection between S3 & S4  
    board.setLineWidth(1.0);
    board.setPenColor(DGtal::Color::Red);
    pair<RealPoint, RealPoint> s1 = getSegment(vecSeg.at(0), ip1.px, d2[0]);
    board.drawLine(s1.first[0], s1.first[1], s1.second[0], s1.second[1]);
    board.setPenColor(DGtal::Color::Cyan);
    pair<RealPoint, RealPoint> s2 = getSegment(vecSeg.at(1), d1[0], ip2.px);
    board.drawLine(s2.first[0], s2.first[1], s2.second[0], s2.second[1]);
    
    board.setPenColor(DGtal::Color::Blue);
    pair<RealPoint, RealPoint> s3 = getSegment(vecSeg.at(2), d1[0], ip2.px);
    board.drawLine(s3.first[0], s3.first[1], s3.second[0], s3.second[1]);
    board.setPenColor(DGtal::Color::Green);
    pair<RealPoint, RealPoint> s4 = getSegment(vecSeg.at(3), ip1.px, d2[0]);
    board.drawLine(s4.first[0], s4.first[1], s4.second[0], s4.second[1]);
    
    //Draw intersection points
    board.setPenColor(DGtal::Color::Red);
    board.drawCircle(ip1.px, ip1.py, 0.1);
    board.setPenColor(DGtal::Color::Blue);
    board.drawCircle(ip2.px, ip2.py, 0.1);
    */
    /*
    //Compute for segments of pixel
    dualPixel dP = getSegmentPixel(p[0], p[1], noise, d1[0], d2[0]);
    //Draw segments
    board.setLineWidth(1.0);
    board.setPenColor(DGtal::Color::Red);
    pair<myRealPoint, myRealPoint> s1 = dP.seg1.first;
    board.drawLine(s1.first.px, s1.first.py, s1.second.px, s1.second.py);
    board.drawArrow(s1.first.px, s1.first.py, s1.first.px, s1.first.py+dP.seg1.second);
    board.setPenColor(DGtal::Color::Cyan);
    pair<myRealPoint, myRealPoint> s2 = dP.seg2.first;
    board.drawLine(s2.first.px, s2.first.py, s2.second.px, s2.second.py);
    board.drawArrow(s2.second.px, s2.second.py, s2.second.px, s2.second.py+dP.seg2.second);
    
    board.setPenColor(DGtal::Color::Blue);
    pair<myRealPoint, myRealPoint> s3 = dP.seg3.first;
    board.drawLine(s3.first.px, s3.first.py, s3.second.px, s3.second.py);
    board.drawArrow(s3.first.px, s3.first.py, s3.first.px, s3.first.py+dP.seg3.second);
    board.setPenColor(DGtal::Color::Green);
    pair<myRealPoint, myRealPoint> s4 = dP.seg4.first;
    board.drawLine(s4.first.px, s4.first.py, s4.second.px, s4.second.py);
    board.drawArrow(s4.second.px, s4.second.py, s4.second.px, s4.second.py+dP.seg4.second);
    */
    dualPixel dP2 = getDualPixel(p[0], p[1], noise, d1[0], d2[0]);
    vector<LibBoard::Point> poly;
    for(size_t it=0; it<dP2.polygon.size(); it++) {
      myRealPoint p = dP2.polygon.at(it);
      poly.push_back(LibBoard::Point(p.px, p.py));
    }
    board.drawClosedPolyline(poly);
  }
  board.saveSVG("DualPixel.svg");
}

void testDrawDualPixelIntersection(dualPixel dP1, dualPixel dP2) {
  Board2D board;
  Point d1(-5,-5);
  Point d2(5,5);
  Domain domain( d1, d2 );
  board << DGtal::SetMode("PointVector", "Both");
  board << domain;// << d1 << d2;
  board.setPenColor(DGtal::Color::Purple);
  board.fillCircle(0,0,0.2);
  
  vector<myRealPoint> poly1 = dP1.polygon;
  vector<myRealPoint> poly2 = dP2.polygon;
  vector<myRealPoint> res = polygonIntersection(poly1, poly2);
  
  vector<LibBoard::Point> poly;
  for(size_t it=0; it<res.size(); it++) {
    myRealPoint p = res.at(it);
    poly.push_back(LibBoard::Point(p.px, p.py));
  }
  board.drawClosedPolyline(poly);
  board.saveSVG("IntersectionDualPixels.svg");
}

void testDrawDualPixelIntersection(vector<Point> vecPts, double noise) {
  Board2D board;
  Point d1(-5,-5);
  Point d2(5,5);
  Domain domain( d1, d2 );
  board << DGtal::SetMode("PointVector", "Both");
  board << domain;// << d1 << d2;
  board.setPenColor(DGtal::Color::Purple);
  board.fillCircle(0,0,0.2);
  
  assert(vecPts.size()>0);
  Point p1 = vecPts.front();
  dualPixel dP1 = getDualPixel(p1[0], p1[1], noise, d1[0], d2[0]);
  vector<myRealPoint> poly1 = dP1.polygon;
  vector<myRealPoint> res = poly1;
  for(size_t it=1; it<vecPts.size(); it++) {
    Point p2 = vecPts.at(it);
    dualPixel dP2 = getDualPixel(p2[0], p2[1], noise, d1[0], d2[0]);
    vector<myRealPoint> poly2 = dP2.polygon;
    res = polygonIntersection(res, poly2);
  }
  
  vector<LibBoard::Point> poly;
  for(size_t it=0; it<res.size(); it++) {
    myRealPoint p = res.at(it);
    poly.push_back(LibBoard::Point(p.px, p.py));
  }
  board.drawClosedPolyline(poly);
  board.saveSVG("IntersectionDualPixels.svg");
}

void testDualSpace() {
  myPoint p1(0,0);
  myPoint p2(1,1);
  myPoint p3(1,2);
  myPoint p4(2,3);
  
  myPoint p5(-1,2);
  vector<double> vecNoise = {1.0, 1.0, 1.0};//, 0.0, 0.0};//
  vector<myPoint> vecPts = {p1, p2, p3};//, p4, p5};//
  dualSpace ds = getDualSpace(vecPts, vecNoise, -5, 5);
  
  
  if(ds.isClosed==true)
    cout<<"Preimage region is closed"<<endl;
  else
    cout<<"Preimage region is open"<<endl;
  if(ds.isEmpty==true)
    cout<<"Preimage region is empty"<<endl;
  else
    cout<<"Preimage region is not empty"<<endl;
  
  if(!ds.isEmpty==true){ //(res.size()!=0)
    Board2D board;
    Point d1(-5,-5);
    Point d2(5,5);
    Domain domain( d1, d2 );
    board << DGtal::SetMode("PointVector", "Both");
    board << domain;// << d1 << d2;
    vector<myRealPoint> res = ds.preImage;
    vector<LibBoard::Point> poly;
    for(size_t it=0; it<res.size(); it++) {
      myRealPoint p = res.at(it);
      poly.push_back(LibBoard::Point(p.px, p.py));
    }
    board.drawClosedPolyline(poly);
    board.saveSVG("DualSpace.svg");
  }
  
  //Add new point
  dualSpace ds_new = addPixel(ds, p4,vecNoise.front());
  if(ds_new.isClosed==true)
    cout<<"New Preimage region is closed"<<endl;
  else
    cout<<"New Preimage region is open"<<endl;
  if(ds_new.isEmpty==true)
    cout<<"New Preimage region is empty"<<endl;
  else
    cout<<"New Preimage region is not empty"<<endl;
  if(!ds_new.isEmpty==true){ //(res.size()!=0)
    Board2D board;
    Point d1(-5,-5);
    Point d2(5,5);
    Domain domain( d1, d2 );
    board << DGtal::SetMode("PointVector", "Both");
    board << domain;// << d1 << d2;
    vector<myRealPoint> res = ds_new.preImage;
    vector<LibBoard::Point> poly;
    for(size_t it=0; it<res.size(); it++) {
      myRealPoint p = res.at(it);
      poly.push_back(LibBoard::Point(p.px, p.py));
    }
    board.drawClosedPolyline(poly);
    board.saveSVG("DualSpace_new.svg");
  }
}

Board2D drawDualPixel(const dualPixel& dp, Board2D& aBoard, Color c = Color(255, 0, 0, 100)) {//Color::Red
  vector<myRealPoint> dpPoly = drawDualPixel(dp);
  vector<LibBoard::Point> poly;
  for(size_t it=0; it<dpPoly.size(); it++) {
    myRealPoint p = dpPoly.at(it);
    poly.push_back(LibBoard::Point(p.px, p.py));
  }
  aBoard.setPenColor(c);
  aBoard.fillPolyline(poly);
  return aBoard;
}

Board2D drawDualPixels(const vector<dualPixel>& dp, Board2D& aBoard, Color c = Color(255, 0, 0, 100)) {
  //Draw the seq of pixels
  for(size_t it=1; it<dp.size(); it++) {
    int p1x = dp.at(it-1).px;
    int p1y = dp.at(it-1).py;
    int p2x = dp.at(it).px;
    int p2y = dp.at(it).py;
    aBoard.drawArrow(p1x, p1y, p2x, p2y);
  }
  //Draw redimensionnement of pixels
  for(size_t it=0; it<dp.size(); it++)
    drawDualPixel(dp.at(it), aBoard, c);
  return aBoard;
}

Board2D drawDualPixels(const vector<dualPixel>& dp, Board2D& aBoard, int idB, int idE, Color c = Color(255, 0, 0, 100)) {
  //Draw the seq of pixels
  for(size_t it=1; it<dp.size(); it++) {
    int p1x = dp.at(it-1).px;
    int p1y = dp.at(it-1).py;
    int p2x = dp.at(it).px;
    int p2y = dp.at(it).py;
    aBoard.drawArrow(p1x, p1y, p2x, p2y);
  }
  //Draw redimensionnement of pixels
  for(size_t it=idB; it<=idE; it++)
    drawDualPixel(dp.at(it), aBoard, c);
  return aBoard;
}

Board2D drawDualSpace (const dualSpace& ds, Board2D& aBoard, Color c = Color(255, 0, 0, 100)) {
  if(ds.isEmpty!=true) {
    aBoard.setPenColor(c);
    vector<myRealPoint> res = ds.preImage;
    vector<LibBoard::Point> poly;
    for(size_t it=0; it<res.size(); it++) {
      myRealPoint p = res.at(it);
      poly.push_back(LibBoard::Point(p.px, p.py));
    }
    aBoard.drawClosedPolyline(poly);
  }
  return aBoard;
}

Board2D drawASolution (myRealPoint line, double xmin, double xmax, Board2D& aBoard, Color c = Color(255, 0, 0, 100), double width=1) {
  //line: y=ax+b
  double ymin = line.px*xmin + line.py;
  double ymax = line.px*xmax + line.py;
  aBoard.setPenColor(c);
  aBoard.setLineWidth(width);
  aBoard.drawLine(xmin, ymin, xmax, ymax);
  return aBoard;
}
  
void testDualSpace2(){
  Point d1(-5,-5);
  Point d2(5,5);
  /*
   myPoint p1(70,26);
   myPoint p2(69, 26);
   myPoint p3(69, 27);
   myPoint p4(70, 27);
   */
  /*
  myPoint p1(0,0);
  myPoint p2(1,1);
  myPoint p3(1,2);
  myPoint p4(2,2);
  */
  myPoint p1(0,0);
  myPoint p2(0,-1);
  myPoint p3(1,-1);
  //myPoint p5(-1,2);
  
  vector<myPoint> vecPt = {p1, p2};//, p3};//, p4};
  vector<double> vecNoise = {1.0, 1.0};//, 1.0};//{1.1, 1.4, 1.3, 1.0};
  dualSpace dualParam = getDualSpace(vecPt, vecNoise, d1[0], d2[0]);
  vector<dualPixel> vecDualPixel = dualParam.vecDualPixel;
  Board2D boardDualPixel;
  Domain domain( d1, d2 );
  boardDualPixel << DGtal::SetMode("PointVector", "Both");
  boardDualPixel << domain;
  drawDualPixels(vecDualPixel, boardDualPixel);
  myRealPoint barryCenter = getBarryCenter(dualParam);
  double xmin = vecPt.front().px-1;
  double xmax = vecPt.back().px+1;
  drawASolution(barryCenter, xmin, xmax, boardDualPixel, Color::Blue);
  
  dualPixel dp3 = getDualPixel(p3.px, p3.py);
  drawDualPixel(dp3, boardDualPixel, Color::Blue);
  dualSpace ds3 = initDualSpace(p3, 1.0);
  
  boardDualPixel.saveSVG("vecDualPixels.svg");
  
  Board2D boardDualSpace;
  boardDualSpace << DGtal::SetMode("PointVector", "Both");
  boardDualSpace << domain;
  drawDualSpace(dualParam, boardDualSpace, Color::Red);
  drawDualSpace(ds3, boardDualSpace, Color::Blue);
  boardDualSpace.setPenColor(Color::Blue);
  boardDualSpace.fillCircle(barryCenter.px, barryCenter.py, 0.1);
  boardDualSpace.saveSVG("vecDualSpace.svg");
}

void testDualSpace3() {  
  myPoint p1(10,10);
  myPoint p2(11,11);
  double xmin=-0.1, xmax=0.1;
  /*
  myPoint p1(0,0);
  myPoint p2(1,1);
  int xmin=-5, xmax=5;
   */
  vector<double> vecNoise = {1.0};//, 0.0, 0.0};//
  vector<myPoint> vecPts = {p1};//, p4, p5};//
  dualSpace ds = getDualSpace(vecPts, vecNoise, xmin, xmax);
  
  if(ds.isClosed==true)
    cout<<"Preimage region is closed"<<endl;
  else
    cout<<"Preimage region is open"<<endl;
  if(ds.isEmpty==true)
    cout<<"Preimage region is empty"<<endl;
  else
    cout<<"Preimage region is not empty"<<endl;
  
  Point d1(xmin,xmin);
  Point d2(xmax,xmax);
  if(!ds.isEmpty==true){ //(res.size()!=0)
    Board2D board;
    Domain domain( d1, d2 );
    board << DGtal::SetMode("PointVector", "Both");
    board << domain;// << d1 << d2;
    vector<myRealPoint> res = ds.preImage;
    vector<LibBoard::Point> poly;
    for(size_t it=0; it<res.size(); it++) {
      myRealPoint p = res.at(it);
      poly.push_back(LibBoard::Point(p.px, p.py));
    }
    board.drawClosedPolyline(poly);
    board.saveSVG("DualSpace.svg");
  }
  
  //Add new point
  dualSpace ds_new = addPixel(ds, p2,vecNoise.front());
  if(ds_new.isClosed==true)
    cout<<"New Preimage region is closed"<<endl;
  else
    cout<<"New Preimage region is open"<<endl;
  if(ds_new.isEmpty==true)
    cout<<"New Preimage region is empty"<<endl;
  else
    cout<<"New Preimage region is not empty"<<endl;
  if(!ds_new.isEmpty==true){ //(res.size()!=0)
    Board2D board;
    Domain domain( d1, d2 );
    board << DGtal::SetMode("PointVector", "Both");
    board << domain;// << d1 << d2;
    vector<myRealPoint> res = ds_new.preImage;
    vector<LibBoard::Point> poly;
    for(size_t it=0; it<res.size(); it++) {
      myRealPoint p = res.at(it);
      poly.push_back(LibBoard::Point(p.px, p.py));
    }
    board.drawClosedPolyline(poly);
    board.saveSVG("DualSpace_new.svg");
  }
  
}

void displayPolygon(const vector<myRealPoint>& vecPts) {
  if(vecPts.size()!=0) {
    cout<<"Polygon [";
    for(size_t it=0; it<vecPts.size()-1; it++)
      cout<<"("<<vecPts.at(it).px<<","<<vecPts.at(it).py<<") -> ";
    cout<<"("<<vecPts.back().px<<","<<vecPts.back().py<<")]"<<endl;
  }
  else
    cout<<"Polygon []"<<endl;
}

Board2D drawPolygon (const vector<myRealPoint>& vecPts, Board2D& aBoard, Color c = Color(255, 0, 0, 100), double width=1) {
  if(vecPts.size()!=0) {
    aBoard.setPenColor(c);
    aBoard.setLineWidth(width);
    vector<LibBoard::Point> poly;
    for(size_t it=0; it<vecPts.size(); it++) {
      myRealPoint p = vecPts.at(it);
      poly.push_back(LibBoard::Point(p.px, p.py));
    }
    aBoard.drawClosedPolyline(poly);
  }
  return aBoard;
}

void testDetection() {
  /********** read data ***************/
  char fileContour[] = "../Data/pentagonNoise4_0.sdp";
  vector<vector<Point>> vecContour = readFile(fileContour);
  assert(vecContour.size()!=0);
  vector<Point> aContour = vecContour.front();
  assert(aContour.size()!=0);
  char noiseLevelMTFile[] = "../Data/pentagonNoise4MeanThicknessOpen_20_10_0.txt";
  vector<double> vecMT=readMeanindfulThicknessFile(noiseLevelMTFile);//*sqrt(2)
  assert(aContour.size()==vecMT.size());
  
  //Get dual pixels
  vector<myPoint> vecPt;
  for(size_t it=0; it<aContour.size(); it++)
    vecPt.push_back(myPoint(aContour.at(it)[0], aContour.at(it)[1]));
  vector<double> vecNoise;
  for(size_t it=0; it<vecMT.size(); it++) {
    vecNoise.push_back(vecMT.at(it));
    //cout<<"vecMT.at(it)="<<vecMT.at(it)<<" vs. vecNoise.at(it)="<<vecNoise.at(it)<<endl;
  }
  pair<myPoint, myPoint> bb=getDomaine(vecPt);
  Point d1(bb.first.px-1, bb.first.py-1);
  Point d2(bb.second.px+1, bb.second.py+1);
  //Point d1(60, 20);
  //Point d2(75, 35);
  vector<dualPixel> vecDualPixel;
  for(size_t it=0; it<vecPt.size(); it++){
    dualPixel dp = getDualPixel(vecPt.at(it).px, vecPt.at(it).py, vecNoise.at(it), -1.0, 1.0);//bb.first.px, bb.second.px);//, d1[0], d2[0]);//
    vecDualPixel.push_back(dp);
  }
  
  //Get dual space
  vector<dualSpace> vecDualSpace;
  bool isDone = true;
  while(isDone) {
    //dualSpace dp_tmp = getDualSpace(myPts, myNoise, bb.first.px, bb.second.px);
    dualSpace dp_tmp = initDualSpace(vecPt.front(), vecNoise.front(), -1.0, 1.0);//bb.first.px, bb.second.px);//, d1[0], d2[0]);//
    cout<<0<<": dp_tmp.isEmpty="<<dp_tmp.isEmpty<<endl;
    //displayPolygon(dp_tmp.preImage);
    int it=1;
    dualSpace dp;
    /*
    Board2D board_tmp;
    drawPolygon(dp_tmp.preImage, board_tmp);
    board_tmp.saveSVG("tmp.svg");
    board_tmp.clear();
    */
    while(it<vecPt.size() && !dp_tmp.isEmpty) {//vecPt.size()
      dp = dp_tmp;
      dp_tmp = addPixel(dp_tmp, vecPt.at(it), vecNoise.at(it));
      /*
      displayPolygon(dp_tmp.preImage);
      drawPolygon(dp_tmp.preImage, board_tmp);
      board_tmp.saveSVG("tmp.svg");
      board_tmp.clear();
       */
      cout<<it<<" ==> dp_tmp.preImage.size="<<dp_tmp.preImage.size()<<" and dp_tmp.isEmpty="<<dp_tmp.isEmpty<<endl;
      it++;
    }
    vecDualSpace.push_back(dp);
    isDone=false;
  }
  cout<<"Final ==> dp_tmp.preImage.size="<<vecDualSpace.front().preImage.size()<<" and dp_tmp.isEmpty="<<vecDualSpace.front().isEmpty<<endl;
  //displayPolygon(vecDualSpace.front().preImage);
  
  //Draw dual pixels with noise
  Board2D boardDualPixel;
  Domain domain( d1, d2 );
  boardDualPixel << DGtal::SetMode("PointVector", "Both");
  boardDualPixel << domain;
  drawDualPixels(vecDualPixel, boardDualPixel);
  myRealPoint barryCenter = getBarryCenter(vecDualSpace.front());
  double xmin = vecDualSpace.front().vecDualPixel.front().px-1;//bb.first.px;//vecPt.front().px-1;
  double xmax = vecDualSpace.front().vecDualPixel.back().px+1;//bb.second.px;//vecPt.back().px+1;
  drawASolution(barryCenter, xmin, xmax, boardDualPixel, Color::Black);
  boardDualPixel.saveSVG("vecContourDualPixels.svg");
  
  //Draw dual space
  Board2D boardDualSpace;
  boardDualSpace << DGtal::SetMode("PointVector", "Both");
  //boardDualSpace << domain;
  drawDualSpace(vecDualSpace.front(), boardDualSpace, Color::Red);
  boardDualSpace.setPenColor(Color::Blue);
  //myRealPoint barryCenter = getBarryCenter(vecDualSpace.front());
  //double xmin = vecPt.front().px-1;
  //double xmax = vecPt.back().px+1;
  //drawASolution(barryCenter, xmin, xmax, boardDualPixel, Color::Blue);
  boardDualSpace.fillCircle(barryCenter.px, barryCenter.py, 0.1);
  boardDualSpace.saveSVG("vecContourDualSpace.svg");
}

vector<pair<RealPoint, RealPoint> > polygonalisation() {
  vector<pair<RealPoint, RealPoint> > vecSeg;
  /********** read data ***************/
  char fileContour[] = "../Data/pentagonNoise4_0.sdp";
  vector<vector<Point>> vecContour = readFile(fileContour);
  assert(vecContour.size()!=0);
  vector<Point> aContour = vecContour.front();
  assert(aContour.size()!=0);
  char noiseLevelMTFile[] = "../Data/pentagonNoise4MeanThicknessOpen_20_10_0.txt";
  vector<double> vecMT=readMeanindfulThicknessFile(noiseLevelMTFile);//*sqrt(2)
  assert(aContour.size()==vecMT.size());
  
  //Get dual pixels
  vector<myPoint> vecPt;
  for(size_t it=0; it<aContour.size(); it++)
    vecPt.push_back(myPoint(aContour.at(it)[0], aContour.at(it)[1]));
  vector<double> vecNoise;
  for(size_t it=0; it<vecMT.size(); it++) {
    vecNoise.push_back(vecMT.at(it));
    //cout<<"vecMT.at(it)="<<vecMT.at(it)<<" vs. vecNoise.at(it)="<<vecNoise.at(it)<<endl;
  }
  pair<myPoint, myPoint> bb=getDomaine(vecPt);
  Point d1(bb.first.px-1, bb.first.py-1);
  Point d2(bb.second.px+1, bb.second.py+1);
  //Point d1(60, 20);
  //Point d2(75, 35);
  vector<dualPixel> vecDualPixel;
  for(size_t it=0; it<vecPt.size(); it++){
    dualPixel dp = getDualPixel(vecPt.at(it).px, vecPt.at(it).py, vecNoise.at(it), -1.0, 1.0);//bb.first.px, bb.second.px);//, d1[0], d2[0]);//
    vecDualPixel.push_back(dp);
  }
  
  int idStart=396;//300;
  vector<pair<int,int> > idSeg;
  vector<dualSpace> vecDualSpace;
  //Get dual spaces
  //while(idStart<vecDualPixel.size()-1) {
  while(vecDualSpace.size()<5) {
    bool isDone = true;
    int idEnd=idStart+1;
    while(isDone && idEnd<vecPt.size()) {
      dualSpace dp_tmp = initDualSpace(vecPt.at(idStart), vecNoise.at(idStart), -1.0, 1.0);
      dualSpace dp;
      while(idEnd<vecPt.size() && !dp_tmp.isEmpty) {
        dp = dp_tmp;
        dp_tmp = addPixel(dp_tmp, vecPt.at(idEnd), vecNoise.at(idEnd));
        idEnd++;
      }
      vecDualSpace.push_back(dp);
      idSeg.push_back(make_pair(idStart, idEnd-1));
      cout<<"idStart="<<idStart<<" --> "<<"idEnd="<<idEnd-1<<endl;
      isDone=false;
      
    }
    idStart = idEnd-1;
  }
  
  //Retrieve segments
  for(size_t it=0; it<vecDualSpace.size(); it++) {
    myRealPoint line = getBarryCenter(vecDualSpace.at(it));
    double xmin = vecPt.at(idSeg.at(it).first).px;//vecDualSpace.at(it).vecDualPixel.front().px-1;
    double xmax = vecPt.at(idSeg.at(it).second).px;//vecDualSpace.at(it).vecDualPixel.back().px+1;
    double ymin = line.px*xmin + line.py;
    double ymax = line.px*xmax + line.py;
    RealPoint p1(xmin, ymin);
    RealPoint p2(xmax, ymax);
    vecSeg.push_back(make_pair(p1, p2));
    cout<<"p1="<<vecPt.at(idSeg.at(it).first)<<" --> p2="<<vecPt.at(idSeg.at(it).second)<<endl;
    cout<<"l:"<<line<<" ==> s1="<<p1<<" --> s2="<<p2<<endl;
  }
  
  //Draw dual pixels with noise + polygon
  Board2D boardDualPixel;
  Domain domain( d1, d2 );
  boardDualPixel << DGtal::SetMode("PointVector", "Both");
  boardDualPixel << domain;
  //drawDualPixels(vecDualPixel, boardDualPixel);
  drawDualPixels(vecDualPixel, boardDualPixel, 396, 398);

  boardDualPixel.setLineStyle(LibBoard::Shape::LineStyle::SolidStyle);
  boardDualPixel.setLineWidth(10.0);
  for(size_t it=0; it<vecSeg.size(); it++) {
    RealPoint p1=vecSeg.at(it).first;
    RealPoint p2=vecSeg.at(it).second;
    //cout<<"====> s1="<<p1<<" --> s2="<<p2<<endl;
    if(it%2==0)
      boardDualPixel.setPenColor(DGtal::Color::Green);
    else
      boardDualPixel.setPenColor(DGtal::Color::Blue);
    boardDualPixel.drawLine(p1[0], p1[1], p2[0], p2[1]);
  }
  boardDualPixel.saveSVG("PolygonDualPixels.svg");
  return vecSeg;
}

int main(int , char**) {
  testDualSpace2();
  return 0;
  /*
  Point p1(10,10);
  Point p2(11,11);
  double noise=1.0;
  double vmin=-0.5, vmax=0.5;
  dualPixel dP1 = getDualPixel(p1[0], p1[1], noise, vmin, vmax);
  dualPixel dP2 = getDualPixel(p2[0], p2[1], noise, vmin, vmax);
  testDrawDualPixelIntersection(dP1, dP2);
  return 0;
  */
  //testDualSpace();
  //testDualSpace2();
  //testDualSpace3();
  //return 0;
  
  polygonalisation();
  return 0;
  /*
  testBoostIntersetion();
  testDualSpace();
  testDualSpace2();
  return 0;
  */
  /*
  double noise=1.0;
  Point p1(0,0);
  Point p2(1,1);
  Point p3(1,2);
  Point p4(-1,2);
  vector<Point> vecPt = {p1, p2, p3};//, p3, p4};
  testDrawDualPixel(vecPt, noise);
  testDrawDualPixelIntersection(vecPt, noise);
  return 0;
  */
  /*
  dualPixel dP1 = getDualPixel(p1[0], p1[1], noise, d1[0], d2[0]);
  dualPixel dP2 = getDualPixel(p2[0], p2[1], noise, d1[0], d2[0]);
  testDrawDualPixelIntersection(dP1, dP2);
  return 0;
  */
}
