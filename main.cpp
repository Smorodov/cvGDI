#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <numeric>
#include "opencv2/opencv.hpp"
#include "spline.h"
#include <cmath>
#include "cvGDI.h"
using namespace std;
using namespace cv;
// 
// http://members.chello.at/easyfilter/bresenham.html
//-----------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------
int main ( int argc, char** argv )
{
	

	Mat canvas=Mat::zeros(512,512,CV_8UC3);
	canvas=cv::Scalar(255,255,255);

	cvGDI gdi(canvas);

	vector<cv::Point> pts;

	pts.push_back( cv::Point(100,100));
	pts.push_back( cv::Point(200,150));
	pts.push_back( cv::Point(200,300));
	pts.push_back( cv::Point(20,350));
	pts.push_back( cv::Point(130,90));
	cvGDIPattern pattern(0x0f,3,0);
	cvGDIPattern pattern_prev(0x0f,3,0);
for(int k=0;k<1000;++k)
{
	
	pattern_prev=pattern;

	canvas=cv::Scalar(255,255,255);

	gdi.CubicSpline(pts,cv::Scalar(0,255,0),5,0.5);
	
	gdi.QuadSpline(pts,cv::Scalar(0,0,255),5,0.5);

	gdi.Line(cv::Point(100,100),cv::Point(120,200),Scalar(255,0,0),1,0.5,&pattern);

	gdi.Line(cv::Point(100,100),cv::Point(220,300),Scalar(255,0,0),2,0.5,&pattern);

	gdi.Rectangle(cv::Point(200,120),cv::Point(400,240),Scalar(255,0,0),1,0.5,&pattern);

	gdi.Circle(cv::Point(200,100),70,Scalar(0,0,0),1,0.5,&pattern);

	gdi.Circle(cv::Point(100,100),70,Scalar(255,255,0),6,0.5,&pattern);

	gdi.Ellipse(cv::Point(100,100),70,30,Scalar(255,0,255),2,0.5,&pattern);

	gdi.Ellipse(cv::Point(100,100),70,30,0.6,Scalar(255,0,255),2,0.5,&pattern);

	pattern=pattern_prev;

	for(int i=0;i<pts.size();++i)
	{
	circle(canvas,pts[i],3,cv::Scalar(0,0,255),-1 ,LINE_AA);
	}

	pattern.advance();

	imshow("canvas",canvas);
	waitKey(20);
}
	waitKey();

}