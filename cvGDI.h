#include <cmath>
#include <vector>
#include <string>
#include "opencv2/opencv.hpp"
#pragma once

class cvGDIPattern
{
private:
	int pattern;
	int pattern_counter;
	int scale;
	int bit_n;
	bool dash;
public:
	cvGDIPattern(int _pattern,int _scale=5,int phase=0)
	{
	pattern=_pattern;
	pattern_counter=0;
	scale=_scale;
	bit_n=phase%(scale);
	dash=1;
	}

	void setPhase(double _phase)
	{
	pattern_counter=cvRound(_phase);
	bit_n=pattern_counter%scale;
	}

	int getScale()
	{
	return scale;
	}
	
	int getPhase()
	{
	return pattern_counter;
	}

	bool advance()
	{
		
		if(pattern_counter%scale==0)
		{
			bit_n++;
			if(bit_n>7){bit_n=0;}		
		}
		dash = ( (1 << bit_n) & pattern )>>bit_n;

		pattern_counter++;

		return dash;
	}

	bool get_dash()
	{
		dash = ( (1 << bit_n) & pattern )>>bit_n;
		return dash;
	}

};

class cvGDI
{
public:
	cvGDI(cv::Mat& canvas);
	~cvGDI(void);

	void Line(cv::Point p1, cv::Point p2, cv::Scalar color, float width=1, double alpha=0, cvGDIPattern* pattern=NULL);
	void Rectangle(cv::Point p1, cv::Point p2, cv::Scalar color, float width=1, double alpha=0, cvGDIPattern* pattern=NULL);
	// ---------------
	// TODO
	// ---------------
	// Pen
	// Brush
	// Clear
	// SetBackground
	// Rounded rectangle
	// Rotated rectangle
	// Rotated rounded rectangle
	// Polygon
	// Polyline
	// Arrow
	// Put image
	// ----------------
	void Circle(cv::Point center, float radius, cv::Scalar color, float width=1, double alpha=0, cvGDIPattern* pattern=NULL);
	void Ellipse(cv::Point center, float a,float b, cv::Scalar color, float width=1, double alpha=0, cvGDIPattern* pattern=NULL);
	void Ellipse(cv::Point center, float a,float b,float angle, cv::Scalar color, float width=1, double alpha=0, cvGDIPattern* pattern=NULL);
	
	void CubicSpline(std::vector<cv::Point>& pts, cv::Scalar color, float width=1, double alpha=0, cvGDIPattern* pattern=NULL);
	void QuadSpline(std::vector<cv::Point>& pts, cv::Scalar color, float width, double alpha, cvGDIPattern* pattern=NULL); // TODO add antialiasing/width support
	
private:
	cv::Mat im_buf;
	cv::Mat canvas;
	//--------------------------------------------------------------------                                                                   *
	//
	//--------------------------------------------------------------------
	void setPixel(cv::Mat& img, int x,int y, cv::Scalar color=cv::Scalar::all(0) , double alpha=0);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotLine(cv::Mat& img,int x0, int y0, int x1, int y1, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotLineAA(cv::Mat& img, int x0, int y0, int x1, int y1, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotLineWidth(cv::Mat& img, double x0,double  y0,double x1,double y1,double th, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotCircle(cv::Mat& img, int xm, int ym, int r, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotCircleAA(cv::Mat& img, int xm, int ym, int r, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotEllipse(cv::Mat& img, int xm, int ym, int a, int b, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotOptimizedEllipse(cv::Mat& img, int xm, int ym, int a, int b, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotEllipseRect(cv::Mat& img,int x0, int y0, int x1, int y1, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotEllipseRectWidth(cv::Mat& img, double  x0,double  y0,double  x1,double  y1,double th, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotQuadBezierSeg(cv::Mat& img, int x0, int y0, int x1, int y1, int x2, int y2, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);	
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotQuadBezier(cv::Mat& img, int x0, int y0, int x1, int y1, int x2, int y2, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//-----------------------------------------------------------------------------------------------
	void plotQuadBezierAA(cv::Mat& img, int x0, int y0, int x1, int y1, int x2, int y2, cv::Scalar color, double alpha, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotQuadRationalBezierSeg(cv::Mat & img, int x0, int y0, int x1, int y1, int x2, int y2, double w, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//-----------------------------------------------------------------------------------------------
	void plotQuadRationalBezierWidthSeg(cv::Mat& img, double x0,double y0,double x1,double y1,double x2,double y2,double w,double th, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotQuadRationalBezier(cv::Mat& img, int x0, int y0, int x1, int y1, int x2, int y2, double w, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//----------------------------------------------------------------------------------------------
	void plotQuadRationalBezierWidth(cv::Mat& img, double x0,double  y0,double  x1,double  y1,double  x2,double  y2,double  w,double  th, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotRotatedEllipse(cv::Mat& img, int x, int y, int a, int b, double angle, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotRotatedEllipseWidth(cv::Mat& img, int x, int y, int a, int b, double angle,double th, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotRotatedEllipseRect(cv::Mat& img, int x0, int y0, int x1, int y1, long zd, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotRotatedEllipseRectWidth(cv::Mat& img, int x0, int y0, int x1, int y1, long zd, double th, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//-----------------------------------------------------------------------------------------------
	void plotCubicBezierSegWidth(cv::Mat& img,double x0,double y0,double  x1,double y1,double  x2,double y2,double  x3,double y3,double  th, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotCubicBezierSeg(cv::Mat& img, int x0, int y0, double x1, double y1, double x2, double y2, int x3, int y3, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotCubicBezier(cv::Mat& img, int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);

	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotCubicBezierWidth(cv::Mat& img, double x0, double  y0, double  x1, double  y1, double  x2, double  y2, double  x3, double  y3,double th, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);

	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotCubicBezierAA(cv::Mat& img, int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotEllipseRectAA(cv::Mat& img, int x0, int y0, int x1, int y1, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotQuadBezierSegAA(cv::Mat& img, int x0, int y0, int x1, int y1, int x2, int y2, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotQuadRationalBezierSegAA(cv::Mat& img, int x0, int y0, int x1, int y1,int x2, int y2, double w, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotCubicBezierSegAA(cv::Mat& img, int x0, int y0, double x1, double y1, double x2, double y2, int x3, int y3, cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotQuadSpline(cv::Mat& img, int n, int x[], int y[],cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotCubicSpline(cv::Mat& img, int n, int x[], int y[],cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotCubicSplineWidth(cv::Mat& img, int n, int x[], int y[],double th,cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//------------------------------------------------------------------------------------------------
	//
	//------------------------------------------------------------------------------------------------
	void plotCubicSplineAA(cv::Mat& img, int n, int x[], int y[],cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//-----------------------------------------------------------------------------------------------------
	//
	//-----------------------------------------------------------------------------------------------------
	void plotQuadSplineAA(cv::Mat& img, int n, int x[], int y[],cv::Scalar color, double alpha=0, cvGDIPattern* pattern=NULL);
	//-----------------------------------------------------------------------------------------------------
	//
	//-----------------------------------------------------------------------------------------------------
};

