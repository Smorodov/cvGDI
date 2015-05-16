#include "cvGDI.h"
using namespace std;
using namespace cv;
//--------------------------------------------------------------------                                                                   
//
//--------------------------------------------------------------------
cvGDI::cvGDI(cv::Mat& _canvas)
{
	im_buf=canvas.clone();
	canvas=_canvas;
}
//--------------------------------------------------------------------                                                                   
//
//--------------------------------------------------------------------
cvGDI::~cvGDI(void)
{
}
//--------------------------------------------------------------------                                                                   
//
//--------------------------------------------------------------------
void cvGDI::Line(cv::Point p1, cv::Point p2, cv::Scalar color, float width, double alpha, cvGDIPattern* pattern)
{
	im_buf=canvas.clone();
	if(width>1)
	{		
		plotLineWidth(canvas,p1.x,p1.y,p2.x,p2.y,width, color,alpha,pattern);
	}
	else
	{
		plotLineAA(canvas,p1.x,p1.y,p2.x,p2.y, color,alpha,pattern);
	}
}
//--------------------------------------------------------------------                                                                   
//
//--------------------------------------------------------------------
void cvGDI::Rectangle(cv::Point p1, cv::Point p2, cv::Scalar color, float width, double alpha, cvGDIPattern* pattern)
{
	im_buf=canvas.clone();
	cv::Point pa(p2.x,p1.y);
	cv::Point pb(p1.x,p2.y);
	if(width>1)
	{		
		plotLineWidth(canvas,p1.x,p1.y+floor(width/2.0),pa.x,pa.y+floor(width/2.0),width, color,alpha,pattern);

		plotLineWidth(canvas,pa.x,pa.y,p2.x,p2.y,width, color,alpha);

		plotLineWidth(canvas,p2.x,p2.y-floor(width/2.0),pb.x,pb.y-floor(width/2.0),width, color,alpha,pattern);

		plotLineWidth(canvas,pb.x,pb.y,p1.x,p1.y,width, color,alpha,pattern);
	}
	else
	{
		plotLineAA(canvas,p1.x,p1.y,pa.x,pa.y,color,alpha,pattern);
		plotLineAA(canvas,pa.x,pa.y,p2.x,p2.y,color,alpha,pattern);
		plotLineAA(canvas,p2.x,p2.y,pb.x,pb.y,color,alpha,pattern);
		plotLineAA(canvas,pb.x,pb.y,p1.x,p1.y,color,alpha,pattern);
	}
}
//--------------------------------------------------------------------                                                                   
//
//-------------------------------------------------------------------
void cvGDI::Circle(cv::Point center, float radius, cv::Scalar color, float width, double alpha, cvGDIPattern* pattern)
{
	im_buf=canvas.clone();
	if(width>1)
	{
		plotEllipseRectWidth(canvas, center.x-radius,center.y-radius,center.x+radius,center.y+radius,width, color, alpha,pattern);
	}else
	{
		plotCircleAA(canvas, center.x, center.y, radius, color,alpha,pattern);
	}

}
//--------------------------------------------------------------------                                                                   
//
//------------------------------------------------------------------
void cvGDI::Ellipse(cv::Point center, float a,float b, cv::Scalar color, float width, double alpha, cvGDIPattern* pattern)
{
	im_buf=canvas.clone();
	if(width>1)
	{
		plotEllipseRectWidth(canvas, center.x-a,center.y-b,center.x+a,center.y+b,width, color, alpha,pattern);
	}else
	{
		plotOptimizedEllipse(canvas, center.x, center.y, a, b, color, alpha,pattern);
	}
}
//--------------------------------------------------------------------                                                                   
//
//------------------------------------------------------------------
void cvGDI::Ellipse(cv::Point center, float a,float b,float angle, cv::Scalar color, float width, double alpha, cvGDIPattern* pattern)
{
	im_buf=canvas.clone();
	if(width>1)
	{
		plotRotatedEllipseWidth(canvas, center.x,center.y,a,b,angle,width, color, alpha,pattern);
	}else
	{
		plotRotatedEllipse(canvas, center.x, center.y, a, b, angle, color, alpha,pattern);
	}
}
//--------------------------------------------------------------------                                                                   
//
//--------------------------------------------------------------------
void cvGDI::CubicSpline(std::vector<cv::Point>& pts, cv::Scalar color, float width, double alpha, cvGDIPattern* pattern)
{
	int* x=new int[pts.size()];
	int* y=new int[pts.size()];
	for(int i=0;i<pts.size();++i)
	{
		x[i]=pts[i].x;
		y[i]=pts[i].y;
	}
	im_buf=canvas.clone();

	if(width>1)
	{
		plotCubicSplineWidth(canvas, pts.size()-1, x, y,width,color,alpha,pattern);
	}else
	{
		plotCubicSplineAA(canvas,pts.size()-1,x,y,color,alpha,pattern);
	}
	delete x;
	delete y;
}
//--------------------------------------------------------------------                                                                   
//
//--------------------------------------------------------------------
void cvGDI::QuadSpline(std::vector<cv::Point>& pts, cv::Scalar color, float width, double alpha, cvGDIPattern* pattern)
{
	int* x=new int[pts.size()];
	int* y=new int[pts.size()];
	for(int i=0;i<pts.size();++i)
	{
		x[i]=pts[i].x;
		y[i]=pts[i].y;
	}
	im_buf=canvas.clone();

	if(width>1)
	{
		plotQuadSplineAA(canvas, pts.size()-1, x, y,color,alpha,pattern);
	}else
	{
		plotQuadSplineAA(canvas,pts.size()-1,x,y,color,alpha,pattern);
	}
	delete x;
	delete y;
}
//--------------------------------------------------------------------                                                                   
//
//--------------------------------------------------------------------
void cvGDI::setPixel(Mat& img, int x,int y, cv::Scalar color , double alpha)
{
	if(alpha > 1.0){alpha=1.0;}
	if(alpha < 0) {alpha=0;}
	if(x>=0 && x<img.cols && y>=0 && y<img.rows)
	{
		Vec3b c=im_buf.at<Vec3b>(y,x);
		cv::Scalar result_color=alpha*cv::Scalar(c[0],c[1],c[2])+(1.0-alpha)*color;
		img.at<Vec3b>(y,x)=Vec3b(result_color[0],result_color[1],result_color[2]);
	}
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotLine(Mat& img,int x0, int y0, int x1, int y1, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	int dx =  abs(x1-x0), sx = x0<x1 ? 1 : -1;
	int dy = -abs(y1-y0), sy = y0<y1 ? 1 : -1;
	int err = dx+dy, e2;                                   // error value e_xy
	bool dash=1;
	for (;;)                                                           // loop
	{
		if(pattern!=NULL)
		{
			dash=pattern->get_dash();
			if(dash)
			{
				setPixel(img,x0,y0, color, alpha);
			}
			pattern->advance();
		}
		e2 = 2*err;
		if (e2 >= dy)                                           // e_xy+e_x > 0
		{
			if (x0 == x1)
			{
				break;
			}
			err += dy;
			x0 += sx;
		}
		if (e2 <= dx)                                           // e_xy+e_y < 0
		{
			if (y0 == y1)
			{
				break;
			}
			err += dx;
			y0 += sy;
		}
	}
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotLineAA(Mat& img, int x0, int y0, int x1, int y1, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// draw a black (0) anti-aliased line on white (255) background
	int sx = x0 < x1 ? 1 : -1, sy = y0 < y1 ? 1 : -1, x2;
	long dx = abs(x1-x0), dy = abs(y1-y0), err = dx*dx+dy*dy;
	long e2 = err == 0 ? 1 : 0xffff7fl/sqrt(err);     // multiplication factor
	double a=0;
	dx *= e2;
	dy *= e2;
	err = dx-dy;                       // error value e_xy

	bool dash=1;

	for ( ; ; )                                                  // pixel loop
	{
		if(pattern!=NULL)
		{
			dash=pattern->get_dash();
		}

		a=double(abs(err-dx+dy)>>16)/255.0;
		a=1.0-((1.0-a)*(1.0-alpha));

		if(dash)
			setPixel(img, x0,y0,color,a);

		e2 = err;
		x2 = x0;
		if (2*e2 >= -dx)                                              // x step
		{
			if (x0 == x1)
			{
				break;
			}
			if (e2+dy < 0xff0000l)
			{
				a=double((e2+dy)>>16)/255.0;
				a=1.0-((1.0-a)*(1.0-alpha));
				if(dash)
					setPixel(img, x0,y0+sy,color, a);
			}
			err -= dy;
			x0 += sx;
		}
		if (2*e2 <= dy)                                               // y step
		{
			if (y0 == y1)
			{
				break;
			}
			if (dx-e2 < 0xff0000l)
			{
				a=double((dx-e2)>>16)/255.0;
				a=1.0-((1.0-a)*(1.0-alpha));
				if(dash)
					setPixel(img, x2+sx,y0,color, a);
			}
			err += dx;
			y0 += sy;
		}

		if(pattern!=NULL)
		{
			pattern->advance();
		}

	}
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotLineWidth(Mat& img, double x0,double  y0,double x1,double y1,double th, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	double dx = abs(x1-x0), sx = x0 < x1 ? 1 : -1;
	double dy = abs(y1-y0), sy = y0 < y1 ? 1 : -1;
	double err, e2 = sqrt(dx*dx+dy*dy);                            // length 
	double acol;
	bool dash=1;
	if (th <= 1 || e2 == 0)
	{
		return plotLineAA(img,x0,y0, x1,y1,color,alpha);    // assert 
	}
	dx *= 255/e2;
	dy *= 255/e2;
	th = 255*(th-1);               // scale values 
	if (dx < dy)                                                 // steep line 
	{
		x1 = floor((e2+th/2)/dy);                          // start offset 
		err = x1*dy-th/2;                  // shift error value to offset width 
		err*=2; // Этого нет в оригинале , но без него хреново
		for (x0 -= x1*sx; ; y0 += sy)
		{
			acol=double( err )/255.0;
			acol=1.0-((1.0-acol)*(1.0-alpha));
			if(pattern!=NULL)
			{
				dash=pattern->get_dash();
				pattern->advance();
			}
			if(dash)
			{
				setPixel(img,x1 = x0, y0,color, acol);                  // aliasing pre-pixel 
			}
			for (e2 = dy-err-th; e2+dy < 255; e2 += dy)
			{
				if(dash)
				{
					setPixel(img,x1 += sx, y0,color,alpha);    // pixel on the line
				}
			}
			acol=double( e2 )/255.0;
			acol=1.0-((1.0-acol)*(1.0-alpha));
			if(dash)
			{
				setPixel(img,x1+sx, y0,color,acol);                    // aliasing post-pixel 
			}
			if (y0 == y1)
			{
				break;
			}
			err += dx;                                                 // y-step 
			if (err > 255)
			{
				err -= dy;    // x-step 
				x0 += sx;
			}
		}
	}
	else                                                          // flat line 
	{
		y1 = floor((e2+th/2)/dx);                          // start offset 
		err = y1*dx-th/2;                  // shift error value to offset width 
		err*=2;		// Этого нет в оригинале , но без него хреново

		for (y0 -= y1*sy; ; x0 += sx)
		{
			if(pattern!=NULL)
			{
				dash=pattern->get_dash();
				pattern->advance();
			}

			acol=double( err )/255.0;
			acol=1.0-((1.0-acol)*(1.0-alpha));
			y1 = y0;
			if(acol<=1)
			{
				if(dash)
				{
					setPixel(img,x0, y1,color, acol);                  // aliasing pre-pixel 
				}
			}

			for (e2 = dx-err-th; e2+dx < 255; e2 += dx)
			{
				if(dash)
				{
					setPixel(img,x0, y1 += sy,color,alpha);    // pixel on the line
				}
			}
			acol=double( e2 )/255.0;
			acol=1.0-((1.0-acol)*(1.0-alpha));
			if(dash)
			{
				setPixel(img,x0, y1+sy, color,acol);                    // aliasing post-pixel 
			}
			if (x0 == x1)
			{
				break;
			}
			err += dy;                                                 // x-step 
			if (err > 255)
			{
				err -= dx;    // y-step 
				y0 += sy;
			}
		}
	}
}


//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotEllipse(Mat& img, int xm, int ym, int a, int b, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	int x = -a, y = 0;           // II. quadrant from bottom left to top right
	long e2 = (long)b*b, err = (long)x*(2*e2+x)+e2;         // error of 1.step
	bool dash=1;
	do
	{
		setPixel(img,xm-x, ym+y,color,alpha);                                 //   I. Quadrant
		setPixel(img,xm+x, ym+y,color,alpha);                                 //  II. Quadrant
		setPixel(img,xm+x, ym-y,color,alpha);                                 // III. Quadrant
		setPixel(img,xm-x, ym-y,color,alpha);                                 //  IV. Quadrant
		e2 = 2*err;
		if (e2 >= (x*2+1)*(long)b*b)                           // e_xy+e_x > 0
		{
			err += (++x*2+1)*(long)b*b;
		}
		if (e2 <= (y*2+1)*(long)a*a)                           // e_xy+e_y < 0
		{
			err += (++y*2+1)*(long)a*a;
		}
	}
	while (x <= 0);
	while (y++ < b)                    // too early stop of flat ellipses a=1,
	{
		setPixel(img,xm, ym+y,color, alpha);                        // -> finish tip of ellipse
		setPixel(img,xm, ym-y,color, alpha);
	}
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotOptimizedEllipse(Mat& img, int xm, int ym, int a, int b, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	long x = -a, y = 0;          // II. quadrant from bottom left to top right
	long e2 = b, dx = (1+2*x)*e2*e2;                       // error increment
	long dy = x*x, err = dx+dy;                             // error of 1.step
	do
	{
		setPixel(img,xm-x, ym+y,color,alpha);                                 //   I. Quadrant
		setPixel(img,xm+x, ym+y,color,alpha);                                 //  II. Quadrant
		setPixel(img,xm+x, ym-y,color,alpha);                                 // III. Quadrant
		setPixel(img,xm-x, ym-y,color,alpha);                                 //  IV. Quadrant
		e2 = 2*err;
		if (e2 >= dx)
		{
			x++;    // x step
			err += dx += 2*(long)b*b;
		}
		if (e2 <= dy)
		{
			y++;    // y step
			err += dy += 2*(long)a*a;
		}
	}
	while (x <= 0);
	while (y++ < b)              // too early stop for flat ellipses with a=1,
	{
		setPixel(img,xm, ym+y,color,alpha);                        // -> finish tip of ellipse
		setPixel(img,xm, ym-y,color,alpha);
	}
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotCircle(Mat& img, int xm, int ym, int r, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	int x = -r, y = 0, err = 2-2*r;                // bottom left to top right
	do
	{
		setPixel(img, xm-x, ym+y, color,alpha);                            //   I. Quadrant +x +y
		setPixel(img, xm-y, ym-x, color,alpha);                            //  II. Quadrant -x +y
		setPixel(img, xm+x, ym-y, color,alpha);                            // III. Quadrant -x -y
		setPixel(img, xm+y, ym+x, color,alpha);                            //  IV. Quadrant +x -y
		r = err;
		if (r <= y)
		{
			err += ++y*2+1;    // e_xy+e_y < 0
		}
		if (r > x || err > y)                  // e_xy+e_x > 0 or no 2nd y-step
		{
			err += ++x*2+1;    // -> x-step now
		}
	}
	while (x < 0);
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotEllipseRect(Mat& img,int x0, int y0, int x1, int y1, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// rectangular parameter enclosing the ellipse
	long a = abs(x1-x0), b = abs(y1-y0), b1 = b&1;                 // diameter
	double dx = 4*(1.0-a)*b*b, dy = 4*(b1+1)*a*a;           // error increment
	double err = dx+dy+b1*a*a, e2;                          // error of 1.step
	if (x0 > x1)
	{
		x0 = x1;    // if called with swapped points
		x1 += a;
	}
	if (y0 > y1)
	{
		y0 = y1;    // .. exchange them
	}
	y0 += (b+1)/2;
	y1 = y0-b1;                               // starting pixel
	a = 8*a*a;
	b1 = 8*b*b;
	do
	{
		setPixel(img,x1, y0,color,alpha);                                      //   I. Quadrant
		setPixel(img,x0, y0,color,alpha);                                      //  II. Quadrant
		setPixel(img,x0, y1,color,alpha);                                      // III. Quadrant
		setPixel(img,x1, y1,color,alpha);                                      //  IV. Quadrant
		e2 = 2*err;
		if (e2 <= dy)
		{
			y0++;    // y step
			y1--;
			err += dy += a;
		}
		if (e2 >= dx || 2*err > dy)
		{
			x0++;    // x step
			x1--;
			err += dx += b1;
		}
	}
	while (x0 <= x1);
	while (y0-y1 <= b)                  // too early stop of flat ellipses a=1
	{
		setPixel(img, x0-1, y0	,color,alpha);                         // -> finish tip of ellipse
		setPixel(img, x1+1, y0++,color,alpha);
		setPixel(img, x0-1, y1	,color,alpha);
		setPixel(img, x1+1, y1--,color,alpha);
	}
}


void cvGDI::plotEllipseRectWidth(Mat& img, double  x0,double  y0,double  x1,double  y1,double th, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// draw anti-aliased ellipse inside rectangle with thick line 
	long a = fabs(x1-x0), b = fabs(y1-y0), b1 = b&1;  // outer diameter 
	double a2 = a-2*th, b2 = b-2*th;                            // inner diameter 
	double dx = 4*(a-1)*b*b, dy = 4*(b1-1)*a*a;                // error increment 
	double i = a+b2, err = b1*a*a, dx2, dy2, e2, ed;
	double acol;
	// thick line correction 
	if (th < 1.5)
	{
		return plotEllipseRectAA(img, x0,y0, x1,y1, color, alpha);
	}
	if ((th-1)*(2*b-th) > a*a)
	{
		b2 = sqrt(a*(b-a)*i*a2)/(a-th);
	}
	if ((th-1)*(2*a-th) > b*b)
	{
		a2 = sqrt(b*(a-b)*i*b2)/(b-th);
		th = (a-a2)/2;
	}
	if (a == 0 || b == 0)
	{
		return plotLineWidth(img,x0,y0, x1,y1,th, color, alpha);
	}
	if (x0 > x1)
	{
		x0 = x1;    // if called with swapped points 
		x1 += a;
	}
	if (y0 > y1)
	{
		y0 = y1;    // .. exchange them 
	}
	if (b2 <= 0)
	{
		th = a;    // filled ellipse 
	}
	e2 = th-floor(th);
	th = x0+th-e2;
	dx2 = 4*(a2+2*e2-1)*b2*b2;
	dy2 = 4*(b1-1)*a2*a2;
	e2 = dx2*e2;
	y0 += (b+1)>>1;
	y1 = y0-b1;                              // starting pixel 
	a = 8*a*a;
	b1 = 8*b*b;
	a2 = 8*a2*a2;
	b2 = 8*b2*b2;
	do
	{
		for (;;)
		{
			if (err < 0 || x0 > x1)
			{
				i = x0;
				break;
			}
			i = min(dx,dy);
			ed = max(dx,dy);
			if (y0 == y1+1 && 2*err > dx && a > b1)
			{
				ed = a/4;    // x-tip 
			}
			else
			{
				ed += 2*ed*i*i/(4*ed*ed+i*i+1)+1;    // approx ed=sqrt(dx*dx+dy*dy) 
			}
			//i = 255*err/ed;                             // outside anti-aliasing 
			acol=double( 255*err/ed )/255.0;
			acol=1.0-((1.0-acol)*(1.0-alpha));
			setPixel(img, x0,y0, color,acol);
			setPixel(img, x0,y1, color,acol);
			setPixel(img, x1,y0, color,acol);
			setPixel(img, x1,y1, color,acol);
			if (err+dy+a < dx)
			{
				i = x0+1;
				break;
			}
			x0++;
			x1--;
			err -= dx;
			dx -= b1;                // x error increment 
		}
		for (; i < th && 2*i <= x0+x1; i++)                  // fill line pixel 
		{
			setPixel(img, i,y0,color,alpha);
			setPixel(img, x0+x1-i,y0,color,alpha);
			setPixel(img, i,y1,color,alpha);
			setPixel(img, x0+x1-i,y1,color,alpha);
		}
		while (e2 > 0 && x0+x1 >= 2*th)                 // inside anti-aliasing 
		{
			i = min(dx2,dy2);
			ed = max(dx2,dy2);
			if (y0 == y1+1 && 2*e2 > dx2 && a2 > b2)
			{
				ed = a2/4;    // x-tip 
			}
			else
			{
				ed += 2*ed*i*i/(4*ed*ed+i*i);    // approximation 
			}
			//i = 255-255*e2/ed;             // get intensity value by pixel error 
			acol=double( 255-255*e2/ed )/255.0;
			acol=1.0-((1.0-acol)*(1.0-alpha));
			setPixel(img, th,y0,color,acol);
			setPixel(img,x0+x1-th,y0,color, acol);
			setPixel(img, th,y1,color,acol);
			setPixel(img,x0+x1-th,y1,color, acol);
			if (e2+dy2+a2 < dx2)
			{
				break;
			}
			th++;
			e2 -= dx2;
			dx2 -= b2;                     // x error increment 
		}
		e2 += dy2 += a2;
		y0++;
		y1--;
		err += dy += a;                                   // y step 
	}
	while (x0 < x1);
	if (y0-y1 <= b)
	{
		if (err > dy+a)
		{
			y0--;
			y1++;
			err -= dy -= a;
		}
		for (; y0-y1 <= b; err += dy += a)   // too early stop of flat ellipses 
		{
			//i = 255*4*err/b1;                        // -> finish tip of ellipse 
			acol=double( 255*4*err/b1 )/255.0;
			acol=1.0-((1.0-acol)*(1.0-alpha));
			setPixel(img,x0,y0, color,acol);
			setPixel(img,x1,y0++, color,acol);
			setPixel(img,x0,y1, color,acol);
			setPixel(img,x1,y1--, color,acol);
		}
	}
}

//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotQuadBezierSeg(Mat& img, int x0, int y0, int x1, int y1, int x2, int y2, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot a limited quadratic Bezier segment
	int sx = x2-x1, sy = y2-y1;
	long xx = x0-x1, yy = y0-y1, xy;              // relative values for checks
	double dx, dy, err, cur = xx*sy-yy*sx;                         // curvature
	assert(xx*sx <= 0 && yy*sy <= 0);       // sign of gradient must not change
	if (sx*(long)sx+sy*(long)sy > xx*xx+yy*yy)        // begin with longer part
	{
		x2 = x0;
		x0 = sx+x1;
		y2 = y0;
		y0 = sy+y1;
		cur = -cur;       // swap P0 P2
	}
	if (cur != 0)                                           // no straight line
	{
		xx += sx;
		xx *= sx = x0 < x2 ? 1 : -1;                // x step direction
		yy += sy;
		yy *= sy = y0 < y2 ? 1 : -1;                // y step direction
		xy = 2*xx*yy;
		xx *= xx;
		yy *= yy;               // differences 2nd degree
		if (cur*sx*sy < 0)                                  // negated curvature?
		{
			xx = -xx;
			yy = -yy;
			xy = -xy;
			cur = -cur;
		}
		dx = 4.0*sy*cur*(x1-x0)+xx-xy;                  // differences 1st degree
		dy = 4.0*sx*cur*(y0-y1)+yy-xy;
		xx += xx;
		yy += yy;
		err = dx+dy+xy;                     // error 1st step
		do
		{
			setPixel(img, x0,y0, color, alpha);                                          // plot curve
			if (x0 == x2 && y0 == y2)
			{
				return;    // last pixel -> curve finished
			}
			y1 = 2*err < dx;                       // save value for test of y step
			if (2*err > dy)
			{
				x0 += sx;    // x step
				dx -= xy;
				err += dy += yy;
			}
			if (    y1    )
			{
				y0 += sy;    // y step
				dy -= xy;
				err += dx += xx;
			}
		}
		while (dy < 0 && dx > 0);          // gradient negates -> algorithm fails
	}
	plotLine(img, x0,y0, x2,y2, color,alpha);                       // plot remaining part to end
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotQuadBezier(Mat& img, int x0, int y0, int x1, int y1, int x2, int y2, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot any quadratic Bezier curve
	int x = x0-x1, y = y0-y1;
	double t = x0-2*x1+x2, r;
	if ((long)x*(x2-x1) > 0)                          // horizontal cut at P4?
	{
		if ((long)y*(y2-y1) > 0)                     // vertical cut at P6 too?
			if (fabs((y0-2*y1+y2)/t*x) > abs(y))                 // which first?
			{
				x0 = x2;
				x2 = x+x1;
				y0 = y2;
				y2 = y+y1;            // swap points
			}                            // now horizontal cut at P4 comes first
			t = (x0-x1)/t;
			r = (1-t)*((1-t)*y0+2.0*t*y1)+t*t*y2;                       // By(t=P4)
			t = (x0*x2-x1*x1)*t/(x0-x1);                       // gradient dP4/dx=0
			x = floor(t+0.5);
			y = floor(r+0.5);
			r = (y1-y0)*(t-x0)/(x1-x0)+y0;                  // intersect P3 | P0 P1
			plotQuadBezierSeg(img, x0,y0, x,floor(r+0.5), x,y,color,alpha);
			r = (y1-y2)*(t-x2)/(x1-x2)+y2;                  // intersect P4 | P1 P2
			x0 = x1 = x;
			y0 = y;
			y1 = floor(r+0.5);             // P0 = P4, P1 = P8
	}
	if ((long)(y0-y1)*(y2-y1) > 0)                      // vertical cut at P6?
	{
		t = y0-2*y1+y2;
		t = (y0-y1)/t;
		r = (1-t)*((1-t)*x0+2.0*t*x1)+t*t*x2;                       // Bx(t=P6)
		t = (y0*y2-y1*y1)*t/(y0-y1);                       // gradient dP6/dy=0
		x = floor(r+0.5);
		y = floor(t+0.5);
		r = (x1-x0)*(t-y0)/(y1-y0)+x0;                  // intersect P6 | P0 P1
		plotQuadBezierSeg(img, x0,y0, floor(r+0.5),y, x,y, color,alpha);
		r = (x1-x2)*(t-y2)/(y1-y2)+x2;                  // intersect P7 | P1 P2
		x0 = x;
		x1 = floor(r+0.5);
		y0 = y1 = y;             // P0 = P6, P1 = P7
	}
	plotQuadBezierSeg(img, x0,y0, x1,y1, x2,y2, color, alpha);                  // remaining part
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotQuadBezierAA(Mat& img, int x0, int y0, int x1, int y1, int x2, int y2, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot any quadratic Bezier curve
	int x = x0-x1, y = y0-y1;
	double t = x0-2*x1+x2, r;
	if ((long)x*(x2-x1) > 0)                          // horizontal cut at P4?
	{
		if ((long)y*(y2-y1) > 0)                     // vertical cut at P6 too?
			if (fabs((y0-2*y1+y2)/t*x) > abs(y))                 // which first?
			{
				x0 = x2;
				x2 = x+x1;
				y0 = y2;
				y2 = y+y1;            // swap points
			}                            // now horizontal cut at P4 comes first
			t = (x0-x1)/t;
			r = (1-t)*((1-t)*y0+2.0*t*y1)+t*t*y2;                       // By(t=P4)
			t = (x0*x2-x1*x1)*t/(x0-x1);                       // gradient dP4/dx=0
			x = floor(t+0.5);
			y = floor(r+0.5);
			r = (y1-y0)*(t-x0)/(x1-x0)+y0;                  // intersect P3 | P0 P1
			plotQuadBezierSegAA(img, x0,y0, x,floor(r+0.5), x,y,color,alpha);
			r = (y1-y2)*(t-x2)/(x1-x2)+y2;                  // intersect P4 | P1 P2
			x0 = x1 = x;
			y0 = y;
			y1 = floor(r+0.5);             // P0 = P4, P1 = P8
	}
	if ((long)(y0-y1)*(y2-y1) > 0)                      // vertical cut at P6?
	{
		t = y0-2*y1+y2;
		t = (y0-y1)/t;
		r = (1-t)*((1-t)*x0+2.0*t*x1)+t*t*x2;                       // Bx(t=P6)
		t = (y0*y2-y1*y1)*t/(y0-y1);                       // gradient dP6/dy=0
		x = floor(r+0.5);
		y = floor(t+0.5);
		r = (x1-x0)*(t-y0)/(y1-y0)+x0;                  // intersect P6 | P0 P1
		plotQuadBezierSegAA(img, x0,y0, floor(r+0.5),y, x,y, color,alpha);
		r = (x1-x2)*(t-y2)/(y1-y2)+x2;                  // intersect P7 | P1 P2
		x0 = x;
		x1 = floor(r+0.5);
		y0 = y1 = y;             // P0 = P6, P1 = P7
	}
	plotQuadBezierSegAA(img, x0,y0, x1,y1, x2,y2, color, alpha);                  // remaining part
}

//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotQuadRationalBezierSeg(Mat & img, int x0, int y0, int x1, int y1, int x2, int y2, double w, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot a limited rational Bezier segment, squared weight
	int sx = x2-x1, sy = y2-y1;                   // relative values for checks
	double dx = x0-x2, dy = y0-y2, xx = x0-x1, yy = y0-y1;
	double xy = xx*sy+yy*sx, cur = xx*sy-yy*sx, err;               // curvature
	assert(xx*sx <= 0.0 && yy*sy <= 0.0);   // sign of gradient must not change
	if (cur != 0.0 && w > 0.0)                              // no straight line
	{
		if (sx*(long)sx+sy*(long)sy > xx*xx+yy*yy)      // begin with longer part
		{
			x2 = x0;
			x0 -= dx;
			y2 = y0;
			y0 -= dy;
			cur = -cur;         // swap P0 P2
		}
		xx = 2.0*(4.0*w*sx*xx+dx*dx);                   // differences 2nd degree
		yy = 2.0*(4.0*w*sy*yy+dy*dy);
		sx = x0 < x2 ? 1 : -1;                                // x step direction
		sy = y0 < y2 ? 1 : -1;                                // y step direction
		xy = -2.0*sx*sy*(2.0*w*xy+dx*dy);
		if (cur*sx*sy < 0.0)                                // negated curvature?
		{
			xx = -xx;
			yy = -yy;
			xy = -xy;
			cur = -cur;
		}
		dx = 4.0*w*(x1-x0)*sy*cur+xx/2.0+xy;            // differences 1st degree
		dy = 4.0*w*(y0-y1)*sx*cur+yy/2.0+xy;
		if (w < 0.5 && (dy > xy || dx < xy))     // flat ellipse, algorithm fails
		{
			cur = (w+1.0)/2.0;
			w = sqrt(w);
			xy = 1.0/(w+1.0);
			sx = floor((x0+2.0*w*x1+x2)*xy/2.0+0.5);    // subdivide curve in half
			sy = floor((y0+2.0*w*y1+y2)*xy/2.0+0.5);
			dx = floor((w*x1+x0)*xy+0.5);
			dy = floor((y1*w+y0)*xy+0.5);
			plotQuadRationalBezierSeg(img,x0,y0, dx,dy, sx,sy, cur,color,alpha);// plot separately
			dx = floor((w*x1+x2)*xy+0.5);
			dy = floor((y1*w+y2)*xy+0.5);
			plotQuadRationalBezierSeg(img,sx,sy, dx,dy, x2,y2, cur,color,alpha);
			return;
		}
		err = dx+dy-xy;                                           // error 1.step
		do
		{
			setPixel(img, x0,y0, color, alpha);                                          // plot curve
			if (x0 == x2 && y0 == y2)
			{
				return;    // last pixel -> curve finished
			}
			x1 = 2*err > dy;
			y1 = 2*(err+yy) < -dy;// save value for test of x step
			if (2*err < dx || y1)
			{
				y0 += sy;    // y step
				dy += xy;
				err += dx += xx;
			}
			if (2*err > dx || x1)
			{
				x0 += sx;    // x step
				dx += xy;
				err += dy += yy;
			}
		}
		while (dy <= xy && dx >= xy);      // gradient negates -> algorithm fails
	}
	plotLine(img, x0,y0, x2,y2, color, alpha);                     // plot remaining needle to end
}
//------------------------------------------------------------------------------------------------
// NOT WORKING
//-----------------------------------------------------------------------------------------------
void cvGDI::plotQuadRationalBezierWidthSeg(Mat& img, double x0,double y0,double x1,double y1,double x2,double y2,double w,double th, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot a limited rational Bezier segment of thickness th, squared weight 
	double sx = x2-x1, sy = y2-y1;                  // relative values for checks 
	double dx = x0-x2, dy = y0-y2, xx = x0-x1, yy = y0-y1;
	double xy = xx*sy+yy*sx, cur = xx*sy-yy*sx, err, e2, ed;         // curvature 
	double acol;
	assert(xx*sx <= 0.0 && yy*sy <= 0.0);  // sign of gradient must not change 
	if (cur != 0.0 && w > 0.0)                             // no straight line 
	{
		if (sx*sx+sy*sy > xx*xx+yy*yy)                // begin with longer part 
		{
			x2 = x0;
			x0 -= dx;
			y2 = y0;
			y0 -= dy;
			cur = -cur;      // swap P0 P2 
		}
		xx = 2.0*(4.0*w*sx*xx+dx*dx);                 // differences 2nd degree 
		yy = 2.0*(4.0*w*sy*yy+dy*dy);
		sx = x0 < x2 ? 1 : -1;                              // x step direction 
		sy = y0 < y2 ? 1 : -1;                              // y step direction 
		xy = -2.0*sx*sy*(2.0*w*xy+dx*dy);
		if (cur*sx*sy < 0)                                // negated curvature? 
		{
			xx = -xx;
			yy = -yy;
			cur = -cur;
			xy = -xy;
		}
		dx = 4.0*w*(x1-x0)*sy*cur+xx/2.0;             // differences 1st degree 
		dy = 4.0*w*(y0-y1)*sx*cur+yy/2.0;
		if (w < 0.5 && (dx+xx <= 0 || dy+yy >= 0))  // flat ellipse, algo fails 
		{
			cur = (w+1.0)/2.0;
			w = sqrt(w);
			xy = 1.0/(w+1.0);
			sx = floor((x0+2.0*w*x1+x2)*xy/2.0+0.5);    // subdivide curve  
			sy = floor((y0+2.0*w*y1+y2)*xy/2.0+0.5);     // plot separately 
			dx = floor((w*x1+x0)*xy+0.5);
			dy = floor((y1*w+y0)*xy+0.5);
			plotQuadRationalBezierWidthSeg(img, x0,y0, dx,dy, sx,sy, cur, th,color, alpha);
			dx = floor((w*x1+x2)*xy+0.5);
			dy = floor((y1*w+y2)*xy+0.5);
			return plotQuadRationalBezierWidthSeg(img, sx,sy, dx,dy, x2,y2, cur, th, color, alpha);
		}
fail:
		for (err = 0; dy+2*yy < 0 && dx+2*xx > 0; ) // loop of steep/flat curve 
			if (dx+dy+xy < 0)                                     // steep curve 
			{
				do
				{
					ed = -dy-2*dy*dx*dx/(4.*dy*dy+dx*dx);      // approximate sqrt 
					w = (th-1)*ed;                             // scale line width 
					x1 = floor((err-ed-w/2)/dy);              // start offset 
					e2 = err-x1*dy-w/2;                   // error value at offset 
					x1 = x0-x1*sx;                                  // start point 
					acol=double( 255*e2/ed )/255.0;
					acol=1.0-((1.0-acol)*(1.0-alpha));
					setPixel(img, x1, y0, color,acol);           // aliasing pre-pixel 
					for (e2 = -w-dy-e2; e2-dy < ed; e2 -= dy)
					{
						setPixel(img, x1 += sx, y0, color, alpha);    // pixel on thick line 
					}
					acol=double( 255*e2/ed )/255.0;
					acol=1.0-((1.0-acol)*(1.0-alpha));
					setPixel(img, x1+sx, y0, color, acol);       // aliasing post-pixel 
					if (y0 == y2)
					{
						return;    // last pixel -> curve finished 
					}
					y0 += sy;
					dy += xy;
					err += dx;
					dx += xx;             // y step 
					if (2*err+dy > 0)                              // e_x+e_xy > 0 
					{
						x0 += sx;
						dx += xy;
						err += dy;
						dy += yy;          // x step 
					}
					if (x0 != x2 && (dx+2*xx <= 0 || dy+2*yy >= 0))
						if (fabs(y2-y0) > fabs(x2-x0))
						{
							break;
						}
						else
						{
							break;    // other curve near 
						}
				}
				while (dx+dy+xy < 0);                    // gradient still steep? 
				// change from steep to flat curve 
				for (cur = err-dy-w/2, y1 = y0; cur < ed; y1 += sy, cur += dx)
				{
					for (e2 = cur, x1 = x0; e2-dy < ed; e2 -= dy)
					{
						setPixel(img, x1 -= sx, y1,color, alpha);    // pixel on thick line 
					}
					acol=double( 255.0*e2/ed )/255.0;
					acol=1.0-((1.0-acol)*(1.0-alpha));
					setPixel(img,x1-sx, y1, color, acol);       // aliasing post-pixel 
				}
			}
			else                                                   // flat curve 
			{
				do
				{
					ed = dx+2*dx*dy*dy/(4.*dx*dx+dy*dy);       // approximate sqrt 
					w = (th-1)*ed;                             // scale line width 
					y1 = floor((err+ed+w/2)/dx);              // start offset 
					e2 = y1*dx-w/2-err;                   // error value at offset 
					y1 = y0-y1*sy;                                  // start point 
					acol=double( 255*e2/ed )/255.0;
					acol=1.0-((1.0-acol)*(1.0-alpha));
					setPixel(img,x0, y1, color,acol);           // aliasing pre-pixel 
					for (e2 = dx-e2-w; e2+dx < ed; e2 += dx)
					{
						setPixel(img, x0, y1 += sy, color, alpha);    // pixel on thick line 
					}
					acol=double( 255*e2/ed )/255.0;
					acol=1.0-((1.0-acol)*(1.0-alpha));
					setPixel(img, x0, y1+sy, color, acol);       // aliasing post-pixel 
					if (x0 == x2)
					{
						return;    // last pixel -> curve finished 
					}
					x0 += sx;
					dx += xy;
					err += dy;
					dy += yy;             // x step 
					if (2*err+dx < 0)                              // e_y+e_xy < 0 
					{
						y0 += sy;
						dy += xy;
						err += dx;
						dx += xx;          // y step 
					}
					if (y0 != y2 && (dx+2*xx <= 0 || dy+2*yy >= 0))
						if (fabs(y2-y0) <= fabs(x2-x0))
						{
							break;
						}
						else
						{
							break;    // other curve near 
						}
				}
				while (dx+dy+xy > 0);                    // gradient still flat? 
				// change from flat to steep curve 
				for (cur = -err+dx-w/2.0, x1 = x0; cur < ed; x1 += sx, cur -= dy)
				{
					for (e2 = cur, y1 = y0; e2+dx < ed; e2 += dx)
					{
						setPixel(img, x1, y1 -= sy,color, alpha);    // pixel on thick line 
					}
					acol=double( 255*e2/ed )/255.0;
					acol=1.0-((1.0-acol)*(1.0-alpha));
					setPixel(img,x1, y1-sy, color,acol);       // aliasing post-pixel 
				}
			}
	}
	plotLineWidth(img,x0,y0, x2,y2, th,color, alpha);            // confusing error values  
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotQuadRationalBezier(Mat& img, int x0, int y0, int x1, int y1, int x2, int y2, double w, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot any quadratic rational Bezier curve
	int x = x0-2*x1+x2, y = y0-2*y1+y2;
	double xx = x0-x1, yy = y0-y1, ww, t, q;
	assert(w >= 0.0);
	if (xx*(x2-x1) > 0)                               // horizontal cut at P4?
	{
		if (yy*(y2-y1) > 0)                          // vertical cut at P6 too?
			if (fabs(xx*y) > fabs(yy*x))                         // which first?
			{
				x0 = x2;
				x2 = xx+x1;
				y0 = y2;
				y2 = yy+y1;          // swap points
			}                            // now horizontal cut at P4 comes first
			if (x0 == x2 || w == 1.0)
			{
				t = (x0-x1)/(double)x;
			}
			else                                   // non-rational or rational case
			{
				q = sqrt(4.0*w*w*(x0-x1)*(x2-x1)+(x2-x0)*(long)(x2-x0));
				if (x1 < x0)
				{
					q = -q;
				}
				t = (2.0*w*(x0-x1)-x0+x2+q)/(2.0*(1.0-w)*(x2-x0));        // t at P4
			}
			q = 1.0/(2.0*t*(1.0-t)*(w-1.0)+1.0);                 // sub-divide at t
			xx = (t*t*(x0-2.0*w*x1+x2)+2.0*t*(w*x1-x0)+x0)*q;               // = P4
			yy = (t*t*(y0-2.0*w*y1+y2)+2.0*t*(w*y1-y0)+y0)*q;
			ww = t*(w-1.0)+1.0;
			ww *= ww*q;                    // squared weight P3
			w = ((1.0-t)*(w-1.0)+1.0)*sqrt(q);                         // weight P8
			x = floor(xx+0.5);
			y = floor(yy+0.5);                             // P4
			yy = (xx-x0)*(y1-y0)/(x1-x0)+y0;                // intersect P3 | P0 P1
			plotQuadRationalBezierSeg(img, x0,y0, x,floor(yy+0.5), x,y, ww, color, alpha);
			yy = (xx-x2)*(y1-y2)/(x1-x2)+y2;                // intersect P4 | P1 P2
			y1 = floor(yy+0.5);
			x0 = x1 = x;
			y0 = y;            // P0 = P4, P1 = P8
	}
	if ((y0-y1)*(long)(y2-y1) > 0)                      // vertical cut at P6?
	{
		if (y0 == y2 || w == 1.0)
		{
			t = (y0-y1)/(y0-2.0*y1+y2);
		}
		else                                   // non-rational or rational case
		{
			q = sqrt(4.0*w*w*(y0-y1)*(y2-y1)+(y2-y0)*(long)(y2-y0));
			if (y1 < y0)
			{
				q = -q;
			}
			t = (2.0*w*(y0-y1)-y0+y2+q)/(2.0*(1.0-w)*(y2-y0));        // t at P6
		}
		q = 1.0/(2.0*t*(1.0-t)*(w-1.0)+1.0);                 // sub-divide at t
		xx = (t*t*(x0-2.0*w*x1+x2)+2.0*t*(w*x1-x0)+x0)*q;               // = P6
		yy = (t*t*(y0-2.0*w*y1+y2)+2.0*t*(w*y1-y0)+y0)*q;
		ww = t*(w-1.0)+1.0;
		ww *= ww*q;                    // squared weight P5
		w = ((1.0-t)*(w-1.0)+1.0)*sqrt(q);                         // weight P7
		x = floor(xx+0.5);
		y = floor(yy+0.5);                             // P6
		xx = (x1-x0)*(yy-y0)/(y1-y0)+x0;                // intersect P6 | P0 P1
		plotQuadRationalBezierSeg(img, x0,y0, floor(xx+0.5),y, x,y, ww, color, alpha);
		xx = (x1-x2)*(yy-y2)/(y1-y2)+x2;                // intersect P7 | P1 P2
		x1 = floor(xx+0.5);
		x0 = x;
		y0 = y1 = y;            // P0 = P6, P1 = P7
	}
	plotQuadRationalBezierSeg(img, x0,y0, x1,y1, x2,y2, w*w, color, alpha);          // remaining
}

//------------------------------------------------------------------------------------------------
//
//----------------------------------------------------------------------------------------------
void cvGDI::plotQuadRationalBezierWidth(Mat& img, double x0,double  y0,double  x1,double  y1,double  x2,double  y2,double  w,double  th, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot any anti-aliased quadratic rational Bezier curve 
	double x = x0-2*x1+x2, y = y0-2*y1+y2;
	double xx = x0-x1, yy = y0-y1, ww, t, q;
	assert(w >= 0.0);
	if (xx*(x2-x1) > 0)                               // horizontal cut at P4? 
	{
		if (yy*(y2-y1) > 0)                          // vertical cut at P6 too? 
			if (fabs(xx*y) > fabs(yy*x))                 // which first? 
			{
				x0 = x2;
				x2 = xx+x1;
				y0 = y2;
				y2 = yy+y1;          // swap points 
			}                            // now horizontal cut at P4 comes first 
			if (x0 == x2 || w == 1.0)
			{
				t = (x0-x1)/x;
			}
			else                                   // non-rational or rational case 
			{
				q = sqrt(4.0*w*w*(x0-x1)*(x2-x1)+(x2-x0)*(x2-x0));
				if (x1 < x0)
				{
					q = -q;
				}
				t = (2.0*w*(x0-x1)-x0+x2+q)/(2.0*(1.0-w)*(x2-x0));        // t at P4 
			}
			q = 1.0/(2.0*t*(1.0-t)*(w-1.0)+1.0);                 // sub-divide at t 
			xx = (t*t*(x0-2.0*w*x1+x2)+2.0*t*(w*x1-x0)+x0)*q;               // = P4 
			yy = (t*t*(y0-2.0*w*y1+y2)+2.0*t*(w*y1-y0)+y0)*q;
			ww = t*(w-1.0)+1.0;
			ww *= ww*q;                    // squared weight P3 
			w = ((1.0-t)*(w-1.0)+1.0)*sqrt(q);                    // weight P8 
			x = floor(xx+0.5);
			y =floor(yy+0.5);                   // P4 
			yy = (xx-x0)*(y1-y0)/(x1-x0)+y0;                // intersect P3 | P0 P1 
			plotQuadRationalBezierWidthSeg(img, x0,y0, x,floor(yy+0.5), x,y, ww, th, color,alpha);
			yy = (xx-x2)*(y1-y2)/(x1-x2)+y2;                // intersect P4 | P1 P2 
			y1 = floor(yy+0.5);
			x0 = x1 = x;
			y0 = y;       // P0 = P4, P1 = P8 
	}
	if ((y0-y1)*(y2-y1) > 0)                            // vertical cut at P6? 
	{
		if (y0 == y2 || w == 1.0)
		{
			t = (y0-y1)/(y0-2.0*y1+y2);
		}
		else                                   // non-rational or rational case 
		{
			q = sqrt(4.0*w*w*(y0-y1)*(y2-y1)+(y2-y0)*(y2-y0));
			if (y1 < y0)
			{
				q = -q;
			}
			t = (2.0*w*(y0-y1)-y0+y2+q)/(2.0*(1.0-w)*(y2-y0));        // t at P6 
		}
		q = 1.0/(2.0*t*(1.0-t)*(w-1.0)+1.0);                 // sub-divide at t 
		xx = (t*t*(x0-2.0*w*x1+x2)+2.0*t*(w*x1-x0)+x0)*q;               // = P6 
		yy = (t*t*(y0-2.0*w*y1+y2)+2.0*t*(w*y1-y0)+y0)*q;
		ww = t*(w-1.0)+1.0;
		ww *= ww*q;                    // squared weight P5 
		w = ((1.0-t)*(w-1.0)+1.0)*sqrt(q);                    // weight P7 
		x = floor(xx+0.5);
		y = floor(yy+0.5);                   // P6 
		xx = (x1-x0)*(yy-y0)/(y1-y0)+x0;                // intersect P6 | P0 P1 
		plotQuadRationalBezierWidthSeg(img,x0,y0, floor(xx+0.5),y, x,y, ww, th, color,alpha);
		xx = (x1-x2)*(yy-y2)/(y1-y2)+x2;                // intersect P7 | P1 P2 
		x1 = floor(xx+0.5);
		x0 = x;
		y0 = y1 = y;       // P0 = P6, P1 = P7 
	}
	plotQuadRationalBezierWidthSeg(img, x0,y0, x1,y1, x2,y2, w*w, th, color,alpha);
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotRotatedEllipse(Mat& img, int x, int y, int a, int b, double angle, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot ellipse rotated by angle (radian)
	double xd = (long)a*a, yd = (long)b*b;
	double s = sin(angle), zd = (xd-yd)*s;                  // ellipse rotation
	xd = sqrt(xd-zd*s), yd = sqrt(yd+zd*s);           // surrounding rectangle
	a = xd+0.5;
	b = yd+0.5;
	zd = zd*a*b/(xd*yd);           // scale to integer
	plotRotatedEllipseRect(img, x-a,y-b, x+a,y+b, (long)(4*zd*cos(angle)), color, alpha);
}

//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotRotatedEllipseWidth(Mat& img, int x, int y, int a, int b, double angle,double th, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot ellipse rotated by angle (radian)
	double xd = (long)a*a, yd = (long)b*b;
	double s = sin(angle), zd = (xd-yd)*s;                  // ellipse rotation
	xd = sqrt(xd-zd*s), yd = sqrt(yd+zd*s);           // surrounding rectangle
	a = xd+0.5;
	b = yd+0.5;
	zd = zd*a*b/(xd*yd);           // scale to integer
	plotRotatedEllipseRectWidth(img, x-a,y-b, x+a,y+b, (long)(4*zd*cos(angle)),th, color, alpha);
}

//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotRotatedEllipseRect(Mat& img, int x0, int y0, int x1, int y1, long zd, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// rectangle enclosing the ellipse, integer rotation angle
	int xd = x1-x0, yd = y1-y0;
	double w = xd*(long)yd;
	if (zd == 0)
	{
		return plotEllipseRect(img, x0,y0, x1,y1, color, alpha);    // looks nicer
	}
	if (w != 0.0)
	{
		w = (w-zd)/(w+w);    // squared weight of P1
	}
	assert(w <= 1.0 && w >= 0.0);                // limit angle to |zd|<=xd*yd
	xd = floor(xd*w+0.5);
	yd = floor(yd*w+0.5);           // snap xe,ye to int
	plotQuadRationalBezierSeg(img, x0,y0+yd, x0,y0, x0+xd,y0, 1.0-w,color, alpha);
	plotQuadRationalBezierSeg(img, x0,y0+yd, x0,y1, x1-xd,y1, w, color, alpha);
	plotQuadRationalBezierSeg(img, x1,y1-yd, x1,y1, x1-xd,y1, 1.0-w, color, alpha);
	plotQuadRationalBezierSeg(img, x1,y1-yd, x1,y0, x0+xd,y0, w, color,alpha);
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotRotatedEllipseRectWidth(Mat& img, int x0, int y0, int x1, int y1, long zd, double th, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// rectangle enclosing the ellipse, integer rotation angle
	int xd = x1-x0, yd = y1-y0;
	double w = xd*(long)yd;
	if (zd == 0)
	{
		return plotEllipseRectWidth(img, x0,y0, x1,y1,th, color, alpha);    // looks nicer
	}
	if (w != 0.0)
	{
		w = (w-zd)/(w+w);    // squared weight of P1
	}
	assert(w <= 1.0 && w >= 0.0);                // limit angle to |zd|<=xd*yd
	xd = floor(xd*w+0.5);
	yd = floor(yd*w+0.5);           // snap xe,ye to int
	plotQuadRationalBezierWidthSeg(img, x0,y0+yd, x0,y0, x0+xd,y0, 1.0-w,th,color, alpha);
	plotQuadRationalBezierWidthSeg(img, x0,y0+yd, x0,y1, x1-xd,y1, w,th, color, alpha);
	plotQuadRationalBezierWidthSeg(img, x1,y1-yd, x1,y1, x1-xd,y1, 1.0-w,th, color, alpha);
	plotQuadRationalBezierWidthSeg(img, x1,y1-yd, x1,y0, x0+xd,y0, w,th, color,alpha);
}
//------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------
void cvGDI::plotCubicBezierSegWidth(Mat& img,double x0,double y0,double  x1,double y1,double  x2,double y2,double  x3,double y3,double  th, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// split cubic Bezier segment in two quadratic segments 
	double x = floor((x0+3.0*x1+3*x2+x3+4.0)/8.0);
	double y = floor((y0+3.0*y1+3*y2+y3+4.0)/8.0);
	plotQuadRationalBezierWidthSeg(img, x0,y0,floor((x0+3*x1+2.0)/4.0),floor((y0+3*y1+2.0)/4.0), x,y, 1,th,color,alpha);
	plotQuadRationalBezierWidthSeg(img, x,y,floor((3*x2+x3+2.0)/4.0),floor((3*y2+y3+2.0)/4.0), x3,y3, 1,th,color,alpha);
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotCubicBezierSeg(Mat& img, int x0, int y0, double x1, double y1, double x2, double y2, int x3, int y3, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot limited cubic Bezier segment
	int f, fx, fy, leg = 1;
	int sx = x0 < x3 ? 1 : -1, sy = y0 < y3 ? 1 : -1;        // step direction
	double xc = -fabs(x0+x1-x2-x3), xa = xc-4*sx*(x1-x2), xb = sx*(x0-x1-x2+x3);
	double yc = -fabs(y0+y1-y2-y3), ya = yc-4*sy*(y1-y2), yb = sy*(y0-y1-y2+y3);
	double ab, ac, bc, cb, xx, xy, yy, dx, dy, ex, *pxy, EP = 0.01;
	// check for curve restrains
	// slope P0-P1 == P2-P3    and  (P0-P3 == P1-P2      or   no slope change)
	assert((x1-x0)*(x2-x3) < EP && ((x3-x0)*(x1-x2) < EP || xb*xb < xa*xc+EP));
	assert((y1-y0)*(y2-y3) < EP && ((y3-y0)*(y1-y2) < EP || yb*yb < ya*yc+EP));
	if (xa == 0 && ya == 0)                                // quadratic Bezier
	{
		sx = floor((3*x1-x0+1)/2);
		sy = floor((3*y1-y0+1)/2);   // new midpoint
		return plotQuadBezierSeg(img, x0,y0, sx,sy, x3,y3, color, alpha);
	}
	x1 = (x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+1;                    // line lengths
	x2 = (x2-x3)*(x2-x3)+(y2-y3)*(y2-y3)+1;
	do                                                  // loop over both ends
	{
		ab = xa*yb-xb*ya;
		ac = xa*yc-xc*ya;
		bc = xb*yc-xc*yb;
		ex = ab*(ab+ac-3*bc)+ac*ac;       // P0 part of self-intersection loop?
		f = ex > 0 ? 1 : sqrt(1+1024/x1);               // calculate resolution
		ab *= f;
		ac *= f;
		bc *= f;
		ex *= f*f;            // increase resolution
		xy = 9*(ab+ac+bc)/8;
		cb = 8*(xa-ya);  // init differences of 1st degree
		dx = 27*(8*ab*(yb*yb-ya*yc)+ex*(ya+2*yb+yc))/64-ya*ya*(xy-ya);
		dy = 27*(8*ab*(xb*xb-xa*xc)-ex*(xa+2*xb+xc))/64-xa*xa*(xy+xa);
		// init differences of 2nd degree
		xx = 3*(3*ab*(3*yb*yb-ya*ya-2*ya*yc)-ya*(3*ac*(ya+yb)+ya*cb))/4;
		yy = 3*(3*ab*(3*xb*xb-xa*xa-2*xa*xc)-xa*(3*ac*(xa+xb)+xa*cb))/4;
		xy = xa*ya*(6*ab+6*ac-3*bc+cb);
		ac = ya*ya;
		cb = xa*xa;
		xy = 3*(xy+9*f*(cb*yb*yc-xb*xc*ac)-18*xb*yb*ab)/8;
		if (ex < 0)           // negate values if inside self-intersection loop
		{
			dx = -dx;
			dy = -dy;
			xx = -xx;
			yy = -yy;
			xy = -xy;
			ac = -ac;
			cb = -cb;
		}                                     // init differences of 3rd degree
		ab = 6*ya*ac;
		ac = -6*xa*ac;
		bc = 6*ya*cb;
		cb = -6*xa*cb;
		dx += xy;
		ex = dx+dy;
		dy += xy;                    // error of 1st step
		for (pxy = &xy, fx = fy = f; x0 != x3 && y0 != y3; )
		{
			setPixel(img, x0,y0, color, alpha);                                       // plot curve
			do                                    // move sub-steps of one pixel
			{
				if (dx > *pxy || dy < *pxy)
				{
					goto exit;    // confusing values
				}
				y1 = 2*ex-dy;                    // save value for test of y step
				if (2*ex >= dx)                                     // x sub-step
				{
					fx--;
					ex += dx += xx;
					dy += xy += ac;
					yy += bc;
					xx += ab;
				}
				if (y1 <= 0)                                        // y sub-step
				{
					fy--;
					ex += dy += yy;
					dx += xy += bc;
					xx += ac;
					yy += cb;
				}
			}
			while (fx > 0 && fy > 0);                         // pixel complete?
			if (2*fx <= f)
			{
				x0 += sx;    // x step
				fx += f;
			}
			if (2*fy <= f)
			{
				y0 += sy;    // y step
				fy += f;
			}
			if (pxy == &xy && dx < 0 && dy > 0)
			{
				pxy = &EP;    // pixel ahead valid
			}
		}
exit:
		xx = x0;
		x0 = x3;
		x3 = xx;
		sx = -sx;
		xb = -xb;             // swap legs
		yy = y0;
		y0 = y3;
		y3 = yy;
		sy = -sy;
		yb = -yb;
		x1 = x2;
	}
	while (leg--);                                            // try other end
	plotLine(img, x0,y0, x3,y3, color, alpha);       // remaining part in case of cusp or crunode
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotCubicBezier(Mat& img, int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot any cubic Bezier curve
	int n = 0, i = 0;
	long xc = x0+x1-x2-x3, xa = xc-4*(x1-x2);
	long xb = x0-x1-x2+x3, xd = xb+4*(x1+x2);
	long yc = y0+y1-y2-y3, ya = yc-4*(y1-y2);
	long yb = y0-y1-y2+y3, yd = yb+4*(y1+y2);
	double fx0 = x0, fx1, fx2, fx3, fy0 = y0, fy1, fy2, fy3;
	double t1 = xb*xb-xa*xc, t2, t[5];
	// sub-divide curve at gradient sign changes
	if (xa == 0)                                                 // horizontal
	{
		if (abs(xc) < 2*abs(xb))
		{
			t[n++] = xc/(2.0*xb);    // one change
		}
	}
	else if (t1 > 0.0)                                          // two changes
	{
		t2 = sqrt(t1);
		t1 = (xb-t2)/xa;
		if (fabs(t1) < 1.0)
		{
			t[n++] = t1;
		}
		t1 = (xb+t2)/xa;
		if (fabs(t1) < 1.0)
		{
			t[n++] = t1;
		}
	}
	t1 = yb*yb-ya*yc;
	if (ya == 0)                                                   // vertical
	{
		if (abs(yc) < 2*abs(yb))
		{
			t[n++] = yc/(2.0*yb);    // one change
		}
	}
	else if (t1 > 0.0)                                          // two changes
	{
		t2 = sqrt(t1);
		t1 = (yb-t2)/ya;
		if (fabs(t1) < 1.0)
		{
			t[n++] = t1;
		}
		t1 = (yb+t2)/ya;
		if (fabs(t1) < 1.0)
		{
			t[n++] = t1;
		}
	}
	for (i = 1; i < n; i++)                         // bubble sort of 4 points
		if ((t1 = t[i-1]) > t[i])
		{
			t[i-1] = t[i];
			t[i] = t1;
			i = 0;
		}
		t1 = -1.0;
		t[n] = 1.0;                                // begin / end point
		for (i = 0; i <= n; i++)                   // plot each segment separately
		{
			t2 = t[i];                                // sub-divide at t[i-1], t[i]
			fx1 = (t1*(t1*xb-2*xc)-t2*(t1*(t1*xa-2*xb)+xc)+xd)/8-fx0;
			fy1 = (t1*(t1*yb-2*yc)-t2*(t1*(t1*ya-2*yb)+yc)+yd)/8-fy0;
			fx2 = (t2*(t2*xb-2*xc)-t1*(t2*(t2*xa-2*xb)+xc)+xd)/8-fx0;
			fy2 = (t2*(t2*yb-2*yc)-t1*(t2*(t2*ya-2*yb)+yc)+yd)/8-fy0;
			fx0 -= fx3 = (t2*(t2*(3*xb-t2*xa)-3*xc)+xd)/8;
			fy0 -= fy3 = (t2*(t2*(3*yb-t2*ya)-3*yc)+yd)/8;
			x3 = floor(fx3+0.5);
			y3 = floor(fy3+0.5);        // scale bounds to int
			if (fx0 != 0.0)
			{
				fx1 *= fx0 = (x0-x3)/fx0;
				fx2 *= fx0;
			}
			if (fy0 != 0.0)
			{
				fy1 *= fy0 = (y0-y3)/fy0;
				fy2 *= fy0;
			}
			if (x0 != x3 || y0 != y3)                            // segment t1 - t2
			{
				plotCubicBezierSeg(img, x0,y0, x0+fx1,y0+fy1, x0+fx2,y0+fy2, x3,y3, color, alpha);
			}
			x0 = x3;
			y0 = y3;
			fx0 = fx3;
			fy0 = fy3;
			t1 = t2;
		}
}


//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotCubicBezierWidth(Mat& img, double x0, double  y0, double  x1, double  y1, double  x2, double  y2, double  x3, double  y3,double th, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot any cubic Bezier curve
	int n = 0;
	int i = 0;
	double xc = x0+x1-x2-x3, xa = xc-4.0*(x1-x2);
	double xb = x0-x1-x2+x3, xd = xb+4.0*(x1+x2);
	double yc = y0+y1-y2-y3, ya = yc-4.0*(y1-y2);
	double yb = y0-y1-y2+y3, yd = yb+4.0*(y1+y2);
	double fx0 = x0, fx1, fx2, fx3, fy0 = y0, fy1, fy2, fy3;
	double t1 = xb*xb-xa*xc, t2, t[5];
	// sub-divide curve at gradient sign changes
	if (xa == 0)                                                 // horizontal
	{
		if (fabs(xc) <= 2.0*fabs(xb))
		{
			t[n++] = xc/(2.0*xb);    // one change
		}
	}
	else if (t1 > 0.0)                                          // two changes
	{
		t2 = sqrt(t1);
		t1 = (xb-t2)/xa;
		if (fabs(t1) < 1.0)
		{
			t[n++] = t1;
		}
		t1 = (xb+t2)/xa;
		if (fabs(t1) < 1.0)
		{
			t[n++] = t1;
		}
	}
	t1 = yb*yb-ya*yc;
	if (ya == 0)                                                   // vertical
	{
		if (abs(yc) < 2.0*fabs(yb))
		{
			t[n++] = yc/(2.0*yb);    // one change
		}
	}
	else if (t1 > 0.0)                                          // two changes
	{
		t2 = sqrt(t1);
		t1 = (yb-t2)/ya;
		if (fabs(t1) < 1.0)
		{
			t[n++] = t1;
		}
		t1 = (yb+t2)/ya;
		if (fabs(t1) < 1.0)
		{
			t[n++] = t1;
		}
	}
	for (i = 1; i < n; i++)                         // bubble sort of 4 points
		if ((t1 = t[i-1]) > t[i])
		{
			t[i-1] = t[i];
			t[i] = t1;
			i = 0;
		}
		t1 = -1.0;
		t[n] = 1.0;                                // begin / end point
		for (i = 0; i <= n; i++)                   // plot each segment separately
		{
			t2 = t[i];                                // sub-divide at t[i-1], t[i]
			fx1 = (t1*(t1*xb-2*xc)-t2*(t1*(t1*xa-2*xb)+xc)+xd)/8.0-fx0;
			fy1 = (t1*(t1*yb-2*yc)-t2*(t1*(t1*ya-2*yb)+yc)+yd)/8.0-fy0;
			fx2 = (t2*(t2*xb-2*xc)-t1*(t2*(t2*xa-2*xb)+xc)+xd)/8.0-fx0;
			fy2 = (t2*(t2*yb-2*yc)-t1*(t2*(t2*ya-2*yb)+yc)+yd)/8.0-fy0;
			fx0 -= fx3 = (t2*(t2*(3*xb-t2*xa)-3*xc)+xd)/8.0;
			fy0 -= fy3 = (t2*(t2*(3*yb-t2*ya)-3*yc)+yd)/8.0;
			x3 = floor(fx3+0.5);
			y3 = floor(fy3+0.5);        // scale bounds to int
			if (fx0 != 0.0)
			{
				fx1 *= fx0 = (x0-x3)/fx0;
				fx2 *= fx0;
			}
			if (fy0 != 0.0)
			{
				fy1 *= fy0 = (y0-y3)/fy0;
				fy2 *= fy0;
			}
			if (x0 != x3 || y0 != y3)                            // segment t1 - t2
			{
				plotCubicBezierSegWidth(img, x0,y0, x0+fx1,y0+fy1, x0+fx2,y0+fy2, x3,y3,th, color, alpha);
			}
			x0 = x3;
			y0 = y3;
			fx0 = fx3;
			fy0 = fy3;
			t1 = t2;
		}
}

//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotCubicBezierAA(Mat& img, int x0, int y0, int x1, int y1, int x2, int y2, int x3, int y3, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot any cubic Bezier curve
	int n = 0, i = 0;
	long xc = x0+x1-x2-x3, xa = xc-4*(x1-x2);
	long xb = x0-x1-x2+x3, xd = xb+4*(x1+x2);
	long yc = y0+y1-y2-y3, ya = yc-4*(y1-y2);
	long yb = y0-y1-y2+y3, yd = yb+4*(y1+y2);
	double fx0 = x0, fx1, fx2, fx3, fy0 = y0, fy1, fy2, fy3;
	double t1 = xb*xb-xa*xc, t2, t[5];
	// sub-divide curve at gradient sign changes
	if (xa == 0)                                                 // horizontal
	{
		if (abs(xc) < 2*abs(xb))
		{
			t[n++] = xc/(2.0*xb);    // one change
		}
	}
	else if (t1 > 0.0)                                          // two changes
	{
		t2 = sqrt(t1);
		t1 = (xb-t2)/xa;
		if (fabs(t1) < 1.0)
		{
			t[n++] = t1;
		}
		t1 = (xb+t2)/xa;
		if (fabs(t1) < 1.0)
		{
			t[n++] = t1;
		}
	}
	t1 = yb*yb-ya*yc;
	if (ya == 0)                                                   // vertical
	{
		if (abs(yc) < 2*abs(yb))
		{
			t[n++] = yc/(2.0*yb);    // one change
		}
	}
	else if (t1 > 0.0)                                          // two changes
	{
		t2 = sqrt(t1);
		t1 = (yb-t2)/ya;
		if (fabs(t1) < 1.0)
		{
			t[n++] = t1;
		}
		t1 = (yb+t2)/ya;
		if (fabs(t1) < 1.0)
		{
			t[n++] = t1;
		}
	}
	for (i = 1; i < n; i++)                         // bubble sort of 4 points
		if ((t1 = t[i-1]) > t[i])
		{
			t[i-1] = t[i];
			t[i] = t1;
			i = 0;
		}
		t1 = -1.0;
		t[n] = 1.0;                                // begin / end point
		for (i = 0; i <= n; i++)                   // plot each segment separately
		{
			t2 = t[i];                                // sub-divide at t[i-1], t[i]
			fx1 = (t1*(t1*xb-2*xc)-t2*(t1*(t1*xa-2*xb)+xc)+xd)/8-fx0;
			fy1 = (t1*(t1*yb-2*yc)-t2*(t1*(t1*ya-2*yb)+yc)+yd)/8-fy0;
			fx2 = (t2*(t2*xb-2*xc)-t1*(t2*(t2*xa-2*xb)+xc)+xd)/8-fx0;
			fy2 = (t2*(t2*yb-2*yc)-t1*(t2*(t2*ya-2*yb)+yc)+yd)/8-fy0;
			fx0 -= fx3 = (t2*(t2*(3*xb-t2*xa)-3*xc)+xd)/8;
			fy0 -= fy3 = (t2*(t2*(3*yb-t2*ya)-3*yc)+yd)/8;
			x3 = floor(fx3+0.5);
			y3 = floor(fy3+0.5);        // scale bounds to int
			if (fx0 != 0.0)
			{
				fx1 *= fx0 = (x0-x3)/fx0;
				fx2 *= fx0;
			}
			if (fy0 != 0.0)
			{
				fy1 *= fy0 = (y0-y3)/fy0;
				fy2 *= fy0;
			}
			if (x0 != x3 || y0 != y3)                            // segment t1 - t2
			{
				plotCubicBezierSegAA(img, x0,y0, x0+fx1,y0+fy1, x0+fx2,y0+fy2, x3,y3, color, alpha);    // << -----------------------------------
			}
			x0 = x3;
			y0 = y3;
			fx0 = fx3;
			fy0 = fy3;
			t1 = t2;
		}
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotCircleAA(Mat& img, int xm, int ym, int r, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// draw a black anti-aliased circle on white background
	int x = -r, y = 0;           // II. quadrant from bottom left to top right
	int x2, e2, err = 2-2*r;                             // error of 1.step
	r = 1-err;
	double a;
	double PI=3.1415926535897932384626433832795;
	double dl=0.5/(double)pattern->getScale();
	double dph=0;
	bool dash=1;
	do
	{
		if(pattern!=NULL)
		{
			dash=pattern->get_dash();
			pattern->advance();
		}
		//i = 255*abs(err-2*(x+y)-2)/r;               // get blend value of pixel
		a=double( 255*abs(err-2*(x+y)-2)/r )/255.0;
		a=1.0-((1.0-a)*(1.0-alpha));
		if(dash)
		{
			setPixel(img, xm-x, ym+y,color,a);                             //   I. Quadrant
		}
		
		if(dash)
		{
			setPixel(img, xm-y, ym-x,color,a);                             //  II. Quadrant
		}

		if(dash)
		{
			setPixel(img, xm+x, ym-y,color,a);                             // III. Quadrant
		}

		if(dash)
		{
			setPixel(img, xm+y, ym+x,color,a);                             //  IV. Quadrant
		}

		e2 = err;
		x2 = x;                                    // remember values
		if (err+y > 0)                                                // x step
		{
			//i = 255*(err-2*x-1)/r;                              // outward pixel
			a=double( 255*(err-2*x-1)/r )/255.0;
			a=1.0-((1.0-a)*(1.0-alpha));
			if (a <= 1)
			{

				if(dash)
				{
					setPixel(img,xm-x, ym+y+1, color,a);
				}

				if(dash)
				{
					setPixel(img,xm-y-1, ym-x, color,a);
				}

				if(dash)
				{
					setPixel(img,xm+x, ym-y-1, color,a);
				}

				if(dash)
				{
					setPixel(img,xm+y+1, ym+x, color,a);
				}
			}
		
		err += ++x*2+1;
	}
	if (e2+x2 <= 0)                                               // y step
	{
		// i = 255*(2*y+3-e2)/r;                                // inward pixel
		a=double( 255*(2*y+3-e2)/r )/255.0;
		a=1.0-((1.0-a)*(1.0-alpha));
		if (a <= 1)
		{

			if(dash)
			{
				setPixel(img,xm-x2-1, ym+y, color,a);
			}

			if(dash)
			{
				setPixel(img,xm-y, ym-x2-1, color,a);
			}

			if(dash)
			{
				setPixel(img,xm+x2+1, ym-y, color,a);
			}

			if(dash)
			{
				setPixel(img,xm+y, ym+x2+1, color,a);

			}
		}
	
	err += ++y*2+1;
}
	}while (x < 0);
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotEllipseRectAA(Mat& img, int x0, int y0, int x1, int y1, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// draw a black anti-aliased rectangular ellipse on white background
	long a = abs(x1-x0), b = abs(y1-y0), b1 = b&1;                 // diameter
	double dx = 4*(a-1.0)*b*b, dy = 4*(b1+1)*a*a;            // error increment
	double i, ed, err = b1*a*a-dx+dy;                        // error of 1.step
	bool f;
	double acol;
	if (a == 0 || b == 0)
	{
		return plotLineAA(img, x0,y0, x1,y1, color, alpha);
	}
	if (x0 > x1)
	{
		x0 = x1;    // if called with swapped points
		x1 += a;
	}
	if (y0 > y1)
	{
		y0 = y1;    // .. exchange them
	}
	y0 += (b+1)/2;
	y1 = y0-b1;                               // starting pixel
	a = 8*a*a;
	b1 = 8*b*b;
	for (;;)                               // approximate ed=sqrt(dx*dx+dy*dy)
	{
		i = min(dx,dy);
		ed = max(dx,dy);
		if (y0 == y1+1 && err > dy && a > b1)
		{
			ed = 255*4./a;    // x-tip
		}
		else
		{
			ed = 255/(ed+2*ed*i*i/(4*ed*ed+i*i));    // approximation
		}
		// i = ed*fabs(err+dx-dy);           // get intensity value by pixel error
		acol=double( ed*fabs(err+dx-dy) )/255.0;
		acol=1.0-((1.0-acol)*(1.0-alpha));
		setPixel(img, x0,y0, color,acol);
		setPixel(img, x0,y1, color,acol);
		setPixel(img, x1,y0, color,acol);
		setPixel(img, x1,y1, color,acol);
		if (f = 2*err+dy >= 0)                    // x step, remember condition
		{
			if (x0 >= x1)
			{
				break;
			}
			//i = ed*(err+dx);
			acol=double( ed*(err+dx) )/255.0;
			acol=1.0-((1.0-acol)*(1.0-alpha));
			if (acol<=1.0)
			{
				setPixel(img, x0,y0+1, color, acol);
				setPixel(img, x0,y1-1, color, acol);
				setPixel(img, x1,y0+1, color, acol);
				setPixel(img, x1,y1-1, color, acol);
			}          // do error increment later since values are still needed
		}
		if (2*err <= dx)                                              // y step
		{
			// i = ed*(dy-err);
			acol=double( ed*(dy-err) )/255.0;
			acol=1.0-((1.0-acol)*(1.0-alpha));
			if (acol <= 1.0)
			{
				setPixel(img, x0+1,y0, color,acol);
				setPixel(img, x1-1,y0, color,acol);
				setPixel(img, x0+1,y1, color,acol);
				setPixel(img, x1-1,y1, color,acol);
			}
			y0++;
			y1--;
			err += dy += a;
		}
		if (f)
		{
			x0++;    // x error increment
			x1--;
			err -= dx -= b1;
		}
	}
	if (--x0 == x1++)                       // too early stop of flat ellipses
		while (y0-y1 < b)
		{
			// i = 255*4*fabs(err+dx)/b1;               // -> finish tip of ellipse
			acol=double( 255*4*fabs(err+dx)/b1 )/255.0;
			acol=1.0-((1.0-acol)*(1.0-alpha));
			setPixel(img, x0,++y0	,color,acol);
			setPixel(img, x1,y0		,color,acol);
			setPixel(img, x0,--y1	,color,acol);
			setPixel(img, x1,y1		,color,acol);
			err += dy += a;
		}
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotQuadBezierSegAA(Mat& img, int x0, int y0, int x1, int y1, int x2, int y2, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// draw an limited anti-aliased quadratic Bezier segment
	int sx = x2-x1, sy = y2-y1;
	long xx = x0-x1, yy = y0-y1, xy;             // relative values for checks
	double dx, dy, err, ed, cur = xx*sy-yy*sx;                    // curvature
	double a;
	assert(xx*sx <= 0 && yy*sy <= 0);      // sign of gradient must not change
	if (sx*(long)sx+sy*(long)sy > xx*xx+yy*yy)       // begin with longer part
	{
		x2 = x0;
		x0 = sx+x1;
		y2 = y0;
		y0 = sy+y1;
		cur = -cur;     // swap P0 P2
	}
	if (cur != 0)
	{
		// no straight line
		xx += sx;
		xx *= sx = x0 < x2 ? 1 : -1;              // x step direction
		yy += sy;
		yy *= sy = y0 < y2 ? 1 : -1;              // y step direction
		xy = 2*xx*yy;
		xx *= xx;
		yy *= yy;             // differences 2nd degree
		if (cur*sx*sy < 0)                                // negated curvature?
		{
			xx = -xx;
			yy = -yy;
			xy = -xy;
			cur = -cur;
		}
		dx = 4.0*sy*(x1-x0)*cur+xx-xy;                // differences 1st degree
		dy = 4.0*sx*(y0-y1)*cur+yy-xy;
		xx += xx;
		yy += yy;
		err = dx+dy+xy;                   // error 1st step
		do
		{
			cur = min(dx+xy,-xy-dy);
			ed = max(dx+xy,-xy-dy);               // approximate error distance
			ed += 2*ed*cur*cur/(4*ed*ed+cur*cur);
			a=double( 255*fabs(err-dx-dy-xy)/ed )/255.0;
			a=1.0-((1.0-a)*(1.0-alpha));
			setPixel(img, x0,y0, color,a);          // plot curve
			if (x0 == x2 || y0 == y2)
			{
				break;    // last pixel -> curve finished
			}
			x1 = x0;
			cur = dx-err;
			y1 = 2*err+dy < 0;
			if (2*err+dx > 0)                                          // x step
			{
				if (err-dy < ed)
				{
					a=double(  255*fabs(err-dy)/ed )/255.0;
					a=1.0-((1.0-a)*(1.0-alpha));
					setPixel(img,x0,y0+sy, color,a);
				}
				x0 += sx;
				dx -= xy;
				err += dy += yy;
			}
			if (y1)                                                    // y step
			{
				if (cur < ed)
				{
					a=double(  255*fabs(cur)/ed )/255.0;
					a=1.0-((1.0-a)*(1.0-alpha));
					setPixel(img,x1+sx,y0, color,a);
				}
				y0 += sy;
				dy -= xy;
				err += dx += xx;
			}
		}
		while (dy < dx);                    // gradient negates -> close curves
	}
	plotLineAA(img,x0,y0, x2,y2,color, alpha);                  // plot remaining needle to end
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotQuadRationalBezierSegAA(Mat& img, int x0, int y0, int x1, int y1,int x2, int y2, double w, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// draw an anti-aliased rational quadratic Bezier segment, squared weight
	int sx = x2-x1, sy = y2-y1;                  // relative values for checks
	double dx = x0-x2, dy = y0-y2, xx = x0-x1, yy = y0-y1;
	double xy = xx*sy+yy*sx, cur = xx*sy-yy*sx, err, ed;          // curvature
	bool f;
	double a;
	assert(xx*sx <= 0.0 && yy*sy <= 0.0);  // sign of gradient must not change
	if (cur != 0.0 && w > 0.0)                             // no straight line
	{
		if (sx*(long)sx+sy*(long)sy > xx*xx+yy*yy)    // begin with longer part
		{
			x2 = x0;
			x0 -= dx;
			y2 = y0;
			y0 -= dy;
			cur = -cur;      // swap P0 P2
		}
		xx = 2.0*(4.0*w*sx*xx+dx*dx);                 // differences 2nd degree
		yy = 2.0*(4.0*w*sy*yy+dy*dy);
		sx = x0 < x2 ? 1 : -1;                              // x step direction
		sy = y0 < y2 ? 1 : -1;                              // y step direction
		xy = -2.0*sx*sy*(2.0*w*xy+dx*dy);
		if (cur*sx*sy < 0)                                // negated curvature?
		{
			xx = -xx;
			yy = -yy;
			cur = -cur;
			xy = -xy;
		}
		dx = 4.0*w*(x1-x0)*sy*cur+xx/2.0+xy;          // differences 1st degree
		dy = 4.0*w*(y0-y1)*sx*cur+yy/2.0+xy;
		if (w < 0.5 && dy > dx)                // flat ellipse, algorithm fails
		{
			cur = (w+1.0)/2.0;
			w = sqrt(w);
			xy = 1.0/(w+1.0);
			sx = floor((x0+2.0*w*x1+x2)*xy/2.0+0.5); // subdivide curve in half
			sy = floor((y0+2.0*w*y1+y2)*xy/2.0+0.5);
			dx = floor((w*x1+x0)*xy+0.5);
			dy = floor((y1*w+y0)*xy+0.5);
			plotQuadRationalBezierSegAA(img,x0,y0, dx,dy, sx,sy, cur,color, alpha); // plot apart
			dx = floor((w*x1+x2)*xy+0.5);
			dy = floor((y1*w+y2)*xy+0.5);
			return plotQuadRationalBezierSegAA(img, sx,sy, dx,dy, x2,y2, cur, color, alpha);
		}
		err = dx+dy-xy;                                       // error 1st step
		do                                                        // pixel loop
		{
			cur = min(dx-xy,xy-dy);
			ed = max(dx-xy,xy-dy);
			ed += 2*ed*cur*cur/(4.*ed*ed+cur*cur); // approximate error distance
			a=double(  255*fabs(err-dx-dy+xy)/ed )/255.0;
			a=1.0-((1.0-a)*(1.0-alpha));
			//x1 = 255*fabs(err-dx-dy+xy)/ed;    // get blend value by pixel error
			if (a<=1)
			{
				setPixel(img,x0,y0, color, a);    // plot curve
			}
			if (f = 2*err+dy < 0)                                      // y step
			{
				if (y0 == y2)
				{
					return;    // last pixel -> curve finished
				}
				if (dx-err < ed)
				{
					a=double(  255*fabs(dx-err)/ed )/255.0;
					a=1.0-((1.0-a)*(1.0-alpha));
					setPixel(img,x0+sx,y0, color,a);
				}
			}
			if (2*err+dx > 0)                                          // x step
			{
				if (x0 == x2)
				{
					return;    // last pixel -> curve finished
				}
				if (err-dy < ed)
				{
					a=double(  255*fabs(err-dy)/ed )/255.0;
					a=1.0-((1.0-a)*(1.0-alpha));
					setPixel(img, x0,y0+sy, color,a);
				}
				x0 += sx;
				dx += xy;
				err += dy += yy;
			}
			if (f)
			{
				y0 += sy;    // y step
				dy += xy;
				err += dx += xx;
			}
		}
		while (dy < dx);                 // gradient negates -> algorithm fails
	}
	plotLineAA(img, x0,y0, x2,y2, color ,alpha);                  // plot remaining needle to end
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotCubicBezierSegAA(Mat& img, int x0, int y0, double x1, double y1, double x2, double y2, int x3, int y3, cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot limited anti-aliased cubic Bezier segment
	int f, fx, fy, leg = 1;
	int sx = x0 < x3 ? 1 : -1, sy = y0 < y3 ? 1 : -1;        // step direction
	double xc = -fabs(x0+x1-x2-x3), xa = xc-4*sx*(x1-x2), xb = sx*(x0-x1-x2+x3);
	double yc = -fabs(y0+y1-y2-y3), ya = yc-4*sy*(y1-y2), yb = sy*(y0-y1-y2+y3);
	double ab, ac, bc, ba, xx, xy, yy, dx, dy, ex, px, py, ed, ip, EP = 0.01;
	double a;
	// check for curve restrains
	// slope P0-P1 == P2-P3     and  (P0-P3 == P1-P2      or  no slope change)
	assert((x1-x0)*(x2-x3) < EP && ((x3-x0)*(x1-x2) < EP || xb*xb < xa*xc+EP));
	assert((y1-y0)*(y2-y3) < EP && ((y3-y0)*(y1-y2) < EP || yb*yb < ya*yc+EP));
	if (xa == 0 && ya == 0)                                // quadratic Bezier
	{
		sx = floor((3*x1-x0+1)/2);
		sy = floor((3*y1-y0+1)/2);   // new midpoint
		return plotQuadBezierSegAA(img, x0,y0, sx,sy, x3,y3, color, alpha);
	}
	x1 = (x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+1;                    // line lengths
	x2 = (x2-x3)*(x2-x3)+(y2-y3)*(y2-y3)+1;
	do                                                  // loop over both ends
	{
		ab = xa*yb-xb*ya;
		ac = xa*yc-xc*ya;
		bc = xb*yc-xc*yb;
		ip = 4*ab*bc-ac*ac;                   // self intersection loop at all?
		ex = ab*(ab+ac-3*bc)+ac*ac;       // P0 part of self-intersection loop?
		f = ex > 0 ? 1 : sqrt(1+1024/x1);               // calculate resolution
		ab *= f;
		ac *= f;
		bc *= f;
		ex *= f*f;            // increase resolution
		xy = 9*(ab+ac+bc)/8;
		ba = 8*(xa-ya);  // init differences of 1st degree
		dx = 27*(8*ab*(yb*yb-ya*yc)+ex*(ya+2*yb+yc))/64-ya*ya*(xy-ya);
		dy = 27*(8*ab*(xb*xb-xa*xc)-ex*(xa+2*xb+xc))/64-xa*xa*(xy+xa);
		// init differences of 2nd degree
		xx = 3*(3*ab*(3*yb*yb-ya*ya-2*ya*yc)-ya*(3*ac*(ya+yb)+ya*ba))/4;
		yy = 3*(3*ab*(3*xb*xb-xa*xa-2*xa*xc)-xa*(3*ac*(xa+xb)+xa*ba))/4;
		xy = xa*ya*(6*ab+6*ac-3*bc+ba);
		ac = ya*ya;
		ba = xa*xa;
		xy = 3*(xy+9*f*(ba*yb*yc-xb*xc*ac)-18*xb*yb*ab)/8;
		if (ex < 0)           // negate values if inside self-intersection loop
		{
			dx = -dx;
			dy = -dy;
			xx = -xx;
			yy = -yy;
			xy = -xy;
			ac = -ac;
			ba = -ba;
		}                                     // init differences of 3rd degree
		ab = 6*ya*ac;
		ac = -6*xa*ac;
		bc = 6*ya*ba;
		ba = -6*xa*ba;
		dx += xy;
		ex = dx+dy;
		dy += xy;                    // error of 1st step
		for (fx = fy = f; x0 != x3 && y0 != y3; )
		{
			y1 = min(fabs(xy-dx),fabs(dy-xy));
			ed = max(fabs(xy-dx),fabs(dy-xy));    // approximate error distance
			ed = f*(ed+2*ed*y1*y1/(4*ed*ed+y1*y1));
			// y1 = 255*fabs(ex-(f-fx+1)*dx-(f-fy+1)*dy+f*xy)/ed;
			a=double(  255*fabs(ex-(f-fx+1)*dx-(f-fy+1)*dy+f*xy)/ed )/255.0;
			a=1.0-((1.0-a)*(1.0-alpha));
			if (a<=1.0)
			{
				setPixel(img,x0, y0,color,a);    // plot curve
			}
			px = fabs(ex-(f-fx+1)*dx+(fy-1)*dy);       // pixel intensity x move
			py = fabs(ex+(fx-1)*dx-(f-fy+1)*dy);       // pixel intensity y move
			y2 = y0;
			do                                    // move sub-steps of one pixel
			{
				if (ip >= -EP)               // intersection possible? -> check..
					if (dx+xx > xy || dy+yy < xy)
					{
						goto exit;    // two x or y steps
					}
					y1 = 2*ex+dx;                    // save value for test of y step
					if (2*ex+dy > 0)                                    // x sub-step
					{
						fx--;
						ex += dx += xx;
						dy += xy += ac;
						yy += bc;
						xx += ab;
					}
					else if (y1 > 0)
					{
						goto exit;    // tiny nearly cusp
					}
					if (y1 <= 0)                                        // y sub-step
					{
						fy--;
						ex += dy += yy;
						dx += xy += bc;
						xx += ac;
						yy += ba;
					}
			}
			while (fx > 0 && fy > 0);                         // pixel complete?
			if (2*fy <= f)                             // x+ anti-aliasing pixel
			{
				if (py < ed)
				{
					a=double(  255*py/ed )/255.0;
					a=1.0-((1.0-a)*(1.0-alpha));
					setPixel(img,x0+sx, y0, color,a);    // plot curve
				}
				y0 += sy;
				fy += f;                                      // y step
			}
			if (2*fx <= f)                             // y+ anti-aliasing pixel
			{
				if (px < ed)
				{
					a=double( 255*px/ed )/255.0;
					a=1.0-((1.0-a)*(1.0-alpha));
					setPixel(img,x0, y2+sy, color,a);    // plot curve
				}
				x0 += sx;
				fx += f;                                      // x step
			}
		}
		break;                                          // finish curve by line
exit:
		if (2*ex < dy && 2*fy <= f+2)           // round x+ approximation pixel
		{
			if (py < ed)
			{
				a=double( 255*py/ed )/255.0;
				a=1.0-((1.0-a)*(1.0-alpha));
				setPixel(img,x0+sx, y0, color,a);    // plot curve
			}
			y0 += sy;
		}
		if (2*ex > dx && 2*fx <= f+2)           // round y+ approximation pixel
		{
			if (px < ed)
			{
				a=double( 255*px/ed )/255.0;
				a=1.0-((1.0-a)*(1.0-alpha));
				setPixel(img,x0, y2+sy, color,a);    // plot curve
			}
			x0 += sx;
		}
		xx = x0;
		x0 = x3;
		x3 = xx;
		sx = -sx;
		xb = -xb;             // swap legs
		yy = y0;
		y0 = y3;
		y3 = yy;
		sy = -sy;
		yb = -yb;
		x1 = x2;
	}
	while (leg--);                                            // try other end
	plotLineAA(img,x0,y0, x3,y3,color,alpha);     // remaining part in case of cusp or crunode
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotQuadSpline(Mat& img, int n, int x[], int y[],cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot quadratic spline, destroys input arrays x,y
#define M_MAX 6
	double mi = 1, m[M_MAX];                    // diagonal constants of matrix
	int i, x0, y0, x1, y1, x2 = x[n], y2 = y[n];
	assert(n > 1);                        //need at least 3 points P[0]..P[n]
	x[1] = x0 = 8*x[1]-2*x[0];                          // first row of matrix
	y[1] = y0 = 8*y[1]-2*y[0];
	for (i = 2; i < n; i++)                                   // forward sweep
	{
		if (i-2 < M_MAX)
		{
			m[i-2] = mi = 1.0/(6.0-mi);
		}
		x[i] = x0 = floor(8*x[i]-x0*mi+0.5);                        // store yi
		y[i] = y0 = floor(8*y[i]-y0*mi+0.5);
	}
	x1 = floor((x0-2*x2)/(5.0-mi)+0.5);                 // correction last row
	y1 = floor((y0-2*y2)/(5.0-mi)+0.5);
	for (i = n-2; i > 0; i--)                             // back substitution
	{
		if (i <= M_MAX)
		{
			mi = m[i-1];
		}
		x0 = floor((x[i]-x1)*mi+0.5);                            // next corner
		y0 = floor((y[i]-y1)*mi+0.5);
		plotQuadBezier(img,(x0+x1)/2.0,(y0+y1)/2.0, x1,y1, x2,y2,color,alpha);
		x2 = (x0+x1)/2.0;
		x1 = x0;
		y2 = (y0+y1)/2.0;
		y1 = y0;
	}
	plotQuadBezier(img,x[0],y[0], x1,y1, x2,y2,color,alpha);
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotQuadSplineAA(Mat& img, int n, int x[], int y[],cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot quadratic spline, destroys input arrays x,y
#define M_MAX 6
	double mi = 1, m[M_MAX];                    // diagonal constants of matrix
	int i, x0, y0, x1, y1, x2 = x[n], y2 = y[n];
	assert(n > 1);                        //need at least 3 points P[0]..P[n]
	x[1] = x0 = 8*x[1]-2*x[0];                          // first row of matrix
	y[1] = y0 = 8*y[1]-2*y[0];
	for (i = 2; i < n; i++)                                   // forward sweep
	{
		if (i-2 < M_MAX)
		{
			m[i-2] = mi = 1.0/(6.0-mi);
		}
		x[i] = x0 = floor(8*x[i]-x0*mi+0.5);                        // store yi
		y[i] = y0 = floor(8*y[i]-y0*mi+0.5);
	}
	x1 = floor((x0-2*x2)/(5.0-mi)+0.5);                 // correction last row
	y1 = floor((y0-2*y2)/(5.0-mi)+0.5);
	for (i = n-2; i > 0; i--)                             // back substitution
	{
		if (i <= M_MAX)
		{
			mi = m[i-1];
		}
		x0 = floor((x[i]-x1)*mi+0.5);                            // next corner
		y0 = floor((y[i]-y1)*mi+0.5);
		plotQuadBezierAA(img,(x0+x1)/2.0,(y0+y1)/2.0, x1,y1, x2,y2,color,alpha);
		x2 = (x0+x1)/2.0;
		x1 = x0;
		y2 = (y0+y1)/2.0;
		y1 = y0;
	}
	plotQuadBezierAA(img,x[0],y[0], x1,y1, x2,y2,color,alpha);
}

//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotCubicSpline(Mat& img, int n, int x[], int y[],cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot cubic spline, destroys input arrays x,y
#define M_MAX 6
	double mi = 0.25, m[M_MAX];                 // diagonal constants of matrix
	int x3 = x[n-1], y3 = y[n-1], x4 = x[n], y4 = y[n];
	int i, x0, y0, x1, y1, x2, y2;
	assert(n > 2);                        // need at least 4 points P[0]..P[n]
	x[1] = x0 = 12*x[1]-3*x[0];                         // first row of matrix
	y[1] = y0 = 12*y[1]-3*y[0];
	for (i = 2; i < n; i++)                                  // foreward sweep
	{
		if (i-2 < M_MAX)
		{
			m[i-2] = mi = 0.25/(2.0-mi);
		}
		x[i] = x0 = floor(12*x[i]-2*x0*mi+0.5);
		y[i] = y0 = floor(12*y[i]-2*y0*mi+0.5);
	}
	x2 = floor((x0-3*x4)/(7-4*mi)+0.5);                    // correct last row
	y2 = floor((y0-3*y4)/(7-4*mi)+0.5);
	plotCubicBezier(img, x3,y3, (x2+x4)/2,(y2+y4)/2, x4,y4, x4,y4, color,alpha);
	if (n-3 < M_MAX)
	{
		mi = m[n-3];
	}
	x1 = floor((x[n-2]-2*x2)*mi+0.5);
	y1 = floor((y[n-2]-2*y2)*mi+0.5);
	for (i = n-3; i > 0; i--)                             // back substitution
	{
		if (i <= M_MAX)
		{
			mi = m[i-1];
		}
		x0 = floor((x[i]-2*x1)*mi+0.5);
		y0 = floor((y[i]-2*y1)*mi+0.5);
		x4 = floor((x0+4*x1+x2+3)/6.0);                     // reconstruct P[i]
		y4 = floor((y0+4*y1+y2+3)/6.0);
		plotCubicBezier(img, x4,y4,
			floor((2*x1+x2)/3+0.5),floor((2*y1+y2)/3+0.5),
			floor((x1+2*x2)/3+0.5),floor((y1+2*y2)/3+0.5),
			x3,y3,color,alpha);
		x3 = x4;
		y3 = y4;
		x2 = x1;
		y2 = y1;
		x1 = x0;
		y1 = y0;
	}
	x0 = x[0];
	x4 = floor((3*x0+7*x1+2*x2+6)/12.0);        // reconstruct P[1]
	y0 = y[0];
	y4 = floor((3*y0+7*y1+2*y2+6)/12.0);
	plotCubicBezier(img, x4,y4, floor((2*x1+x2)/3+0.5),floor((2*y1+y2)/3+0.5), floor((x1+2*x2)/3+0.5),floor((y1+2*y2)/3+0.5), x3,y3, color,alpha);
	plotCubicBezier(img, x0,y0, x0,y0, (x0+x1)/2,(y0+y1)/2, x4,y4, color,alpha);
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotCubicSplineWidth(Mat& img, int n, int x[], int y[],double th,cv::Scalar color, double alpha, cvGDIPattern* pattern)
{

	int x_0=x[0];
	int x_1=x[1];
	int x_2=x[2];
	int x_3=x[3];

	// plot cubic spline, destroys input arrays x,y
#define M_MAX 6
	double mi = 0.25, m[M_MAX];                 // diagonal constants of matrix
	int x3 = x[n-1], y3 = y[n-1], x4 = x[n], y4 = y[n];
	int i, x0, y0, x1, y1, x2, y2;
	assert(n > 2);                        // need at least 4 points P[0]..P[n]
	x[1] = x0 = 12.0*x[1]-3*x[0];                         // first row of matrix
	y[1] = y0 = 12.0*y[1]-3*y[0];
	for (i = 2; i < n; i++)                                  // foreward sweep
	{
		if (i-2 < M_MAX)
		{
			m[i-2] = mi = 0.25/(2.0-mi);
		}
		x[i] = x0 = floor(12*x[i]-2*x0*mi+0.5);
		y[i] = y0 = floor(12*y[i]-2*y0*mi+0.5);
	}
	x2 = floor((x0-3*x4)/(7.0-4.0*mi)+0.5);                    // correct last row
	y2 = floor((y0-3*y4)/(7.0-4.0*mi)+0.5);
	plotCubicBezierWidth(img, x3,y3, (x2+x4)*0.5,(y2+y4)*0.5, x4,y4, x4,y4,th, color,alpha);
	if (n-3 < M_MAX)
	{
		mi = m[n-3];
	}
	x1 = floor((x[n-2]-2.0*x2)*mi+0.5);
	y1 = floor((y[n-2]-2.0*y2)*mi+0.5);
	for (i = n-3; i > 0; i--)                             // back substitution
	{
		if (i <= M_MAX)
		{
			mi = m[i-1];
		}
		x0 = floor((x[i]-2.0*x1)*mi+0.5);
		y0 = floor((y[i]-2.0*y1)*mi+0.5);
		x4 = floor((x0+4.0*x1+x2+3.0)/6.0);                     // reconstruct P[i]
		y4 = floor((y0+4.0*y1+y2+3.0)/6.0);
		plotCubicBezierWidth(img, x4,y4,floor((2.0*x1+x2)/3.0+0.5),floor((2.0*y1+y2)/3.0+0.5),floor((x1+2.0*x2)/3.0+0.5),floor((y1+2*y2)/3.0+0.5),x3,y3,th,color,alpha);
		x3 = x4;
		y3 = y4;
		x2 = x1;
		y2 = y1;
		x1 = x0;
		y1 = y0;
	}
	x0 = x[0];
	x4 = floor((3*x0+7*x1+2*x2+6)/12.0);        // reconstruct P[1]
	y0 = y[0];
	y4 = floor((3*y0+7*y1+2*y2+6)/12.0);
	plotCubicBezierWidth(img, x4,y4, floor((2*x1+x2)/3.0+0.5),floor((2*y1+y2)/3.0+0.5), floor((x1+2*x2)/3.0+0.5),floor((y1+2*y2)/3.0+0.5), x3,y3,th, color,alpha);
	plotCubicBezierWidth(img, x0,y0, x0,y0, (x0+x1)*0.5,(y0+y1)*0.5, x4,y4,th, color,alpha);
}
//------------------------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------------------------
void cvGDI::plotCubicSplineAA(Mat& img, int n, int x[], int y[],cv::Scalar color, double alpha, cvGDIPattern* pattern)
{
	// plot cubic spline, destroys input arrays x,y
#define M_MAX 6
	double mi = 0.25, m[M_MAX];                 // diagonal constants of matrix
	int x3 = x[n-1], y3 = y[n-1], x4 = x[n], y4 = y[n];
	int i, x0, y0, x1, y1, x2, y2;
	assert(n > 2);                        // need at least 4 points P[0]..P[n]
	x[1] = x0 = 12*x[1]-3*x[0];                         // first row of matrix
	y[1] = y0 = 12*y[1]-3*y[0];
	for (i = 2; i < n; i++)                                  // foreward sweep
	{
		if (i-2 < M_MAX)
		{
			m[i-2] = mi = 0.25/(2.0-mi);
		}
		x[i] = x0 = floor(12*x[i]-2*x0*mi+0.5);
		y[i] = y0 = floor(12*y[i]-2*y0*mi+0.5);
	}
	x2 = floor((x0-3*x4)/(7-4*mi)+0.5);                    // correct last row
	y2 = floor((y0-3*y4)/(7-4*mi)+0.5);
	plotCubicBezierAA(img, x3,y3, (x2+x4)/2,(y2+y4)/2, x4,y4, x4,y4, color,alpha,pattern);
	if (n-3 < M_MAX)
	{
		mi = m[n-3];
	}
	x1 = floor((x[n-2]-2*x2)*mi+0.5);
	y1 = floor((y[n-2]-2*y2)*mi+0.5);
	for (i = n-3; i > 0; i--)                             // back substitution
	{
		if (i <= M_MAX)
		{
			mi = m[i-1];
		}
		x0 = floor((x[i]-2*x1)*mi+0.5);
		y0 = floor((y[i]-2*y1)*mi+0.5);
		x4 = floor((x0+4*x1+x2+3)/6.0);                     // reconstruct P[i]
		y4 = floor((y0+4*y1+y2+3)/6.0);
		plotCubicBezierAA(img,x4,y4,
			floor((2*x1+x2)/3+0.5),floor((2*y1+y2)/3+0.5),
			floor((x1+2*x2)/3+0.5),floor((y1+2*y2)/3+0.5),
			x3,y3, color,alpha,pattern);
		x3 = x4;
		y3 = y4;
		x2 = x1;
		y2 = y1;
		x1 = x0;
		y1 = y0;
	}
	x0 = x[0];
	x4 = floor((3*x0+7*x1+2*x2+6)/12.0);        // reconstruct P[1]
	y0 = y[0];
	y4 = floor((3*y0+7*y1+2*y2+6)/12.0);
	plotCubicBezierAA(img,x4,y4, floor((2*x1+x2)/3+0.5),floor((2*y1+y2)/3+0.5), floor((x1+2*x2)/3+0.5),floor((y1+2*y2)/3+0.5), x3,y3, color,alpha,pattern);
	plotCubicBezierAA(img,x0,y0, x0,y0, (x0+x1)/2,(y0+y1)/2, x4,y4, color,alpha,pattern);
}
//-----------------------------------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------------------------------
