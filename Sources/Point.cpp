/*
 *  Point.cpp
 *
 *  Created by Silvanus Alt on 3/25/12.
 *  Copyright 2012 Max Planck Society
 *
 */

#include "Point.h"
#define PI 3.1415926535897932384

Point::Point()
{
	x=0;
    y=0;
    z=0;
}

Point::Point(const double &px, const double &py, const double &pz) : x(px), y(py), z(pz)
{
}

double Point::Distance(const Point &p) const
{
	return sqrt(pow(x-p.x,2)+pow(y-p.y,2)+pow(z-p.z,2));
}

// vector product: this x p
Point Point::VectorProduct(const Point &p) const
{
	return(Point(y*p.z-z*p.y,
                 z*p.x-x*p.z,
                 x*p.y-y*p.x));
}


Point Point::operator * (const Point &p) const	
{
	return Point(y*p.z-z*p.y,z*p.x-x*p.z,x*p.y-y*p.x);
}


// multiplication with scalar, has been newly introduced
Point Point::operator*(const double &l) const
{
	return(Point(l*x,
                 l*y,
                 l*z));
}

// addition of two points
Point Point::operator + (const Point &p) const
{
	return Point(x+p.x,
				 y+p.y,
				 z+p.z);
};

// subtraction of a point from another, has been changed
Point Point::operator - (const Point &p) const
{
	return Point(x-p.x,
				 y-p.y,
				 z-p.z);
};

// comparison of two points - they are equal, iff all the coordinates are the same; has been changed
bool Point::operator==(const Point &p) const
{
	return((x==p.x)&
           (y==p.y)&
           (z==p.z));
	
}

// division by double
Point Point::operator/(const double &div) const
{
	if(div==0)
    {
        std::cout << "Attempted division of a point by scalar zero." << std::endl;
        //QApplication::notify();
    }
	return Point(x/div,
                 y/div,
                 z/div);
}		   

bool Point::isZero()
{
	return (x==0 & y==0 & z==0);
}


void Point::Print()
{
	std::cout<<"("<<x<<","<<y<<","<< std::setprecision(10) << z<<")\n";
}

void Point::Print_in_line()
{
	std::cout<<"("<<x<<","<<y<<","<<z<<")";
}

// project point on ellipsoid with axis ex, ey, ez
Point Point::projectOnEllipse(double ex, double ez)
{
    bool pseudoProject = false;
    if(pseudoProject) { //pseudo projection, where the point is moved towards the origin until it crosses the ellipse
    
        double lambda = 1./(sqrt( pow(x/ex,2) + pow(y/ex,2) + pow(z/ez,2))  );
        
        // move the point back on the ellipse
        return Point(x,y,z)*lambda;
    }
    else { //proper projection, where the point is moved towards the origin until it crosses the ellipse
        double t = 0;         // the initial guess close and < to the actual root
        double tol = 0.0001;    // 0.0001 is the error level we wish
        double told;
        
        int iteration = 0;
        
        do
        {
            iteration++;
            told = t;
            t = newtonRaphson(t,ex,ez,x,y,z);
            
            //std::cout << std::abs(told - t) << std::endl;
            
            if(iteration>1000) {
                
                std::cout << "Problematic point: (" << x << " , " <<  y << " , " << z << ")\n";
            }
        }
        while (std::abs(told - t) > tol);     // stop criterion
        
        //std::cout << "t: " << t << std::endl;
        
        double x_p = x/(1+t/std::pow(ex,2));
        double y_p = y/(1+t/std::pow(ex,2));
        double z_p = z/(1+t/std::pow(ez,2));
        
        // test algorithm by calculating the normal, adding it to the projected point and see if it gives the right result
        bool testAlgorithm = false;
        if(testAlgorithm) {
            // calculate the normal on the ellipsoid
            Point Normal = Point(2*x_p/std::pow(ex,2),  2*y_p/std::pow(ex,2), 2*z_p/std::pow(ez ,2));
            
            double dist = Norm(Point(x_p,y_p,z_p)+Normal*(t/2) - *this);
            
            std::cout << "The error is " << dist << " and the result was reached after " << iteration << " iterations." << std::endl;
            
        }
        
        
        return Point(x_p,y_p,z_p);
    }
    
}

// pseudo distance, where the point is moved towards the origin until it crosses the ellipse
double Point::pseudoDistanceFromEllipse(double ex, double ez)
{
    double lambda = 1./(sqrt( pow(x/ex,2) + pow(y/ex,2) + pow(z/ez,2))  );
    
    // move the point back on the ellipse
    return Norm(Point(x,y,z)) * (1-lambda);
}


Point Point::projectOnYZCylinder(double x_radius)
{
    double l = std::sqrt(y*y + z*z);
    return Point(x,x_radius*y/l,x_radius*z/l);
}