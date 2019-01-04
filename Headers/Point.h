#ifndef _Point_h_
#define _Point_h_
/// most basic geometrical class of the system that contains the information about a point/vector

#include <cmath>
#include <math.h>
#include <iostream>
#include <iomanip> // std::setpricision()
#include <algorithm>

#define MIN_LENGTH 10e-10
#define MIN_AREA 10e-10

template <class T> const T& min (const T& a, const T& b) 
{
	return !(b<a)?a:b;     // or: return !comp(b,a)?a:b; for version (2)
}

inline double mod(const double &p, const double &d)
{ return p - std::floor(p/d)*d ;}

class Point{
public:	
	/// x coordinate
	double x;
	/// y coordinate
	double y;
	/// z coordinate
	double z;

	/// std constructor
	Point();
	
	/// constructor that takes the (x,y,z)-coordinates of the point
	Point(const double &px, const double &py, const double &pz);
	
	/// the 2-norm of the vector between the point and the origin
	static double Norm(const Point &p)
	{
		return(sqrt(pow(p.x,2)+pow(p.y,2)+pow(p.z,2)));
	};

	/// 2-norm of the vector between this point and point p
	double Distance(const Point &p) const;
	
	
	/// scalar product of this point with point p
	double dot(const Point &p) const {return x*p.x + y*p.y + z * p.z;};
	
	/// returns the vector (cross) product of this point with point p
	Point VectorProduct(const Point &p) const;
	
	/// equality overloaded: returns 1 if this point has the same coordinates as point p and 0 otherwise
	bool operator == (const Point &p) const;
	
	/// increment overloaded
	Point operator += (const Point &p) const;
		
	/// multiplication with scalar overloaded: returns point with coordinates (l*x,l*y,l*z)
	Point operator * (const double &l) const;
	
	/// division by scalar overloaded: returns point with coordinates (x/l,y/l,z/l)
	Point operator / (const double &div) const;
	
	/// multiplication with point overloaded: returns the vector (cross) product of the two vectors
	Point operator * (const Point &p) const;
	
	/// cross product
	Point cross(const Point &p) const {return Point(y*p.z-z*p.y,z*p.x-x*p.z,x*p.y-y*p.x);};
	
	
	/// addition of points overloaded: returns point that is the sum of the two points
	Point operator + (const Point &p) const;
		
	/// subtraction of points overloaded: returns point that is the difference of the two points
	Point operator - (const Point &p) const;
	
	//Point& operator += (const Quadrant &q) {x += q.x; y += q.y; return *this;};
	Point& operator += (const Point &p) {x += p.x; y += p.y; z += p.z; return *this;};

	//Point& operator += (const Quadrant &q) {x += q.x; y += q.y; return *this;};
	Point& operator -= (const Point &p) {x -= p.x; y -= p.y; z -= p.z; return *this;};
	
	//Point operator + (const Quadrant &q) const { return Point(x+q.x,y+q.y,z); }

	//Quadrant quadrant_a, quadrant_b;
	
	/// the information if the vector is the zero-vector
	bool isZero();

	void Print();
	void Print_in_line();

	/// area of the triangle that is spanned between the three points p1,p2,p3
	static double TriangleArea(const Point &p1, const Point &p2, const Point &p3)
	{
		return(Norm((p2-p1)*(p3-p1))/2.);
	}
	
	/// area of the triangle that is spanned between the three points p1,p2 and the origin
	static double TriangleArea(const Point &p1, const Point &p2)
	{
		//sqrt(pow(p1.x,2) + p2.y**2 )
		
		return(Norm(p1*p2)/2);
	}
	
	/// normal of the triangle that is spanned between the three points p1,p2 and the origin
	static Point TriangleNormal(const Point &p1, const Point &p2)
	{
		Point CP = p1*p2;
		return(CP/Point::Norm(CP));
	}
	
	
	/// the elementwise modulus of the point p
	static Point Abs(const Point &p)
	{
		return Point(std::abs(p.x), std::abs(p.y), std::abs(p.z));
	}
	
	/// the derivation of the area of triangle (pi,pj,p0) with respect to the position of pi
	static Point TriangleAreaChange(const Point &pi, const Point &pj, const Point &p0) //derivation of the Area of the triangle (pi,pj,p0) with respect to vertex i
	{
		Point u = p0-pj;
		
        if(Norm((pi-pj)*(p0-pj)) < MIN_AREA) return Point(0,0,0); // avoid numerical instabilities due to division by small number
        
        Point y = pj + u*(u.dot(pi-pj)/u.dot(u));
        
//        Point cp =  (pi-p0)*(pj-p0);
//        Point w = pj - p0;
//        Point test = Point(-cp.y*w.z+cp.z*w.y,cp.x*w.z-cp.z*w.x,-cp.x*w.y+cp.y*w.x)/(2*Norm(cp));
        
        return (pi-y)*(Norm(p0-pj)/Norm(pi-y)*0.5);;
	}
	
	static Point TriangleAreaChange2(const Point &v1, const Point &v2, const Point &b, const double &n)
    {
        Point cp = (v1-b)*(v2-b);
        Point w = v2*((n-1)/n)-(b-v1/n);
        if(Norm(cp) < MIN_AREA) return Point(0,0,0); // avoid numerical instabilities due to division by small number
        else return Point(-cp.y*w.z+cp.z*w.y,cp.x*w.z-cp.z*w.x,-cp.x*w.y+cp.y*w.x)/(2*Norm(cp));
    };
	
    static Point TriangleAreaChange3(const Point &v1, const Point &v2, const Point &b, const double &n)
    {
        Point cp = (v1-b)*(v2-b);
        Point w = v2*((n-1)/n)-(b-v1/n);
        //w = Point(1,1,1);
        if(Norm(cp) < MIN_AREA) return Point(0,0,0); // avoid numerical instabilities due to division by small number
        else return Point(-cp.y*w.z+cp.z*w.y,cp.x*w.z-cp.z*w.x,-cp.x*w.y+cp.y*w.x)/(2*Norm(cp)) - TriangleAreaChange(b,v1,v2)/n;
    };
    
	/// the (signed) volume of the tetrahedron that is spanned between p1,p2,p3 and the origin, given by <(p3-p1)x(p3-p2),p3>/6
	static inline double TetrahedronSignedVolume(const Point &p1, const Point &p2, const Point &p3)
	{
		return(p3.dot(p1*p2)/6);
	}
	
	
//	static Point dV_dp(const Point &pj, const Point &pl)
//	{
//		return (pl*pj)/6;
//	}
//
    /// the derivation of the volume of the tetrahedron that is spanned between p1,p2,p3 and the origin with respect to p3 - is independent of p1!
    static inline Point dV_dp(const Point &p2, const Point &p3)
    {
        return (p2*p3)/6;
    }
    static inline Point dV_dp(const Point &p1, const Point &p2, const Point &p3)
    {
        return (p2*p3)/6;
    }

    
    static inline Point dV_dp_numerical(Point &p1, Point &p2, Point &p3)
    {
        double V_plusx, V_minusx, V_plusy, V_minusy, V_plusz, V_minusz;
        double d = 10e-5, d_2 = d/2.;
        Point P_x_d_2(d_2,0,0),P_x_d(-d,0,0),P_y_d_2(0,d_2,0),P_y_d(0,-d,0),P_z_d_2(0,0,d_2),P_z_d(0,0,-d),P_null(0,0,0);
        
        V_plusx = TetrahedronSignedVolume(p1+P_x_d_2,p2,p3);
        V_minusx = TetrahedronSignedVolume(p1-P_x_d_2,p2,p3);
        
        V_plusy = TetrahedronSignedVolume(p1+P_y_d_2,p2,p3);
        V_minusy = TetrahedronSignedVolume(p1-P_y_d_2,p2,p3);
        
        V_plusz = TetrahedronSignedVolume(p1+P_z_d_2,p2,p3);
        V_minusz = TetrahedronSignedVolume(p1-P_z_d_2,p2,p3);
        
        return Point(V_plusx-V_minusx, V_plusy-V_minusy, V_plusz-V_minusz)/d;
    }

    static inline Point dV_dp_b(Point &p1, Point &p2, Point &b, double n)
    {
        Point cp = (p2-p1)*(b-p1);
        Point w = p2*(1.-1./n)-(b - p1/n);
        
        Point ret;
        ret.x = ( cp.x + p1.dot(Point(1,0,0)*w)) / 6.;
        ret.y = ( cp.y + p1.dot(Point(0,1,0)*w)) / 6.;
        ret.z = ( cp.z + p1.dot(Point(0,0,1)*w)) / 6.;
        
        return ret;
    }
    
    static inline Point dV_dp_b_1(Point &p1, Point &p2, Point &b, double n)
    {        
        Point cp = (p2-p1)*(b-p1);
        Point w = p2*(1.-1./n)-(b - p1/n);
        
        Point ret;
        ret.x = ( cp.x + p1.dot(Point(1,0,0)*w)) / 6.;
        ret.y = ( cp.y + p1.dot(Point(0,1,0)*w)) / 6.;
        ret.z = ( cp.z + p1.dot(Point(0,0,1)*w)) / 6.;
        
        // subtract contribution from the same triangle that has been added without barycenter contribution
        ret -= dV_dp(p1,p2)/n;
        
        return ret;
    }
    
    static inline Point dV_dp_b_numerical(Point &p1, Point &p2, Point &b, int n)
    {
        double V_plusx, V_minusx, V_plusy, V_minusy, V_plusz, V_minusz;
        double d = 10e-5, d_2 = d/2.;
        Point P_x_d_2(d_2,0,0),P_x_d(-d,0,0),P_y_d_2(0,d_2,0),P_y_d(0,-d,0),P_z_d_2(0,0,d_2),P_z_d(0,0,-d),P_null(0,0,0);

        V_plusx = TetrahedronSignedVolume(p1+P_x_d_2,p2,b+P_x_d_2/n);
        V_minusx = TetrahedronSignedVolume(p1-P_x_d_2,p2,b-P_x_d_2/n);
        
        V_plusy = TetrahedronSignedVolume(p1+P_y_d_2,p2,b+P_y_d_2/n);
        V_minusy = TetrahedronSignedVolume(p1-P_y_d_2,p2,b-P_y_d_2/n);
        
        V_plusz = TetrahedronSignedVolume(p1+P_z_d_2,p2,b+P_z_d_2/n);
        V_minusz = TetrahedronSignedVolume(p1-P_z_d_2,p2,b-P_z_d_2/n);
        
        return Point(V_plusx-V_minusx, V_plusy-V_minusy, V_plusz-V_minusz)/d;
    }
        
    /// the derivation of the distance between p0 and p1 with respect to p1
	static Point LineChange(const Point &p0, const Point &p1) 
	{		
		double dx, dy, dz;
		
		dx=p1.x-p0.x;
		dy=p1.y-p0.y;
		dz=p1.z-p0.z;
		
		return(Point(dx,dy,dz)/(sqrt(pow(dx, 2)+pow(dy,2)+pow(dz,2))));
	}
	
	/// the elementwise product of the two points
	static Point ElementwiseProduct(const Point &p0, const Point &p1)
	{
		return Point(p0.x*p1.x,p0.y*p1.y,p0.z*p1.z);
	}
	
	/// the areas of the projection of the triangle between p1, p2 and the origin on the three planes x - (y,z), y - (x,z), z - (x,y)
	static Point Projection3D(const Point &p1, const Point &p2)
	{
		// calculates the projection of the triangles given by p1 and p2 on the 3 axes
		// get the cross product of the triangle, which contains information about the normal and the surface area
		Point CrossProduct=p1*p2;
		
		double Area=Point::TriangleArea(p1,p2);
		
		// normal to the triangle
		Point Normal=CrossProduct/Area*2;
		
		// the force per unit exerted in the x(,y,z) direction is given by (1-abs(Normal_x)), which then has to be multiplied by the area to get the integrated tension		
		return (Point(1-std::abs(Normal.x),1-std::abs(Normal.y),1-std::abs(Normal.z))*Area);
	}
	
	static Point surfaceStress(const Point &R1,const  Point &R2,const  Point &R3)
	{
		// calculate D
		Point D = (R1-R2)*(R1-R3);
		if(Point::Norm(D)<10e-14) return Point(0,0,0); // avoid numerical instabilities by division by zero
        else return Point((pow(D.y,2)+pow(D.z,2)),(pow(D.x,2)+pow(D.z,2)),0)/(2*Point::Norm(D));
	}
	
	static Point lineStress(const Point &R1,const  Point &R2)
	{
		Point d = R1-R2;
        
        return Point(pow(d.x,2), pow(d.y,2),0)/Point::Norm(d);
	}
	
	/// the elementwise arccos
	static Point Arccos(const Point &p)
	{
		return Point(acos(p.x),acos(p.y),acos(p.z));
	}
	
	static bool getXsectionWithLine(const Point &P1, const Point &P2, const double &xCross, Point &section)
	{
		double CrossCheck = (P1.x-xCross)*(P2.x-xCross);
		
		if(CrossCheck<0)
		{
			double ratio = (xCross - P1.x)/(P2.x-P1.x);
			section = Point(xCross, P1.y + ratio * (P2.y-P1.y), P1.z + ratio * (P2.z-P1.z));
			return true;			
		}
		
		return false;
	}	
	
	static bool getYsectionWithLine(const Point &P1, const Point &P2, const double &yCross, Point &section)
	{
		double CrossCheck = (P1.y-yCross)*(P2.y-yCross);
		
		if(CrossCheck<0)
		{
			double ratio = (yCross - P1.y)/(P2.y-P1.y);
			section = Point(P1.x + ratio * (P2.x-P1.x), yCross,  P1.z + ratio * (P2.z-P1.z));
			return true;			
		}
		
		return false;
	}
	
	// calculates the intersection of triangle (A,B,C) with the plane given by its normal and the position vector
	static bool TrianglePlaneSection(const Point& A , const Point& B , const Point& C , const Point& NormalVector_plane , const Point& PositionVector_plane, Point& P1, Point& P2, bool& liesInPlane, bool& crossesAB)
	{
		// shift the triangle by the negative position vector
		Point A1 = A - PositionVector_plane;
		Point B1 = B - PositionVector_plane;
		Point C1 = C - PositionVector_plane;
		
		double d;
		
		// bool to indicate if an intersection has taken place
		bool found = false;
		liesInPlane=false; // also use this to indicate if AB lies inside the plane
		crossesAB = false;
		
		// check if the line AB crosses the plane
		// l = direction vector of the edge
		Point l = A1 - B1;
		
		double enumerator = B1.dot(NormalVector_plane);
		double denominator = l.dot(NormalVector_plane);
		
		if(denominator==0 && enumerator!=0)
		{
			// line is parallel to the plane but not inside the plane
		}
		else if(denominator==0 && enumerator==0)
		{
			// line lies inside the plane
			liesInPlane = true;
			P1 = A;
			P2 = B;
		}
		else 
		{
			d = -enumerator/denominator;
			if(d>=0 && d<=1)
			{
				P1 = B + l*d;
				found = true;
				crossesAB = true;
			}
		}
		
		// check if the line BC crosses the plane
		// l = direction vector of the edge
		l = B1 - C1;
		
		enumerator = C1.dot(NormalVector_plane);
		denominator = l.dot(NormalVector_plane);
		
		if(denominator==0 && enumerator!=0)
		{
			// line is parallel to the plane but not inside the plane
		}
		else if(denominator==0 && enumerator==0)
		{
			// line lies inside the plane
			// if already AB lied in the plane, then the whole triangle lies in the plane
			if(liesInPlane)
			{
				return true;
			}
			else 
			{
				P1 = B;
				P2 = C;
				return true;
			}
			
			
		}
		else 
		{
			d = -enumerator/denominator;
			if(d>=0 && d<=1)
			{
				if(found)
				{
					P2 = C+ l*d;
					return true;
				}
				
				P1 = C + l*d;
				found = true;
			} 
		}
		
		// check if the plane intersects with line CA
		// l = direction vector of the edge
		l = C1 - A1;
		
		enumerator = A1.dot(NormalVector_plane);
		denominator = l.dot(NormalVector_plane);
		
		if(denominator==0 && enumerator!=0)
		{
			// line is parallel to the plane but not inside the plane
		}
		else if(denominator==0 && enumerator==0)
		{
			P1 = C;
			P2 = A;
		
			return true;
			
		}
		else 
		{
			d = -enumerator/denominator;
			if(d>=0 && d<=1)
			{
				P2 = A + l*d;
				return true;
			}
		}
		
		return false;
		
		
	}
    
    void rotate(double phi){ // rotate vector v around the z axis by angle phi
        double x_old(x);
        x=std::cos(phi)*x-std::sin(phi)*y;
        y=std::sin(phi)*x_old + std::cos(phi)*y;
    }
    
    /// projection of this point on the ellipse given by X^2/a^2 + y^2/b^2 + z^2/c^2=1
    double pseudoDistanceFromEllipse(double ex, double ez);
    Point projectOnEllipse(double ex, double ez);
    double inline newtonRaphson(double t_n, double ex, double ez, double x1, double y1, double z1)
    {
        return t_n + 0.5* (  std::pow(ex,2)*(std::pow(x1,2)+std::pow(y1,2))/std::pow(std::pow(ex,2)+t_n,2) +  std::pow(ez*z1,2)/std::pow(std::pow(ez,2)+t_n,2) - 1  ) / (      std::pow(ex,2)*(std::pow(x1,2)+std::pow(y1,2))/std::pow(std::pow(ex,2)+t_n,3)     +        std::pow(ez*z1,2)/std::pow(std::pow(ez,2)+t_n,3)              );
    };
    
	Point projectOnYZCylinder(double x_radius);
};

// the quadrant class contains the information about the position of the vertex if the tissue is periodic
class Quadrant
{
public:	
	// constructors
	Quadrant(const double &x,const double &y): x(x), y(y){};
	
	Quadrant() : x(0), y(0){};
	
	// operators
	Quadrant operator*(const double &l) const {return Quadrant(l*x,l*y);}
	
	Quadrant operator*(const Quadrant &q) const {return Quadrant(q.x*x,q.y*y);}
	
	Quadrant operator+(const Quadrant &q) const {return Quadrant(q.x+x,q.y+y);}
	Quadrant operator-(const Quadrant &q) const {return Quadrant(x-q.x,y-q.y);}

	void operator+=(const Quadrant &q) {x+=q.x;y+=q.y;}
	void operator-=(const Quadrant &q) {x-=q.x;y-=q.y;}
	
	bool operator==(const Quadrant &q) const {return x==q.x && y==q.y;}
	
	void operator=(const Quadrant &q) {x = q.x; y=q.y;}
	
	bool isZero() const {return x==0 && y==0;}
	
	//Point operator+(const Point &p) const {return Point(p.x+x, p.y+y, p.z);}
	
	// members
	double x, y;
	
	/// point p to the modulus mod_p, i.e. it returns (p.x % mod_p.x , p.y 5 mod_p.y, p.z % mod_p.z)
	static Point Modulus(const Point &p, const Quadrant &mod_p)
	{
		return Point(mod(p.x,mod_p.x),mod(p.y,mod_p.y),p.z);
	}
	
	/// the 2-norm of the vector between the point and the origin for a periodic tissue with size SystemSize
	static double Norm(const Point &p,const Quadrant &SystemSize)
	{
		double dx = std::min(pow(p.x,2),pow(p.x-SystemSize.x,2));
		double dy = std::min(pow(p.y,2),pow(p.y-SystemSize.y,2));
		return(sqrt(dx+dy+pow(p.z,2)));
	};
    
    void Print()
    {
        std::cout<<"("<<x<<","<<y<<")\n";
    };
    
    void rotate(double phi){ // rotate vector v around the z axis by angle phi
        double x_old(x);
        x=std::cos(phi)*x-std::sin(phi)*y;
        y=std::sin(phi)*x_old + std::cos(phi)*y;
    };
	
	
};



static Point shift(const Point &p, const Quadrant &q) {return Point(p.x+q.x, p.y+q.y,p.z);}


#endif

