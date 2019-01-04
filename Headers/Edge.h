#include<math.h>
#include "Point.h"
#include "Tissue.h"

#include <list>
#ifndef _Edge_h_
#define _Edge_h_

class Vertex;
class Cell;
class Tissue;

/// each object of that class contains two vertices and hence defines a vertical interface of the tissue
class Edge{
//private:
public:
	/// pointer to vertex 1
	Vertex *v1;
	/// pointer to vertex 2
	Vertex *v2;
	
	/// adjacent cell, which where the apical side of the edge is clockwise
	Cell *c1;
	
	/// adjacent cell, which where the apical side of the edge is counterclockwise
	Cell *c2;
	
	/// (optional) type of the edge - is set 0 by default 
	int Type;
	
	/// number that enables the unique identification of the edge
	int Number; 
	
	/// pointer to the tissue that contains the edge
	Tissue* T; 
	
	/// counter for the number of created edges
	static int NumberOfCreatedEdges;
	
    /// counter for the number of deleted edges
	static int NumberOfDeletedEdges;
	
	//static void resetCounter(){NumberOfCreatedEdges=0;NumberOfDeletedEdges=0;};
	
	/// line tension constant of the apical edge
	double G_a;
    
    /// line tension constant of the basal edge
	double G_b;
    
    /// anisotropic apical line tension along the z axis
    double G_a_ai;
    
    /// line tension decreases to zero from l0_a on
    double l0_a;
    
    /// line tension decreases to zero from l0_b on
    double l0_b;


    /// return the current line tension, which can be a function of the current length
    double Lambda_a();

    /// return the current line tension, which can be a function of the current length
    double Lambda_b();
    
    /// tension of the face inbetween the apical and the basal edge
	double T_l;
	
    /// lengths of the apical and basal bond
	double la, lb;
    
    /// update the currents lengths
    void UpdateLengths();
    
    double height();
    
public:
	/// std constructor
	Edge(Tissue *T);
	/// constructor that gets the two vertices 
	Edge(Tissue *T, Vertex *v1, Vertex *v2);
	/// constructor that gets the two vertices, the line tensions of the apical and basal connections and the surface tension of the interface
	Edge(Tissue *T, Vertex *v1, Vertex *v2, double G_a, double G_b, double T_l);
	/// complete constructor that gets the two vertices, the line tensions of the apical and basal connections, the surface tension of the interface and pointers to the connecting cells
	Edge(Tissue *T, Vertex *v1, Vertex *v2, Cell* c1, Cell* c2, double G_a, double G_b, double T_l) ;
	/// std destructor
	~Edge(){static int NumberOfDeletedEdges(0);NumberOfDeletedEdges++;}
	
    //static void resetCounter(){NumberOfCreatedEdges=0;}
	
	/// getter for pointer to vertex 1
	Vertex* getVertex1(){return v1;}
	/// getter for pointer to vertex 2
	Vertex* getVertex2(){return v2;}
	
	/// setter for pointer to vertex 1
	void setVertex1(Vertex* newVertex1){v1=newVertex1;}
	/// setter for pointer to vertex 2
	void setVertex2(Vertex* newVertex2){v2=newVertex2;}
	
	
	/// getter for pointer to tissue
	Tissue* getTissue(){return T;}
	/// setter for pointer to tissue
	void setTissue(Tissue* newTissue){T=newTissue;}
	
	/// setter for the number of the edge
	void setNumber(int newNumber){Number=newNumber;}
	/// getter for the number of the edge
	int getNumber(){return Number;}
	
	/// function returns if the face under the edge somewhere crosses the boundary
	bool CrossBoundary();
	
	
	/// getter for pointer of cell 1
    Cell* getCell1(){return c1;}
	/// setter for pointer of cell 1
	void setCell1(Cell* newCell1){c1=newCell1;}	
	/// getter for pointer of cell 2
	Cell* getCell2(){return c2;}
	/// setter for pointer of cell 2
	void setCell2(Cell* newCell2){c2=newCell2;}
	
	/// getter for stress
	Point getStress(){return Stress;}
	/// getter for energy
	//double getEnergy(){return LineEnergy+SurfaceEnergy;}	
	/// getter for surface energy
	//double getSurfaceEnergy(){return SurfaceEnergy;}	
	/// getter for line energy
	//double getLineEnergy(){return LineEnergy;}	
	
	/// equality operator overloaded
	bool operator==(const Edge &e);
	
	
	/// calculates the length of the apical edge
	double l_a();
	
	/// calculates the length of the basal edge
	double l_b();
	
	/// update the current line forces
	Point ApicalLineForce();
	Point BasalLineForce();

	
	/// update the present barycenter
	bool UpdateBarycenter();
	void UpdateForcesFromBarycenters();

    //void UpdateForcesAndEnergy();
    
    void update_forces();
    //double update_energy();
    
    double W();
    
    //void UpdateForcesAndEnergy_minArea();
    //void update_energy_minArea();

    
	double CalculateAreaInbetween();
	
	// the sum of the areas of the six triangles (four inside, one apical, one basal)
	double CalculateContributingArea();
	
	// returns 1 if v is v1, 2 if v is v2 and 0 else
	int VertexIsPart(Vertex* v);
	
	// print relevant information
	void Print();
	
	
	// in which quadrant does the second vertex of the edge lie? if periodicity in both directions:
	// (-1, 1)  (0, 1)  (1, 1)
	// (-1, 0)  (0, 0)  (1, 0)
	// (-1,-1)  (0,-1)  (1,-1)
	// default: 5
	// if only periodicity in x direction: (-1,0)  (0,0)  (1,0) = -1  0  1
	

	Point MP_a(); // returns the apical middle point
	Point MP_a(Cell* c); // return apical middle point with respect to the apical barycenter of the edge c
    Point MP_b(); // ... and the basal middle point
	Point MP_b(Cell* c); // return basal middle point with respect to the apical barycenter of the edge c
	
    double area(); // returns the surface area of the lateral edge
    
	//Point getApicalLineForce(){return ApicalLineForce;}
	//Point getBasalLineForce(){return BasalLineForce;}
	
	Point getBarycenter(){return BC_l;}
	
	/// getter of ForceFromBarycenter_Area
	Point getForceFromBarycenter_Area(){return ForceFromBarycenter_Area;}
	/// getter of dV_dBC
	Point getdV_dBC(){return dV_dBC;}
	
	/// removes the edge from the tissue and returns the vertex
	Vertex* T1_remove();
	Vertex* T1_remove_v2();

	/// checks if the edge should be removed from the tissue and removes it if neccessary
	bool CheckAndRemoveFromTissue();
	
	
	/// reset number of created edges
	static void resetNumberOfCreatedEdges(){NumberOfCreatedEdges=0;}
	
	//void getXcrossing(double xCross, Point* ApicalCross, bool* doesApicalCross, Point* BasalCross, bool* doesBasalCross, Point* SecondApicalCross, bool* isApicallyAligned, Point* SecondBasalCross, bool* isBasallyAligned);
	
	//void getYcrossing(double yCross, Point* ApicalCross, bool* doesApicalCross, Point* BasalCross, bool* doesBasalCross, Point* SecondApicalCross, bool* isApicallyAligned, Point* SecondBasalCross, bool* isBasallyAligned);
	
	//void updatePlaneSection();//Point positionVectorPlane, Point normalVectorPlane);
	
	//Point planeSection_apical, planeSection_basal;
	
	//bool planeCrosses_apical, planeCrosses_basal;
	
	
	// calculate the force coming from the pressure in cell c that acts on the T2 transition of Vertex v 
	Point Force_PressureT2Transition(Cell *c, Vertex *v);
	
	// calculate the force coming from Cell c that acts toward a T2 transition of Vertex v 
	Point Force_SurfaceTensionT2Transition(Cell *c, Vertex *v, Point shift);
	
//private:
	/// does the edge cross the boundary?
	bool crossBoundary;
	
    /// is the edge an outer edge of the tissue? (i.e. should it be drawn?)
    bool isAtBoundary();
    
	/// the stress that is imposed by the edge
	Point Stress;
	
	/// the energy contribution of the surface area to the energy of the system
	double Energy;

	/// coordinates of the barycenter of the interface
	Point BC_l;	
	Quadrant q_BC_l;
	
	Quadrant q_v2;
	
	/// current quadrant of the apical barycenter of c1/c2 with respect to the vertex, that is first in the direction of the cell c1/c2
	Quadrant q_c1;
	Quadrant q_c2;
		
	/// current force that is acting on v1.coord (towards v2.coord) due to line tension
	//Point ApicalLineForce;
	/// current force that is acting on v1.coord_b (towards v2.coord_b) due to line tension
	//Point BasalLineForce;

	
	/**	derivation of the area of the edge with respect to the barycenter divided by four times the surface tension 
	-> dA(R,B(R))/dR * T_l=(dA/dR + dA/dB*dB/dR)*T_l=(dA/dR + 1/4 * dA/dB)*T_l -> here we calculate 1/4 * dA/dB*T_l as this force contribution is the same for all points */
	Point ForceFromBarycenter_Area;
    Point ForceFromBarycenter_Area_numerical;
	
	/** derivation of the volume contribution of the edge with respect to the barycenter divided by four - dV(R,B(R))/dR=dV/dR + dV/dB*dB/dR=dV/dR + 1/4 * dV/dB 
	here we calculate 1/4 * dV/dB since this is the same for all the points of the cell and hence is a property of the cell */
	Point dV_dBC;
	Point dVc1_dBC;
	Point dVc2_dBC;
	
	Point Fa_v1, Fa_v2, Fb_v1, Fb_v2, b, v1a, v1b, v2a, v2b, F_barycenter_area, dV_db;
	Point b_v1a,b_v1b,b_v2a,b_v2b;
    
    std::list<Edge*>* neighbouringEdges();
    
    // return the surface and line tension energy contributions of the lateral interface and the apical and basal contributions to cell c
    double energy(Cell* c);

    // return the surface and line tension energy contributions of all cells that are attached to the edge
    double SurroundingEnergy();
    
    void setMechanicalProperties();
};
#endif
