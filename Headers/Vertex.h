#ifndef _Vertex_h_
#define _Vertex_h_
#include <list>
#include <vector>
#include <deque>
#include <map>
#include <iostream>
#include "Point.h"


class Tissue;
class Edge;
class Cell;

struct ForceFromCell
{
	Point ForceFromCell_apical;
	Point ForceFromCell_basal;
};

/// represents two topological identifiable points in the tissue, one on the apical and one on the basal side, and also the information about the connection of those two
class Vertex{
public:		
	/// coordinates of the apical point
	Point coord;
	
	/// coordinates of the basal point
	Point coord_b;
	
	/// the line tension of the lateral bond between the two points
	double G_l;
	
	Quadrant q_b;
	
	/// (optional) information about the type of the vertex - is set 0 by default. could e.g. set the line tension 
	int Type;
	
	/// unique number of the vertex, that allows for numbering of in- and output 
	int Number;
	
	/// counter for created vertices
	static int NumberOfCreatedVertices;
	
	/// counter for deleted vertices
	static int NumberOfDeletedVertices;
	
	//static void resetCounter(){NumberOfCreatedVertices=0;NumberOfDeletedVertices=0;};
	
	/// pointer to the tissue that contains the vertex
	Tissue* T;
	
	
//	/// list of pointers to the edges that are connected to the vertex, i.e. the edges of which the vertex is part
//	std::list<Edge*> ListEdges;

	
	
	/// the stress that is caused by the edge inbetween due to the line tension
	Point Stress;
	
	/// the energy contribution of the edge inbetween to the system's energy
	double Energy;
	
	/// vector that contains the force currently acting on the apical point of the vertex
	Point dW_dva;
	Point dW_dva_numerical;
    
	/// vector that contains the force currently acting on the basal point of the vertexw
	Point dW_dvb;
	Point dW_dvb_numerical;
    
    /// externally applied apical force
    Point ExternalForce_a;
    
    /// externally applied basal force
    Point ExternalForce_b;
    
	/// double ended queue where the coordinates of the previous timesteps are saved
	std::deque<Point> ApicalHistory, BasalHistory;
	
    /// current length of the lateral bond
    double l_l;
    void UpdateLength(){l_l = getLength();}
    
public:
	/// list of pointers to the edges that are connected to the vertex -> does not contain crucial information about the tissue but is kept to save calculation time
	std::list<Edge*> NeighbouringEdges_sorted;
    
    /// updates the neighbouring edges -> only needs to be done at the first appearance of the tissue/vertex and after T1 transitions, apoptosis and divisions
	void updateNeighbouringEdges();

	/// std constructor
	Vertex(Tissue *T);
	
	/// constructor that takes the coordinates of the apical and basal point
	Vertex(Tissue *T, Point coord, Point coord_b);
	
	/// constructor that receives the coordinates and the list of of pointers of the connected edges
	Vertex(Tissue *T, Point coord_a, Point coord_b, std::list<Edge*> NeighbouringEdges_sorted);
	
	/// std destructor
	~Vertex(){static int NumberOfDeletedVertices(0);NumberOfDeletedVertices++;}
	
	/// setter for the apical coordinates
	void setCoord(Point newApicalPoint){coord=newApicalPoint;}
	
	/// setter for the basal coordinates
	void setBasalCoord(Point newBasalPoint){coord_b=newBasalPoint;}
	
	/// getter for the apical coordinates
	Point getCoord() const {return coord;}
	
	/// getter for the basal coordinates
	Point getBasalCoord() const {return coord_b;}
	
	void setdW_dva(Point f_apical){dW_dva=f_apical;}
	void setdW_dvb(Point f_basal){dW_dvb=f_basal;}
	void resetForces(){dW_dva=Point();dW_dvb=Point();}
	
	void update_forces();
    double W();

	double getLength(){return Point::Norm(coord - coord_b);}
	
	/// set quadrant of the basal point
	void setBasalQuadrant(Quadrant newBasalQuadrant){q_b=newBasalQuadrant;}
	
	/// get the point that describes the quadrant of the basal point relative to the apical point
	Quadrant getBasalQuadrant(){return q_b;}

	
	/// set the type of the vertex
	void setType(int newType){Type=newType;}
	
	/// get the type of the vertex
	int getType(){return Type;}
	
	/// get the number of the vertex
	int getNumber(){return Number;}
	
	/// get the pointer to the tissue
	Tissue* getTissue(){return T;}

	/// set the pointer to the tissue
	void setTissue(Tissue* newTissue){ T=newTissue; }
	
	/// set the number of the vertex
	void setNumber(int newNumber){Number=newNumber;}
	
	/// equality operator overloaded: returns if the other vertices has the same position - seems redundant, if no threshold is included
	bool operator==(const Vertex& v);
	
	/// the information if the apical point can be moved in x-direction (==0) or if it is fixed (==1)
	bool xafixed;
	
	/// the information if the apical point can be moved in y-direction (==0) or if it is fixed (==1)
	bool yafixed;
	
	/// the information if the apical point can be moved in z-direction (==0) or if it is fixed (==1)
	bool zafixed;
	
	/// the information if the basal point can be moved in x-direction (==0) or if it is fixed (==1)
	bool xbfixed;
	
	/// the information if the basal point can be moved in y-direction (==0) or if it is fixed (==1)
	bool ybfixed;	
	
	/// the information if the basal point can be moved in z-direction (==0) or if it is fixed (==1)
	bool zbfixed;
	
	/// function returns if the line under the vertex somewhere crosses the boundary
	bool CrossBoundary();
	
	/// get the stress on the line inbetween
	Point getStress(){return Stress;};
	
	/// get the energy contribution of the edge
	double getEnergy(){return Energy;}
	
	
	/// returns the force acting on the basal vertex, due to the line tension between the apical and the basal vertices
	Point lineForceOnApicalPoint();
		
	/// updates the stress that is currently caused by the vertex
	Point UpdateLineStress();
	
	/// updates the energy contribution of the vertex to the energy of the system
	//void update_energy();
	
	/// print the properties of the vertex in the console
	void Print();
	
	/// return the current force that acts on the apical vertex
	Point getdW_dva(){return dW_dva;}
	
	/// return the current force that acts on the basal vertex
	Point getdW_dvb(){return dW_dvb;}
	
	/// iterates the vertex with the step size mobility in the periodic case
	double movePoints(double mobility, double &Lx_new, double &Ly_new);
	
	/// iterates the vertex with the step size mobility
	void movePoints(double mobility){coord=coord+dW_dva*mobility;coord_b=coord_b+dW_dvb*mobility;}
	
    /// change the vertex' apical and basal positions according to the movement, taling into accoutn possible boundary restrictions
    Quadrant moveVertex(Point apicalMovement, Point basalMovement);
	
	/// add a change to the two points, e.g. during the iteration of the system
	void movePoints(Point movementApicalPoint, Point movementBasalPoint){coord=coord+movementApicalPoint;coord_b=coord_b+movementBasalPoint;}
	
	/// reset number of created vertices
	static void resetNumberOfCreatedVertices(){NumberOfCreatedVertices=0;}
	
	void clearHistory(){ApicalHistory.clear();BasalHistory.clear();}
	
	/// takes an edge, draws a plane with normal direction through the vertex, checks which edges lie on which side and then creates a new edge at (coord-direction/2, coord + direction/2)	
	int UnpopEdge(Point direction, int* flag, Edge* e);

	Point planeSection;
	double planeCrosses;
	
	/// a map that saves for every neighbouring cell of the vertex the projected force that the cell exerts on the vertex in the plane orthogonal to the lateral edge
	std::map<Cell*,Point> ForcesFromCells;
	
    // split and compare energies
	bool checkAndPop_energy(double edgeLength);
	// split and check for separation forces
    bool checkAndPop_forces_original(double edgeLength);
    bool checkAndPop_forces(double edgeLength);
    
	int ForcesFromPop(Cell* c1, Cell* c2, Point* popVector, std::list<Cell*> *CellGroup1,std::list<Cell*> *CellGroup2, double SurfaceTensionNewEdge, double G_lNewEdge);

	int pop(Cell* c1, Cell* c2, std::list<Cell*> CellGroup1, std::list<Cell*> CellGroup2, Point popVector, Edge* e);
	
	// algorithm that returns if the vertex is at the boundary
	bool get_isAtBoundary();
	bool isAtBoundary;
    
    bool inBoundaries(Quadrant Q);
	
	// values that are used at every iteration step, to save the time of the initialization
	//Point LineForce_apical, LineForce_basal, SurfaceForce_apical, SurfaceForce_basal, VolumeForce_apical, VolumeForce_basal, ApicalVolume_Barycenter, BasalVolume_Barycenter, ForceCell1, ForceCell2;
	//Point P1, P2, P, P_b, B, P0, P0_b;
	//Point apicalDer_x, apicalDer_y, apicalDer_z, basalDer_x, basalDer_y, basalDer_z;
	Point dP_dBa,dP_dBb, ForceFromBarycenter_a,ForceFromBarycenter_b;
	Point normal;
    
    std::list<Cell*> getNeighbouringCells();
    
    
    /// a function that returns true, iff all the neighbouring cells are of type 1
    bool allNeighboursType1();
    
    double SurroundingEnergy();

    bool neighbouringEdgeCrossesBoundary();
    
   // std::list<Edge*>* neighbouringEdges();
    
};
#include "Edge.h"
#endif