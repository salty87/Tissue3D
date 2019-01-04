#ifndef _Cell_h_
#define _Cell_h_
#include <list>
#include <map>
#include <utility>
#include <math.h>
#include "Point.h"


struct CellType
{
	double V0;
	double T_a;
	double T_b;
	double K;
    double A0_a;
    double A0_b;
};

class Edge;
class Vertex;
class Tissue;

class Cell{
	
public:
		
	/* constructors */
	Cell();
	Cell(std::list<Edge*> ListEdges,double K, double Ta, double Tb, double V0);
	
	// Destructor
	~Cell(){ListEdges.clear();}
	
	
	/// the number of elements of class Cell that has been created
	static int NumberOfCreatedCells;
	/// the number of elements of class Cell that has been created
	static int NumberOfDeletedCells;
	
	//static void resetCounter() {NumberOfCreatedCells=0;NumberOfDeletedCells=0;};
		
	//std::list<Vertex*> ListVertex; //list of vertices forming the cell

	std::list<Edge*> ListEdges; //an ordered list that contains all the edges on the apical side of the cell (ordered means, that connections are as follows e1->e2->...->en->e1)
	
	
	// Tissue which contains the cell
	Tissue* T;
	
	// number of the cell, which is used to identify it in the tissue
	int Number;
	
	double K; //the elastic modulus of the cell regarding volume compression
	double V0;//preferred volume of the cell;
	double Volume; // current volume of the cell
    // stresses contains the stresses in x, y and z direction respectively - in cells just the stress coming from pressure plays a role
	Point VolumeStress;
    double ElasticEnergy;
	double Pressure;

	int Type; //number of the type of the cell
	
	// barycenters of all the faces
	Point BC_a;// barycenter of the apical face
	Point BC_b;// barycenter of the basal face
	Quadrant q_BC_b; // how is the basal barycenter located relative to the apical barycenter
	
	double Ta;   //surface tension of the apical face
	double Tb; //surface tension of the basal face
	
    // term that governs from which area the cells do not behave like under an active tension any longer, but rather exhibit elasticity around A0
    double A0_a;
    // term that governs from which area the cells do not behave like under an active tension any longer, but rather exhibit elasticity around A0
    double A0_b;

    /// tension at zero in terms of elasticity
    double T0_a;
    double T0_b;
    
    /// return the current apical surface tension, that might be dependent on the current apical surface area
    double T_a();
    /// return the current basal surface tension, that might be dependent on the current basal surface area
    double T_b();
    
    /// if Ta is a function of Aa, the derivative takes a somewhat more complicated form than only Ta
    double dW_dAa();
    /// if Tb is a function of Ab, the derivative takes a somewhat more complicated form than only Tb
    double dW_dAb();
    
	bool crossBoundary;
	bool CrossBoundary();
	
	double Perimeter_a;
	double Perimeter_b;
	
	//std::list<Vertex*> getIncludedVertices();
	
	Point dAa_dBCa;
    Point dAa_dBCa_numerical;
	Point ForceFromApicalBarycenter_Volume;
	
    Point dV_dBCa, dV_dBCb;
    
	Point dAb_dBCb;
	Point ForceFromBasalBarycenter_Volume;
	
	bool operator==(Cell c);
		
	/// get the direction of an edge 
	int DirectionOfEdge(Edge* e);
	
	/// calculate the contribution on the forces coming from the change in the position of the barycenter
	void UpdateForcesFromBarycenters();
	
	/// update the perimeter values that are stored in the object
	void updatePerimeter();

	/// calculate the apical and basal center of the cell using the edge 
	void UpdateCenter();

	/// calculate the current volume of the cell
	double UpdateVolume();
	//void UpdateVolume_test();
	
	/// returns the whole surface of the cell
	double CalculateSurfaceArea();
	
	/// returns the area of the apical surface
	double CalculateApicalSurfaceArea();
	
	/// calculate the area of the lateral side of the cells
	double CalculateLateralSurfaceArea();
	
	/// calculate the height of the cells, defined by the distance of the barycenters
    double height();
	
    /// angle between the cell and (0,0,1)
    double angleWithZPlane();
    
	/// calculate the apical perimeter of the cells
	double CalculateApicalPerimeter();
	
	/// calculate the basal perimeter of the cells
	double CalculateBasalPerimeter();
	
	/// returns the area of the basal surface
	double CalculateBasalSurfaceArea();
	
	/// calculates the pressure based on the current volume, the preferred volume and the bulk modulus K
	void UpdateElasticEnergy();
	
	/// update the energy contribution that stems from the cell 
	//bool UpdateVolumeEnergy();
	
	/// writes the value from CellType CT into the cell
	void assignCellType(int typeNo);
	
	void Print();
	

	
	//void addVertex(Vertex* v){if(find(ListVertex.begin(),ListVertex.end(),v)==ListVertex.end()) ListVertex.push_front(v);}
	void addEdge(Edge* e){if(find(ListEdges.begin(),ListEdges.end(),e)==ListEdges.end()) ListEdges.push_front(e);}
	
	
	static void resetNumberOfCreatedCells(){NumberOfCreatedCells=0;};

	//CellIntersection getXintersection(double x_sect);
	
	/// let cell undergo apoptosis
	Vertex* Apoptosis();
    
    /// cell divides
    void Division(Edge* e1, Edge* e2, Edge* e_new, double ratio_e1, double ratio_e2, Cell *c_new);

    /// cell divides randomly (calls Division with some random parameters)
    void RandomDivision();
    
	/// reorient all edges in the cell such that they are all oriented in counter clock wise direction to the cell (needed for simplified cell division)
    void reorientAllEdges();
    void orderEdges();

    
    /// returns 0 if the two cells dont have a common edge or 1 and the first common edge that has been found
	bool hasCommonEdgeWith(Cell * c, Edge ** e);
	    
    // returns sum of all energy contributions from volume elasticity of the cell and line and surface tensions of the surrounding edges
    double energy();
    
    /// contribution from this cell to the total energy: bulk elasticity, apical and basal surface contributions
    double W();
    
    void setMechanicalProperties();
    
    bool positionFixed;
};
#include "Edge.h"
#endif
