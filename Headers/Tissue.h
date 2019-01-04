#ifndef _Tissue_h_
#define _Tissue_h_
#include <math.h>
#include "Vertex.h"
#include "Edge.h"
#include "Cell.h"

#ifdef _USE_QT_
#include <QtGui>
#endif // _USE_QT_

#include <QString>
#include <list>
#include <map>
#include <vector>

#define NO_CELL_TYPES 4
//#define VERTEX_NUMBER 72
//#define CELL_NUMBER 36

enum BarycenterType {PointCenter, ContourCenter};

struct CellProperties
{
	double T_a;
    double A0_a;
	double T0_a;
    double T_b;
    double A0_b;
	double T0_b;
	double V0;
	double K;
};

struct TissueData
{
//	std::list< std::vector<double> > GeometricalData_type1;//[8]; // height - apical surface area - basal surface area - lateral surface area - volume - pressure - apical perimeter - basal perimeter
//	
//	std::list< std::vector<double> > Cell_averages_Types;//[10][8]; // 0. column is overall average
//	std::list< std::vector<double> > Cell_variances_Types;//[10][8]; // 0. column is overall variance

	// positions of apical and basal midpoints
    std::vector<double> BC_a_z[NO_CELL_TYPES+1]; // z values of the apical midpoints of the cells of type 1/2
    double mean_BC_a_z[NO_CELL_TYPES+1], max_BC_a_z[NO_CELL_TYPES+1], min_BC_a_z[NO_CELL_TYPES+1];
    std::vector<double> BC_b_z[NO_CELL_TYPES+1]; // z values of the basal midpoints of the cells of type 1/2
    double mean_BC_b_z[NO_CELL_TYPES+1], max_BC_b_z[NO_CELL_TYPES+1], min_BC_b_z[NO_CELL_TYPES+1];
    Point BC_a_midCell, BC_b_midCell;
    
    // heights of the cells
    std::vector<double> height[NO_CELL_TYPES+1]; // distance of the apical and basal midpoints of the cells of type 1/2
    double mean_height[NO_CELL_TYPES+1], max_height[NO_CELL_TYPES+1], min_height[NO_CELL_TYPES+1];
    
    // apical and basal areas
    std::vector<double> apicalArea[NO_CELL_TYPES+1]; // apical areas of the cells ordered by type
    double mean_apicalArea[NO_CELL_TYPES+1], max_apicalArea[NO_CELL_TYPES+1], min_apicalArea[NO_CELL_TYPES+1], var_apicalArea[NO_CELL_TYPES+1];
    std::vector<double> basalArea[NO_CELL_TYPES+1]; // basal areas of the cells of type 1/2
    double mean_basalArea[NO_CELL_TYPES+1], max_basalArea[NO_CELL_TYPES+1], min_basalArea[NO_CELL_TYPES+1];

    // angles of the cells
    std::vector<double> angle[NO_CELL_TYPES+1]; // angles between the cells and (0,0,1)
    double mean_angle[NO_CELL_TYPES+1], max_angle[NO_CELL_TYPES+1], min_angle[NO_CELL_TYPES+1];
    
    // number of cells
    unsigned numberOfCells[NO_CELL_TYPES+1];

    // volumes
    std::vector<double> volume[NO_CELL_TYPES+1];
    double mean_volume[NO_CELL_TYPES+1], min_volume[NO_CELL_TYPES+1], max_volume[NO_CELL_TYPES+1];
    
    // how many edges do the cell types have in common
    unsigned numberOfCommonEdges[NO_CELL_TYPES+1][NO_CELL_TYPES+1];

    // how big is the common interface (apical/basal line, lateral surface)
    double apicalCircumference[NO_CELL_TYPES+1][NO_CELL_TYPES+1]; // apical length of the common interface between cell types 1 and 2
    double basalCircumference[NO_CELL_TYPES+1][NO_CELL_TYPES+1]; // basal length of the common interface between cell types 1 and 2
    double lateralInterface[NO_CELL_TYPES+1][NO_CELL_TYPES+1]; // surface area of the interface between cell types 1 and 2    double max_height[NO_CELL_TYPES+1]; // highest z difference between apical and basal barycenter for cells of type 1

    // average z-positions of the circumferences between different cell types
    double averageApicalZheight[NO_CELL_TYPES+1][NO_CELL_TYPES+1];
    double averageBasalZheight[NO_CELL_TYPES+1][NO_CELL_TYPES+1];
    
    // 
    double apicalIndentation;
    
    double apicalIndentation_max[NO_CELL_TYPES+1]; // maximal z-distance of the apical boundary from the initial flat state
    double apicalIndentation_min[NO_CELL_TYPES+1]; // maximal z-distance of the apical boundary from the initial flat state

    double basalIndentation_max[NO_CELL_TYPES+1]; // maximal z-distance of the basal boundary from the initial flat state
    double basalIndentation_min[NO_CELL_TYPES+1]; // maximal z-distance of the basal boundary from the initial flat state
    
    double max_R[NO_CELL_TYPES+1];
    double min_R[NO_CELL_TYPES+1];
    
    Point midpoint;
    double Radius;
};


struct MechanicalTissueProperties
{
	CellProperties CellPropVector[NO_CELL_TYPES + 1];
	double Intercell_SurfaceTension[NO_CELL_TYPES + 1][NO_CELL_TYPES + 1];
	double Intercell_ApicalLineTension[NO_CELL_TYPES + 1][NO_CELL_TYPES + 1];
	double Intercell_ApicalLineTension_anisotropic[NO_CELL_TYPES + 1][NO_CELL_TYPES + 1];
	double Intercell_BasalLineTension[NO_CELL_TYPES + 1][NO_CELL_TYPES + 1];
	double Vertex_LineTension;
};


enum CurrentTab {VERTEX,EDGE,CELL,GENERAL};

#ifdef _USE_QT_
class Tissue : public QObject
{
    Q_OBJECT
#else
class Tissue
{
#endif
    
public:
    Tissue(); /// constructor
    //~Tissue(); /// destructor
    // Tissue operator=(const Tissue &T);  /// copy operator
	void clear(); /// delete all elements of the tissue
    
    BarycenterType BaryType; /// kind of barycenter calculation (PointCenter or ContourCenter)
    
    /*   TISSUE STRUCTURE INFORMATION   */
  	
    std::list<Vertex> ListVertex; /// list of all vertices
	std::list<Edge> ListEdge; /// list of all edges
	std::list<Cell> ListCell; /// list of all cells
	
    std::list<Cell*> pointerCellList(); /// return a list of pointers on the cells

    /*  MECHANICAL PROPERTIES */
	
    std::vector<CellType> CellTypes; /// vector of all existing defined cell types
	int number_of_celltypes; /// number of cell types
    //double IntercellSurfaceTensions[MAX_NO_CELL_TYPES][MAX_NO_CELL_TYPES+1]; /// surface tensions inbetween cells of different types
    //double IntercellApicalLineTensions[MAX_NO_CELL_TYPES][MAX_NO_CELL_TYPES+1]; /// apical line tensions inbetween cells of different types
    //double IntercellBasalLineTensions[MAX_NO_CELL_TYPES][MAX_NO_CELL_TYPES+1]; /// basal line tensions inbetween cells of different types
    bool IntercellSurfaceTensions_isset, IntercellApicalLineTensions_isset, IntercellBasalLineTensions_isset;
    MechanicalTissueProperties MechProp; /// struct that contains all mechanical properties of the tissue
    bool anisotropicLineTension;
    
    
    double springConstant_basementMembrane;
    
    /* VARIABLE MECHANICAL INFORMATION ON THE TISSUE LEVEL */
    
    Point dW_dL; /// derivative of the potential with respect to the system sizes
    //Point ext_dW_dL; /// externally applied force on the system size
    double T_ext;
    Point dW_dL_numerical; /// numerically obtained force on the system size (should equal dW_dL!)
    bool Lx_fix; /// fixation of the size in x direction
    bool Ly_fix; /// fixation of the size in x direction
	// Point Stress;
	double Energy; /// effective potential of the tissue
	double Energy_Surface; /// potential contribution from all surfaces of the tissue
    double Energy_Volume; /// potential contribution from the bulk modulus of the cells
    double Energy_Line; /// potential contribution from all edges in the system
    
    /* GEOMETRICAL INFORMATION ABOUT THE TISSUE */
    
    bool isPeriodic; /// periodicity of the tissue in x- or y-direction
	bool boundary_isCylinder; /// tissue is attached to a cylinder (vs apical and basal vertices can be attached to different circles)
	bool boundary_canDetach; /// tissue is bounded by a cylinder but can detach

	Quadrant SystemSize; /// system size for the case of periodic boundary conditions
    bool boundary_isSet; /// if the system is not periodic: is there a boundary around the tissue?
	double boundary_externaldW_dva; /// force applied to the apical 
	double boundary_externaldW_dvb; /// 
	double boundary_apicalRadius;
	double boundary_basalRadius;
    int midCell_number, midVertex_number;

    double BasalMembrane_SpringConstant;
    bool STIFF_BASAL_MEMBRANE;
    
    // ECM
    double flatECM_apicalPosition;
    double flatECM_basalPosition;
    bool flatECM_basalPosition_bothways;
    double flatECM_stiffness;
    // pre-factor for vertices that are not completely surrounded by cells of type 1
    double flatECM_stiffness_alpha_nonCell1;
    
    bool apical_wall, basal_wall; // used in set_pos to make sure that the apical/basal vertices do not cross the basal membrane at apical_wall_position/basal_wall_position
    double apical_wall_position, basal_wall_position;
    
    // tissue is confined into an ellipsoid x^2/a^2+y^2/b^2+z^2/c^2=1
    double ellipsoid_ex, ellipsoid_ey, ellipsoid_ez, cylindricalConfinement, ellipsoid_K;
    
    // save computational time if only apical sides of the cells are taken into account
    bool onlyApicalContributions;
    
    /*  FUNCTIONS THAT ARE NEEDED TO UPDATE THE TISSUE AND RELATED QUANTITIES */
    void UpdateLengths();
    void UpdateBarycenters();
    void UpdateVolumeOfCells();
	void UpdateElasticEnergyOfCells();
	void UpdateForcesFromBarycenters();
	void ResetVertexForces();
	void update();
    //bool update_wo_pressure();
    bool TestUpdate();
    //void UpdateForcesAndEnergy();
    void update_forces();
    //void UpdateForcesAndEnergy_wo_pressure();
    void update_forces_numeric();
    //double update_energy();
    double W();
    double W(bool update_BCs);
    //void UpdateForces_nonPeriodic_numeric()
	void ChangeLxAndLy(double newLx, double newLy);  // after adapting the system sizes, the coordinates of the vertices are scaled accordingly
	void reestablishPeriodicity(); 	/// after a movement of the vertices in the periodic tissue, the periodicity is checked and if necessary reestablished
	bool Changed; // carries the information if the tissue has been changed since the last update:
    void NumericalDifferentiation();
    
    /*  FUNCTIONS TO CHANGE THE TISSUE STRUCTURE */
	bool PopEdge(Edge *e);
   	void AddVertex(Vertex *v);
	void RemoveVertex(Vertex * v);
	void AddEdge(Edge* e){ListEdge.push_front(*e);ListEdge.front().setTissue(this);}
	void RemoveEdge(Edge *e);
	void AddCell(Cell *c);
	void RemoveCell(Cell *c);
	void RemoveShortEdges();
    Vertex* Apoptosis(Cell *c);
	//bool T1Transition(Edge* e);
	bool CheckForPop;
    bool PopAndUnpopAllEdges();
    double pop_minLength, pop_initLength;
    void T1_remove(Edge* e);
    void randomCellDivisions(unsigned N);
    void allCellsDivide();
    void divideCells(int cellType, int maxDivisions = -1);
    
    /* INITIALIZATION FUNCTIONS */
	bool createHexagonalTissue(int rows, int columns, double height, double length, double width);
	bool createPeriodicHexagonalTissue(int rows, int columns, double height, double length, double width);
	bool createRectangularTissue(int rows, int columns, double height, double length, double width);
	int setMechanicalProperties();
	bool SetTissueStructure();
	bool UnpopVertices();
	int testTissue(); // function that checks to some extents if the tissue is well defined
    void correct_orientations();
    
    /* FUNCTIONS TO READ IN AND READ OUT THE TISSUE DATA */
    Vertex* SimilarOrNewVertex(Point ApicalCoord, Point BasalCoord);
	Edge* SimilarOrNewEdge(Vertex* Ver1, Vertex* Ver2, bool* sameDirection, bool* newEdge);
	Edge* SimilarOrNewEdge(Vertex* Ver1, Vertex* Ver2, bool* sameDirection, bool* newEdge, Quadrant &quadrant_v2);
	void PrintEdges();
	void PrintCells();

    int loadFileFromTxt(const QString &fileName);
    
    int writeFileToTxt(const char* fileName);
    
    /* READ OUT ADDITIONAL INFORMATION FROM THE TISSUE */
    void updateData();
	void writeTissueData(const char *fileName);
    void extractDataFromDir(std::string dirName);
	TissueData tissueData;
    void UpdateExtPositions(); // in the non-periodic case this saves the outmost points of the tissue
    Point extPosition_bottomLeft, extPosition_topRight;
    
	/*  FUNCTIONS TO OBTAIN AND TO CHANGE THE CURRENTLY CHOSEN CELLS/VERTICES/EDGES */
    std::list<Vertex*> currentVertices;
    std::list<Cell*> currentCells;
    std::list<Edge*> currentEdges;
    Cell* closestCell(const Point &ClickedPoint, std::list<Cell*> cellList);
    Edge* closestEdge(const Point &ClickedPoint);
    Vertex* closestVertex(const Point &ClickedPoint);
    
    /* INCLUSION OF LUMEN PRESSURE */
    double V_Lumen_a, V0_Lumen_a, K_Lumen_a, P0_Lumen_a, P_Lumen_a;
    double V_Lumen_b, V0_Lumen_b, K_Lumen_b, P0_Lumen_b, P_Lumen_b;
    void update_V_Lumen_a();
    void update_V_Lumen_b();
    void update_P_Lumen_a();
    void update_P_Lumen_b();
    
    /* ESPECIALLY FOR GASTRULATION: AVERAGE BULK AND SHEAR MODULUS AND KAPPA, AS WELL AS PREFERRED AREA OF CELLS    */
    double bulkModulus_2D, shearModulus_2D, kappa_2D, preferredArea_2D, areaElasticity_2D, lambda_2D;
    
    
    /* EXTRACT AND WRITE TISSUE INFORMATION IN LIST FORM */
    void getCoordinatesAndForces(std::list<Point> *ApicalCoordinates, std::list<Point> *BasalCoordinates, std::list<Point> *dW_dvas, std::list<Point> *dW_dvbs);
    void getForces(std::list<Point> *dW_dvas, std::list<Point> *dW_dvbs);
    void setCoordinatesAndForces(std::list<Point> *ApicalCoordinates, std::list<Point> *BasalCoordinates, std::list<Point> *dW_dvas, std::list<Point> *dW_dvbs);
    void setCoordinates(std::list<Point> *ApicalCoordinates, std::list<Point> *BasalCoordinates);
    
    
    /* TEST THE TISSUE */
    bool TEST_FORCES;
    void applyPositionalNoise(double stdDev, unsigned distribution);
    double random_normal(); // create random number from N(0,1)
    void moveVertex(Vertex *v, Point apicalMovement, Point basalMovement);
    Point* VolumeDerivativesApical_numerically; // VolumeDerivatives_numerically[cell][vertex]
    Point* VolumeDerivativesBasal_numerically;

    Point* VolumeDerivativesApical_analytically; // VolumeDerivatives_numerically[cell][vertex]
    Point* VolumeDerivativesBasal_analytically;
    
    void ResetVolumeDerivatives();
    void UpdateVolumeDerivatives_numerical();
    
    Cell* CellFromNumber(int no);
    
    /* FUNCTIONS FOR THE GUI CONTROL */
    int mainWindow_tab;
    
    void addCurrentCell(double tx, double ty);
    void currentCell(double tx, double ty);
    void addCurrentVertex(double tx, double ty);
    void currentVertex(double tx, double ty);
    void addCurrentEdge(double tx, double ty);
    void currentEdge(double tx, double ty);
    
    void updateMainWindow();
        
    void update_pos(bool including_L);
    void update_der(bool including_L);
    void set_pos(bool including_L);
    void update_pos_and_der(bool including_L);
    
    std::vector<double> pos;
    std::vector<double> der;
    
    void change_size(double new_Lx, double new_Ly);

    Cell* get_random_cell();
    Cell* get_biggest_cell();
    
    static void vector_statistics(const std::vector<double> &vec, double &mean, double &var, double &max, double &min)
    {
        if(vec.size()==0) {mean = 0; var = 0; min = 0; max = 0; return;}
        
        // initiate min, mean and max
        max = vec[0];
        min = vec[0];
        mean = 0;
        for(unsigned i=0; i<vec.size(); i++)
        {
            if(vec[i]>max) max = vec[i];
            if(vec[i]<min) min = vec[i];
            mean+=vec[i];
        }
        
        mean /= vec.size();
        
        // second step: calculate the variance
        var=0;
        for(unsigned i=0; i<vec.size(); i++) var+=std::pow(vec[i]-mean,2);
        var/=vec.size();
    }
    
    // return the biggest number of Vertex/Edge/Cell of the tissue, in order to find a unique name for the next Vertex
    int maxVertexNumber();
    int maxEdgeNumber();
    int maxCellNumber();
    
    /// return the cell in the middle of the clone of type, but all the distance are Euclidean and the periodicity is not taken into account. hence the method makes mostly sense, when the clone is connected
    Cell* cell_in_middle(int type);
    
#ifdef _USE_QT_
public slots:
#endif //_USE_QT_
    /* THE CURRENT LIST<> HAS BEEN CHANGED */
    void slot_currentPoint(double, double);
    void slot_addCurrentPoint(double, double);
    // slots for signals from the main window, where some values have been changed
    void slot_setFirstCellNumber(int c_no){if(!currentCells.empty()) (*currentCells.begin())->Number=c_no; }
    void slot_setCellsPreferredVolume(double c_V0){if(!currentCells.empty()) (*currentCells.begin())->V0=c_V0; }
    void slot_setCellsK(double c_K){if(!currentCells.empty()) (*currentCells.begin())->K=c_K; }
    void slot_setCellsType(int c_type);
    void slot_firstCellApoptosis(){if(!currentCells.empty()) Apoptosis(*currentCells.begin()); currentCells.pop_front();}
    void slot_firstCellRandomDivision(){if(!currentCells.empty()) (*currentCells.begin())->RandomDivision(); }
    void slot_setFirstEdge_Number(int e_no){if(!currentEdges.empty()) (*currentEdges.begin())->Number=e_no; }
    void slot_setFirstEdge_Tl(double e_Tl){if(!currentEdges.empty()) (*currentEdges.begin())->T_l=e_Tl; }
    void slot_setFirstEdge_Ga(double e_Ga){if(!currentEdges.empty()) (*currentEdges.begin())->G_a=e_Ga; }
    void slot_setFirstEdge_Gb(double e_Gb){if(!currentEdges.empty()) (*currentEdges.begin())->G_b=e_Gb;}
    void slot_setFirstEdge_T1(){if(!currentEdges.empty()) T1_remove(*currentEdges.begin()); currentEdges.pop_front();}
    void slot_change_lumenV0_a(double new_V0_Lumen_b){V0_Lumen_a = new_V0_Lumen_b;}
    void slot_change_lumenK_a(double new_K_Lumen){K_Lumen_a = new_K_Lumen;}
    void slot_change_lumenP0_a(double new_P0_Lumen){P0_Lumen_a = new_P0_Lumen;}
    void slot_change_lumenV0_b(double new_V0_Lumen_b){V0_Lumen_b = new_V0_Lumen_b;}
    void slot_change_lumenK_b(double new_K_Lumen){K_Lumen_b = new_K_Lumen;}
    void slot_change_lumenP0_b(double new_P0_Lumen){P0_Lumen_b = new_P0_Lumen;}
    void slot_change_ellipsoid_ex(double new_ellipsoid_ex) {ellipsoid_ex = new_ellipsoid_ex;}// cylindricalConfinement = new_ellipsoid_ex;}
    void slot_change_ellipsoid_ey(double new_ellipsoid_ey) {}
    void slot_change_ellipsoid_ez(double new_ellipsoid_ez) {ellipsoid_ez = new_ellipsoid_ez;}
    void slot_change_ellipsoid_K(double new_ellipsoid_K) {ellipsoid_K = new_ellipsoid_K;}
    
    void slot_updateGastrulationProperties();
    
    void slot_setText(double T){T_ext=T;}
    //void slot_setExtFLy(double Fly){ext_dW_dL.y=Fly;}
    
    void slot_setMainWindowTab(int newTab){mainWindow_tab=newTab;}
    void slot_setFixLx(int x);
    void slot_setFixLy(int y);
    
    void slot_Ta(double f){MechProp.CellPropVector[1].T_a = f; MechProp.CellPropVector[2].T_a = f; setMechanicalProperties();}
    void slot_Tb(double f){MechProp.CellPropVector[1].T_b = f; MechProp.CellPropVector[2].T_b = f; setMechanicalProperties();}
    void slot_K(double f){MechProp.CellPropVector[1].K = f; MechProp.CellPropVector[2].K = f; setMechanicalProperties();}
    // void slot_V0(double f){MechProp.CellPropVector[1].V0 = f; MechProp.CellPropVector[2].V0 = f; setMechanicalProperties();}
    
    void slot_Ga(double f){MechProp.Intercell_ApicalLineTension[1][1] = f;MechProp.Intercell_ApicalLineTension[2][2] = f;setMechanicalProperties();}
    void slot_Gb(double f){MechProp.Intercell_BasalLineTension[1][1] = f;MechProp.Intercell_BasalLineTension[2][2] = f;setMechanicalProperties();}
    void slot_Tl(double f){MechProp.Intercell_SurfaceTension[1][1] = f;MechProp.Intercell_SurfaceTension[2][2] = f;setMechanicalProperties();}
    
    void slot_Ga_factor(double f){MechProp.Intercell_ApicalLineTension[2][1] = f*MechProp.Intercell_ApicalLineTension[1][1];setMechanicalProperties();}
    void slot_Gb_factor(double f){MechProp.Intercell_BasalLineTension[2][1] = f*MechProp.Intercell_BasalLineTension[1][1];setMechanicalProperties();}
    void slot_Tl_factor(double f){MechProp.Intercell_SurfaceTension[2][1] = f*MechProp.Intercell_SurfaceTension[1][1];setMechanicalProperties();}
    
    //void slot_A0(double f){cells_A0=f;}
    
    void slot_writeFileToTxt(const char* fileName){writeFileToTxt(fileName);}
    
    /// changes all number_of_closest_Cells closest cells to the tissues midpoint to type 2 & all other to type 1
    void mutate_cells_in_middle(int number_of_closest_Cells);
    /// changes all cells in radius to the midpoint of the tissue to type 2 & all other to type 1
    void mutate_cells_in_radius(double radius); // all cells whose apical BC is in radius to the middlepoint of the system become type 2, all others are type 1
    void mutate_cells_in_annulus(double r1, double r2); // all cells whose apical BC is on the annulus around the center with radii r1 and r2 will become cell type newType
    /// changes all number_of_closest_Cells closest cells to the current cells midpoint to type 2 & all other to type 1
    void mutate_cells_in_circle(int number_of_closest_Cells);
    /// all cells that are in distance to the midline of the tissue x = Lx/2 become type 1, others become type 2
    void mutate_cells_in_stripe(double distance); // all cell;s whose apical BC is close to the line x = Lx/2 are mutated
    /// all cells whose apical BC is close to the line x = Lx/3 become type 2, all cells whose apical BC is close to the line x = 2*Lx/3 are type 3
    void mutate_cells_in_two_stripes(double width); // all cells whose apical BC is close to the line x = Lx/3 are type 2, all cells whose apical BC is close to the line x = 2*Lx/3 are type 3
    
    /// for non-periodic tissues only: all cells whose midpoints (B_a+B_b)/2 is in the box are assigned the cell_type
    void mutate_cells_in_box(int cell_type, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max);

    bool isFlat();
    
    void createRoundHexagonalTissue(int radiusInCells, double height, double length, double MaximalRadius);
    
    void setBasalTensionGradient(double basalTensionLeft, double basalTensionRight);

    std::list<Cell*> extractCellsInPlane(Point positionVector_plane, Point normalVector_plane);
    
    // void slot_setCurrentZForce(double f){if(!currentVertices.empty()) {for(std::list<Vertex*>::iterator itv = currentVertices.begin(); itv!=currentVertices.end(); itv++) (*itv)->zForce = f;}}
    
    //void slot_currentVertex(int no);
    //void slot_addCurrentEdge(int no);
    //void slot_currentEdge(int no);
    
#ifdef _USE_QT_
signals:
    void signal_currentCell(int c_no, double c_V0, double c_V, int c_type, double c_K, double c_P, double c_T_a, double c_T_b, double c_A_a, double c_AR);
    void signal_currentVertex(int v_no, Point &v_Ra, Point &v_Rb, Point &v_Fa, Point &v_Fb, double G_l);
    void signal_currentEdge(int e_no, int c1, int c2, double e_T_l, double e_G_a, double e_G_b, double e_length_a, double e_length_b);
	void signal_currentTissue(int t_NoCells, double Lx, double Ly, double FLx, double FLy, double T_ext, bool fixLx, bool fix_Ly, double Energy);
    void signal_currentLumen_a(double new_V_Lumen);
    void signal_currentLumen_b(double new_V_Lumen);
    void signal_changeFlatECM(double flatECM_apicalPosition, double flatECM_basalPosition, double flatECM_stiffness);
    void signal_changeLumen_a(double V0, double K, double P0);
    void signal_changeLumen_b(double V0, double K, double P0);
    void signal_changeEllipsoid(double x, double y, double z, double K);
    void signal_updateGastrulationProperties(double preferredArea_2D,double areaElasticity_2D,double bulkModulus_2D,double shearModulus_2D,double kappa_2D);
    //void signal_update_cyst_data(double,double,double,double,double,double,double,double,double,double,double);
    void updateMainWindowData();
    
#else //_USE_QT_
    void signal_currentCell(int c_no, double c_V0, double c_V, int c_type, double c_K, double c_P, double c_T_a, double c_T_b, double c_A_a, double c_AR){};
    void signal_currentVertex(int v_no, Point &v_Ra, Point &v_Rb, Point &v_Fa, Point &v_Fb, double G_l){};
    void signal_currentEdge(int e_no, int c1, int c2, double e_T_l, double e_G_a, double e_G_b, double e_length_a, double e_length_b){};
	void signal_currentTissue(int t_NoCells, double Lx, double Ly, double FLx, double FLy, double T_ext, bool fixLx, bool fix_Ly, double Energy){};
    void signal_currentLumen(double new_V_Lumen){};
    //void signal_update_cyst_data(double,double,double,double,double,double,double,double,double,double){};
#endif

};

#include "Integrator.h"

#endif
