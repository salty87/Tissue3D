#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <time.h>
//#include <random>
#include "Tissue.h"
#include "Edge.h"
#include <time.h>
#include <math.h>
#include <QFile>
#include <QTextStream>
#include <cstdlib>

//#include "QsLog.h"
//#include "QsLogDest.h"
#define PI 3.1415926535897932384

Tissue::Tissue()
{
	CheckForPop=true;
	
	Changed=true;
	isPeriodic=false;
    
	IntercellSurfaceTensions_isset = false;
	IntercellApicalLineTensions_isset = false;
	IntercellBasalLineTensions_isset = false;
    onlyApicalContributions = false;
    
    // initialize the list of selected cells, vertices, edges
    currentCells=std::list<Cell*>();
    currentVertices=std::list<Vertex*>();
    currentEdges=std::list<Edge*>();
    
    // initialize the lists
    ListVertex=std::list<Vertex>();
    ListEdge=std::list<Edge>();
    ListCell=std::list<Cell>();
    
    // set the way, the barycenter is calculated
    BaryType = PointCenter;

    // allow the system sizes to change
    Lx_fix = false;
    Ly_fix = false;
        
    // no external force on the system size
    T_ext = 0;
    
    // no spring like behaviour of the basal vertices around z = 0
    springConstant_basementMembrane = 0;
    
    // ECM which forbids any movement of basal vertices in z direction
    STIFF_BASAL_MEMBRANE = false;
        
    // basal wall which forbids vertices to move through the plane z = basal_wall_position
    basal_wall = false;
    basal_wall_position = 0;
    apical_wall = false;
    apical_wall_position = 0;
    
    anisotropicLineTension = false;
    
    flatECM_basalPosition_bothways = false;
    
    // initiate the tissue's mechanical properties
    MechProp.CellPropVector[1].V0=      0;  //V0=0 -> is calculated to equal the initial volume of the cells
    MechProp.CellPropVector[1].K=       0.01;
    MechProp.CellPropVector[1].T_a =    0.5;
    MechProp.CellPropVector[1].T_b =    2;
    MechProp.CellPropVector[1].T0_a =   -100;
    MechProp.CellPropVector[1].T0_b =   -100;
    MechProp.CellPropVector[1].A0_a =   50;
    MechProp.CellPropVector[1].A0_b =   50;
    
    for(unsigned i = 2; i<=NO_CELL_TYPES+1; i++) MechProp.CellPropVector[i] = MechProp.CellPropVector[1];
    
    for(unsigned i = 0; i<=NO_CELL_TYPES; i++) {
        for(unsigned j = 0; j<=NO_CELL_TYPES; j++) {
            if(i>=j) {
                MechProp.Intercell_SurfaceTension[i][j] = 1;
                MechProp.Intercell_ApicalLineTension[i][j] = 50;
                MechProp.Intercell_ApicalLineTension_anisotropic[i][j] = 0;
                MechProp.Intercell_BasalLineTension[i][j] = 50;
            } else {
                MechProp.Intercell_SurfaceTension[i][j] = 0;
                MechProp.Intercell_ApicalLineTension[i][j] = 0;
                MechProp.Intercell_ApicalLineTension_anisotropic[i][j] = 0;
                MechProp.Intercell_BasalLineTension[i][j] = 0;
            }
        }
    }

    MechProp.Vertex_LineTension = 0;

    flatECM_stiffness = 0;
    flatECM_apicalPosition = 0;
    flatECM_basalPosition = 0;
    flatECM_stiffness_alpha_nonCell1 = 1;
    
    T_ext = 0;
    
    // initialize the lumen pressure constants
    V0_Lumen_a = 0;
    K_Lumen_a  = 0;
    P0_Lumen_a = 0;
    V0_Lumen_b = 0;
    K_Lumen_b  = 0;
    P0_Lumen_b = 0;
    
    pop_minLength = 1;
    pop_initLength = 1;
}

// take the current state of the system (i.e. positions of the vertices) and calculate the values, that are neccessary for the calculation of the evolution
void Tissue::update()
{
    UpdateLengths();
	UpdateBarycenters();
	if(!onlyApicalContributions) UpdateElasticEnergyOfCells();
    if(K_Lumen_a!=0 || P0_Lumen_a!=0) update_P_Lumen_a();
    else P_Lumen_a = 0;
    
    if(!onlyApicalContributions && (K_Lumen_b!=0 || P0_Lumen_b!=0)) update_P_Lumen_b();
    else P_Lumen_b = 0;
    
    W();
    update_forces();
    
    
    //TestUpdate();

	Changed = false;
}

void Tissue::slot_updateGastrulationProperties()
{
    // extract properties of cells' apical sides
    double T_a = MechProp.CellPropVector[1].T_a;
    double A0_a = MechProp.CellPropVector[1].A0_a;
    double T0_a = MechProp.CellPropVector[1].T0_a;
    
    // extract cell densities
    updateData();
    double rho = 1.0/tissueData.mean_apicalArea[1];
    
    // calculate preferred cell apical surface A0 and cell surface elasticity K2D
    preferredArea_2D = A0_a*T0_a/(T0_a-T_a);
    areaElasticity_2D = (T_a-T0_a)/A0_a;
    lambda_2D = MechProp.Intercell_ApicalLineTension[1][1];
    
    // calculate bulk & shear moduli
    bulkModulus_2D = areaElasticity_2D/rho - 0.4653*lambda_2D*std::sqrt(rho);
    shearModulus_2D = 0.4653*lambda_2D*std::sqrt(rho);
    
    // and last but not least calculate kappa
    kappa_2D = bulkModulus_2D/shearModulus_2D;
    
    signal_updateGastrulationProperties(preferredArea_2D,areaElasticity_2D,bulkModulus_2D,shearModulus_2D,kappa_2D);
}

bool Tissue::SetTissueStructure()
{
	// add the sorted lists of edges for every vertex
	for(std::list<Vertex>::iterator it_v = ListVertex.begin(); it_v != ListVertex.end() ; ++it_v)
    {
        it_v->updateNeighbouringEdges();
        
    }
	
	// if the tissue is periodic inform all the cells if they cross the boundary
	if(isPeriodic)
	{
		for(std::list<Cell>::iterator itCell = ListCell.begin(); itCell!=ListCell.end();itCell++)
		{
			itCell->CrossBoundary();
		}
	}
    
    for(std::list<Cell>::iterator itCell = ListCell.begin(); itCell!=ListCell.end();itCell++)
    {
        itCell->orderEdges();
    }
    
	return (true);
}

//bool Tissue::update_wo_pressure()
//{
//    UpdateLengths();
//	UpdateBarycenters();
//	UpdateElasticEnergyOfCells();
//	UpdateForcesFromBarycenters();
//    
//	UpdateForcesAndEnergy();
//    
//	Changed = false;
//}


bool Tissue::TestUpdate()
{
    UpdateLengths();
	UpdateBarycenters();
	UpdateElasticEnergyOfCells();
    update_forces();
    W();

//    UpdateLengths();
//	UpdateBarycenters();
//	UpdateElasticEnergyOfCells();
//    update_forces();
//    W();
    
//    bool TEST_VOLUME_DER = false;
//    
//    // test if the volumes are derived properly
//    if(TEST_VOLUME_DER)
//    {
//        UpdateVolumeDerivatives_numerical();
//        
//        for(unsigned cell = 0; cell<CELL_NUMBER; cell++)
//        {
//            
//            for(unsigned vertex = 0; vertex<VERTEX_NUMBER; vertex++)
//            {
//                
//                if(Point::Norm(VolumeDerivativesApical_numerically[cell*VERTEX_NUMBER + vertex]-VolumeDerivativesApical_analytically[cell*VERTEX_NUMBER + vertex])>10e-5 ||Point::Norm(VolumeDerivativesBasal_numerically[cell*VERTEX_NUMBER + vertex]-VolumeDerivativesBasal_analytically[cell*VERTEX_NUMBER + vertex])>10e-5)
//                {
//                    std::cout << "cell: " << cell << "  Crosses Boundary?  " <<  CellFromNumber(cell)->CrossBoundary() << std::endl;
//                    std::cout << "      vertex: " << vertex << std::endl << "analytical    ";
//                    
//                    VolumeDerivativesApical_numerically[cell*VERTEX_NUMBER + vertex].Print();
//                    std::cout << "numerical    ";
//                    VolumeDerivativesApical_analytically[cell*VERTEX_NUMBER + vertex].Print();
//                    std::cout << "analytical   ";
//                    VolumeDerivativesBasal_numerically[cell*VERTEX_NUMBER + vertex].Print();
//                    std::cout << "numerical    ";
//                    VolumeDerivativesBasal_analytically[cell*VERTEX_NUMBER + vertex].Print();
//                }
//            }
//        }
//    }
//       
//    
    TEST_FORCES = true;

    // test if all forces have been calculated exactly
    if(TEST_FORCES)
    {
        // do the numerical differentiation
        update_forces_numeric();
        
        bool ProbVertexFound = false;
        
        for(std::list<Vertex>::iterator it=ListVertex.begin();it!=ListVertex.end();it++)
        {
            
            double difference = std::max(Point::Norm(it->dW_dva_numerical-it->dW_dva)/std::max(Point::Norm(it->dW_dva_numerical),Point::Norm(it->dW_dva)),Point::Norm(it->dW_dvb_numerical-it->dW_dvb)/std::max(Point::Norm(it->dW_dvb_numerical),Point::Norm(it->dW_dvb)));
            
            //std::cout << difference << std::endl;
            if((std::max(Point::Norm(it->dW_dva_numerical),Point::Norm(it->dW_dva))>0.0001 || std::max(Point::Norm(it->dW_dvb_numerical),Point::Norm(it->dW_dvb)) > 0.0001) && difference>10e-4)
            {
                std::cout << "difference = " << difference << std::endl;
                
                if(!ProbVertexFound)
                {
                    std::cout << "Forces vary for the vertices: "<< std::endl;ProbVertexFound=true;
                }
                
                std::cout << "Vertex no: " << it->Number << "\nRelative Force difference = " << Point::Norm(it->dW_dva_numerical-it->dW_dva)/Point::Norm(it->dW_dva_numerical) << " , " << Point::Norm(it->dW_dvb_numerical-it->dW_dvb)/Point::Norm(it->dW_dvb_numerical) << std::endl;
                
                std::cout << "Apical Force analytically = ";
                it->dW_dva.Print();
                std::cout << "Apical force numerically = ";
                it->dW_dva_numerical.Print();
                
                std::cout << "Basal Force analytically = ";
                it->dW_dvb.Print();
                std::cout << "Basal force numerically = ";
                it->dW_dvb_numerical.Print();
                if(!it->q_b.isZero()) std::cout << "lateral edge crosses boundary!";
                
                // check the neighbouring cells
                std::list<Cell*> cellList = it->getNeighbouringCells();
                
                if(isPeriodic){
                    for(std::list<Cell*>::iterator it_cell = cellList.begin(); it_cell!=cellList.end(); it_cell++)
                    {
                        if(!(*it_cell)->q_BC_b.isZero()) std::cout << "in neighbouring cell q_BC_b!=0"<< std::endl;
                    }
                }
            }
        }
        
        if(!ProbVertexFound) {
            std::cout << "Differentiation wrt vertex position works" << std::endl;
        }
        else {
            std::cout << "Differentiation wrt vertex problematic" << std::endl;
            return 0;
        }
        //
        if(isPeriodic && Point::Norm(dW_dL-dW_dL_numerical)>10e-3)
        {
            std::cout<<"Force on system size deviates:" << std::endl;
            dW_dL.Print();
            dW_dL_numerical.Print();
            return 0;
        }
    }
    
    return 1;
}


void Tissue::UpdateLengths()
{
    for(std::list<Vertex>::iterator it_vertex=ListVertex.begin();it_vertex!=ListVertex.end();it_vertex++) it_vertex->UpdateLength();
    for(std::list<Edge>::iterator it_edge=ListEdge.begin();it_edge!=ListEdge.end();it_edge++) it_edge->UpdateLengths();
}


void Tissue::UpdateBarycenters()
{
	for(std::list<Cell>::iterator it_cell=ListCell.begin();it_cell!=ListCell.end();it_cell++) it_cell->UpdateCenter();
    if(onlyApicalContributions) return; // no need to update the lateral barycenters
    for(std::list<Edge>::iterator it_edge=ListEdge.begin();it_edge!=ListEdge.end();it_edge++) it_edge->UpdateBarycenter();
}

// take the current state of the system (i.e. positions of the vertices) and calculate the values, that are neccessary for the calculation of the evolution
//void Tissue::UpdateForcesAndEnergy()
//{
//	// SET THE FORCES AND THE ENERGY TO ZERO
//	for(std::list<Vertex>::iterator it_v = ListVertex.begin();it_v!=ListVertex.end();it_v++)
//	{
//		(*it_v).resetForces();
//	}
//    
//    dW_dL = Point();
//	
//    Energy_Volume = 0;
//	Energy_Line = 0;
//	Energy_Surface = 0;
//    
//    // CALCULATE FORCE CONTRIBUTIONS FROM THE FORCES ACTING ON THE BARYCENTERS
//    UpdateForcesFromBarycenters();
//    
//    // ADD THE FORCE CONTRIBUTIONS FROM THE EDGES TO THE CONNECTED EDGES
//    for(std::list<Edge>::iterator it_edge = ListEdge.begin(); it_edge!=ListEdge.end(); it_edge++)   it_edge->UpdateForcesAndEnergy();
//    
//	// LATERAL LINE TENSIONS
//	for(std::list<Vertex>::iterator it_v = ListVertex.begin();it_v!=ListVertex.end();it_v++) it_v->UpdateForcesAndEnergy();
//	    
//	double VolumeStress=0;
//	
//	// VOLUME PRESSURE FROM THE CELLS
//	for(std::list<Cell>::iterator it_cell=ListCell.begin();it_cell!=ListCell.end();it_cell++) {
//        // energy contribution
//		Energy_Volume += it_cell->ElasticEnergy;
//		// force on system size
//		VolumeStress += it_cell->Pressure*it_cell->Volume;
//	}
//	
//    Energy = Energy_Volume+Energy_Line+Energy_Surface;
//    
//    if(isPeriodic) { // in the periodic case calculate the forces on the periodic box and include the external tension to the energy
//        // add the volume contributions
//        dW_dL.x+=VolumeStress;
//        dW_dL.y+=VolumeStress;
//        
//        // TO OBTAIN THE FORCE ACTING ON THE SYSTEM SIZE THE OBTAINED DERIVATIVE HAS TO BE DIVIDED BY THE SYSTEM SIZES (CF. NOTES.PDF)
//        dW_dL.x /= SystemSize.x;
//        dW_dL.y /= SystemSize.y;
//        
//        // ADD CONTRIBUTIONS FROM EXTERNAL TENSION
//        dW_dL.x -= T_ext*SystemSize.y;
//        dW_dL.y -= T_ext*SystemSize.x;
//        
//        if(Lx_fix) dW_dL.x = 0;
//        if(Ly_fix) dW_dL.y = 0;
//        
//        Energy += T_ext*SystemSize.x*SystemSize.y;
//        
//        if(STIFF_BASAL_MEMBRANE) { // if we assume stiff ECM reset all the forces of the basal vertices in z direction to zero
//            for(std::list<Vertex>::iterator it_v = ListVertex.begin(); it_v!=ListVertex.end(); it_v++) it_v->dW_dvb.z = 0;
//        }
//    }
//    
//    //
////    for(std::list<Vertex>::iterator it_v = ListVertex.begin(); it_v!=ListVertex.end(); it_v++) {
////        
////        if( std::abs(Point::Norm(it_v->dW_dva)-Point::Norm(it_v->dW_dvb))>1e-4)
////        {
////            std::cout << "Difference in apical and basal force at vertex " << it_v->Number << ", the difference being " << Point::Norm(it_v->dW_dva)-Point::Norm(it_v->dW_dvb);
////            if(it_v->neighbouringEdgeCrossesBoundary()) {std::cout << ", crosses = 1!\n";}
////            else { std::cout << ", crosses = 0!\n";}
////            
////            it_v->dW_dva.Print();
////            it_v->dW_dvb.Print();
////            std::cout << "SystemSize = ( " << SystemSize.x << " , " << SystemSize.y << " )" << std::endl;
////        }
////        
////    }
//    
//}

void Tissue::update_forces()
{
	// SET THE FORCES AND THE ENERGY TO ZERO
	for(std::list<Vertex>::iterator it_v = ListVertex.begin();it_v!=ListVertex.end();it_v++)
	{
		(*it_v).resetForces();
	}
    
    dW_dL = Point();
	
    // CALCULATE FORCE CONTRIBUTIONS FROM THE FORCES ACTING ON THE BARYCENTERS
    UpdateForcesFromBarycenters();
    
    // ADD THE FORCE CONTRIBUTIONS FROM THE EDGES TO THE CONNECTED EDGES
    for(std::list<Edge>::iterator it_edge = ListEdge.begin(); it_edge!=ListEdge.end(); it_edge++)   it_edge->update_forces();
    
    bool apicalPositionsFixed = false;
    bool moveBasalLikeApical = false;

	// LATERAL LINE TENSIONS AND ELASTICITY OF BASEMENT MEMBRANE
	for(std::list<Vertex>::iterator it_v = ListVertex.begin();it_v!=ListVertex.end();it_v++)
    {
        it_v->update_forces();
        if(apicalPositionsFixed) it_v->dW_dva = Point();
        if(moveBasalLikeApical) it_v->dW_dvb = Point(it_v->dW_dva.x,it_v->dW_dva.y,-it_v->dW_dva.z);
    }
    
    
    if(isPeriodic) { // in the periodic case calculate the forces on the periodic box and include the external tension to the energy

       	double VolumeStress=0;
        // VOLUME PRESSURE FROM THE CELLS
        for(std::list<Cell>::iterator it_cell=ListCell.begin();it_cell!=ListCell.end();it_cell++) {
            // force on system size
            VolumeStress += it_cell->Pressure*it_cell->Volume;
        }

        
        // add the volume contributions
        dW_dL.x+=VolumeStress;
        dW_dL.y+=VolumeStress;
        
        // TO OBTAIN THE FORCE ACTING ON THE SYSTEM SIZE THE OBTAINED DERIVATIVE HAS TO BE DIVIDED BY THE SYSTEM SIZES (CF. NOTES.PDF)
        dW_dL.x /= SystemSize.x;
        dW_dL.y /= SystemSize.y;
        
        // ADD CONTRIBUTIONS FROM EXTERNAL TENSION
        dW_dL.x -= T_ext*SystemSize.y;
        dW_dL.y -= T_ext*SystemSize.x;
        
        if(Lx_fix) dW_dL.x = 0;
        if(Ly_fix) dW_dL.y = 0;
        
        if(STIFF_BASAL_MEMBRANE) { // if we assume stiff ECM reset all the forces of the basal vertices in z direction to zero
            for(std::list<Vertex>::iterator it_v = ListVertex.begin(); it_v!=ListVertex.end(); it_v++) it_v->dW_dvb.z = 0;
        }
    }
    
    bool gastrulation_checkForMove = true;
    
    if(gastrulation_checkForMove)
    {
        // CHECK IF VERTEX IS ALLOWED TO MOVE DUE TO CELL CONSTRAINTS
        for(std::list<Vertex>::iterator it_v = ListVertex.begin();it_v!=ListVertex.end();it_v++)
        {
            std::list<Cell*> CellList = it_v->getNeighbouringCells();
            
            // if any of the cells attached to the vertex is fixed, set the forces on the vertex to zero
            for(std::list<Cell*>::iterator it_c = CellList.begin(); it_c!=CellList.end();it_c++)
            {
                if((*it_c)->positionFixed)
                {
                    it_v->dW_dva = Point();
                    it_v->dW_dvb = Point();
                    break;
                }
            }
        }
    }

    
}

//void Tissue::UpdateForcesAndEnergy_wo_pressure()
//{
//	// SET THE FORCES AND THE ENERGY TO ZERO
//	//ResetVertexForces();
//    
//    
//	dW_dL = Point();
//	Energy_Volume = 0;
//	Energy_Line = 0;
//	Energy_Surface = 0;
//    
//    // set pressure to zero
//    for(std::list<Cell>::iterator it_cell=ListCell.begin();it_cell!=ListCell.end();it_cell++) {
//        Cell *c = &(*it_cell);
//		c->Pressure=0;
//    }
//    
//    // ADD THE FORCE CONTRIBUTIONS FROM THE EDGES
//    for(std::list<Edge>::iterator it_edge = ListEdge.begin(); it_edge!=ListEdge.end(); it_edge++)
//	{
//		it_edge->UpdateForcesAndEnergy();
//	}
//    
//	// LATERAL LINE TENSIONS
//	for(std::list<Vertex>::iterator it_v = ListVertex.begin();it_v!=ListVertex.end();it_v++)
//	{
//        it_v->UpdateForcesAndEnergy();
//	}
//    
//	double VolumeStress=0;
//	
//	// VOLUME PRESSURE FROM THE CELLS
//	for(std::list<Cell>::iterator it_cell=ListCell.begin();it_cell!=ListCell.end();it_cell++) {
//        Cell *c = &(*it_cell);
//		c->UpdateElasticEnergy();
//        // energy contribution
//		Energy_Volume += c->ElasticEnergy;
//		// force on system size
//		VolumeStress += c->Pressure*c->Volume;
//	}
//	
//    Energy = Energy_Volume+Energy_Line+Energy_Surface;
//    
//    if(isPeriodic) { // in the periodic case calculate the forces on the periodic box and include the external tension to the energy
//        // add the volume contributions
//        dW_dL.x+=VolumeStress;
//        dW_dL.y+=VolumeStress;
//        
//        // TO OBTAIN THE FORCE ACTING ON THE SYSTEM SIZE THE OBTAINED DERIVATIVE HAS TO BE DIVIDED BY THE SYSTEM SIZES (CF. NOTES.PDF)
//        dW_dL.x /= SystemSize.x;
//        dW_dL.y /= SystemSize.y;
//        
//        // ADD CONTRIBUTIONS FROM EXTERNAL TENSION
//        dW_dL.x -= T_ext*SystemSize.y;
//        dW_dL.y -= T_ext*SystemSize.x;
//        
//        Energy += T_ext*SystemSize.x*SystemSize.y;
//        
//        if(STIFF_BASAL_MEMBRANE) { // if we assume stiff ECM reset all the forces of the basal vertices in z direction to zero
//            for(std::list<Vertex>::iterator it_v = ListVertex.begin(); it_v!=ListVertex.end(); it_v++) it_v->dW_dvb.z = 0;
//        }
//    }
////}
//
//double Tissue::update_energy()
//{
//    // UPDATE THE BARYCENTER AND THE VOLUME OF THE CELLS
//    UpdateBarycenters();
//    
//	// SET THE FORCES AND THE ENERGY TO ZERO
//	Energy_Volume = 0;
//	Energy_Line = 0;
//	Energy_Surface = 0;
//	
//    // ADD THE ENERGY CONTRIBUTIONS FROM THE EDGES
//    for(std::list<Edge>::iterator it_edge = ListEdge.begin(); it_edge!=ListEdge.end(); it_edge++) it_edge->update_energy();
//
//	// LATERAL LINE TENSIONS
//	for(std::list<Vertex>::iterator it_v = ListVertex.begin();it_v!=ListVertex.end();it_v++) {
//        Energy_Line += it_v->getLength()*it_v->G_l;
//	}
//    	
//	// VOLUME PRESSURE FROM THE CELLS
//	for(std::list<Cell>::iterator it_cell=ListCell.begin();it_cell!=ListCell.end();it_cell++) {
//        Energy_Volume += it_cell->K/2.*std::pow(it_cell->UpdateVolume()-it_cell->V0,2);
//    }
//	
//	Energy = Energy_Volume+Energy_Line+Energy_Surface;
//    
//    if(isPeriodic) Energy += T_ext*SystemSize.x*SystemSize.y;
//    
//    return Energy;
//}

double Tissue::W()
{
    // UPDATE THE BARYCENTER AND THE VOLUME OF THE CELLS
    UpdateBarycenters();
    
	// SET THE FORCES AND THE ENERGY TO ZERO
    Energy = 0;
    
    // ADD THE ENERGY CONTRIBUTIONS FROM THE EDGES (lateral surface tension, apical and basal line tensions)
    for(std::list<Edge>::iterator it_edge = ListEdge.begin(); it_edge!=ListEdge.end(); it_edge++) Energy += it_edge->W();
    
	// LATERAL LINE TENSIONS & ENERGETICAL CONTRIBUTIONS FROM THE ELLIPSOIDAL POTENTIAL
	for(std::list<Vertex>::iterator it_v = ListVertex.begin();it_v!=ListVertex.end();it_v++) {
        Energy += it_v->W();
	}
    
	// VOLUME PRESSURE FROM THE CELLS, APICAL AND BASAL SURFACE TENSIONS
	for(std::list<Cell>::iterator it_cell=ListCell.begin();it_cell!=ListCell.end();it_cell++) {
        Energy += it_cell->W();
    }
    
    // VOLUME PRESSURE FROM THE INTERNAL LUMEN
    bool Lumen_a = (K_Lumen_a!=0 || P0_Lumen_a!=0); // updating the volume of a already updates the volume of b
    
    if(Lumen_a) {
        update_V_Lumen_a();
        Energy += K_Lumen_a/2.*std::pow(V_Lumen_a-V0_Lumen_a,2) - P0_Lumen_a*V_Lumen_a;
    }

    if(K_Lumen_b!=0 || P0_Lumen_b!=0) {
        update_V_Lumen_b();
        if(!Lumen_a) update_V_Lumen_b(); // updating the volume only if neccessary
        Energy += K_Lumen_b/2.*std::pow(V_Lumen_b-V0_Lumen_b,2) - P0_Lumen_b*V_Lumen_b;
    }

    
    if(isPeriodic) Energy += T_ext*SystemSize.x*SystemSize.y;
    
    return Energy;
}

double Tissue::W(bool update_BCs)
{
    // UPDATE THE BARYCENTER AND THE VOLUME OF THE CELLS
    if(update_BCs) UpdateBarycenters();
    
	// SET THE FORCES AND THE ENERGY TO ZERO
    Energy = 0;
    
    // ADD THE ENERGY CONTRIBUTIONS FROM THE EDGES (lateral surface tension, apical and basal line tensions)
    for(std::list<Edge>::iterator it_edge = ListEdge.begin(); it_edge!=ListEdge.end(); it_edge++) Energy += it_edge->W();
    
	// LATERAL LINE TENSIONS & ENERGETICAL CONTRIBUTIONS FROM THE ELLIPSOIDAL POTENTIAL
	for(std::list<Vertex>::iterator it_v = ListVertex.begin();it_v!=ListVertex.end();it_v++) {
        Energy += it_v->W();
	}
    
	// VOLUME PRESSURE FROM THE CELLS & APICAL AND BASAL SURFACE TENSIONS
	for(std::list<Cell>::iterator it_cell=ListCell.begin();it_cell!=ListCell.end();it_cell++) {
        Energy += it_cell->W();
    }
    
    // VOLUME PRESSURE FROM THE INTERNAL LUMEN
    bool Lumen_a = (K_Lumen_a!=0 || P0_Lumen_a!=0); // updating the volume of a already updates the volume of b
    
    if(Lumen_a) {
        update_V_Lumen_a();
        Energy += K_Lumen_a/2.*std::pow(V_Lumen_a-V0_Lumen_a,2) - P0_Lumen_a*V_Lumen_a;
    }
    
    if(K_Lumen_b!=0 || P0_Lumen_b!=0) {
        if(!Lumen_a) update_V_Lumen_b(); // updating the volume only if neccessary
        Energy += K_Lumen_b/2.*std::pow(V_Lumen_b-V0_Lumen_b,2) - P0_Lumen_b*V_Lumen_b;
    }
    
    if(isPeriodic) Energy += T_ext*SystemSize.x*SystemSize.y;
    
    return Energy;
}


/// use central difference method to numerically calculate the forces in the system
void Tissue::update_forces_numeric()
{
    // 1) choose small dLx, change Lx->Lx+dLx/2 & update the vertex positions
    // 2) calulate the energy E(Lx+dLx/2, Ly):=dEx_plus
    // 3) decrease Lx->Lx-dLx/2 & update the vertex positions
    // 4) calulate the energy E(Lx-dLx/2, Ly):=dEx_minus
    // 5) force on Lx := (dEx_plus-dEx_minus)/dLx
    // 6)-10) repeat the steps 1)-5) for Ly and dLy
    // 11) change the system back to its original state
    
    bool update_BCs = true;
    
    if(isPeriodic) { // in the periodic case also calculate the derivative of the energy w.r.t. the system size
        
        
        double E_x_plus, E_x_minus, E_y_plus, E_y_minus;
        
        double dLx = 1e-6;
        /* 1 */
        ChangeLxAndLy(SystemSize.x+dLx/2, SystemSize.y);
        /* 2 */
        E_x_plus = W(update_BCs);;
        /* 3 */
        ChangeLxAndLy(SystemSize.x-dLx, SystemSize.y);
        /* 4 */
        E_x_minus = W(update_BCs);
        
        /* 5 */
        dW_dL_numerical.x = (E_x_minus-E_x_plus)/dLx;
        
        double dLy = 1e-6;
        /* 6 */
        ChangeLxAndLy(SystemSize.x+dLx/2, SystemSize.y+dLy/2);
        /* 7 */
        E_y_plus = W(update_BCs);
        /* 8 */
        ChangeLxAndLy(SystemSize.x, SystemSize.y-dLy);
        /* 9 */
        E_y_minus = W(update_BCs);
        
        /* 10 */
        dW_dL_numerical.y = (E_y_minus-E_y_plus)/dLy;
        
        /* 11 */
        ChangeLxAndLy(SystemSize.x, SystemSize.y+dLy/2);
    }
    
    // iterate through all vertices and calculate the forces on the vertices numerically
    for(std::list<Vertex>::iterator it_vertex = ListVertex.begin(); it_vertex != ListVertex.end(); it_vertex++)
    {
        Vertex* v = &(*it_vertex);
        
        double d = 1e-4;
        double d_2 = d/2.;
        
        double E_x_plus, E_x_minus,E_y_plus, E_y_minus,E_z_plus, E_z_minus;
        
        Point P_x_d_2(d_2,0,0),P_x_d(-d,0,0),P_y_d_2(0,d_2,0),P_y_d(0,-d,0),P_z_d_2(0,0,d_2),P_z_d(0,0,-d),P_null(0,0,0);
        
        Quadrant qpx = v->moveVertex(P_x_d_2,P_null);
        E_x_plus = W(update_BCs);
        Quadrant qmx = v->moveVertex(P_x_d,P_null);
        E_x_minus = W(update_BCs);
        v->moveVertex(P_x_d_2,P_null);
        
        Quadrant qpy = v->moveVertex(P_y_d_2,P_null);
        UpdateBarycenters();
        E_y_plus = W(update_BCs);
        Quadrant qmy = v->moveVertex(P_y_d,P_null);
        UpdateBarycenters();
        E_y_minus = W(update_BCs);
        v->moveVertex(P_y_d_2,P_null);
        
        Quadrant qpz = v->moveVertex(P_z_d_2,P_null);
        E_z_plus = W(update_BCs);
        Quadrant qmz = v->moveVertex(P_z_d,P_null);
        E_z_minus = W(update_BCs);
        v->moveVertex(P_z_d_2,P_null);
        
        v->dW_dva_numerical.x = (E_x_minus - E_x_plus)/d;
        v->dW_dva_numerical.y = (E_y_minus - E_y_plus)/d;
        v->dW_dva_numerical.z = (E_z_minus - E_z_plus)/d;
        
        qpx=v->moveVertex(P_null,P_x_d_2);
        E_x_plus = W(update_BCs);
        qmx=v->moveVertex(P_null,P_x_d);
        E_x_minus = W(update_BCs);
        v->moveVertex(P_null,P_x_d_2);
        
        qpy=v->moveVertex(P_null,P_y_d_2);
        E_y_plus = W(update_BCs);
        qmy=v->moveVertex(P_null,P_y_d);
        E_y_minus = W(update_BCs);
        v->moveVertex(P_null,P_y_d_2);
        
        qpz=v->moveVertex(P_null,P_z_d_2);
        E_z_plus = W(update_BCs);
        qmz=v->moveVertex(P_null,P_z_d);
        E_z_minus = W(update_BCs);
        v->moveVertex(P_null,P_z_d_2);
        
        v->dW_dvb_numerical.x = (E_x_minus - E_x_plus)/d;
        v->dW_dvb_numerical.y = (E_y_minus - E_y_plus)/d;
        v->dW_dvb_numerical.z = (E_z_minus - E_z_plus)/d;
    }
}


///// use central difference method to numerically calculate the forces in the system
//void Tissue::UpdateForces_numeric()
//{
//    // 1) choose small dLx, change Lx->Lx+dLx/2 & update the vertex positions
//    // 2) calulate the energy E(Lx+dLx/2, Ly):=dEx_plus
//    // 3) decrease Lx->Lx-dLx/2 & update the vertex positions
//    // 4) calulate the energy E(Lx-dLx/2, Ly):=dEx_minus
//    // 5) force on Lx := (dEx_plus-dEx_minus)/dLx
//    // 6)-10) repeat the steps 1)-5) for Ly and dLy
//    // 11) change the system back to its original state
//   
//    if(isPeriodic) { // in the periodic case also calculate the derivative of the energy w.r.t. the system size
//        
//        
//        double E_x_plus, E_x_minus, E_y_plus, E_y_minus;
//        
//        double dLx = 1e-5;
//        /* 1 */
//        ChangeLxAndLy(SystemSize.x+dLx/2, SystemSize.y);
//        /* 2 */
//        update_energy();
//        E_x_plus = Energy;
//        /* 3 */
//        ChangeLxAndLy(SystemSize.x-dLx, SystemSize.y);
//        /* 4 */
//        update_energy();
//        E_x_minus = Energy;
//        
//        /* 5 */
//        dW_dL_numerical.x = (E_x_minus-E_x_plus)/dLx;
//        
//        double dLy = 1e-5;
//        /* 6 */
//        ChangeLxAndLy(SystemSize.x+dLx/2, SystemSize.y+dLy/2);
//        /* 7 */
//        update_energy();
//        E_y_plus = Energy;
//        /* 8 */
//        ChangeLxAndLy(SystemSize.x, SystemSize.y-dLy);
//        /* 9 */
//        update_energy();
//        E_y_minus = Energy;
//
//        /* 10 */
//        dW_dL_numerical.y = (E_y_minus-E_y_plus)/dLy;
//        
//        /* 11 */
//        ChangeLxAndLy(SystemSize.x, SystemSize.y+dLy/2);
//    }
//    
//    // iterate through all vertices and calculate the forces on the vertices numerically
//    for(std::list<Vertex>::iterator it_vertex = ListVertex.begin(); it_vertex != ListVertex.end(); it_vertex++)
//    {
//        Vertex* v = &(*it_vertex);
//        
//        double d = 1e-5;
//        double d_2 = d/2.;
//        
//        double E_x_plus, E_x_minus,E_y_plus, E_y_minus,E_z_plus, E_z_minus;
//                
//        Point P_x_d_2(d_2,0,0),P_x_d(-d,0,0),P_y_d_2(0,d_2,0),P_y_d(0,-d,0),P_z_d_2(0,0,d_2),P_z_d(0,0,-d),P_null(0,0,0);
//        
//        Quadrant qpx = v->moveVertex(P_x_d_2,P_null);
//        update_energy();
//        E_x_plus = Energy;
//        Quadrant qmx = v->moveVertex(P_x_d,P_null);
//        update_energy();
//        E_x_minus = Energy;
//        v->moveVertex(P_x_d_2,P_null);
//        
//        Quadrant qpy = v->moveVertex(P_y_d_2,P_null);
//        update_energy();
//        E_y_plus = Energy;
//        Quadrant qmy = v->moveVertex(P_y_d,P_null);
//        update_energy();
//        E_y_minus = Energy;
//        v->moveVertex(P_y_d_2,P_null);
//        
//        Quadrant qpz = v->moveVertex(P_z_d_2,P_null);
//        update_energy();
//        E_z_plus = Energy;
//        Quadrant qmz = v->moveVertex(P_z_d,P_null);
//        update_energy();
//        E_z_minus = Energy;
//        v->moveVertex(P_z_d_2,P_null);
//        
//        v->dW_dva_numerical.x = (E_x_minus - E_x_plus)/d;
//        v->dW_dva_numerical.y = (E_y_minus - E_y_plus)/d;
//        v->dW_dva_numerical.z = (E_z_minus - E_z_plus)/d;
//        
//        qpx=v->moveVertex(P_null,P_x_d_2);
//        update_energy();
//        E_x_plus = Energy;
//        qmx=v->moveVertex(P_null,P_x_d);
//        update_energy();
//        E_x_minus = Energy;
//        v->moveVertex(P_null,P_x_d_2);
//        
//        qpy=v->moveVertex(P_null,P_y_d_2);
//        update_energy();
//        E_y_plus = Energy;
//        qmy=v->moveVertex(P_null,P_y_d);
//        update_energy();
//        E_y_minus = Energy;
//        v->moveVertex(P_null,P_y_d_2);
//        
//        qpz=v->moveVertex(P_null,P_z_d_2);
//        update_energy();
//        E_z_plus = Energy;
//        qmz=v->moveVertex(P_null,P_z_d);
//        update_energy();
//        E_z_minus = Energy;
//        v->moveVertex(P_null,P_z_d_2);
//        update_energy();
//        
//        v->dW_dvb_numerical.x = (E_x_minus - E_x_plus)/d;
//        v->dW_dvb_numerical.y = (E_y_minus - E_y_plus)/d;
//        v->dW_dvb_numerical.z = (E_z_minus - E_z_plus)/d;
//    }
//}
//    // iterate through all edges and calculate the forces on the lateral barycenter numerically
//    for(std::list<Edge>::iterator it_edge = ListEdge.begin(); it_edge != ListEdge.end(); it_edge++)
//    {
//        UpdateBarycenters();
//        
//        Edge* e = &(*it_edge);
//        
//        double d = 1e-5;
//        double d_2 = d/2.;
//        
//        double E_x_plus, E_x_minus,E_y_plus, E_y_minus,E_z_plus, E_z_minus;
//        
//        //        if(v->coord.x+d_2>=SystemSize.x || v->coord.x-d_2<0 ||v->coord.y+d_2>=SystemSize.y || v->coord.y-d_2<0 || v->coord_b.x+d_2>=SystemSize.x || v->coord_b.x-d_2<0 ||v->coord_b.y+d_2>=SystemSize.y || v->coord_b.y-d_2<0)
//        //        {
//        //            std::cout << "in Tissue::UpdateForces_numeric() crossed the boundary!" << std::endl;
//        //        }
//        
//        Point P_x_d_2(d_2/4,0,0),P_x_d(-d/4,0,0),P_y_d_2(0,d_2/4,0),P_y_d(0,-d/4,0),P_z_d_2(0,0,d_2/4),P_z_d(0,0,-d/4),P_null(0,0,0);
//        
//        e->BC_l += P_x_d_2;
//        update_energy();
//        E_x_plus = Energy;
//        e->BC_l += P_x_d;
//        update_energy();
//        E_x_minus = Energy;
//        e->BC_l += P_x_d_2;
//
//        e->BC_l += P_y_d_2;
//        update_energy();
//        E_y_plus = Energy;
//        e->BC_l += P_y_d;
//        update_energy();
//        E_y_minus = Energy;
//        e->BC_l += P_y_d_2;
//
//        e->BC_l += P_z_d_2;
//        update_energy();
//        E_z_plus = Energy;
//        e->BC_l += P_z_d;
//        update_energy();
//        E_z_minus = Energy;
//        e->BC_l += P_z_d_2;
//        
//        e->ForceFromBarycenter_Area_numerical.x = (E_x_minus - E_x_plus)/d;
//        e->ForceFromBarycenter_Area_numerical.y = (E_y_minus - E_y_plus)/d;
//        e->ForceFromBarycenter_Area_numerical.z = (E_z_minus - E_z_plus)/d;
//        
//        update_energy();
//    }
    
    // iterate through all edges and calculate the forces on the apical barycenter numerically
//    for(std::list<Cell>::iterator it_cell = ListCell.begin(); it_cell != ListCell.end(); it_cell++)
//    {
//        UpdateBarycenters();
//        
//        Cell* c = &(*it_cell);
//        
//        double d = 1e-5;
//        double d_2 = d/2.;
//        
//        double E_x_plus, E_x_minus,E_y_plus, E_y_minus,E_z_plus, E_z_minus;
//        //        //        if(v->coord.x+d_2>=SystemSize.x || v->coord.x-d_2<0 ||v->coord.y+d_2>=SystemSize.y || v->coord.y-d_2<0 || v->coord_b.x+d_2>=SystemSize.x || v->coord_b.x-d_2<0 ||v->coord_b.y+d_2>=SystemSize.y || v->coord_b.y-d_2<0)
//        //        {
//        //            std::cout << "in Tissue::UpdateForces_numeric() crossed the boundary!" << std::endl;
//        //        }
//        
//        int no = c->ListEdges.size();
//        
//        Point P_x_d_2(d_2,0,0),P_x_d(-d,0,0),P_y_d_2(0,d_2,0),P_y_d(0,-d,0),P_z_d_2(0,0,d_2),P_z_d(0,0,-d),P_null(0,0,0);
//        
//        c->BC_a += P_x_d_2;
//        update_energy(false);
//        E_x_plus = Energy;
//        c->BC_a += P_x_d;
//        update_energy(false);
//        E_x_minus = Energy;
//        c->BC_a += P_x_d_2;
//        
//        c->BC_a += P_y_d_2;
//        update_energy(false);
//        E_y_plus = Energy;
//        c->BC_a += P_y_d;
//        update_energy(false);
//        E_y_minus = Energy;
//        c->BC_a += P_y_d_2;
//        
//        c->BC_a += P_z_d_2;
//        update_energy(false);
//        E_z_plus = Energy;
//        c->BC_a += P_z_d;
//        update_energy(false);
//        E_z_minus = Energy;
//        c->BC_a += P_z_d_2;
//        
//        c->dAa_dBCa_numerical.x = (E_x_plus - E_x_minus)/d;
//        c->dAa_dBCa_numerical.y = (E_y_plus - E_y_minus)/d;
//        c->dAa_dBCa_numerical.z = (E_z_plus - E_z_minus)/d;
//        
//        update_energy();
//    }


void Tissue::getCoordinatesAndForces(std::list<Point> *ApicalCoordinates, std::list<Point> *BasalCoordinates, std::list<Point> *dW_dvas, std::list<Point> *dW_dvbs)
{
	// clear the input lists
	ApicalCoordinates->clear();
	BasalCoordinates->clear();
	dW_dvas->clear();
	dW_dvbs->clear();
	
	
	std::list<Vertex>::iterator itv = ListVertex.begin();
	for(;itv!=ListVertex.end();itv++)
	{
		ApicalCoordinates->push_back((*itv).getCoord());
		ApicalCoordinates->push_back((*itv).getBasalCoord());
		dW_dvas->push_back((*itv).getdW_dva());
		dW_dvbs->push_back((*itv).getdW_dvb());
	}
}

void Tissue::getForces(std::list<Point> *dW_dvas, std::list<Point> *dW_dvbs)
{
	// clear the input lists
	dW_dvas->clear();
	dW_dvbs->clear();
	
	
	std::list<Vertex>::iterator itv = ListVertex.begin();
	for(;itv!=ListVertex.end();itv++)
	{
		dW_dvas->push_back((*itv).getdW_dva());
		dW_dvbs->push_back((*itv).getdW_dvb());
	}
}

void Tissue::setCoordinatesAndForces(std::list<Point> *ApicalCoordinates, std::list<Point> *BasalCoordinates, std::list<Point> *dW_dvas, std::list<Point> *dW_dvbs)
{	
	std::list<Vertex>::iterator itv = ListVertex.begin();
	std::list<Point>::iterator it_apicalCoordinates = ApicalCoordinates->begin();
	std::list<Point>::iterator it_basalCoordinates = BasalCoordinates->begin();
	std::list<Point>::iterator it_dW_dvas = dW_dvas->begin();
	std::list<Point>::iterator it_dW_dvbs = dW_dvbs->begin();
	
	for(;itv!=ListVertex.end();itv++)
	{
		it_apicalCoordinates++;it_basalCoordinates++;it_dW_dvas++;it_dW_dvbs++;
		
		(*itv).setCoord(*it_apicalCoordinates);
		(*itv).setBasalCoord(*it_basalCoordinates);
		(*itv).setdW_dva(*it_dW_dvas);
		(*itv).setdW_dvb(*it_dW_dvbs);
	}
	
}


void Tissue::setCoordinates(std::list<Point> *ApicalCoordinates, std::list<Point> *BasalCoordinates)
{
	//if(ApicalCoordinates->size()!=ListVertex.size()) return -1;
	
	std::list<Vertex>::iterator itv = ListVertex.begin();
	std::list<Point>::iterator it_apicalCoordinates = ApicalCoordinates->begin();
	std::list<Point>::iterator it_basalCoordinates = BasalCoordinates->begin();
	
	for(;itv!=ListVertex.end();itv++)
	{
		it_apicalCoordinates++;it_basalCoordinates++;
		
		(*itv).setCoord(*it_apicalCoordinates);
		(*itv).setBasalCoord(*it_basalCoordinates);
	}
}


void Tissue::UpdateVolumeOfCells()
{
	std::list<Cell>::iterator it_cell;
	for(it_cell=ListCell.begin();it_cell!=ListCell.end();it_cell++)
	{
		it_cell->UpdateVolume();
	}
}

void Tissue::UpdateElasticEnergyOfCells()
{
	for(std::list<Cell>::iterator it_cell=ListCell.begin();it_cell!=ListCell.end();it_cell++)
	{
		it_cell->UpdateElasticEnergy();
	}
}

void Tissue::UpdateForcesFromBarycenters()
{
	// set all the forces in the cells to be zero as they are added up in the edge iteration
	for(std::list<Cell>::iterator it_cell=ListCell.begin();it_cell!=ListCell.end();it_cell++)
	{
		it_cell->dAa_dBCa=Point();
		it_cell->dAb_dBCb=Point();
		it_cell->dV_dBCa=Point();
		it_cell->dV_dBCb=Point();
	}
    
	std::list<Edge>::iterator it_edge;
	for(it_edge=ListEdge.begin();it_edge!=ListEdge.end();it_edge++)
	{
		it_edge->UpdateForcesFromBarycenters();
	}
	
	// divide the contributions in every cell by the number of its edges to obtain the derivative
	for(std::list<Cell>::iterator it_cell=ListCell.begin();it_cell!=ListCell.end();it_cell++)
	{
		double cell_noEdges = it_cell->ListEdges.size();
		it_cell->ForceFromApicalBarycenter_Volume=it_cell->ForceFromApicalBarycenter_Volume*((-1)*it_cell->Pressure/cell_noEdges);
		if(!onlyApicalContributions) it_cell->ForceFromBasalBarycenter_Volume=it_cell->ForceFromBasalBarycenter_Volume*((-1)*it_cell->Pressure/cell_noEdges);
 	}
}

//bool Tissue::UpdateForcesFromLineTensions()
//{
//	std::list<Edge>::iterator it_edge;
//	std::list<Vertex>::iterator it_vertex;
//	int size = ListEdge.size();
//
//	for(it_edge=ListEdge.begin();it_edge!=ListEdge.end();it_edge++)	(*it_edge)->UpdateLineForces(); // forces due to line tension inbetween the two vertices
//	for(it_vertex=ListVertex.begin();it_vertex!=ListVertex.end();it_vertex++)	(*it_vertex)->UpdateLineForces(); // forces due to line tension inbetween the two vertices
//
//}

void Tissue::RemoveShortEdges()
{
	// check if in the end any change has been done to the tissue to know if its neccessary to update the tissue
	bool anyEdgeRemoved = false;
	
	std::list<Edge>::iterator it_edge = ListEdge.begin();
	
	double EdgeThreshold_apical = 1;
	double EdgeThreshold_basal = 1;
	
	while(it_edge!=ListEdge.end())
	{
		// if the edge fulfills the shortness criterion calculate the forces
		if((*it_edge).l_a()<EdgeThreshold_apical || (*it_edge).l_b()<EdgeThreshold_basal)
		{
			
			if((*it_edge++).CheckAndRemoveFromTissue())
			{
				anyEdgeRemoved = true;
			}
            
		}
		else
		{
			it_edge++;
		}
        
		
	}
    
	if(anyEdgeRemoved) Changed = true;
}

void Tissue::AddVertex(Vertex *v)
{
	ListVertex.push_front(*v);
	ListVertex.front().setTissue(this);
}


void Tissue::AddCell(Cell *c)
{
	
	//Point center;
    //	center=c->Centroid;
	ListCell.push_front(*c);
	ListCell.front().T=this;
	//c->T=this;
    //	ListCell.front().ListEdges.clear();//on enleve tout car on va tout remplacer par des pointeurs qui vont vers des objets du tissu
    //	//il faut egalement actualiser les pointeurs, pour qu'ils pointent vers les objets qui sont dans la liste du tissu
    //	std::list<Edge*>::iterator Celledge;
    //	for(Celledge=((*c).ListEdges).begin(); Celledge!=((*c).ListEdges).end();Celledge++)
    //	{//pour chaque cote on regarde s'il est dans la liste du tissu
    //		//std::cout<<"edge recherche : \n";
    //		//(**Celledge).Print();
    //		//std::cout<<"etat de la liste des edges :\n";
    //		//std::list<Edge>::iterator it;
    //		//for(it=ListEdge.begin(); it!=ListEdge.end(); it++)
    //		//{(*it).Print();
    //		//}
    //		std::list<Edge>::iterator findCelledge;
    //		findCelledge = std::find(ListEdge.begin(), ListEdge.end(), **Celledge);
    //
    //		if( findCelledge == ListEdge.end() )
    //		{//dans ce cas l'edge n'est pas deja dans la liste
    //			//std::cout<<"en ajoutant une cellule, un nouveau cote a ete ajoute\n";
    //			AddEdge(*Celledge);//copie l'edge dans le tissu s'il n'y etait pas
    //			((ListCell.front())).ListEdges.push_front(&((ListEdge.front())));//ajoute le cote a la liste de la cellule qui est dans le tissu;
    //			//ici je suppose que les cotes sont deja au courant de quel cellule ils cotoient, car c'est ce que j'ai mis dans le constructeur de cellule, mais peut etre il faudrait rajouter qqch ici
    //			//en fait je m'assure quand meme que le cote qu'on vient d'ajouter sait bien ou est sa cellule
    //			//en fait c'est probablement absolument necessaire pour etre sur, encore une fois, qu'on pointe bien dans la cellule qui est rangee dans le tissu.
    //
    //			//if(ListEdge.front().Director().ProduitVectoriel(ListEdge.front().v1->coord.Vector(center))>0)//teste a l'aide d'un produit vectoriel dans quel sens est oriente l'edge
    //			if(!(c->Isc2(&ListEdge.front())))
    //			{(ListEdge.front()).c1=&(ListCell.front());}
    //			else {(ListEdge.front()).c2=&(ListCell.front());}
    //		}
    //		else
    //		{//si on a trouve un edge qui est deja dans la liste, on l'ajoute a la cellule (nouvelle version dans la liste ListCell)
    //
    //			ListCell.front().ListEdges.push_front(&(*findCelledge));//ajoute le cote a la liste de la cellule qui est dans le tissu;
    //			//et j'ajoute la cellule a l'edge en question
    //			//std::cout<<"coordonnee du centroide : ("<<center.x<<","<<center.y<<")\n";
    //			//
    //			//if((*findCelledge).Director().ProduitVectoriel((*findCelledge).v1->coord.Vector(center))>0 )//teste a l'aide d'un produit vectoriel dans quel sens est oriente l'edge
    //			if(!(c->Isc2(&(*findCelledge))))
    //			{(*findCelledge).c1=&((ListCell.front()));
    //				//std::cout<<"la cellule a ete ajoute en position 1"<<"\n";
    //			}
    //			else
    //			{
    //				(*findCelledge).c2=&((ListCell.front()));
    //				//std::cout<<"la cellule a ete ajoute en position 2"<<"\n";
    //			}
    //			//std::cout<<"edge retrouve appartenant a la cellule "<<ListCell.front().Number<<":\n";
    //			//(*findCelledge).Print();
    //		}
    //
    //	}
	
	
}

void Tissue::RemoveCell(Cell *c)
{
	ListCell.remove(*c);
}

Vertex* Tissue::Apoptosis(Cell *c)
{
    Vertex* v = c->Apoptosis();
    ListCell.remove(*c);
    return v;
}

void Tissue::clear()
{
	ListVertex.clear();
	ListEdge.clear();
	ListCell.clear();
	
	currentVertices = std::list<Vertex*>();
	currentEdges = std::list<Edge*>();
	currentCells = std::list<Cell*>();
}

void Tissue::PrintEdges()
{
	std::list<Edge>::iterator it;
	for(it=ListEdge.begin();it!=ListEdge.end();it++)
	{
		Edge e = (*it);
		it->Print();
	}
}

void Tissue::PrintCells()
{
	std::list<Cell>::iterator it;
	for(it=ListCell.begin();it!=ListCell.end();it++)
	{
		it->Print();
	}
}

Cell* Tissue::closestCell(const Point &ClickedPoint, std::list<Cell*> cellList)
{
	// take the distance to the projected barycenter of the apical side of the cell as reference!
	
	std::list<Cell*>::iterator it=cellList.begin();
	
	Cell* currentCell=(*it);
	
    Point BC_a = (*it)->BC_a;
    
	double distancemin=ClickedPoint.Distance(BC_a);
	it++;
	
	for(;it!=cellList.end();it++)
	{
        BC_a = (*it)->BC_a;
		double distance=ClickedPoint.Distance(BC_a);
		if(distance<distancemin) {currentCell=(*it);distancemin=distance;}
	}
    
	return currentCell;
}

Edge* Tissue::closestEdge(const Point &ClickedPoint)
{
    Edge* CurrentClosestEdge = &ListEdge.back();
    
    Point proj_MP_a = CurrentClosestEdge->MP_a();
    proj_MP_a.z = 0;
    
    double min_distance=ClickedPoint.Distance(proj_MP_a);
    
    for(std::list<Edge>::iterator it=++ListEdge.begin();it!=ListEdge.end();it++)
    {
        Point proj_MP_a = it->MP_a();
        proj_MP_a.z = 0;
        
        double new_distance=ClickedPoint.Distance(proj_MP_a);
        
        if(new_distance<min_distance)
        {
            min_distance=new_distance;
            CurrentClosestEdge=&(*it);
        }
    }
    
    return CurrentClosestEdge;
}

Vertex* Tissue::closestVertex(const Point &ClickedPoint)
{
    
    Vertex* CurrentClosestVertex = &ListVertex.back();
    
    Point proj_MP_a = CurrentClosestVertex->coord;
    proj_MP_a.z = 0;
    
    double min_distance=ClickedPoint.Distance(proj_MP_a);

    for(std::list<Vertex>::iterator itvertex=ListVertex.begin();itvertex!=ListVertex.end();itvertex++)
    {
        proj_MP_a = itvertex->coord;
        proj_MP_a.z = 0;

        if(ClickedPoint.Distance(proj_MP_a)<min_distance)
        {
            min_distance=ClickedPoint.Distance(proj_MP_a);
            CurrentClosestVertex=&(*itvertex);
        }
    }
    
    return CurrentClosestVertex;
}

// check!
int Tissue::testTissue()
{
	// the tissue structure has to be tested if it is consistent with the idea of the cells
	//std::list<Cell>::iterator itcell;
    //
    //	for(itcell=ListCell.begin();itcell!=ListCell.end();itcell++)
    //	{
    //
    //		// what should be tested?
    //		std::list<Edge*>::iterator it;
    //
    //		for(it=(*itcell).ListEdges.begin();it!=(*itcell).ListEdges.end();it++)
    //		{
    //			// are the addresses of first vertices (apical and basal) of the next edge (i.e. v1) the last of that one (v2)?
    //			if(!(**it).v2==(**(it+1)).v1) return (-1);
    //		}
    //
    //		//
    //		//if(!((**(*itcell).ListEdges.begin()).v1)==((**((*itcell).ListEdges.end())).v2)) return (-1);
    //	}
    
	// ...
	// does any edge contain twice the same cells as reference?
	//for(std::<Edge>::iterator it=it
	
	return 1;
}

int Tissue::loadFileFromTxt(const QString &fileName)
{
	bool isNumber=true;
	boundary_isSet = false;
	
    STIFF_BASAL_MEMBRANE = false;
    apical_wall = false;
    basal_wall = false;
    flatECM_basalPosition_bothways = false;
    onlyApicalContributions = false;
    
    int numberOfCellTypes = 4;
	
	// open file and if thats not possible throw error
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text))
	{
        //QMessageBox::warning(this, tr("Application"),tr("Cannot read file %1:\n%2.").arg(fileName).arg(file.errorString()));
        
        std::cout << "Cannot read file " << fileName.toStdString() << "!!" << std::endl;
        return -1;
    }
	
	// erase existing tissue
  	clear();
    
    // set anisotropic line tensions to zero
    anisotropicLineTension = false;
    for(unsigned i = 0; i<=NO_CELL_TYPES; i++)
        for(unsigned j = 0; j<=NO_CELL_TYPES; j++)
                MechProp.Intercell_ApicalLineTension_anisotropic[i][j] = 0;
  	
    P0_Lumen_a = 0;
    K_Lumen_a = 0;
    P0_Lumen_b = 0;
    K_Lumen_b = 0;
    ellipsoid_K = 0;
    
	// set default tissue to be periodic
	SystemSize=Quadrant();
	isPeriodic=false;
	
	//begin to read the file
	QTextStream in(&file);
	QString line;
	line=in.readLine();
	
	bool vertices_isset(false), edges_isset(false), cells_isset(false),celltype_isset(false);
	
	// iterate through file and look out for key words
	
	do
	{
        
		/*********************************************** SystemSize *************************************************/
		// if in the first passage the SystemSize is defined
		// be aware: periodicity just in y direction is not allowed, could be done later though
		if((QString::compare(line, "periodicity", Qt::CaseInsensitive))==0)
		{
			line=in.readLine();
			QTextStream lineStream(&line);
			QString QSproperty;
			lineStream>>QSproperty;
			if(QSproperty=="Lx")
			{
				QString QSvalue;
				lineStream>>QSvalue;
				bool isDouble;
				double value=QSvalue.toDouble(&isDouble);
				
				if(isDouble && value>0)
				{
					SystemSize.x=value;
					isPeriodic=true;
					
					line=in.readLine();
					
					QTextStream lineStream2(&line);
					lineStream2>>QSproperty;
					
					if(QSproperty=="Ly")
					{
						lineStream2>>QSvalue;
						value=QSvalue.toDouble(&isDouble);
						if(isDouble && value>0)	SystemSize.y=value;
						else SystemSize.y=-1;
						line=in.readLine();
					}
					
				}
				else // no periodicity in the tissue
				{
					SystemSize=Quadrant();
					isPeriodic=false;
				}
			}
            
            QTextStream lineStream3(&line);
            lineStream3>>QSproperty;
            
            if(QSproperty=="T_ext")
			{
				QString QSvalue;

                lineStream3>>QSvalue;
				bool isDouble;
				double value=QSvalue.toDouble(&isDouble);
				
				if(isDouble){
					T_ext = value;
					line=in.readLine();
				}
                                        
			}
		}

        if((QString::compare(line, "lumen_a", Qt::CaseInsensitive))==0)
		{
			line=in.readLine();
			QTextStream lineStream(&line);
			QString QSproperty;
			lineStream>>QSproperty;
			if(QSproperty=="V0")
			{
				QString QSvalue;
				lineStream>>QSvalue;
				bool isDouble;
				double value=QSvalue.toDouble(&isDouble);
				
				if(isDouble)
                {
                    V0_Lumen_a=value;
					line=in.readLine();
                }
			}
            
            QTextStream lineStream2(&line);
            lineStream2>>QSproperty;
            
            if(QSproperty=="K")
			{
				QString QSvalue;
                
                lineStream2>>QSvalue;
				bool isDouble;
				double value=QSvalue.toDouble(&isDouble);
				
				if(isDouble){
					K_Lumen_a = value;
					line=in.readLine();
				}
                
			}
            
            QTextStream lineStream3(&line);
            lineStream3>>QSproperty;
            
            if(QSproperty=="P0")
			{
				QString QSvalue;
                
                lineStream3>>QSvalue;
				bool isDouble;
				double value=QSvalue.toDouble(&isDouble);
				
				if(isDouble){
					P0_Lumen_a = value;
					line=in.readLine();
				}

			}
            
            signal_changeLumen_a(V0_Lumen_a,K_Lumen_a, P0_Lumen_a);
		}

        if((QString::compare(line, "lumen", Qt::CaseInsensitive))==0 || (QString::compare(line, "lumen_b", Qt::CaseInsensitive))==0)
		{
			line=in.readLine();
			QTextStream lineStream(&line);
			QString QSproperty;
			lineStream>>QSproperty;
			if(QSproperty=="V0")
			{
				QString QSvalue;
				lineStream>>QSvalue;
				bool isDouble;
				double value=QSvalue.toDouble(&isDouble);
				
				if(isDouble)
                {
                    V0_Lumen_b=value;
					line=in.readLine();
                }
			}
            
            QTextStream lineStream2(&line);
            lineStream2>>QSproperty;
            
            if(QSproperty=="K")
			{
				QString QSvalue;
                
                lineStream2>>QSvalue;
				bool isDouble;
				double value=QSvalue.toDouble(&isDouble);
				
				if(isDouble){
					K_Lumen_b = value;
					line=in.readLine();
				}
                
			}
            
            QTextStream lineStream3(&line);
            lineStream3>>QSproperty;
            
            if(QSproperty=="P0")
			{
				QString QSvalue;
                
                lineStream3>>QSvalue;
				bool isDouble;
				double value=QSvalue.toDouble(&isDouble);
				
				if(isDouble){
					P0_Lumen_b = value;
					line=in.readLine();
				}
                
			}

            
            signal_changeLumen_b(V0_Lumen_b,K_Lumen_b, P0_Lumen_b);
		}
        
        /*********************************************** Ellipsoidal Confinement *************************************************/
        if((QString::compare(line, "Ellipsoid", Qt::CaseInsensitive))==0)
		{
			line=in.readLine();
			QTextStream lineStream(&line);
			QString QSproperty;
			lineStream>>QSproperty;
			if(QSproperty=="x")
			{
				QString QSvalue;
				lineStream>>QSvalue;
				bool isDouble;
				double value=QSvalue.toDouble(&isDouble);
				
				if(isDouble)
                {
                    ellipsoid_ex=value;
					line=in.readLine();
                }
			}

            
            QTextStream lineStream2(&line);
            lineStream2>>QSproperty;
			if(QSproperty=="y")
			{
				QString QSvalue;
				lineStream2>>QSvalue;
				bool isDouble;
				double value=QSvalue.toDouble(&isDouble);
				
				if(isDouble)
                {
                    ellipsoid_ey=value;
					line=in.readLine();
                }
			}

            QTextStream lineStream3(&line);
            lineStream3>>QSproperty;
			if(QSproperty=="z")
			{
				QString QSvalue;
				lineStream3>>QSvalue;
				bool isDouble;
				double value=QSvalue.toDouble(&isDouble);
				
				if(isDouble)
                {
                    ellipsoid_ez=value;
					line=in.readLine();
                }
			}
            
            QTextStream lineStream4(&line);
            lineStream4>>QSproperty;
			if(QSproperty=="K")
			{
				QString QSvalue;
				lineStream4>>QSvalue;
				bool isDouble;
				double value=QSvalue.toDouble(&isDouble);
				
				if(isDouble)
                {
                    ellipsoid_K=value;
					line=in.readLine();
                }
			}
            signal_changeEllipsoid(ellipsoid_ex, ellipsoid_ey, ellipsoid_ez, ellipsoid_K);
		}

		/*********************************************** FLAT ECM CONFINEMENT *************************************************/
        if((QString::compare(line, "flatecm", Qt::CaseInsensitive))==0)
		{
			line=in.readLine();
			QTextStream lineStream(&line);
			QString QSproperty;
			lineStream>>QSproperty;
			if(QSproperty=="a")
			{
				QString QSvalue;
				lineStream>>QSvalue;
				bool isDouble;
				double value=QSvalue.toDouble(&isDouble);
				
				if(isDouble)
                {
                    flatECM_apicalPosition=value;
					line=in.readLine();
                }
			} else {
                std::cout << "Error in reading flatECM properties";
                return -1;
            }
            
            
            QTextStream lineStream2(&line);
            lineStream2>>QSproperty;
			if(QSproperty=="b")
			{
				QString QSvalue;
				lineStream2>>QSvalue;
				bool isDouble;
				double value=QSvalue.toDouble(&isDouble);
				
				if(isDouble)
                {
                    flatECM_basalPosition=value;
					line=in.readLine();
                }
			} else {
                std::cout << "Error in reading flatECM properties";
                return -1;
            }
            
            QTextStream lineStream3(&line);
            lineStream3>>QSproperty;
			if(QSproperty=="K")
			{
				QString QSvalue;
				lineStream3>>QSvalue;
				bool isDouble;
				double value=QSvalue.toDouble(&isDouble);
				
				if(isDouble)
                {
                    flatECM_stiffness=value;
					line=in.readLine();
                }
			} else {
                std::cout << "Error in reading flatECM properties";
                return -1;
            }
            
            signal_changeFlatECM(flatECM_apicalPosition, flatECM_basalPosition, flatECM_stiffness);
        }
        
        else if((QString::compare(line, "flatECM_basalPosition_bothways", Qt::CaseInsensitive))==0)
        {
            flatECM_basalPosition_bothways = true;
            line=in.readLine();
        }

        else if((QString::compare(line, "flatECM_stiffness_alpha_nonCell1", Qt::CaseInsensitive))==0)
        {
            line=in.readLine();
            QTextStream lineStream(&line);
            lineStream>>flatECM_stiffness_alpha_nonCell1;
            line=in.readLine();
        }
        
        else if((QString::compare(line, "onlyApicalContributions", Qt::CaseInsensitive))==0)
        {
            onlyApicalContributions = true;
            line=in.readLine();
        }
        
        
		/*********************************************** Cell Types *************************************************/
		else if((QString::compare(line, "Cell Types", Qt::CaseInsensitive))==0)
		{
            celltype_isset=true;
            
            // read properties of cell types 1 and 2            
            QString propertyValue;
            double Ta, Tb, K, V0, A0_a, A0_b, T0_a, T0_b;
            bool isDouble;
            
            line = in.readLine();
            for(int i=1;i<numberOfCellTypes+1;i++){ // do for every cell type
                QTextStream lineStream(&line);

                lineStream>>propertyValue;
                Ta = propertyValue.toDouble(&isDouble);
                MechProp.CellPropVector[i].T_a = Ta;
                
                lineStream>>propertyValue;
                Tb = propertyValue.toDouble(&isDouble);
                MechProp.CellPropVector[i].T_b = Tb;

                lineStream>>propertyValue;
                V0 = propertyValue.toDouble(&isDouble);
                MechProp.CellPropVector[i].V0 = V0;

                lineStream>>propertyValue;
                K = propertyValue.toDouble(&isDouble);
                MechProp.CellPropVector[i].K = K;
                
                // check if lineStream reached the end of the line and if not read in the parameters for the area elasticity
                if(!lineStream.atEnd())
                {
                    lineStream>>propertyValue;
                    A0_a = propertyValue.toDouble(&isDouble);
                    MechProp.CellPropVector[i].A0_a = A0_a;
                    lineStream>>propertyValue;
                    T0_a = propertyValue.toDouble(&isDouble);
                    MechProp.CellPropVector[i].T0_a = T0_a;
                    lineStream>>propertyValue;
                    A0_b = propertyValue.toDouble(&isDouble);
                    MechProp.CellPropVector[i].A0_b = A0_b;
                    lineStream>>propertyValue;
                    T0_b = propertyValue.toDouble(&isDouble);
                    MechProp.CellPropVector[i].T0_b = T0_b;
                } else {
                    MechProp.CellPropVector[i].A0_a = 0;
                    MechProp.CellPropVector[i].T0_a = Ta;
                    MechProp.CellPropVector[i].A0_b = 0;
                    MechProp.CellPropVector[i].T0_b = Tb;
                }
                
                
                    line = in.readLine();
            }
            
            // read in the intercell surface tensions
            for(int i=1;i<numberOfCellTypes+1;i++){ // do for every cell type
                QTextStream lineStream(&line);
                
                for(int j=0;j<=i;j++) {
                    lineStream>>propertyValue;
                    double T_l = propertyValue.toDouble(&isDouble);
                    MechProp.Intercell_SurfaceTension[i][j]=T_l;
                }
                
                line = in.readLine();
            }
        
            // read in the apical line tensions
            for(int i=1;i<numberOfCellTypes+1;i++){ // do so for every cell type
                QTextStream lineStream(&line);
                
                for(int j=0;j<=i;j++) { // up to the same cell type
                    lineStream>>propertyValue;
                    double G_a = propertyValue.toDouble(&isDouble);
                    MechProp.Intercell_ApicalLineTension[i][j]=G_a;
                }
                
                line = in.readLine();
            }
            
            // read in the basal line tensions
            for(int i=1;i<numberOfCellTypes+1;i++){ // do so for every cell type
                QTextStream lineStream(&line);
                
                for(int j=0;j<=i;j++) { // up to the same cell type
                    lineStream>>propertyValue;
                    double G_b = propertyValue.toDouble(&isDouble);
                    MechProp.Intercell_BasalLineTension[i][j]=G_b;
                }
                
                line = in.readLine();
            }
        }
        
        /*********************************************** anisotropic line tensions *************************************************/
        else if((QString::compare(line, "anisotropicLineTension", Qt::CaseInsensitive))==0)
        {
            anisotropicLineTension = true;
            
            line = in.readLine();
            QString propertyValue;
            bool isDouble;
            // read in the intercell surface tensions
            for(int i=1;i<numberOfCellTypes+1;i++){ // do for every cell type
                QTextStream lineStream(&line);
                
                for(int j=0;j<=i;j++) {
                    lineStream>>propertyValue;
                    double La_ai = propertyValue.toDouble(&isDouble);
                    MechProp.Intercell_ApicalLineTension_anisotropic[i][j]=La_ai;
                }
                
                line = in.readLine();
            }
        }
        
		/*********************************************** Basement Membrane *************************************************/
        else if((QString::compare(line, "ECM", Qt::CaseInsensitive))==0)
        {
            line=in.readLine();
            QTextStream lineStream(&line);
            double hasECM;
            lineStream>>hasECM;
            STIFF_BASAL_MEMBRANE = hasECM==1;
            if(STIFF_BASAL_MEMBRANE) {
                basal_wall = false;
                apical_wall = false;
            }
            line=in.readLine();
        }

        /*********************************************** Basal Wall *************************************************/
        else if((QString::compare(line, "basalWall", Qt::CaseInsensitive))==0)
        {
            line=in.readLine();
            basal_wall = true;
            STIFF_BASAL_MEMBRANE = false;
            QTextStream lineStream(&line);
            lineStream>>basal_wall_position;
            line=in.readLine();
        }
        
        /*********************************************** Apical Wall *************************************************/
        else if((QString::compare(line, "apicalWall", Qt::CaseInsensitive))==0)
        {
            line=in.readLine();
            apical_wall = true;
            QTextStream lineStream(&line);
            lineStream>>apical_wall_position;
            line=in.readLine();
        }
        
		/*********************************************** Vertices *************************************************/
		else if((QString::compare(line, "Vertices", Qt::CaseInsensitive))==0)
		{
			vertices_isset=true;
            
			isNumber=true;
			while(isNumber)
			{
				// get next line of the file
				line=in.readLine();
				
				// define dummies for parameters
				double coordx;
				double coordy;
				double coordz;
				double coordx_b;
				double coordy_b;
				double coordz_b;
				QString dummy;
				//double LineTension;
				
				QString word;
				
				// read in the next line
				QTextStream lineStream(&line);//cree un stream pour la ligne en question
				lineStream>>word;
				lineStream>>coordx;
				lineStream>>coordy;
				lineStream>>coordz;
				lineStream>>coordx_b;
				lineStream>>coordy_b;
				lineStream>>coordz_b;
				
				if(SystemSize.x>0)
				{
					// coordinates must not be greater than the periodic limits
					if(coordx>SystemSize.x | coordx_b>SystemSize.x)	return -3;
				}
				if(SystemSize.y>0)
				{
					// coordinates must not be greater than the periodic limits
					if(coordy>SystemSize.y | coordy_b>SystemSize.y)
						return -4;
				}
				
				int Number=word.toInt(&isNumber,10);
				
				if(!isNumber) break;
				
				// define the vertex using the streamed coordinates
				Vertex* v=new Vertex(this,Point(coordx,coordy,coordz),Point(coordx_b,coordy_b,coordz_b));
				
				// set number of the vertex
				v->Number = Number;
                
                //std::cout << Number << std::endl;
				
				// now read in the remaining properties
				
				QString property;
				QString propertyValue;
				double propertyValue2double;
				bool isDouble;
				
				do
				{
					if(lineStream.atEnd()) break;

					lineStream>>property;
                    
					if(property=="fixed")
					{
						v->xafixed=true;
						v->xbfixed=true;
						v->yafixed=true;
						v->ybfixed=true;
						v->zafixed=true;
						v->zbfixed=true;
					}
					if(property=="xfixed")
					{
						v->xafixed=true;
						v->xbfixed=true;
						continue;
					}
					if(property=="yfixed")
					{
						v->yafixed=true;
						v->ybfixed=true;
						continue;
					}
					if(property=="zfixed")
					{
						v->zafixed=true;
						v->zbfixed=true;
						continue;
					}
					
					if(property=="xafixed")
					{
						v->xafixed=true;
						continue;
					}
					if(property=="yafixed")
					{
						v->yafixed=true;
						continue;
					}
					if(property=="zafixed")
					{
						v->zafixed=true;
						continue;
					}
					if(property=="xbfixed")
					{
						v->xbfixed=true;
						continue;
					}
					if(property=="ybfixed")
					{
						v->ybfixed=true;
						continue;
					}
					if(property=="zbfixed")
					{
						v->zbfixed=true;
						continue;
					}
					if(isPeriodic)
					{
						if(property=="cb_left" | property=="left" | property=="l")
						{
							v->q_b.x=-1;
							continue;
						}
						if(property=="cb_right" | property=="right" | property=="r")
						{
							v->q_b.x=1;
							continue;
						}
						if(property=="cb_top" | property=="top" | property=="t")
						{
							v->q_b.y=1;
							continue;
						}
						if(property=="cb_bottom" | property=="bottom" | property=="b")
						{
							v->q_b.y=-1;
							continue;
						}
                    }
					
					if(lineStream.atEnd())
					{
						// unknown property, that is not followed by a value
						return -20;
						//break;
					}
					
					
					lineStream>>propertyValue;
					
					propertyValue2double=propertyValue.toDouble(&isDouble);
            
					// check if the string after the property is a double
					if(!isDouble)
					{
						//QMessageBox::warning(this, tr("Application"),tr("Property value cannot be identified"));
                        return -5;
					}
					
					if(property=="LineTension") v->G_l = propertyValue2double;
					else
					{
						//QMessageBox::warning(this, tr("Application"),tr("Unknown Property Name"));
						return -6;
					}
					
				}
				while(1);
				
				
				// ... and push the new vertex into the list of vertices
				AddVertex(v);
			}
		}
        
        
		/*********************************************** Edges *************************************************/
		else if((QString::compare(line, "Edges", Qt::CaseInsensitive))==0)
		{
			if(!vertices_isset) return -7; //vertices have to be set before edges can be defined
			
			edges_isset=true;
			
			Vertex* vertex1;
			Vertex* vertex2;
			QString word;
			
			//bool TensionsWritten=false;
			while(true)
			{
				// read next line
				line=in.readLine();
				
				Edge* e=new Edge(this);
				
				// initialise variables
				QString NumberVertex1_s;
				QString NumberVertex2_s;
				
				double NumberVertex1;
				double NumberVertex2;
				
				int numberEdge;
				
				// create stream
				QTextStream lineStream(&line);
				
				// first numbers are : number, vertex1, vertex2
				
				lineStream>>word;
				numberEdge=word.toInt(&isNumber,10);
				
				e->Number = numberEdge;
				
				if(!isNumber) break;
				
				lineStream>>NumberVertex1_s;
				lineStream>>NumberVertex2_s;
				
				NumberVertex1=NumberVertex1_s.toDouble(&isNumber);
				if(!isNumber)
				{
					//QMessageBox::warning(this, tr("Application"),tr("Number of vertex cannot be identified"));
					return -8;
				}
				
				NumberVertex2=NumberVertex2_s.toDouble(&isNumber);
				if(!isNumber)
				{
					//QMessageBox::warning(this, tr("Application"),tr("Number of vertex cannot be identified"));
					return -9;
				}
				
				QString property;
				QString propertyValue;
				double propertyValue2double;
				bool isDouble;
				
				e->crossBoundary = false;
				
				lineStream>>property;
				
				do
				{
					if(isPeriodic)
					{
						if(property=="cb_left" | property=="left" | property=="l")
						{
							e->q_v2.x=-1;
							lineStream>>property;
							e->crossBoundary = true;
							if(!lineStream.atEnd())	continue;
						}
						if(property=="cb_right" | property=="right" | property=="r")
						{
							e->q_v2.x=1;
							lineStream>>property;
							e->crossBoundary = true;
							if(!lineStream.atEnd()) continue;
						}
						if(property=="cb_top" | property=="top" | property=="t")
						{
							e->q_v2.y=1;
							lineStream>>property;
							e->crossBoundary = true;
							if(!lineStream.atEnd())	continue;
						}
						if(property=="cb_bottom" | property=="bottom" | property=="b")
						{
							e->q_v2.y=-1;
							lineStream>>property;
							e->crossBoundary = true;
							if(!lineStream.atEnd())	continue;
						}
//                        if(property=="q_c1")
//                        {
//                            lineStream>>property;
//                            e->q_c1.x = property.toInt(&isNumber);
//                            lineStream>>property;
//                            e->q_c1.y = property.toInt(&isNumber);
//                            lineStream>>property;
//                            if(!lineStream.atEnd())	continue;
//                        }
//                        if(property=="q_c2")
//                        {
//                            lineStream>>property;
//                            e->q_c2.x = property.toInt(&isNumber);
//                            lineStream>>property;
//                            e->q_c2.y = property.toInt(&isNumber);
//                            lineStream>>property;
//                            if(!lineStream.atEnd())	continue;
//                        }
//                        if(property=="q_BC_l")
//                        {
//                            lineStream>>property;
//                            e->q_BC_l.x = property.toInt(&isNumber);
//                            lineStream>>property;
//                            e->q_BC_l.y = property.toInt(&isNumber);
//                            lineStream>>property;
//                            if(!lineStream.atEnd())	continue;
//                        }
                        if(property=="c1r")
                        {
                            e->q_c1.x = 1;
                            e->crossBoundary = true;
                            lineStream>>property;
                            if(!lineStream.atEnd())	continue;
                        }
                        if(property=="c1l")
                        {
                            e->q_c1.x = -1;
                            e->crossBoundary = true;
                            lineStream>>property;
                            if(!lineStream.atEnd())	continue;
                        }
                        if(property=="c1t")
                        {
                            e->q_c1.y = 1;
                            e->crossBoundary = true;
                            lineStream>>property;
                            if(!lineStream.atEnd())	continue;
                        }
                        if(property=="c1b")
                        {
                            e->q_c1.y = -1;
                            e->crossBoundary = true;
                            lineStream>>property;
                            if(!lineStream.atEnd())	continue;
                        }
                        if(property=="c2r")
                        {
                            e->q_c2.x = 1;
                            e->crossBoundary = true;
                            lineStream>>property;
                            if(!lineStream.atEnd())	continue;
                        }
                        if(property=="c2l")
                        {
                            e->q_c2.x = -1;
                            e->crossBoundary = true;
                            lineStream>>property;
                            if(!lineStream.atEnd())	continue;
                        }
                        if(property=="c2t")
                        {
                            e->q_c2.y = 1;
                            e->crossBoundary = true;
                            lineStream>>property;
                            if(!lineStream.atEnd())	continue;
                        }
                        if(property=="c2b")
                        {
                            e->q_c2.y = -1;
                            e->crossBoundary = true;
                            lineStream>>property;
                            if(!lineStream.atEnd())	continue;
                        }
                    }
					
					if(lineStream.atEnd()) break;
					
					
					//if(lineStream.atEnd())
                    //					{
                    //						// unknown property at the end of the line, that is not followed by a value
                    //						return -20;
                    //					}
                    //
					lineStream>>propertyValue;
					propertyValue2double=propertyValue.toDouble(&isDouble);
					
					if(!isDouble)
					{
						// unknown property without a following value
						return -10;
					}
					
					if(property=="Type") ;
					else if	(property=="LineTensionApical")
					{
						e->G_a = propertyValue2double;
					}
					else if(property=="LineTensionBasal")
					{
						e->G_b = propertyValue2double;
					}
					else if(property=="SurfaceTension")
					{
						e->T_l = propertyValue2double;
					}
					else
					{
						//QMessageBox::warning(this, tr("Application"),tr("Unknown Property Name"));
						return -21;
					}
					
					lineStream>>property;
					
				}
				while(1);
				
				
				
				// bools that make sure, that both vertices are in the list
				bool v1found=false;
				bool v2found=false;
				
				std::list<Vertex>::iterator it;
				
				// check if vertices are in the list
				for(it=ListVertex.begin();it!=ListVertex.end();++it)
				{
					if (it->getNumber()==NumberVertex1) {vertex1=&(*it);/*it->ListEdges.push_back(e);*/ v1found=true;}
					
					else if (it->getNumber()==NumberVertex2) {vertex2=&(*it);/*it->ListEdges.push_back(e);*/v2found=true;}
					
					if(v1found && v2found) break;
				}
				
				// edge refers to a nonexisting vertex!
				if(!(v1found&&v2found))
				{
					//QMessageBox::warning(this, tr("Application"),tr("Unknown Vertex"));
					return -10;
				}
				
				// set number of edge
				e->setNumber(numberEdge);
				
				// construct edge from Vertices
				e->setVertex1(vertex1);
				e->setVertex2(vertex2);
				
				e->CrossBoundary();
				
				//e->Print();
				
				AddEdge(e);
			}
		}
		
		/*********************************************** Cells *************************************************/
		else if((QString::compare(line, "Cells", Qt::CaseInsensitive))==0)
		{
			if(!edges_isset) return -8; // the definition of the cells requires the definition of the edges
			cells_isset=true;
			QString dummy;
			
			while(true)
			{
				line=in.readLine();
				
				QTextStream lineStream2(&line);
				
				QString word;
				lineStream2>>word;
				int numberCell=word.toInt(&isNumber,10);
				
				if(!isNumber) break;
				
				// construct cell
				Cell cinter;
				// add cell cinter to tissue T
				AddCell(&cinter);
				Cell* c=&(ListCell.front());
				
				// set cell number
				c->Number=numberCell;
				
				c->Type = 1;
				
				// set default elastic mudulus
				c->K=1;//tissueWindow->globalK;
				
				// read in the numbers of the corresponding edges
				while(!(lineStream2.atEnd()))
				{
					// read next word
					lineStream2>>dummy;
					
					// check if the next number is an int
					int numberEdge=dummy.toInt(&isNumber,10);
					
					if(!isNumber) break;
					
                    //					if(numberCell==1)
                    //					{
                    //						int a=1+1;
                    //					}
					
					//iterate through all found edges and check, if an edge has the same number
					std::list<Edge>::iterator it;
					// has the edge been found?
					bool found=false;
					for(it=ListEdge.begin();it!=ListEdge.end();it++)
					{
						// if absolute value of number equals the number of the edge, then the edge is part of the cell
						//std::cout<<it->getNumber();
						if (it->getNumber()==abs(numberEdge))
						{
							// assign edge to cell
							c->ListEdges.push_back(&(*it));
							found=true;
							
							if(it->CrossBoundary())
							{
								c->crossBoundary=true;
							}
							
							// assign cell to edge
							if(numberEdge>0)
							{
								// assign cell to to the edge
								it->setCell1(c);
								
								//std::cout << it->c1->Number << std::endl;
								break;
							}
							else
							{
								it->setCell2(c);
								break;
							}
						}
					}
					
					if(!found)
					{
						//QMessageBox::warning(this, tr("Application"),tr("Could not find the edge!").arg(fileName).arg(file.errorString()));
						return -12;
					}
				}
				
				
				if(c->ListEdges.size()<3)
				{
					//QMessageBox::warning(this, tr("Application"),tr("Every Edge has to contain at least 3 edges.").arg(fileName).arg(file.errorString()));
					return -13;
				}
				
				// define bool, that saves if the reading in has worked
				//bool hasWorked(true);
				// contains the value of the parameter
				//double value;
				//int typeNumber;
				
				// now the edges are over and the other properties are read in
				// first set them to "standard values", for the case that they are not provided in the file
				c->V0=0;
				//c->Type=1;
				c->K=1;
				c->Ta=1;
				c->Tb=1;
				
				QString property=dummy;
				QString propertyValue;
				double propertyValue2double;
				bool isDouble;
				
				while(!lineStream2.atEnd())
				{
					lineStream2>>propertyValue;
					
					propertyValue2double=propertyValue.toDouble(&isDouble);
					
					// check if the string after the property name is a double
					if(!isDouble)
					{
						//QMessageBox::warning(this, tr("Application"),tr("Property value cannot be identified"));
						return -14;
					}
					bool V0set=false;
					if(property=="V0")
					{
						V0set = true;
						c->V0=propertyValue2double;
					}
					//else if(property=="Type") c->Type=(int)propertyValue2double;
					else if(property=="K") c->K=propertyValue2double;
					else if(property=="SurfaceTension" || property=="SurfaceTensionApical") c->Ta=propertyValue2double;
					else if(property=="SurfaceTension_b"|| property=="SurfaceTensionBasal") c->Tb=propertyValue2double;
                    else if(property=="positionFixed") c->positionFixed = propertyValue2double>0;
					else if(property=="Type")
					{
						int typeNo=propertyValue2double;
						
						if(typeNo<1 | typeNo>numberOfCellTypes)
						{
							// type number is not available
							return -35;
						}
						
						c->assignCellType(typeNo);
						
						c->Type = typeNo;
						
					}
					else
					{
						//QMessageBox::warning(this, tr("Application"),tr("Unknown Property Name for cell /numberCell"));
						return -15;
					}
					
					lineStream2>>property;
				}
			}
		}
		
		else
		{
			// keyword not found
            
			return -22;
		}
		
	}
	while(!in.atEnd());
	
	//saveFile("/Users/silvanus/MPI/3Dshape/Codes/Log Files/Tissue.rtf");
	
	// print all the cells of the single edges
//	for(std::list<Edge>::iterator it_e = ListEdge.begin();it_e!=ListEdge.end();it_e++)
//	{
//		std::cout << (*it_e).Number << ": ";
//        
//		if((*it_e).c1)
//		{
//			std::cout << "1="<< (*it_e).c1->Number << " ";
//		}
//		
//		if((*it_e).c2)
//		{
//			std::cout << "2="<< (*it_e).c2->Number;
//		}
//		std::cout << std::endl;
//	}
	
	
	SetTissueStructure();
	
	// print all the cells of the single edges
    //	for(std::list<Edge>::iterator it_e = ListEdge.begin();it_e!=ListEdge.end();it_e++)
    //	{
    //		//std::cout << (*it_e).Number << ": ";
    //		if((*it_e).c1)
    //		{
    //			std::cout << (*it_e).c1->Number << " ";
    //		}
    //
    //		if((*it_e).c2)
    //		{
    //			std::cout << (*it_e).c2->Number;
    //		}
    //		//std::cout << std::endl;
    //	}
    //
    //
    //	std::list<Cell>::iterator it_cell = ListCell.begin();
    //
    //	for(;it_cell!=ListCell.end();it_cell++)
    //	{
    //
    //		//std::cout << (*it_cell).ListVertex.size() << std::endl;
    //		if((*it_cell).ListVertex.size() == 0)
    //		{
    //			std::list<Edge*>::iterator it_e = (*it_cell).ListEdges.begin();
    //			std::cout << std::endl << (*it_cell).Number;
    //
    //			for(;it_e!=(*it_cell).ListEdges.end();it_e++)
    //			{
    //				std::cout << std::endl << (*it_e)->Number;
    //			}
    //
    //		}
    //	}
	
	UpdateBarycenters();
	
	for(std::list<Cell>::iterator it_cell = ListCell.begin();it_cell!=ListCell.end();it_cell++)
	{
		//std::cout << (*it_cell).ListVertex.size();
		
		if((*it_cell).V0==0)
		{
            (*it_cell).UpdateVolume();
			(*it_cell).V0 = (*it_cell).Volume;
		}
	}
	
	setMechanicalProperties();
	
	//if(IntercellSurfaceTensions_isset) setIntercellSurfaceTensions();
	//if(IntercellApicalLineTensions_isset) setIntercellApicalLineTensions();
	//if(IntercellBasalLineTensions_isset) setIntercellBasalLineTensions();
	
	
	update();
	
	
	// give every vertex the information if at lies at the boundary
	for(std::list<Vertex>::iterator it_vertex=ListVertex.begin();it_vertex!=ListVertex.end();it_vertex++)
	{
		(*it_vertex).isAtBoundary = (*it_vertex).get_isAtBoundary();
	}
	
//	for(std::list<Edge>::iterator it_edge=ListEdge.begin();it_edge!=ListEdge.end();it_edge++)
//	{
//		Edge *e = &(*it_edge);
//		
//		//std::cout << (*it_edge).G_a << "  ";
//	}
    //update();
    
	return 1;
}


// set the properties of the cells and the edges depending on the type of the cell (apical & basal surface tension, pressure)
int Tissue::setMechanicalProperties()
{
	UpdateBarycenters();
	
	for(std::list<Cell>::iterator it_cell = ListCell.begin();it_cell!=ListCell.end();it_cell++)
	{
		it_cell->Ta = MechProp.CellPropVector[it_cell->Type].T_a;
		it_cell->Tb = MechProp.CellPropVector[it_cell->Type].T_b;
		it_cell->T0_a = MechProp.CellPropVector[it_cell->Type].T0_a;
		it_cell->T0_b = MechProp.CellPropVector[it_cell->Type].T0_b;
		it_cell->A0_a = MechProp.CellPropVector[it_cell->Type].A0_a;
		it_cell->A0_b = MechProp.CellPropVector[it_cell->Type].A0_b;
		it_cell->K = MechProp.CellPropVector[it_cell->Type].K;
		
		if(MechProp.CellPropVector[(*it_cell).Type].V0<=0)
		{
			(*it_cell).UpdateVolume();
			(*it_cell).V0 = (*it_cell).Volume;
		}
		else
		{
			(*it_cell).V0 = MechProp.CellPropVector[(*it_cell).Type].V0;
		}
	}
	
	bool cell1isset, cell2isset, cell1hastype, cell2hastype;
	
	for(std::list<Edge>::iterator it_edge = ListEdge.begin();it_edge!=ListEdge.end();it_edge++)
	{
		
		cell1isset = (it_edge->c1!=NULL);
		cell2isset = (it_edge->c2!=NULL);
		cell1hastype = cell1isset && (it_edge->c1->Type>0);
		cell2hastype = cell2isset && (it_edge->c2->Type>0);
		
		// if both cells have a type
		if(cell1hastype && cell2hastype)
		{
			//double a = IntercellSurfaceTensions[std::min((*it_edge).c1->Type, (*it_edge).c2->Type)][std::max((*it_edge).c1->Type, (*it_edge).c2->Type)+1];
			it_edge->T_l = MechProp.Intercell_SurfaceTension[std::max(it_edge->c1->Type, it_edge->c2->Type)][std::min(it_edge->c1->Type, it_edge->c2->Type)];
			it_edge->G_a = MechProp.Intercell_ApicalLineTension[std::max(it_edge->c1->Type, it_edge->c2->Type)][std::min(it_edge->c1->Type, it_edge->c2->Type)];
			it_edge->G_a_ai = MechProp.Intercell_ApicalLineTension_anisotropic[std::max(it_edge->c1->Type, it_edge->c2->Type)][std::min(it_edge->c1->Type, it_edge->c2->Type)];
            it_edge->G_b = MechProp.Intercell_BasalLineTension[std::max(it_edge->c1->Type, it_edge->c2->Type)][std::min(it_edge->c1->Type, it_edge->c2->Type)];
		}
		// if the edge only is adjacent to cell1 and cell1 has a type defined
		else if(cell1hastype && !cell2isset)
		{
			it_edge->T_l = MechProp.Intercell_SurfaceTension[it_edge->c1->Type][0];
			it_edge->G_a = MechProp.Intercell_ApicalLineTension[it_edge->c1->Type][0];
            it_edge->G_a_ai = MechProp.Intercell_ApicalLineTension_anisotropic[std::max(it_edge->c1->Type, it_edge->c2->Type)][std::min(it_edge->c1->Type, it_edge->c2->Type)];
			it_edge->G_b = MechProp.Intercell_BasalLineTension[it_edge->c1->Type][0];
		}
		// if the edge only is adjacent to cell2 and cell2 has a type defined
		else if(cell2hastype && !cell1isset)
		{
			it_edge->T_l = MechProp.Intercell_SurfaceTension[it_edge->c2->Type][0];
			it_edge->G_a = MechProp.Intercell_ApicalLineTension[it_edge->c2->Type][0];
			it_edge->G_b = MechProp.Intercell_BasalLineTension[it_edge->c2->Type][0];
            it_edge->G_a_ai = MechProp.Intercell_ApicalLineTension_anisotropic[it_edge->c2->Type][0];
		}
	}
	
	
	for(std::list<Vertex>::iterator it_vertex = ListVertex.begin();it_vertex!=ListVertex.end();it_vertex++)
	{
		it_vertex->G_l=MechProp.Vertex_LineTension;
	}
    
    updateMainWindowData();
    
    return 1;
	
}

// set the properties of the cells and the edges depending on the type of the cell (apical & basal surface tension, pressure)
//int Tissue::setMechanicalPropertiesWithoutV0()
//{
//	std::list<Cell>::iterator it_cell = ListCell.begin();
//	
//	UpdateBarycenters();
//	
//	for(;it_cell!=ListCell.end();it_cell++)
//	{
//		it_cell->Ta = MechProp.CellPropVector[it_cell->Type].T_a;
//		it_cell->Tb = MechProp.CellPropVector[it_cell->Type].T_b;
//        it_cell->T0_a = MechProp.CellPropVector[it_cell->Type].T0_a;
//		it_cell->T0_b = MechProp.CellPropVector[it_cell->Type].T0_b;
//        it_cell->A0_a = MechProp.CellPropVector[it_cell->Type].A0_a;
//		it_cell->A0_b = MechProp.CellPropVector[it_cell->Type].A0_b;
//
//		it_cell->K =   MechProp.CellPropVector[it_cell->Type].K;
//		
//	}
//	
//	std::list<Edge>::iterator it_edge = ListEdge.begin();
//	
//	bool cell1isset, cell2isset, cell1hastype, cell2hastype;
//	
//	for(;it_edge!=ListEdge.end();it_edge++)
//	{
//		
//		cell1isset = (it_edge->c1!=NULL);
//		cell2isset = (it_edge->c2!=NULL);
//		cell1hastype = cell1isset && (it_edge->c1->Type>0);
//		cell2hastype = cell2isset && (it_edge->c2->Type>0);
//		
//		// if both cells have a type
//		if(cell1hastype && cell2hastype)
//		{
//			//double a = IntercellSurfaceTensions[std::min((*it_edge).c1->Type, (*it_edge).c2->Type)][std::max((*it_edge).c1->Type, (*it_edge).c2->Type)+1];
//			it_edge->T_l = MechProp.Intercell_SurfaceTension[std::max(it_edge->c1->Type, it_edge->c2->Type)][std::min(it_edge->c1->Type, it_edge->c2->Type)];
//			it_edge->G_a = MechProp.Intercell_ApicalLineTension[std::max(it_edge->c1->Type, it_edge->c2->Type)][std::min(it_edge->c1->Type, it_edge->c2->Type)];
//			it_edge->G_b = MechProp.Intercell_BasalLineTension[std::max(it_edge->c1->Type, it_edge->c2->Type)][std::min(it_edge->c1->Type, it_edge->c2->Type)];
//		}
//		// if the edge only is adjacent to cell1 and cell1 has a type defined
//		else if(cell1hastype && !cell2isset)
//		{
//			it_edge->T_l = MechProp.Intercell_SurfaceTension[it_edge->c1->Type][0];
//			it_edge->G_a = MechProp.Intercell_ApicalLineTension[it_edge->c1->Type][0];
//			it_edge->G_b = MechProp.Intercell_BasalLineTension[it_edge->c1->Type][0];
//		}
//		// if the edge only is adjacent to cell2 and cell2 has a type defined
//		else if(cell2hastype && !cell1isset)
//		{
//			it_edge->T_l = MechProp.Intercell_SurfaceTension[it_edge->c2->Type][0];
//			it_edge->G_a = MechProp.Intercell_ApicalLineTension[it_edge->c2->Type][0];
//			it_edge->G_b = MechProp.Intercell_BasalLineTension[it_edge->c2->Type][0];
//		}
//	}
//	
//	
//	std::list<Vertex>::iterator it_vertex = ListVertex.begin();
//	for(;it_vertex!=ListVertex.end();it_vertex++)
//	{
//		it_vertex->G_l = 0;
//	}
//    
//    //signal_update_cyst_data(MechProp.CellPropVector[1].T_a,MechProp.CellPropVector[1].T_b,MechProp.Intercell_SurfaceTension[1][1],MechProp.CellPropVector[1].V0,MechProp.CellPropVector[1].K,MechProp.Intercell_ApicalLineTension[1][1],MechProp.Intercell_BasalLineTension[1][1],MechProp.Intercell_SurfaceTension[2][1],MechProp.Intercell_ApicalLineTension[2][1],MechProp.Intercell_BasalLineTension[2][1]);
//}
// function writes the current tissue into the file named fileName
int Tissue::writeFileToTxt(const char* fileName)
{
	//Open File;
	std::ofstream output;
	output.open(fileName);// initial position is set to be at the end of the file
	
    int numberOfCellTypes = 4;
    
	// Write into file
	if(output.is_open())
	{
        // Periodic Boundary Conditions
        if(isPeriodic)
        {
            output << "Periodicity\nLx " << std::setprecision(9) << SystemSize.x << std::endl << "Ly " << std::setprecision(9) << SystemSize.y << std::endl;
            output << "T_ext " << T_ext << std::endl;
        }
        
        // lumen elasticities
        if(K_Lumen_a>0 || P0_Lumen_a>0)
        {
            output << "Lumen_a\nV0 " << std::setprecision(10) << V0_Lumen_a << std::endl << "K " << K_Lumen_a << std::endl << "P0 " << P0_Lumen_a << std::endl;
        }
        if(K_Lumen_b>0 || P0_Lumen_b>0)
        {
            output << "Lumen_b\nV0 " << std::setprecision(10) << V0_Lumen_b << std::endl << "K " << K_Lumen_b << std::endl << "P0 " << P0_Lumen_b << std::endl;
        }

        // ellipsoidal confinement
        if(ellipsoid_ex>0 && ellipsoid_ez>0 && ellipsoid_K>0)
        {
            output << "Ellipsoid\nx " << std::setprecision(10) << ellipsoid_ex << std::endl << "y " << ellipsoid_ey << std::endl << "z " << ellipsoid_ez << "\nK " << ellipsoid_K << std::endl;
        }

        // flat ECM
        if(flatECM_stiffness>0)
        {
            output << "flatECM\na " << std::setprecision(10) << flatECM_apicalPosition << "\nb " << flatECM_basalPosition << "\nK " << flatECM_stiffness << std::endl;
            if(flatECM_basalPosition_bothways) output << "flatECM_basalPosition_bothways\n";
        }
        
        if(flatECM_stiffness_alpha_nonCell1!=1) {
            output << "flatECM_stiffness_alpha_nonCell1\n" << flatECM_stiffness_alpha_nonCell1 << std::endl;
        }        
        
        // only apical contributions
        if(onlyApicalContributions) {
            output << "onlyApicalContributions\n";
        }
        
        // definition of the cell types and the edges inbetween
        output << "Cell Types\n";// << numberOfCellTypes << std::endl;
        
        for(int i = 1; i<=numberOfCellTypes; i++) { // properties of cell types 1 and 2
            output<<MechProp.CellPropVector[i].T_a << " " << MechProp.CellPropVector[i].T_b << " " << MechProp.CellPropVector[i].V0 << " " << MechProp.CellPropVector[i].K << " " << MechProp.CellPropVector[i].A0_a << " " << MechProp.CellPropVector[i].T0_a << " " << MechProp.CellPropVector[i].A0_b << " " << MechProp.CellPropVector[i].T0_b << std::endl;
        }
        
        // properties of interface between cell types 1 and 2
        for(int i = 1; i<=numberOfCellTypes; i++) { // lateral surface tension
            for(int j=0; j<i; j++) {
                output<<MechProp.Intercell_SurfaceTension[i][j]<<" ";
            }
            output<<MechProp.Intercell_SurfaceTension[i][i]<<std::endl;
        }
        
        for(int i = 1; i<=numberOfCellTypes; i++) { // apical line tension
            for(int j=0; j<i; j++) {
                output<<MechProp.Intercell_ApicalLineTension[i][j]<<" ";
            }
            output<<MechProp.Intercell_ApicalLineTension[i][i]<<std::endl;
        }
        
        for(int i = 1; i<=numberOfCellTypes; i++) { // basal line tension
            for(int j=0; j<i; j++) {
                output<<MechProp.Intercell_BasalLineTension[i][j]<<" ";
            }
            output<<MechProp.Intercell_BasalLineTension[i][i]<<std::endl;
        }
        
        if(anisotropicLineTension) {
           output <<"anisotropicLineTension\n";
        
            for(int i = 1; i<=numberOfCellTypes; i++) { // basal line tension
                for(int j=0; j<i; j++) {
                    output<<MechProp.Intercell_ApicalLineTension_anisotropic[i][j]<<" ";
                }
                output<<MechProp.Intercell_ApicalLineTension_anisotropic[i][i]<<std::endl;
            }
        }
        
		output << "Vertices\n";
		
		std::list<Vertex>::iterator it;
		for(it=ListVertex.begin();it!=ListVertex.end();it++)
		{
			output << it->getNumber() << " " << std::setprecision(9) << it->getCoord().x << " "<< it->getCoord().y << " "<< it->getCoord().z << " "<< it->getBasalCoord().x << " "<< it->getBasalCoord().y << " "<< it->getBasalCoord().z << " LineTension " << it->G_l;
            if(it->q_b.x == -1) output << " l";
            else if(it->q_b.x == 1) output << " r";
            if(it->q_b.y == -1) output << " b";
            else if(it->q_b.y == 1) output << " t";
            output << std::endl;
		}
		
		output << "Edges\n";
		
		std::list<Edge>::iterator it_e;
		for(it_e=ListEdge.begin();it_e!=ListEdge.end();it_e++)
		{
			output << it_e->getNumber() << "  " << it_e->getVertex1()->getNumber() << " " << it_e->getVertex2()->getNumber();// << " LineTensionApical "<< it_e->G_a << " LineTensionBasal " << it_e->G_b << " SurfaceTension " << it_e->T_l;
            if(!isPeriodic)
            {
                output << std::endl;
            }
            else
            {
                if(it_e->q_v2.x==1) output << " r";
                else if(it_e->q_v2.x==-1) output << " l";
                if(it_e->q_v2.y==1) output << " t";
                else if(it_e->q_v2.y==-1) output << " b";
                
                // quadrants of cells wrt the edge
                if(it_e->q_c1.x==1) output << " c1l";
                if(it_e->q_c1.x==-1) output << " c1r";
                if(it_e->q_c1.y==1) output << " c1t";
                if(it_e->q_c1.y==-1) output << " c1b";
                if(it_e->q_c2.x==1) output << " c2l";
                if(it_e->q_c2.x==-1) output << " c2r";
                if(it_e->q_c2.y==1) output << " c2t";
                if(it_e->q_c2.y==-1) output << " c2b";
                
                output << std::endl;
            }
		}
		
		
		output << "Cells\n";
		std::list<Cell>::iterator it_c;
		for(it_c=ListCell.begin();it_c!=ListCell.end();it_c++)
		{
			// print number
			output << it_c->Number << " ";
			// print the edges with the corresponding signatures
			for(std::list<Edge*>::iterator it_edge = (it_c->ListEdges.begin());it_edge!=it_c->ListEdges.end();it_edge++)
			{
				// get direction of the edge with respect to the cell
				if((*it_edge)->c1 != &(*it_c))
					output << "-";
				
				output << (*it_edge)->getNumber() << " ";
			}
            
            if(it_c->positionFixed) output << "positionFixed 1 ";
			
			if(it_c->Type!=0) output << "Type " << it_c->Type << std::endl;
            else  output << "V0 " << it_c->V0 << " K " << it_c->K << " SurfaceTensionApical " << it_c->T_a() << " SurfaceTensionBasal " << it_c->T_b() << std::endl;
			
		}
		
        // output << "END OF FILE";
        
		output.close();
		
	}
	
	else
	{
		//QMessageBox::warning(this, tr("Application"),tr("Cannot write to file %1:\n%2.").arg(fileName).arg(file.errorString()));
	}
    
    return 1;
    
}


// function that removes an edge from the tissue - up to now just for the nonperiodic case!!
bool Tissue::PopAndUnpopAllEdges()
{
    double minLength = pop_minLength;
    double lengthNewEdge = pop_initLength;
    
    
    bool unpopTookPlace , min_1_unpopTookPlace (false);
    
    // iterate through all vertices and check if they should be unpopped
    // first create a list of pointers on the vertices in the tissue
    std::list<Vertex*> ptr_ListVertex;
    for(std::list<Vertex>::iterator it_v = ListVertex.begin(); it_v!=ListVertex.end(); it_v++)
    {
        ptr_ListVertex.push_back(&(*it_v));
    }
    
    for(std::list<Vertex*>::iterator it_v = ptr_ListVertex.begin(); it_v!=ptr_ListVertex.end(); it_v++)
    {
        min_1_unpopTookPlace  = min_1_unpopTookPlace  || (*it_v)->checkAndPop_forces(lengthNewEdge);
        //        min_1_popTookPlace  = min_1_popTookPlace  || (*it_v)->checkAndPop_energy(1);
    }
    
    do
    {
        unpopTookPlace  = false;
        
        // iterate through all vertices and check if they should be unpopped
        // first create a list of pointers on the vertices in the tissue
        std::list<Vertex*> ptr_ListVertex;
        for(std::list<Vertex>::iterator it_v = ListVertex.begin(); it_v!=ListVertex.end(); it_v++)
        {
            ptr_ListVertex.push_back(&(*it_v));
        }
        
        for(std::list<Vertex*>::iterator it_v = ptr_ListVertex.begin(); it_v!=ptr_ListVertex.end(); it_v++)
        {
            if((*it_v)->checkAndPop_forces(lengthNewEdge))
            {
                min_1_unpopTookPlace = true;
                unpopTookPlace = true;
                break;
            }
            
        }
    } while(unpopTookPlace);

    
    std::cout << "Checking for T1s!\n";
    // iterate through all the edges and check if they have to be removed because they are too small
    bool popTookPlace , min_1_popTookPlace (false);
    do
    {
        popTookPlace  = false;
        
        for(std::list<Edge>::iterator it_edge = ListEdge.begin();it_edge!=ListEdge.end();it_edge++) {
            
            Point vector_a = it_edge->v1->coord-it_edge->v2->coord;
            Point vector_b = it_edge->v1->coord_b-it_edge->v2->coord_b;
            
            if(it_edge->CrossBoundary()) continue;

            if(Point::Norm(vector_a) < minLength)// && Point::Norm(vector_b) < minLength )
            {
//                T1_remove(&(*it_edge));
//                currentEdges.remove(&(*it_edge));
//                popTookPlace  = true;
//
//                break;
                
//                Point force_a = it_edge->v2->dW_dva-it_edge->v1->dW_dva;
//                Point force_b = it_edge->v2->dW_dvb-it_edge->v1->dW_dvb;
//                
//                // if the edge will be further shortened apically and basally
//                if(vector_a.dot(force_a)>0 && vector_b.dot(force_b)>0) {
                    //double TINY = 1e-10;
                    //double NormBefore = Point::Norm(vector_a);
                    //double NormAfter = Point::Norm(vector_a+force_a*TINY);
                    //bool popTest = (NormBefore-NormAfter)>0;

                min_1_popTookPlace  = true;
                std::cout << "Edge " << it_edge->Number << " has been removed. With " << ListEdge.size() << " remaining.\n";

                if(it_edge->c1) it_edge->c1->ListEdges.remove(&(*it_edge));
                if(it_edge->c2) it_edge->c2->ListEdges.remove(&(*it_edge));
                
                T1_remove(&(*it_edge));
                
                currentEdges.remove(&(*it_edge));
                popTookPlace  = true;
                break;
            }
        }
    } while(popTookPlace);
    
//    // check if all edges are well defined!
//    for(std::list<Edge>::iterator ite = ListEdge.begin(); ite!=ListEdge.end();)
//    {
//        if(!ite->T)
//        {
//            ListEdge.erase(ite);
//        }
//        else
//        {
//            ite++;
//        }
//    }
//    
//    // check if all edges are well defined!
//    for(std::list<Cell>::iterator itc = ListCell.begin(); itc!=ListCell.end();)
//    {
//        if(!itc->T)
//        {
//            ListCell.erase(itc);
//        }
//        else
//        {
//            itc++;
//        }
//    }
//
	update();
    
    //update_pos(minPar.minimization_including_SystemSize);
    
    return min_1_popTookPlace || min_1_unpopTookPlace ;
}


void Tissue::T1_remove(Edge* e)
{
    e->T1_remove();
    ListEdge.remove(*e);
}

// check if there are any edges that have to be unpopped and unpop then
//bool Tissue::CheckForPop(double minLength)
//{
//	// iterate through all the cells and check if they have to be unpopped
//	std::list<Edge>::iterator it_edge;
//	if((*it_edge).ApicalLength() < minLength )
//	{
//        Edge* e = &(*it_edge);
//        
//        // if the edge will be further shortened
//        if((e->v1->ApicalForce-e->v2->ApicalForce).dot(e->v2->coord-e->v1->coord)>0)
//        {
//            PopEdge(&(*it_edge));
//            
//        }
//	}
//
//	return 1;
//}


//double Tissue::update_energy()
//{
//
//	Energy=0;
//
//	// iterate through the cells and add those contributions
//	std::list<Cell>::iterator it_cell;
//	for(it_cell=ListCell.begin();it_cell!=ListCell.end();it_cell++)
//	{
//		Energy+= (*it_cell).ElasticEnergy;
//	}
//
//	// iterate through the edges and add those contributions
//	std::list<Edge>::iterator it_edge;
//	for(it_edge=ListEdge.begin();it_edge!=ListEdge.end();it_edge++)
//	{
//		Energy += (*it_edge).update_energy();
//	}
//
//	// iterate through the cells and add those contributions
//	std::list<Vertex>::iterator it_vertex;
//	for(it_vertex=ListVertex.begin();it_vertex!=ListVertex.end();it_vertex++)
//	{
//		Energy += (*it_vertex).getLength()*(*it_vertex).LineTension;
//	}
//
//	return Energy;
//}

Vertex* Tissue::SimilarOrNewVertex(Point ApicalCoord, Point BasalCoord)
{
	if(!isPeriodic)
	{
		double epsilon = 10e-3;
		// iterate through all the vertices in the system and check if the vertex has the same coordinates
		std::list<Vertex>::iterator it = ListVertex.begin();
		for(;it!=ListVertex.end();it++)
		{
			Vertex* v = &(*it);
			if(Point::Norm(v->getCoord()-ApicalCoord)<epsilon && Point::Norm(v->getBasalCoord()-BasalCoord)<epsilon)
			{
				return v;
			}
		}
		
		
		// if the there doesnt exist a vertex with similar coordinates them create it
		Vertex* v=new Vertex(this,ApicalCoord,BasalCoord);
		
		AddVertex(v);
	}
	else
	{
		double epsilon = 10e-3;
		// iterate through all the vertices in the system and check if the vertex has the same coordinates also over the box boundaries
		std::list<Vertex>::iterator it = ListVertex.begin();
		for(;it!=ListVertex.end();it++)
		{
			Vertex* v = &(*it);
			if(Quadrant::Norm(v->getCoord()-ApicalCoord, SystemSize)<epsilon && Quadrant::Norm(v->getBasalCoord()-BasalCoord,SystemSize)<epsilon)
			{
				return v;
			}
		}
		
		
		// if the there doesnt exist a vertex with similar coordinates them create it
		Vertex* v=new Vertex(this,ApicalCoord,BasalCoord);
		
		// set the basal quadrant
		v->q_b = Quadrant();
		
		AddVertex(v);
	}
	
	return &ListVertex.front();
	
}


Edge* Tissue::SimilarOrNewEdge(Vertex* Ver1, Vertex* Ver2, bool* sameDirection, bool* newEdge)
{
	// iterate through all the edges in the system and check if another Edge has the same vertices
	std::list<Edge>::iterator it = ListEdge.begin();
	for(;it!=ListEdge.end();it++)
	{
		Edge* e = &(*it);
		if(e->getVertex1() == Ver1 && e->getVertex2() == Ver2)
		{
			*sameDirection=true;
			*newEdge=false;
			return e;
		}
		if(e->getVertex1() == Ver2 && e->getVertex2() == Ver1)
		{
			*sameDirection=false;
			*newEdge=false;
			return e;
		}
		
	}
	
	
	// if the there doesnt exist a edge with similar vertices then create it
	Edge* e=new Edge(this, Ver1, Ver2);
	*newEdge=true;
	*sameDirection=true;
	AddEdge(e);
	
	return &ListEdge.front();
	
}

Edge* Tissue::SimilarOrNewEdge(Vertex* Ver1, Vertex* Ver2, bool* sameDirection, bool* newEdge, Quadrant &quadrant_v2)
{
	// iterate through all the edges in the system and check if another Edge has the same vertices and orientation
	std::list<Edge>::iterator it = ListEdge.begin();
	for(;it!=ListEdge.end();it++)
	{
		Edge* e = &(*it);
		if(e->getVertex1() == Ver1 && e->getVertex2() == Ver2 && e->q_v2==quadrant_v2)
		{
			*sameDirection=true;
			*newEdge=false;
			return e;
		}
		if(e->getVertex1() == Ver2 && e->getVertex2() == Ver1 && e->q_v2==(quadrant_v2*(-1)))
		{
			*sameDirection=false;
			*newEdge=false;
			return e;
		}
		
	}
	
	
	// if the there doesnt exist an edge with similar vertices then create it
	Edge* e=new Edge(this, Ver1, Ver2);
	// set the quadrant of v2
	e->q_v2 = quadrant_v2;
	*newEdge=true;
	*sameDirection=true;
	AddEdge(e);
	
	return &ListEdge.front();
	
}

bool Tissue::createRectangularTissue(int rows, int columns, double height, double length, double width)
{
	// check the parameters:
	if(rows<1 || columns<1 || height<0 || length<0 || width<0)
		return false;
	
	// clear the tissue
	clear();
	
	//Cell::resetCounter();
	//Edge::resetCounter();
	//Vertex::resetCounter();
	
	Point basalDistance = Point(0,0,height);
	
	Point V1, V2, V3, V4, V5, V6, V1b, V2b, V3b, V4b, V5b, V6b;
	Vertex *Ver1, *Ver2, *Ver3, *Ver4;
	Edge *Edge1, *Edge2, *Edge3, *Edge4;
	bool sameDirection, sameEdge;
	
	
	// create all the cells
	for(int c_row = 1; c_row<=rows; c_row++)
	{
		for(int c_col = 1; c_col<=columns; c_col++)
		{
            // create the cell and add it to the tissue
			Cell* c = new Cell();
			AddCell(c);
			c = &(ListCell.front());
            
			c->Type = 1;
			
			// the square in the middle becomes different properties
			if(rows == 12 && columns==12)
			{
				if(c_row>4 && c_row<9 && c_col>4 && c_col<9)
					c->Type = 1;
				else
					c->Type = 2;
			}
			// create 4 vertices
			V1 = Point((c_col-1)*width,(c_row-1)*length,0);
			V2 = Point((c_col)*width,(c_row-1)*length,0);
			V3 = Point((c_col)*width,(c_row)*length,0);
			V4 = Point((c_col-1)*width,(c_row)*length,0);
			
			V1b = V1 - basalDistance;
			V2b = V2 - basalDistance;
			V3b = V3 - basalDistance;
			V4b = V4 - basalDistance;
			
			Ver1 = SimilarOrNewVertex(V1, V1b);
			Ver2 = SimilarOrNewVertex(V2, V2b);
			Ver3 = SimilarOrNewVertex(V3, V3b);
			Ver4 = SimilarOrNewVertex(V4, V4b);
			
			Edge1 = SimilarOrNewEdge(Ver1, Ver2, &sameDirection, &sameEdge);
			if(sameDirection) Edge1->setCell1(c);
			else Edge1->setCell2(c);
			c->addEdge(Edge1);
			
			Edge2 = SimilarOrNewEdge(Ver2, Ver3, &sameDirection, &sameEdge);
			if(sameDirection) Edge2->setCell1(c);
			else Edge2->setCell2(c);
			c->addEdge(Edge2);
			
			Edge3 = SimilarOrNewEdge(Ver3, Ver4, &sameDirection, &sameEdge);
			if(sameDirection) Edge3->setCell1(c);
			else Edge3->setCell2(c);
			c->addEdge(Edge3);
			
			Edge4 = SimilarOrNewEdge(Ver4, Ver1, &sameDirection, &sameEdge);
			if(sameDirection) Edge4->setCell1(c);
			else Edge4->setCell2(c);
			c->addEdge(Edge4);
		}
		
	}
	
	// only the edges know their connections - set the complete structure!
	SetTissueStructure();
	
	// the mechanical properties of the cells remain to be set (all cells, vertices, apical edges and basal edges are equal!!) by MechProp to have a completely defined system
	
	
	setMechanicalProperties();
	
	int size = ListEdge.size();
	size = ListVertex.size();
	update();
	
	return true;
	
	
}

bool Tissue::createHexagonalTissue(int rows, int columns, double height, double length, double width)
{
    isPeriodic = false;
    
 	// check the parameters:
	if(rows<1 || columns<1 || height<0 || length<0 || width<0)
		return false;
	
	// clear the tissue
	clear();
	
	//Cell::resetCounter();
	//Edge::resetCounter();
	//Vertex::resetCounter();
    
	// calculate the basic vectors for a single hexagon
	double SQRT3OVER2 = sqrt(3)/2;
	Point v1 = Point(width,0,0);
	Point v2 = Point(length * 0.5, length * SQRT3OVER2,0);
	Point v3 = Point( -length * 0.5, length * SQRT3OVER2,0);
	Point v4 = Point(-width,0,0);
	Point v5 = v2 * (-1);
	
	Point ApicalToBasal = Point(0,0,-height);
	
	// calculate the vectors to get from one cell to another cell
	Point IteratorInColumn = Point(0,2*length*SQRT3OVER2,0);
	Point IteratorInOddRow = v1 + v2; // if you are in an odd column you can get to the next cell in the same row by
	Point IteratorInEvenRow = v1 - v3;
	
	// set up the geometry
	Point startingPoint;
	Point baseStartingPoint = Point(0,0,0);
	
	Point V1, V2, V3, V4, V5, V6, V1b, V2b, V3b, V4b, V5b, V6b;
	Vertex *Ver1, *Ver2, *Ver3, *Ver4, *Ver5, *Ver6;
	Edge *Edge1, *Edge2, *Edge3, *Edge4, *Edge5, *Edge6;
	
	
	// Test
	//Vertex v = Vertex(Point(1,1,1),Point(1,1,0));
	
	//SimilarOrNewVertex(, )
	
	for(int c = 1;c<=columns;c++)
	{
		V1 = baseStartingPoint;
		
		for(int r=1;r<=rows;r++)
		{
			// get starting point for the next cell in the same column
			
			// create the cell and add it to the tissue
			Cell* c = new Cell();
			AddCell(c);
			c = &(ListCell.front());
			
			c->Type = 1;
			
			// the points of the cell
			V1b = V1 + ApicalToBasal;
			V2 = V1 + v1; V2b = V2 + ApicalToBasal;
			V3 = V2 + v2; V3b = V3 + ApicalToBasal;
			V4 = V3 + v3; V4b = V4 + ApicalToBasal;
			V5 = V4 + v4; V5b = V5 + ApicalToBasal;
			V6 = V5 + v5; V6b = V6 + ApicalToBasal;
			
			// create or find all the vertices that belong to the cell
			Ver1 = SimilarOrNewVertex(V1, V1b);
			Ver2 = SimilarOrNewVertex(V2, V2b);
			Ver3 = SimilarOrNewVertex(V3, V3b);
			Ver4 = SimilarOrNewVertex(V4, V4b);
			Ver5 = SimilarOrNewVertex(V5, V5b);
			Ver6 = SimilarOrNewVertex(V6, V6b);
			
			
			
			
			// create or find all the edges that belong to the cell and then set cell1 or cell2, depending on the direction
			bool rightDirection, newEdge;
			// edge 1
			Edge1 = SimilarOrNewEdge(Ver1, Ver2, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge1->setCell1(c);
			}
			else
			{
				Edge1->setCell2(c);
			}
			
			// edge 2
			Edge2 = SimilarOrNewEdge(Ver2, Ver3, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge2->setCell1(c);
			}
			else
			{
				Edge2->setCell2(c);
			}
			// edge 3
			Edge3 = SimilarOrNewEdge(Ver3, Ver4, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge3->setCell1(c);
			}
			else
			{
				Edge3->setCell2(c);
			}
			// edge 4
			Edge4 = SimilarOrNewEdge(Ver4, Ver5, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge4->setCell1(c);
			}
			else
			{
				Edge4->setCell2(c);
			}
			// edge 5
			Edge5 = SimilarOrNewEdge(Ver5, Ver6, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge5->setCell1(c);
			}
			else
			{
				Edge5->setCell2(c);
			}
			// edge 6
			Edge6 = SimilarOrNewEdge(Ver6, Ver1, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge6->setCell1(c);
			}
			else
			{
				Edge6->setCell2(c);
			}
			
			c->addEdge(Edge1);
			c->addEdge(Edge2);
			c->addEdge(Edge3);
			c->addEdge(Edge4);
			c->addEdge(Edge5);
			c->addEdge(Edge6);
			
			V1 = V1 + IteratorInColumn;
			
		}
		
		// get starting point for new column
		if(c % 2 == 1) //currently in odd row
		{
			baseStartingPoint = baseStartingPoint + IteratorInOddRow;
		}
		else // currently in even row
		{
			baseStartingPoint = baseStartingPoint + IteratorInEvenRow;
		}
		
	}
	
	// only the edges know their connections - set the complete structure!
	SetTissueStructure();
	setMechanicalProperties();	
    update();

	// the mechanical properties of the cells remain to be set (all cells, vertices, apical edges and basal edges are equal!!) by MechProp to have a completely defined system
	
    //std::list<Edge>::iterator it_edge = ListEdge.begin();
	//for(;it_edge!=ListEdge.end();it_edge++)
    //	{
    //		(*it_edge).setSurfaceTension(MechProp.Edge_SurfaceTension);
    //		(*it_edge).setLineTensionApical(MechProp.Edge_ApicalLineTension);
    //		(*it_edge).setLineTensionBasal(MechProp.Edge_BasalLineTension);
    //
    //		// if the edge is at the boundary fix the z-value
    //		if((*it_edge).c1()==NULL || (*it_edge).c2==NULL)
    //		{
    //			//(*it_edge).getVertex1()->zafixed=true;
    ////			(*it_edge).getVertex1()->zbfixed=true;
    ////			(*it_edge).getVertex2()->zafixed=true;
    ////			(*it_edge).getVertex2()->zbfixed=true;
    ////			(*it_edge).getVertex1()->xafixed=true;
    ////			(*it_edge).getVertex1()->xbfixed=true;
    ////			(*it_edge).getVertex2()->xafixed=true;
    ////			(*it_edge).getVertex2()->xbfixed=true;
    ////			(*it_edge).getVertex1()->yafixed=true;
    ////			(*it_edge).getVertex1()->ybfixed=true;
    ////			(*it_edge).getVertex2()->yafixed=true;
    ////			(*it_edge).getVertex2()->ybfixed=true;
    //		}
    //	}
	
	
//	int size = ListEdge.size();
//	size = ListVertex.size();
//
//	// give every vertex the information if at lies at the boundary
//	for(std::list<Vertex>::iterator it_vertex=ListVertex.begin();it_vertex!=ListVertex.end();it_vertex++)
//	{
//		(*it_vertex).isAtBoundary = (*it_vertex).get_isAtBoundary();
//	}
//	
	return true;
	
}

bool Tissue::createPeriodicHexagonalTissue(int rows, int columns, double height, double length, double width)
{
    // check the parameters:
	if(rows<1 || columns<1 || height<0 || length<0 || width<0 || columns%2==1)
		return false;
	
	// delete the existing tissue
	clear();
    
    // no basal boundary constraint
    apical_wall = false;
    basal_wall = false;
    STIFF_BASAL_MEMBRANE = false;
    onlyApicalContributions = false;
    
    // the tissue is periodic!
	isPeriodic = true;
        
	// calculate the basic vectors for a single hexagon
	double SQRT3OVER2 = sqrt(3)/2;
	Point v1 = Point(width,0,0);
	Point v2 = Point(length * 0.5, length * SQRT3OVER2,0);
	Point v3 = Point( -length * 0.5, length * SQRT3OVER2,0);
	Point v4 = v1 * (-1);
	Point v5 = v2 * (-1);
	
	Point ApicalToBasal = Point(0,0,-height);
	
    double rotationAngle = 0;
    
    v1.rotate(rotationAngle);
    v2.rotate(rotationAngle);
    v3.rotate(rotationAngle);
    v4.rotate(rotationAngle);
    v5.rotate(rotationAngle);

    //std::cout << "TEST = " << std::cos(2*PI);
    //v1.Print();
    //(v1+v2+v3+v4+v5-v3).Print();
    
	// calculate the vectors to get from one cell to another cell
	Point IteratorInColumn = Point(0,2*length*SQRT3OVER2,0);
	Point IteratorInOddRow = v1 + v2; // if you are in an odd column you can get to the next cell in the same row by
	Point IteratorInEvenRow = v1 - v3;
	
    IteratorInColumn.rotate(rotationAngle);
    IteratorInOddRow.rotate(rotationAngle);
    IteratorInEvenRow.rotate(rotationAngle);
    
	Point V1, V2, V3, V4, V5, V6, V1b, V2b, V3b, V4b, V5b, V6b;
	Vertex *Ver1, *Ver2, *Ver3, *Ver4, *Ver5, *Ver6;
	Edge *Edge1, *Edge2, *Edge3, *Edge4, *Edge5, *Edge6;
	
	// calculate the system size
	double Ly = sqrt(3)*length * rows;
	double Lx = (length/2.0 + width) * columns;
	
	// set up the geometry
    //Point startingPoint = Point(Lx/2, Ly/2,0);
	//Point baseStartingPoint = Point(Lx/2, Ly/2, -height);
	
    Point startingPoint = Point(60, 60,height);
	
	SystemSize.x = Lx;
	SystemSize.y = Ly;
	//SystemSize.rotate(rotationAngle);

    //Lx = SystemSize.x;
	//Ly = SystemSize.y;

    
	// Test
	//Vertex v = Vertex(Point(1,1,1),Point(1,1,0));
	
	//SimilarOrNewVertex(, )
	
	for(int c = 1;c<=columns;++c)
	{
		V1 = startingPoint;
		
		for(int r=1;r<=rows;++r)
		{
			// create the cell and add it to the tissue
			Cell* c = new Cell();
			AddCell(c);
			c = &(ListCell.front());
			
			c->Type = 1;
			
			// the points of the cell
			V2 = V1 + v1;
			V3 = V2 + v2;
			V4 = V3 + v3;
			V5 = V4 + v4;
			V6 = V5 + v5;
			
			// apply the periodicity and save the quadrants
			V1.x = V1.x/Lx; int V1_quadx = floor(V1.x); V1.x = (V1.x-V1_quadx)*Lx;
			V1.y = V1.y/Ly; int V1_quady = floor(V1.y); V1.y = (V1.y-V1_quady)*Ly;
			V2.x = V2.x/Lx; int V2_quadx = floor(V2.x); V2.x = (V2.x-V2_quadx)*Lx;
			V2.y = V2.y/Ly; int V2_quady = floor(V2.y); V2.y = (V2.y-V2_quady)*Ly;
			V3.x = V3.x/Lx; int V3_quadx = floor(V3.x); V3.x = (V3.x-V3_quadx)*Lx;
			V3.y = V3.y/Ly; int V3_quady = floor(V3.y); V3.y = (V3.y-V3_quady)*Ly;
			V4.x = V4.x/Lx; int V4_quadx = floor(V4.x); V4.x = (V4.x-V4_quadx)*Lx;
			V4.y = V4.y/Ly; int V4_quady = floor(V4.y); V4.y = (V4.y-V4_quady)*Ly;
			V5.x = V5.x/Lx; int V5_quadx = floor(V5.x); V5.x = (V5.x-V5_quadx)*Lx;
			V5.y = V5.y/Ly; int V5_quady = floor(V5.y); V5.y = (V5.y-V5_quady)*Ly;
			V6.x = V6.x/Lx; int V6_quadx = floor(V6.x); V6.x = (V6.x-V6_quadx)*Lx;
			V6.y = V6.y/Ly; int V6_quady = floor(V6.y); V6.y = (V6.y-V6_quady)*Ly;
			
			// create the basal points
			V1b = V1 + ApicalToBasal;
			V2b = V2 + ApicalToBasal;
			V3b = V3 + ApicalToBasal;
			V4b = V4 + ApicalToBasal;
			V5b = V5 + ApicalToBasal;
			V6b = V6 + ApicalToBasal;
			
			// create the quadrants of the edges
			Quadrant quad_e12 = Quadrant(V2_quadx-V1_quadx,V2_quady-V1_quady);
			Quadrant quad_e23 = Quadrant(V3_quadx-V2_quadx,V3_quady-V2_quady);
			Quadrant quad_e34 = Quadrant(V4_quadx-V3_quadx,V4_quady-V3_quady);
			Quadrant quad_e45 = Quadrant(V5_quadx-V4_quadx,V5_quady-V4_quady);
			Quadrant quad_e56 = Quadrant(V6_quadx-V5_quadx,V6_quady-V5_quady);
			Quadrant quad_e61 = Quadrant(V1_quadx-V6_quadx,V1_quady-V6_quady);
			
			// create or find all the vertices that belong to the cell
			Ver1 = SimilarOrNewVertex(V1, V1b);
			Ver2 = SimilarOrNewVertex(V2, V2b);
			Ver3 = SimilarOrNewVertex(V3, V3b);
			Ver4 = SimilarOrNewVertex(V4, V4b);
			Ver5 = SimilarOrNewVertex(V5, V5b);
			Ver6 = SimilarOrNewVertex(V6, V6b);
			
			// create or find all the edges that belong to the cell and then set cell1 or cell2, depending on the direction
			bool rightDirection, newEdge;
			// edge 1
			Edge1 = SimilarOrNewEdge(Ver1, Ver2, &rightDirection, &newEdge, quad_e12);
			if(rightDirection)
			{
				Edge1->setCell1(c);
			}
			else
			{
				Edge1->setCell2(c);
			}
			
			// edge 2
			Edge2 = SimilarOrNewEdge(Ver2, Ver3, &rightDirection, &newEdge, quad_e23);
			if(rightDirection)
			{
				Edge2->setCell1(c);
			}
			else
			{
				Edge2->setCell2(c);
			}
			// edge 3
			Edge3 = SimilarOrNewEdge(Ver3, Ver4, &rightDirection, &newEdge, quad_e34);
			if(rightDirection)
			{
				Edge3->setCell1(c);
			}
			else
			{
				Edge3->setCell2(c);
			}
			// edge 4
			Edge4 = SimilarOrNewEdge(Ver4, Ver5, &rightDirection, &newEdge, quad_e45);
			if(rightDirection)
			{
				Edge4->setCell1(c);
			}
			else
			{
				Edge4->setCell2(c);
			}
			// edge 5
			Edge5 = SimilarOrNewEdge(Ver5, Ver6, &rightDirection, &newEdge, quad_e56);
			if(rightDirection)
			{
				Edge5->setCell1(c);
			}
			else
			{
				Edge5->setCell2(c);
			}
			// edge 6
			Edge6 = SimilarOrNewEdge(Ver6, Ver1, &rightDirection, &newEdge, quad_e61);
			if(rightDirection)
			{
				Edge6->setCell1(c);
			}
			else
			{
				Edge6->setCell2(c);
			}
			
            c->addEdge(Edge6);
			c->addEdge(Edge5);
			c->addEdge(Edge4);
			c->addEdge(Edge3);
			c->addEdge(Edge2);
			c->addEdge(Edge1);
			
			V1 = V1 + IteratorInColumn;
			
		}
		
		// get starting point for new column
		if(c % 2 == 1) //currently in odd row
		{
			startingPoint = startingPoint + IteratorInOddRow;
		}
		else // currently in even row
		{
			startingPoint = startingPoint + IteratorInEvenRow;
		}
		
	}
	
	
	// only the edges know their connections - set the complete structure!
	SetTissueStructure();
	
	setMechanicalProperties();
    
    update();

	return true;
	
}

void Tissue::createRoundHexagonalTissue(int radiusInCells, double height, double length, double MaximalRadius)
{
    // set periodicity
    isPeriodic = false;
    
	// check the basic parameters:
	if(radiusInCells<1 || height<0 || length<0)
		return;
	
    // clear the tissue
	clear();
	
	//Cell::resetCounter();
    //	Edge::resetCounter();
    //	Vertex::resetCounter();
    //	// calculate the basic vectors for a single hexagon
	double SQRT3OVER2 = sqrt(3)/2;
	Point v1 = Point(length,0,0);
	Point v2 = Point(length * 0.5, length * SQRT3OVER2,0);
	Point v3 = Point( -length * 0.5, length * SQRT3OVER2,0);
	Point v4 = Point(-length,0,0);
	Point v5 = v2 * (-1);
	
	// calculate the vectors that allow the navigation through the tissue
	Point U1 = Point(3.0/2*length,-SQRT3OVER2*length,0);
	Point U2 = Point(0,-SQRT3OVER2*length*2,0);
	Point U3 = Point(-3.0/2*length,-SQRT3OVER2*length,0);
	
	// vector from the middle point to the first point of the hexagon
	Point R = Point(-length/2,-SQRT3OVER2*length,0);
	
	Point ApicalToBasal = Point(0,0,-height);
	
	
	// set up the geometry
	Point startingPoint;
	
	Point V1, V2, V3, V4, V5, V6, V1b, V2b, V3b, V4b, V5b, V6b, M;
	Vertex *Ver1, *Ver2, *Ver3, *Ver4, *Ver5, *Ver6;
	Edge *Edge1, *Edge2, *Edge3, *Edge4, *Edge5, *Edge6;
	
	// define the cylindrical boundary
	double radius = (radiusInCells+0.33)*2*SQRT3OVER2*length;
	boundary_isSet=true;
	boundary_apicalRadius = radius;
	boundary_basalRadius = radius;
	boundary_isCylinder=false; // apical and basal restrictions are independent!
	boundary_externaldW_dva = 0;
	boundary_externaldW_dvb = 0;
	boundary_canDetach=false;
	
	V1 = Point(0,0,0);
	
	
	for(int circle = 0;circle<=radiusInCells;circle++)
	{
		// middle point of the first cell of the circle - all aligned line above (0,0,0)
		M = U2*(-circle);
        
		
		// now iterate through the middle points and for every point create a hexagon around
		for(int i = 0;i<=circle;i++)
		{
			
			// construct first point of the cell (bottom left) by first going to middle point and from there to the vertex
			if(i!=0) M = M+U1;
			
			if(MaximalRadius>0 && Point::Norm(M)>MaximalRadius) continue;
			
			// create the cell and add it to the tissue
			Cell* c = new Cell();
			AddCell(c);
            
			c = &(ListCell.front());
			
			V1 = M + R;
			
			// the points of the cell
			V1b = V1 + ApicalToBasal;
			V2 = V1 + v1; V2b = V2 + ApicalToBasal;
			V3 = V2 + v2; V3b = V3 + ApicalToBasal;
			V4 = V3 + v3; V4b = V4 + ApicalToBasal;
			V5 = V4 + v4; V5b = V5 + ApicalToBasal;
			V6 = V5 + v5; V6b = V6 + ApicalToBasal;
			
			// create or find all the vertices that belong to the cell
			Ver1 = SimilarOrNewVertex(V1, V1b);
			Ver2 = SimilarOrNewVertex(V2, V2b);
			Ver3 = SimilarOrNewVertex(V3, V3b);
			Ver4 = SimilarOrNewVertex(V4, V4b);
			Ver5 = SimilarOrNewVertex(V5, V5b);
			Ver6 = SimilarOrNewVertex(V6, V6b);
			
			
			// create or find all the edges that belong to the cell and then set cell1 or cell2, depending on the direction
			bool rightDirection, newEdge;
			// edge 1
			Edge1 = SimilarOrNewEdge(Ver1, Ver2, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge1->setCell1(c);
			}
			else
			{
				Edge1->setCell2(c);
			}
			
			// edge 2
			Edge2 = SimilarOrNewEdge(Ver2, Ver3, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge2->setCell1(c);
			}
			else
			{
				Edge2->setCell2(c);
			}
			// edge 3
			Edge3 = SimilarOrNewEdge(Ver3, Ver4, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge3->setCell1(c);
			}
			else
			{
				Edge3->setCell2(c);
			}
			// edge 4
			Edge4 = SimilarOrNewEdge(Ver4, Ver5, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge4->setCell1(c);
			}
			else
			{
				Edge4->setCell2(c);
			}
			// edge 5
			Edge5 = SimilarOrNewEdge(Ver5, Ver6, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge5->setCell1(c);
			}
			else
			{
				Edge5->setCell2(c);
			}
			// edge 6
			Edge6 = SimilarOrNewEdge(Ver6, Ver1, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge6->setCell1(c);
			}
			else
			{
				Edge6->setCell2(c);
			}
			
			c->addEdge(Edge1);
			c->addEdge(Edge2);
			c->addEdge(Edge3);
			c->addEdge(Edge4);
			c->addEdge(Edge5);
			c->addEdge(Edge6);
            
		}
        
		
		// now iterate through the middle points and for every point create a hexagon around
		for(int i = 1;i<=circle;i++)
		{
			M = M + U2;
			if(!MaximalRadius<=0 && Point::Norm(M)>MaximalRadius) continue;
            
			// create the cell and add it to the tissue
			Cell* c = new Cell();
			AddCell(c);
			c = &(ListCell.front());
			
			// construct first point of the cell (bottom left) by first going to middle point and from there to the vertex
			V1 = M + R;
			
			// the points of the cell
			V1b = V1 + ApicalToBasal;
			V2 = V1 + v1; V2b = V2 + ApicalToBasal;
			V3 = V2 + v2; V3b = V3 + ApicalToBasal;
			V4 = V3 + v3; V4b = V4 + ApicalToBasal;
			V5 = V4 + v4; V5b = V5 + ApicalToBasal;
			V6 = V5 + v5; V6b = V6 + ApicalToBasal;
			
			// create or find all the vertices that belong to the cell
			Ver1 = SimilarOrNewVertex(V1, V1b);
			Ver2 = SimilarOrNewVertex(V2, V2b);
			Ver3 = SimilarOrNewVertex(V3, V3b);
			Ver4 = SimilarOrNewVertex(V4, V4b);
			Ver5 = SimilarOrNewVertex(V5, V5b);
			Ver6 = SimilarOrNewVertex(V6, V6b);
			
			// create or find all the edges that belong to the cell and then set cell1 or cell2, depending on the direction
			bool rightDirection, newEdge;
			// edge 1
			Edge1 = SimilarOrNewEdge(Ver1, Ver2, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge1->setCell1(c);
			}
			else
			{
				Edge1->setCell2(c);
			}
			
			// edge 2
			Edge2 = SimilarOrNewEdge(Ver2, Ver3, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge2->setCell1(c);
			}
			else
			{
				Edge2->setCell2(c);
			}
			// edge 3
			Edge3 = SimilarOrNewEdge(Ver3, Ver4, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge3->setCell1(c);
			}
			else
			{
				Edge3->setCell2(c);
			}
			// edge 4
			Edge4 = SimilarOrNewEdge(Ver4, Ver5, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge4->setCell1(c);
			}
			else
			{
				Edge4->setCell2(c);
			}
			// edge 5
			Edge5 = SimilarOrNewEdge(Ver5, Ver6, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge5->setCell1(c);
			}
			else
			{
				Edge5->setCell2(c);
			}
			// edge 6
			Edge6 = SimilarOrNewEdge(Ver6, Ver1, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge6->setCell1(c);
			}
			else
			{
				Edge6->setCell2(c);
			}
			
			c->addEdge(Edge1);
			c->addEdge(Edge2);
			c->addEdge(Edge3);
			c->addEdge(Edge4);
			c->addEdge(Edge5);
			c->addEdge(Edge6);
			
		}
        
		// now iterate through the middle points and for every point create a hexagon around
		for(int i = 1;i<=circle;i++)
		{
			M = M + U3;
			if(!MaximalRadius<=0 && Point::Norm(M)>MaximalRadius) continue;
			
			// create the cell and add it to the tissue
			Cell* c = new Cell();
			AddCell(c);
			c = &(ListCell.front());
			
			// construct first point of the cell (bottom left) by first going to middle point and from there to the vertex
			V1 = M + R;
			
			// the points of the cell
			V1b = V1 + ApicalToBasal;
			V2 = V1 + v1; V2b = V2 + ApicalToBasal;
			V3 = V2 + v2; V3b = V3 + ApicalToBasal;
			V4 = V3 + v3; V4b = V4 + ApicalToBasal;
			V5 = V4 + v4; V5b = V5 + ApicalToBasal;
			V6 = V5 + v5; V6b = V6 + ApicalToBasal;
			
			// create or find all the vertices that belong to the cell
			Ver1 = SimilarOrNewVertex(V1, V1b);
			Ver2 = SimilarOrNewVertex(V2, V2b);
			Ver3 = SimilarOrNewVertex(V3, V3b);
			Ver4 = SimilarOrNewVertex(V4, V4b);
			Ver5 = SimilarOrNewVertex(V5, V5b);
			Ver6 = SimilarOrNewVertex(V6, V6b);
			
			// create or find all the edges that belong to the cell and then set cell1 or cell2, depending on the direction
			bool rightDirection, newEdge;
			// edge 1
			Edge1 = SimilarOrNewEdge(Ver1, Ver2, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge1->setCell1(c);
			}
			else
			{
				Edge1->setCell2(c);
			}
			
			// edge 2
			Edge2 = SimilarOrNewEdge(Ver2, Ver3, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge2->setCell1(c);
			}
			else
			{
				Edge2->setCell2(c);
			}
			// edge 3
			Edge3 = SimilarOrNewEdge(Ver3, Ver4, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge3->setCell1(c);
			}
			else
			{
				Edge3->setCell2(c);
			}
			// edge 4
			Edge4 = SimilarOrNewEdge(Ver4, Ver5, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge4->setCell1(c);
			}
			else
			{
				Edge4->setCell2(c);
			}
			// edge 5
			Edge5 = SimilarOrNewEdge(Ver5, Ver6, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge5->setCell1(c);
			}
			else
			{
				Edge5->setCell2(c);
			}
			// edge 6
			Edge6 = SimilarOrNewEdge(Ver6, Ver1, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge6->setCell1(c);
			}
			else
			{
				Edge6->setCell2(c);
			}
			
			c->addEdge(Edge1);
			c->addEdge(Edge2);
			c->addEdge(Edge3);
			c->addEdge(Edge4);
			c->addEdge(Edge5);
			c->addEdge(Edge6);
			
		}
		
		// now iterate through the middle points and for every point create a hexagon around
		for(int i = 1;i<=circle;i++)
		{
			M = M - U1;
			if(!MaximalRadius<=0 && Point::Norm(M)>MaximalRadius) continue;
			
			// create the cell and add it to the tissue
			Cell* c = new Cell();
			AddCell(c);
			c = &(ListCell.front());
			
			// construct first point of the cell (bottom left) by first going to middle point and from there to the vertex
			V1 = M + R;
			
			// the points of the cell
			V1b = V1 + ApicalToBasal;
			V2 = V1 + v1; V2b = V2 + ApicalToBasal;
			V3 = V2 + v2; V3b = V3 + ApicalToBasal;
			V4 = V3 + v3; V4b = V4 + ApicalToBasal;
			V5 = V4 + v4; V5b = V5 + ApicalToBasal;
			V6 = V5 + v5; V6b = V6 + ApicalToBasal;
			
			// create or find all the vertices that belong to the cell
			Ver1 = SimilarOrNewVertex(V1, V1b);
			Ver2 = SimilarOrNewVertex(V2, V2b);
			Ver3 = SimilarOrNewVertex(V3, V3b);
			Ver4 = SimilarOrNewVertex(V4, V4b);
			Ver5 = SimilarOrNewVertex(V5, V5b);
			Ver6 = SimilarOrNewVertex(V6, V6b);
			
			// create or find all the edges that belong to the cell and then set cell1 or cell2, depending on the direction
			bool rightDirection, newEdge;
			// edge 1
			Edge1 = SimilarOrNewEdge(Ver1, Ver2, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge1->setCell1(c);
			}
			else
			{
				Edge1->setCell2(c);
			}
			
			// edge 2
			Edge2 = SimilarOrNewEdge(Ver2, Ver3, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge2->setCell1(c);
			}
			else
			{
				Edge2->setCell2(c);
			}
			// edge 3
			Edge3 = SimilarOrNewEdge(Ver3, Ver4, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge3->setCell1(c);
			}
			else
			{
				Edge3->setCell2(c);
			}
			// edge 4
			Edge4 = SimilarOrNewEdge(Ver4, Ver5, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge4->setCell1(c);
			}
			else
			{
				Edge4->setCell2(c);
			}
			// edge 5
			Edge5 = SimilarOrNewEdge(Ver5, Ver6, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge5->setCell1(c);
			}
			else
			{
				Edge5->setCell2(c);
			}
			// edge 6
			Edge6 = SimilarOrNewEdge(Ver6, Ver1, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge6->setCell1(c);
			}
			else
			{
				Edge6->setCell2(c);
			}
			
			c->addEdge(Edge1);
			c->addEdge(Edge2);
			c->addEdge(Edge3);
			c->addEdge(Edge4);
			c->addEdge(Edge5);
			c->addEdge(Edge6);
			
		}
		
		// now iterate through the middle points and for every point create a hexagon around
		for(int i = 1;i<=circle;i++)
		{
			M = M - U2;
			if(!MaximalRadius<=0 && Point::Norm(M)>MaximalRadius) continue;
			
			// create the cell and add it to the tissue
			Cell* c = new Cell();
			AddCell(c);
			c = &(ListCell.front());
			
			// construct first point of the cell (bottom left) by first going to middle point and from there to the vertex
			V1 = M + R;
			
			// the points of the cell
			V1b = V1 + ApicalToBasal;
			V2 = V1 + v1; V2b = V2 + ApicalToBasal;
			V3 = V2 + v2; V3b = V3 + ApicalToBasal;
			V4 = V3 + v3; V4b = V4 + ApicalToBasal;
			V5 = V4 + v4; V5b = V5 + ApicalToBasal;
			V6 = V5 + v5; V6b = V6 + ApicalToBasal;
			
			// create or find all the vertices that belong to the cell
			Ver1 = SimilarOrNewVertex(V1, V1b);
			Ver2 = SimilarOrNewVertex(V2, V2b);
			Ver3 = SimilarOrNewVertex(V3, V3b);
			Ver4 = SimilarOrNewVertex(V4, V4b);
			Ver5 = SimilarOrNewVertex(V5, V5b);
			Ver6 = SimilarOrNewVertex(V6, V6b);
			
			// create or find all the edges that belong to the cell and then set cell1 or cell2, depending on the direction
			bool rightDirection, newEdge;
			// edge 1
			Edge1 = SimilarOrNewEdge(Ver1, Ver2, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge1->setCell1(c);
			}
			else
			{
				Edge1->setCell2(c);
			}
			
			// edge 2
			Edge2 = SimilarOrNewEdge(Ver2, Ver3, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge2->setCell1(c);
			}
			else
			{
				Edge2->setCell2(c);
			}
			// edge 3
			Edge3 = SimilarOrNewEdge(Ver3, Ver4, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge3->setCell1(c);
			}
			else
			{
				Edge3->setCell2(c);
			}
			// edge 4
			Edge4 = SimilarOrNewEdge(Ver4, Ver5, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge4->setCell1(c);
			}
			else
			{
				Edge4->setCell2(c);
			}
			// edge 5
			Edge5 = SimilarOrNewEdge(Ver5, Ver6, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge5->setCell1(c);
			}
			else
			{
				Edge5->setCell2(c);
			}
			// edge 6
			Edge6 = SimilarOrNewEdge(Ver6, Ver1, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge6->setCell1(c);
			}
			else
			{
				Edge6->setCell2(c);
			}
			
			c->addEdge(Edge1);
			c->addEdge(Edge2);
			c->addEdge(Edge3);
			c->addEdge(Edge4);
			c->addEdge(Edge5);
			c->addEdge(Edge6);
			
		}
		
		// now iterate through the middle points and for every point create a hexagon around
		for(int i = 1;i<circle;i++)
		{
			M = M - U3;
			if(!MaximalRadius<=0 && Point::Norm(M)>MaximalRadius) continue;
			
			// create the cell and add it to the tissue
			Cell* c = new Cell();
			AddCell(c);
			c = &(ListCell.front());
			
			// construct first point of the cell (bottom left) by first going to middle point and from there to the vertex
			V1 = M + R;
			
			// the points of the cell
			V1b = V1 + ApicalToBasal;
			V2 = V1 + v1; V2b = V2 + ApicalToBasal;
			V3 = V2 + v2; V3b = V3 + ApicalToBasal;
			V4 = V3 + v3; V4b = V4 + ApicalToBasal;
			V5 = V4 + v4; V5b = V5 + ApicalToBasal;
			V6 = V5 + v5; V6b = V6 + ApicalToBasal;
			
			// create or find all the vertices that belong to the cell
			Ver1 = SimilarOrNewVertex(V1, V1b);
			Ver2 = SimilarOrNewVertex(V2, V2b);
			Ver3 = SimilarOrNewVertex(V3, V3b);
			Ver4 = SimilarOrNewVertex(V4, V4b);
			Ver5 = SimilarOrNewVertex(V5, V5b);
			Ver6 = SimilarOrNewVertex(V6, V6b);
			
			// create or find all the edges that belong to the cell and then set cell1 or cell2, depending on the direction
			bool rightDirection, newEdge;
			// edge 1
			Edge1 = SimilarOrNewEdge(Ver1, Ver2, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge1->setCell1(c);
			}
			else
			{
				Edge1->setCell2(c);
			}
			
			// edge 2
			Edge2 = SimilarOrNewEdge(Ver2, Ver3, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge2->setCell1(c);
			}
			else
			{
				Edge2->setCell2(c);
			}
			// edge 3
			Edge3 = SimilarOrNewEdge(Ver3, Ver4, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge3->setCell1(c);
			}
			else
			{
				Edge3->setCell2(c);
			}
			// edge 4
			Edge4 = SimilarOrNewEdge(Ver4, Ver5, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge4->setCell1(c);
			}
			else
			{
				Edge4->setCell2(c);
			}
			// edge 5
			Edge5 = SimilarOrNewEdge(Ver5, Ver6, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge5->setCell1(c);
			}
			else
			{
				Edge5->setCell2(c);
			}
			// edge 6
			Edge6 = SimilarOrNewEdge(Ver6, Ver1, &rightDirection, &newEdge);
			if(rightDirection)
			{
				Edge6->setCell1(c);
			}
			else
			{
				Edge6->setCell2(c);
			}
            
			c->addEdge(Edge1);
			c->addEdge(Edge2);
			c->addEdge(Edge3);
			c->addEdge(Edge4);
			c->addEdge(Edge5);
			c->addEdge(Edge6);
			
		}
	}
	
	// only the edges know their connections - set the complete structure!
	SetTissueStructure();
	
	// the mechanical properties of the cells remain to be set (all cells, vertices, apical edges and basal edges are equal!!) by MechProp to have a completely defined system
	
	for(std::list<Cell>::iterator itc=ListCell.begin();itc!=ListCell.end();itc++)
		(*itc).Type=1;
	
	int size = ListEdge.size();
	size = ListVertex.size();
	update();
	
	// give every vertex the information if at lies at the boundary
	for(std::list<Vertex>::iterator it_vertex=ListVertex.begin();it_vertex!=ListVertex.end();it_vertex++)
	{
		Vertex* v=&(*it_vertex);
        
		v->isAtBoundary = v->get_isAtBoundary();
        
		
		// project all the points that are supposed to lie on the exterior of the cylinder back on the cylinder along a radius
		if(v->isAtBoundary)
		{
			// how much lies the vertex outside the cylinder??
			double excess = Point::Norm(v->getCoord())-boundary_apicalRadius;
			double scaling = boundary_apicalRadius/(boundary_apicalRadius+excess);
			v->setCoord(Point(v->getCoord().x*scaling,v->getCoord().y*scaling,v->getCoord().z));
			//excess = Point::Norm(v->getBasalCoord())-boundary_basalRadius;
			//scaling = boundary_basalRadius/(boundary_basalRadius+excess);
			v->setBasalCoord(Point(v->getBasalCoord().x*scaling,v->getBasalCoord().y*scaling,v->getBasalCoord().z));
		}
	}
	
	setMechanicalProperties();
}

void Tissue::updateData()
{
	// remove all the information from tissue data
	tissueData = TissueData();
    
    // cell properties
	for(std::list<Cell>::iterator it_cell = ListCell.begin();it_cell != ListCell.end();it_cell++) {
        if(it_cell->Type>=0 && it_cell->Type<=NO_CELL_TYPES){
            tissueData.BC_a_z[it_cell->Type].push_back(it_cell->BC_a.z);
            tissueData.BC_b_z[it_cell->Type].push_back(it_cell->BC_b.z);
            tissueData.height[it_cell->Type].push_back(it_cell->height());
            tissueData.apicalArea[it_cell->Type].push_back(it_cell->CalculateApicalSurfaceArea());
            tissueData.basalArea[it_cell->Type].push_back(it_cell->CalculateBasalSurfaceArea());
            tissueData.volume[it_cell->Type].push_back(it_cell->Volume);
            tissueData.angle[it_cell->Type].push_back(it_cell->angleWithZPlane());
        }
    }
    
    for(unsigned i=0;i<=NO_CELL_TYPES;i++) {

        double var;
        
        vector_statistics(tissueData.height[i], tissueData.mean_height[i], var, tissueData.max_height[i], tissueData.min_height[i]);
        vector_statistics(tissueData.BC_a_z[i], tissueData.mean_BC_a_z[i], var, tissueData.max_BC_a_z[i], tissueData.min_BC_a_z[i]);
        vector_statistics(tissueData.BC_b_z[i], tissueData.mean_BC_b_z[i], var, tissueData.max_BC_b_z[i], tissueData.min_BC_b_z[i]);
        vector_statistics(tissueData.apicalArea[i], tissueData.mean_apicalArea[i], tissueData.var_apicalArea[i], tissueData.max_apicalArea[i], tissueData.min_apicalArea[i]);
        vector_statistics(tissueData.basalArea[i], tissueData.mean_basalArea[i], var, tissueData.max_basalArea[i], tissueData.min_basalArea[i]);
        vector_statistics(tissueData.angle[i], tissueData.mean_angle[i], var, tissueData.max_angle[i], tissueData.min_angle[i]);
        vector_statistics(tissueData.volume[i], tissueData.mean_volume[i], var, tissueData.max_volume[i], tissueData.min_volume[i]);
        
        tissueData.numberOfCells[i] = tissueData.height[i].size();
    }
	
    // initialize interfaces / edge properties
    for(unsigned i=0; i<=NO_CELL_TYPES; i++) {
        for(unsigned j=0; j<=NO_CELL_TYPES; j++) {
            tissueData.apicalCircumference[i][j]=0;
            tissueData.basalCircumference[i][j]=0;
            tissueData.lateralInterface[i][j]=0;
            tissueData.numberOfCommonEdges[i][j]=0;
        }
    }
    
    for(std::list<Edge>::iterator it_e = ListEdge.begin(); it_e!=ListEdge.end(); it_e++)
    {
        unsigned b = std::max(it_e->c1->Type,it_e->c2->Type);
        unsigned s = std::min(it_e->c1->Type,it_e->c2->Type);
        
        tissueData.apicalCircumference[b][s] += it_e->l_a();
        tissueData.basalCircumference[b][s] += it_e->l_b();
        tissueData.lateralInterface[b][s] += it_e->area();
        tissueData.numberOfCommonEdges[b][s]++;
        tissueData.averageApicalZheight[b][s]+= (it_e->v1->coord.z + it_e->v2->coord.z)/2;
        tissueData.averageBasalZheight[b][s]+= (it_e->v1->coord_b.z + it_e->v2->coord_b.z)/2;
    }
    
    // calculate the averages
    for(unsigned i = 0;i<=NO_CELL_TYPES; i++) {
        
        for(unsigned j = 0;j<=i; j++) {
            
            if(tissueData.numberOfCommonEdges[i][j]==0) continue; // avoid division by zero...
            
            tissueData.averageApicalZheight[i][j] = tissueData.averageApicalZheight[i][j]/tissueData.numberOfCommonEdges[i][j];
            tissueData.averageBasalZheight[i][j] = tissueData.averageBasalZheight[i][j]/tissueData.numberOfCommonEdges[i][j];
        }
    }
    
    // calculate the indentations of the cyst
    for(unsigned i=0; i<=NO_CELL_TYPES; i++) {
        // distinguish 3 cases:
        //
        // 1. ring < min < max (apical ring z-position is smaller than both min and maximal z-positions of the cells in the clone: mush room/bulge out)
        // 2. min < ring < max (cyst goes down but comes up again
        
        tissueData.apicalIndentation_min[i] =  tissueData.min_BC_a_z[i] -  tissueData.max_BC_a_z[1];
        tissueData.apicalIndentation_max[i] =  tissueData.max_BC_a_z[1] -  tissueData.min_BC_a_z[i];

        tissueData.basalIndentation_min[i] =  tissueData.min_BC_b_z[i] - tissueData.min_BC_b_z[1];
        tissueData.basalIndentation_max[i] =  tissueData.max_BC_b_z[i] - tissueData.min_BC_b_z[1];
    }
    
    // calculate the indentations of the cyst
    double ring = tissueData.averageApicalZheight[2][1];
    double maximumZposition = tissueData.max_BC_a_z[2];
    double minimumZposition = tissueData.min_BC_a_z[2];
    
    if(ring>maximumZposition) tissueData.apicalIndentation = minimumZposition - ring;
    else tissueData.apicalIndentation = maximumZposition - ring;
    

    // for now the midcell is known to be number 143
    //midCell_number = 1343;
    //Cell* midCell = CellFromNumber(midCell_number);
    
    // cell in the middle of the clone
    Cell* midCell = cell_in_middle(2);//CellFromNumber(midCell_number);
    
    if(!midCell) midCell = &(*ListCell.begin());
    
    tissueData.BC_a_midCell = midCell->BC_a;
    tissueData.BC_b_midCell = midCell->BC_b;
    
    
    // calculate midpoint of apical sides of all cells
    tissueData.midpoint = Point();
    for(std::list<Cell>::iterator itc = ListCell.begin(); itc!=ListCell.end(); itc++)
    {
        tissueData.midpoint += itc->BC_a;
    }
    tissueData.midpoint = tissueData.midpoint/ListCell.size();
    
    tissueData.Radius = 0;
    for(std::list<Cell>::iterator itc = ListCell.begin(); itc!=ListCell.end(); itc++)
    {
        tissueData.Radius += Point::Norm(itc->BC_a-tissueData.midpoint);
    }
    
    tissueData.Radius = tissueData.Radius/ListCell.size();
    
    // maximum distance from tissue midpoint of all cell types
    for(unsigned i=0; i<=NO_CELL_TYPES; i++) {
        tissueData.max_R[i] = 0;
        tissueData.min_R[i] = 9999999999999;
    }
    
    for(std::list<Cell>::iterator itc = ListCell.begin(); itc!=ListCell.end(); itc++)
    {
        double R = Point::Norm(itc->BC_a-tissueData.midpoint);
//        double R = Point::Norm(itc->BC_a-Point(0,0,0));

        if(tissueData.max_R[itc->Type]<R) tissueData.max_R[itc->Type]=R ;
        if(tissueData.min_R[itc->Type]>R) tissueData.min_R[itc->Type]=R ;
        
    }
}


void Tissue::writeTissueData(const char *fileName)
{
	updateData();
	
	//Open File;
//	std::ofstream output;
//	output.open(fileName);// initial position is set to be at the end of the file
//	
//	std::list< std::vector<double> >::iterator it = tissueData.GeometricalData.begin();
//	std::list<Cell>::iterator it_c = ListCell.begin();
//	
//	output << "cell number/height/apical surface area/basal surface area/lateral surface area/volume/pressure/apical perimeter/basal perimeter" << std::endl;
//	
//	for(;it!=tissueData.GeometricalData.end();it++)
//	{
//		output << (*it_c).Number << " / " << (*it)[0] << " / " << (*it)[1] << " / " << (*it)[2] << " / " << (*it)[3] << " / " << (*it)[4] << " / " << (*it)[5] << " / " << (*it)[6] << " / " << (*it)[7] << std::endl;
//		it_c++;
//	}
//	
//	
	
	
}

void Tissue::extractDataFromDir(std::string dirName)
{
    // of stream
    std::ofstream output;
    
    // open the file named fileName where the extracted data shall be written to
    std::string outputName = dirName+"/results.3dvm_ext";
    const char * c = outputName.c_str();

    // us the Qt file system in order to work with the directory and extract all files from it
    QDir dir = QDir(QString::fromStdString(dirName));
    
    // list of all files in the selected directory
    QStringList listOfFiles = dir.entryList(QDir::Files);
    
    bool atLeastOneReadableFile = false;
    
    for(int i = 0; i<listOfFiles.size(); i++) { // iterate through the files
        
        // open the file
        QString fileName = listOfFiles.at(i);
        
        //std::cout << fileName.remove(".txt").append(".pdf").toStdString() << " in folder " << dirName.toStdString() << std::endl;
        
        if(fileName.endsWith(".3dvm")) {
            int flag = loadFileFromTxt(QString::fromStdString(dirName)+"/"+fileName);
            
            // get all the geometrical data of the current tissue
            updateData();
            
            
            
            if(flag==1) {
                
                enum DataMode {Gastrulation, StabilityCheck, CystFormation};
                
                //DataMode dataMode = Gastrulation;
                DataMode dataMode = Gastrulation;
                
                if(dataMode == StabilityCheck)
                {
                    if(!atLeastOneReadableFile) {
                        output.open(c, std::ofstream::trunc);
                        // write directory name to the top of the file
                        output << "# " << dirName << std::endl;
                        output << "fileName | systemSize | Ts | rho | isFlat" << std::endl;
                        atLeastOneReadableFile = !atLeastOneReadableFile;
                    }
                    
                    double rho = 1./tissueData.mean_apicalArea[1]; // average cell density
                    
                    // if the file has been opened properly
                    output << fileName.toStdString()
                    << " | "  << tissueData.numberOfCells[1]            // 32
                    << " | "  << MechProp.CellPropVector[1].T_a + MechProp.CellPropVector[1].T_b            // 32
                    << " | "  << rho
                    << " | "  << isFlat()
                    << std::endl;
                }
                else if(dataMode == Gastrulation) {
                    if(!atLeastOneReadableFile) {
                        output.open(c, std::ofstream::trunc);
                        // write directory name to the top of the file
                        output << "# " << dirName << std::endl;
                        output << "fileName | cellNumber | V_Lumen_a | V0_Lumen_a | K_Lumen_a | K_cell1 | A0_cell1 | Lambda1 | rho1 | bulkModulus1 | shearModulus1 | kappa1 | K_cell2 | A0_cell2 | Lambda2  | rho2 | bulkModulus2 | shearModulus2 | kappa2 | R | R_test | zeta1 | zeta2 | zeta_laplace | deltaZeta21_over_zeta1 |  max_R_1 | min_R_1 | max_R_2 | min_R_2" << std::endl;
//                        | midpoint_x | midpoint_y | midpoint_z" << std::endl;
//                        mean_apical_area_1 | mean_apical_area_2 |
                        atLeastOneReadableFile = !atLeastOneReadableFile;
                    }
                    
                    // calculate the effective 2D bulk and shear modulus of tissue with type 1
                    double Lambda1 = MechProp.Intercell_ApicalLineTension[1][1];
                    double Lambda2 = MechProp.Intercell_ApicalLineTension[2][2];
                    
                    // calculate effective tensions for tissue of type 1 and tissue of type 2
                    double K_cell1 = (MechProp.CellPropVector[1].T_a-MechProp.CellPropVector[1].T0_a)/MechProp.CellPropVector[1].A0_a;
                    double A0_cell1 = MechProp.CellPropVector[1].A0_a*MechProp.CellPropVector[1].T0_a/(MechProp.CellPropVector[1].T0_a-MechProp.CellPropVector[1].T_a);
                    double rho1 = 1./tissueData.mean_apicalArea[1]; // average cell density
                    double bulkModulus1 = K_cell1/rho1 - std::pow(3,1./4.)/std::pow(2,3./2.)*Lambda1*std::sqrt(rho1);
                    double shearModulus1 = std::pow(3,1./4.)/std::pow(2,3./2.)*Lambda1*std::sqrt(rho1);
                    double zeta1 = K_cell1*(1./rho1-A0_cell1) + std::pow(3.,0.25)/std::pow(2.,0.5)*Lambda1*std::sqrt(rho1);
                    
                    double K_cell2 = (MechProp.CellPropVector[2].T_a-MechProp.CellPropVector[2].T0_a)/MechProp.CellPropVector[2].A0_a;
                    double A0_cell2 = MechProp.CellPropVector[2].A0_a*MechProp.CellPropVector[2].T0_a/(MechProp.CellPropVector[2].T0_a-MechProp.CellPropVector[2].T_a);
                    double rho2 = 1./tissueData.mean_apicalArea[2]; // average cell density
                    double bulkModulus2 = K_cell2/rho2 - std::pow(3,1./4.)/std::pow(2,3./2.)*Lambda2*std::sqrt(rho2);
                    double shearModulus2 = std::pow(3,1./4.)/std::pow(2,3./2.)*Lambda2*std::sqrt(rho2);
                    double zeta2 = K_cell2*(1./rho1-A0_cell2) + std::pow(3.,0.25)/std::pow(2.,0.5)*Lambda2*std::sqrt(rho1);
                    
                    update_P_Lumen_a();
                    double R = std::pow(3./(4.*PI)*V_Lumen_a,1./3.);
                    double zeta_laplace = 0.5*R*P_Lumen_a;
                    
                    double deltaZeta21_over_zeta1 = (zeta2-zeta1)/zeta1;                    
                    
                    // if the file has been opened properly
                    output << fileName.toStdString()
                    << " | "  << tissueData.numberOfCells[1]+tissueData.numberOfCells[2]            // 32
                    << " | "  << V_Lumen_a
                    << " | "  << V0_Lumen_a
                    << " | "  << K_Lumen_a
                    << " | "  << K_cell1
                    << " | "  << A0_cell1
                    << " | "  << Lambda1
                    << " | "  << rho1
                    << " | "  << bulkModulus1
                    << " | "  << shearModulus1
                    << " | "  << bulkModulus1/shearModulus1 // kappa1
                    << " | "  << K_cell2
                    << " | "  << A0_cell2
                    << " | "  << Lambda2
                    << " | "  << rho2
                    << " | "  << bulkModulus2
                    << " | "  << shearModulus2
                    << " | "  << bulkModulus2/shearModulus2 // kappa2
                    << " | "  << R
                    << " | "  << tissueData.Radius
                    << " | "  << zeta1
                    << " | "  << zeta2
                    << " | "  << zeta_laplace
                    << " | "  << deltaZeta21_over_zeta1
//                  << " | "  << tissueData.mean_apicalArea[1]
//                  << " | "  << tissueData.mean_apicalArea[2]
                    << " | "  << tissueData.max_R[1]
                    << " | "  << tissueData.min_R[1]
                    << " | "  << tissueData.max_R[2]
                    << " | "  << tissueData.min_R[2]
//                  << " | "  << tissueData.midpoint.x
//                  << " | "  << tissueData.midpoint.y
//                  << " | "  << tissueData.midpoint.z
                    << std::endl;
                }
                else if(dataMode == CystFormation) {
                    
                    if(!atLeastOneReadableFile) {
                        output.open(c, std::ofstream::trunc);
                        // write directory name to the top of the file
                        output << "# " << dirName << std::endl;
                        output << "fileName | cloneCellNumber | perimeter_a | perimeter_b | z_ring_a | z_ring_b | height_wt_max | height_mut_max | height_mut_min | V_Lumen_a | mean_apical_area_1 | mean_apical_area_2 | var_apical_area_1 | var_apical_area_2 | V0 | volume_mut_mean | z_wt_max_a | z_wt_min_a | z_wt_max_b | z_wt_min_b | z_mut_max_a | z_mut_min_a | z_mut_max_b | z_mut_min_b | BC_a_midCell_z | BC_b_midCell_z | z_wt_mean_a | z_wt_mean_b | z_mid_a | z_mid_b | z_interface_a | z_interface_b | width_a | width_b | max_R_1 | min_R_1 | max_R_2 | min_R_2 | midpoint_x | midpoint_y | midpoint_y | radius" << std::endl;
                        atLeastOneReadableFile = !atLeastOneReadableFile;
                    }

                    // if the file has been opened properly
                    output << fileName.toStdString()
                    << " | "  << tissueData.numberOfCells[2]            // 32
                    << " | "  << tissueData.apicalCircumference[2][1]   // 27
                    << " | "  << tissueData.basalCircumference[2][1]    // 28
                    << " | "  << tissueData.averageApicalZheight[2][1]
                    << " | "  << tissueData.averageBasalZheight[2][1]
                    << " | "  << tissueData.max_height[1]               // 16
                    << " | "  << tissueData.max_height[2]               // 17
                    << " | "  << tissueData.min_height[2]               // 17
                    << " | "  << V_Lumen_a
                    << " | "  << tissueData.mean_apicalArea[1]
                    << " | "  << tissueData.mean_apicalArea[2]
                    << " | "  << tissueData.var_apicalArea[1]
                    << " | "  << tissueData.var_apicalArea[2]
                    << " | "  << MechProp.CellPropVector[1].V0
                    << " | "  << tissueData.mean_volume[2]               // 25
                    << " | "  << tissueData.max_BC_a_z[1]               // 6
                    << " | "  << tissueData.min_BC_a_z[1]               // 7
                    << " | "  << tissueData.max_BC_b_z[1]               // 10
                    << " | "  << tissueData.min_BC_b_z[1]               // 11
                    << " | "  << tissueData.max_BC_a_z[2]               // 8
                    << " | "  << tissueData.min_BC_a_z[2]               // 9
                    << " | "  << tissueData.max_BC_b_z[2]               // 12
                    << " | "  << tissueData.min_BC_b_z[2]               // 13
                    << " | "  << tissueData.BC_a_midCell.z
                    << " | "  << tissueData.BC_b_midCell.z
                    << " | "  << tissueData.mean_BC_a_z[1]
                    << " | "  << tissueData.mean_BC_b_z[1]
                    << " | "  << tissueData.BC_a_midCell.z
                    << " | "  << tissueData.BC_b_midCell.z
                    << " | "  << tissueData.averageApicalZheight[2][1]
                    << " | "  << tissueData.averageBasalZheight[2][1]
                    << " | "  << tissueData.apicalCircumference[2][1]/PI
                    << " | "  << tissueData.basalCircumference[2][1]/PI
                    << " | "  << tissueData.max_R[1]
                    << " | "  << tissueData.min_R[1]
                    << " | "  << tissueData.max_R[2]
                    << " | "  << tissueData.min_R[2]
                    << " | "  << tissueData.midpoint.x
                    << " | "  << tissueData.midpoint.y
                    << " | "  << tissueData.midpoint.z
                    << " | "  << tissueData.Radius
                    << std::endl;

                    
                    
                } else  {
                    
                    if(!atLeastOneReadableFile) {
                        output.open(c, std::ofstream::trunc);
                        // write directory name to the top of the file
                        output << "# " << dirName << std::endl;
                        output << "fileName | max_apicalIndentation | min_apicalInde ntation | max_basalIndentation | min_basalIndentation | max_BC_a_z1 | min_BC_a_z1 | max_BC_a_z2 | min_BC_a_z2 | max_BC_b_z1 | min_BC_b_z1 | max_BC_b_z2 | min_BC_b_z2 | mean_height1 | mean_height2 | max_height1 | max_height2 | min_height1 | min_height2 | mean_apicalArea1 | mean_apicalArea2 | mean_basalArea1 | mean_basalArea2 | min_volume1 | min_volume2 | lateralSurface12 | apicalCirc12 | basalCirc12 | max_angle1 | max_angle2 | #cells1 | #cells2" << std::endl;
                        atLeastOneReadableFile = !atLeastOneReadableFile;
                    }

                    // if the file has been opened properly
                    output << fileName.toStdString()
                        << " | "  << tissueData.apicalIndentation_max[2]    // 2
                        << " | "  << tissueData.apicalIndentation_min[2]    // 3
                        << " | "  << tissueData.basalIndentation_max[2]     // 4
                        << " | "  << tissueData.basalIndentation_min[2]     // 5
                        << " | "  << tissueData.max_BC_a_z[1]               // 6
                        << " | "  << tissueData.min_BC_a_z[1]               // 7
                        << " | "  << tissueData.max_BC_a_z[2]               // 8
                        << " | "  << tissueData.min_BC_a_z[2]               // 9
                        << " | "  << tissueData.max_BC_b_z[1]               // 10
                        << " | "  << tissueData.min_BC_b_z[1]               // 11
                        << " | "  << tissueData.max_BC_b_z[2]               // 12
                        << " | "  << tissueData.min_BC_b_z[2]               // 13
                        << " | "  << tissueData.mean_height[1]              // 14
                        << " | "  << tissueData.mean_height[2]              // 15
                        << " | "  << tissueData.max_height[1]               // 16
                        << " | "  << tissueData.max_height[2]               // 17
                        << " | "  << tissueData.min_height[1]               // 18
                        << " | "  << tissueData.min_height[2]               // 19
                        << " | "  << tissueData.mean_apicalArea[1]          // 20
                        << " | "  << tissueData.mean_apicalArea[2]          // 21
                        << " | "  << tissueData.mean_basalArea[1]           // 22
                        << " | "  << tissueData.mean_basalArea[2]           // 23
                        << " | "  << tissueData.min_volume[1]               // 24
                        << " | "  << tissueData.min_volume[2]               // 25
                        << " | "  << tissueData.lateralInterface[2][1]      // 26
                        << " | "  << tissueData.apicalCircumference[2][1]   // 27
                        << " | "  << tissueData.basalCircumference[2][1]    // 28
                        << " | "  << tissueData.max_angle[1]                // 29
                        << " | "  << tissueData.max_angle[2]                // 30
                        << " | "  << tissueData.numberOfCells[1]            // 31
                        << " | "  << tissueData.numberOfCells[2]            // 32
                        << std::endl;
                }
            }
        }
        
    }

    if(output.is_open()) output.close();
}

// after every step the tissue needs to be tested, if any vertex left the periodic boundaries and the periodicity information in the edges have to be adjusted accordingly
void Tissue::reestablishPeriodicity()
{
	std::list<Vertex>::iterator it=ListVertex.begin();
	
	for(;it!=ListVertex.end();it++)
	{
		Vertex *v = &(*it);
		
		if(v->coord.x >= SystemSize.x)
		{
			// bring the coordinate back into [0,1)
			v->coord.x = v->coord.x-SystemSize.x;
			// update the information about the basal quadrant
			v->q_b.x += 1;
			
			// update the adjacent edges
			for(std::list<Edge*>::iterator ite = v->NeighbouringEdges_sorted.begin();ite!=v->NeighbouringEdges_sorted.end();ite++)
			{
				if(v == (*ite)->getVertex1())
					(*ite)->q_v2.x -= 1;
				else
					(*ite)->q_v2.x +=1;
			}
			
		}
		else if(v->coord.x < 0)
		{
			// bring the coordinate back into [0,1)
			v->coord.x = v->coord.x+SystemSize.x;
			// update the information about the basal quadrant
			v->q_b.x -= 1;
			
			// update the adjacent edges
			for(std::list<Edge*>::iterator ite = v->NeighbouringEdges_sorted.begin();ite!=v->NeighbouringEdges_sorted.end();ite++)
			{
				if(v == (*ite)->getVertex1())
					(*ite)->q_v2.x += 1;
				else
					(*ite)->q_v2.x -=1;
			}
			
		}
		
		
		if(v->coord.y >= SystemSize.y)
		{
			// bring the coordinate back into [0,1)
			v->coord.y = v->coord.y-SystemSize.y;
			// update the information about the basal quadrant
			v->q_b.y += 1;
			
			// update the adjacent edges
			for(std::list<Edge*>::iterator ite = v->NeighbouringEdges_sorted.begin();ite!=v->NeighbouringEdges_sorted.end();ite++)
			{
				if(v == (*ite)->getVertex1())
					(*ite)->q_v2.y -= 1;
				else
					(*ite)->q_v2.y +=1;
			}
			
		}
		else if(v->coord.y < 0)
		{
			// bring the coordinate back into [0,1)
			v->coord.y = v->coord.y+SystemSize.y;
			// update the information about the basal quadrant
			v->q_b.y -= 1;
			
			// update the adjacent edges
			for(std::list<Edge*>::iterator ite = v->NeighbouringEdges_sorted.begin();ite!=v->NeighbouringEdges_sorted.end();ite++)
			{
				if(v == (*ite)->getVertex1())
					(*ite)->q_v2.y += 1;
				else
					(*ite)->q_v2.y -=1;
			}
			
		}
		
	}
	
}

void Tissue::ChangeLxAndLy(double newLx, double newLy)
{
	// calculate the ratios - all positions have to be rescaled using those
	double ratio_Lx = newLx/SystemSize.x;
	double ratio_Ly = newLy/SystemSize.y;
	
	// replace the old system sizes
	SystemSize.x = newLx;
	SystemSize.y = newLy;
	
	// change all the coordinates
	for(std::list<Vertex>::iterator it=ListVertex.begin();it!=ListVertex.end();it++)
	{
		(*it).coord.x *= ratio_Lx;
		(*it).coord.y *= ratio_Ly;
		(*it).coord_b.x *= ratio_Lx;
		(*it).coord_b.y *= ratio_Ly;
	}
}

std::list<Cell*> Tissue::extractCellsInPlane(Point positionVector_plane, Point normalVector_plane)
{
    // list of all cells that are cut through by the current plane
    std::list<Cell*> cellsInPlane;
    
    // obtain the normal vector from the two angles alpha and beta
    //Point normalVector_plane = Point(sin(sectionPlane_beta)*cos(sectionPlane_alpha),sin(sectionPlane_beta)*sin(sectionPlane_alpha),cos(sectionPlane_beta));
    
    // allocate local variables (vertex and barycenter positions)
    Point v1a, v1b, v2a, v2b, b_e, b_a, b_b;
    
    // allocate pointers for easier handling
    Edge *e;
    Vertex *v1, *v2;
    
    Point P1, P2, dummy;
    bool liesInPlane, crossesAB;
    
    // iterate over all cells of the tissue
    for(std::list<Cell>::iterator it_cell = ListCell.begin();it_cell != ListCell.end();++it_cell)
    {
        // name cell c
        Cell* c = &(*it_cell);
        
        // iterate over all the edges of the cells
        for(std::list<Edge*>::iterator it_edge = it_cell->ListEdges.begin();it_edge != it_cell->ListEdges.end();it_edge++)
        {
            // name edge and vertices
            e = *it_edge;
            v1 = e->getVertex1();
            v2 = e->getVertex2();
            
            // apical barycenter is only fixed point of the cell
            b_a = c->BC_a;
            
            // calculate the positions of the vertices such that they are all in the same quadrant as the barycenter of the cell
            b_b = shift( c->BC_b , SystemSize*(c->q_BC_b));
            
            if(c == e->c1)
            {
                v1a = shift( v1->getCoord() ,	   SystemSize*(e->q_c1)*(-1) ) ;
                v1b = shift( v1->getBasalCoord() , SystemSize*(e->q_c1 - v1->q_b)*(-1)  ) ;
                v2a = shift( v2->getCoord() ,	   SystemSize*(e->q_c1 - e->q_v2)*(-1) );
                v2b = shift( v2->getBasalCoord() , SystemSize*(e->q_c1 - e->q_v2 - v2->q_b)*(-1) );
                b_e = shift( e->BC_l, SystemSize*(e->q_c1 - e->q_BC_l)*(-1) );
            }
            else
            {
                v2a = shift( v2->getCoord() ,	   SystemSize*(e->q_c2)*(-1) ) ;
                v2b = shift( v2->getBasalCoord() , SystemSize*(e->q_c2 - v2->q_b)*(-1)  ) ;
                v1a = shift( v1->getCoord() ,	   SystemSize*(e->q_c2 + e->q_v2)*(-1) );
                v1b = shift( v1->getBasalCoord() , SystemSize*(e->q_c2 + e->q_v2 - v1->q_b)*(-1) );
                b_e = shift( e->BC_l, SystemSize*(e->q_c2 + e->q_v2 - e->q_BC_l)*(-1) );
            }
            
            // check if the plane cuts through the cell at any triangle
            // triangle towards v1
            if(Point::TrianglePlaneSection(v1a,v2a,b_e, normalVector_plane, positionVector_plane,  P1,  P2,  liesInPlane,  crossesAB)) {cellsInPlane.push_back(c);break;}
            else if(Point::TrianglePlaneSection( v1a, v2a,b_a, normalVector_plane, positionVector_plane,  P1,  P2,  liesInPlane,  crossesAB)) {cellsInPlane.push_back(c);break;}
        }
        
    }
    
    return cellsInPlane;
}


//void Tissue::updateXprofile(double xCross)
//{
//	// empty the profiles
//	xProfileApical.clear();
//	xProfileBasal.clear();
//
//	std::list<Edge>::iterator it_edge = ListEdge.begin();
//
//	Point ApicalCrossing, SecondApicalCrossing, BasalCrossing, SecondBasalCrossing;
//	bool doesApicalCrossing, doesBasalCrossing, isApicallyAligned, isBasallyAligned;
//
//	for(;it_edge!=ListEdge.end();it_edge++)
//	{
//		// check for the section of the edge with the xCross line
//		(*it_edge).getXcrossing(xCross, &ApicalCrossing, &doesApicalCrossing, &BasalCrossing, &doesBasalCrossing, &SecondApicalCrossing, &isApicallyAligned, &SecondBasalCrossing, &isBasallyAligned );
//		// if both points of the apical side of the  edge are situated in the plane
//		if(isApicallyAligned)
//		{
//			xProfileApical.insert(std::pair<double,Point>(ApicalCrossing.y,ApicalCrossing));
//			xProfileApical.insert(std::pair<double,Point>(SecondApicalCrossing.y,SecondApicalCrossing));
//		}
//		// if the edge crosses the plane
//		else if(doesApicalCrossing)
//		{
//			xProfileApical.insert(std::pair<double,Point>(ApicalCrossing.y,ApicalCrossing));
//		}
//
//		// if both points of the basal side of the  edge are situated in the plane
//		if(isBasallyAligned)
//		{
//			xProfileBasal.insert(std::pair<double,Point>(BasalCrossing.y,BasalCrossing));
//			xProfileBasal.insert(std::pair<double,Point>(SecondBasalCrossing.y,SecondBasalCrossing));
//		}
//		// if the edge crosses the plane
//		else if(doesBasalCrossing)
//		{
//			xProfileBasal.insert(std::pair<double,Point>(BasalCrossing.y,BasalCrossing));
//		}
//
//	}
//}
//
//void Tissue::updateYprofile(double yCross)
//{
//	// empty the profiles
//	yProfileApical.clear();
//	yProfileBasal.clear();
//
//	std::list<Edge>::iterator it_edge = ListEdge.begin();
//
//	Point ApicalCrossing, SecondApicalCrossing, BasalCrossing, SecondBasalCrossing;
//	bool doesApicalCrossing, doesBasalCrossing, isApicallyAligned, isBasallyAligned;
//
//	for(;it_edge!=ListEdge.end();it_edge++)
//	{
//		// check for the section of the edge with the xCross line
//		(*it_edge).getYcrossing(yCross, &ApicalCrossing, &doesApicalCrossing, &BasalCrossing, &doesBasalCrossing, &SecondApicalCrossing, &isApicallyAligned, &SecondBasalCrossing, &isBasallyAligned );
//		// if both points of the apical side of the  edge are situated in the plane
//		if(isApicallyAligned)
//		{
//			yProfileApical.insert(std::pair<double,Point>(ApicalCrossing.x,ApicalCrossing));
//			yProfileApical.insert(std::pair<double,Point>(SecondApicalCrossing.x,SecondApicalCrossing));
//		}
//		// if the edge crosses the plane
//		else if(doesApicalCrossing)
//		{
//			yProfileApical.insert(std::pair<double,Point>(ApicalCrossing.x,ApicalCrossing));
//		}
//
//		// if both points of the basal side of the  edge are situated in the plane
//		if(isBasallyAligned)
//		{
//			yProfileBasal.insert(std::pair<double,Point>(BasalCrossing.x,BasalCrossing));
//			yProfileBasal.insert(std::pair<double,Point>(SecondBasalCrossing.x,SecondBasalCrossing));
//		}
//		// if the edge crosses the plane
//		else if(doesBasalCrossing)
//		{
//			yProfileBasal.insert(std::pair<double,Point>(BasalCrossing.x,BasalCrossing));
//		}
//
//	}
//
//
//}
//
//void Tissue::updateXprofile2(double xCross)
//{
//	xIntersection.clear();
//
//	std::list<Cell>::iterator it_cell = ListCell.begin();
//
//	for(;it_cell!=ListCell.end();it_cell++)
//	{
//		CellIntersection cellXcross=(*it_cell).getXintersection(xCross);
//		if(cellXcross.size()!=0)
//			xIntersection.push_back(cellXcross);
//	}
//
//}

//int Tissue::UpdatePlaneSection(Point currentPosition)
//{
//
//	SectionPlane_PositionVector.x = currentPosition.x;
//	SectionPlane_PositionVector.y = currentPosition.y;
//
//	// obtain the normal vector from the two angles alpha and beta
//	SectionPlane_NormalVector = Point(sin(SectionPlane_beta)*cos(SectionPlane_alpha),sin(SectionPlane_beta)*sin(SectionPlane_alpha),cos(SectionPlane_beta));
//
//	Point e1, e2;
//
//	if(SectionPlane_NormalVector.x==0 && SectionPlane_NormalVector.y==0)
//	{
//		e1 = Point(1,0,0);
//		e2 = Point(0,1,0);
//	}
//
//	else
//	{
//		e1 = SectionPlane_NormalVector * Point(0,0,1);
//		e1 = e1 / Point::Norm(e1);
//
//		e2 = SectionPlane_NormalVector * e1;
//	}
//
//	// update sections with the edges
//	std::list<Edge>::iterator it_edge = ListEdge.begin();
//	for(;it_edge!=ListEdge.end();it_edge++)
//	{
//		(*it_edge).updatePlaneSection();
//		// if there is an apical section transform the coordinates
//		if((*it_edge).planeCrosses_apical)
//		{
//			Point dummy = (*it_edge).planeSection_apical;
//			(*it_edge).planeSection_apical.x = dummy.dot(e1);
//			(*it_edge).planeSection_apical.y = dummy.dot(e2);
//		}
//
//		// if there is an basal section transform the coordinates
//		if((*it_edge).planeCrosses_basal)
//		{
//			Point dummy = (*it_edge).planeSection_basal;
//			(*it_edge).planeSection_basal.x = dummy.dot(e1);
//			(*it_edge).planeSection_basal.y = dummy.dot(e2);
//		}
//	}
//
//	// update the sections with the lateral edges inbetween the vertices
//	std::list<Vertex>::iterator it_vertex = ListVertex.begin();
//	for(;it_vertex!=ListVertex.end();it_vertex++)
//	{
//		(*it_vertex).updatePlaneSection();
//
//		// if there is a section transform the coordinates
//		if((*it_vertex).planeCrosses)
//		{
//			Point dummy = (*it_vertex).planeSection;
//			(*it_vertex).planeSection.x = dummy.dot(e1);
//			(*it_vertex).planeSection.y = dummy.dot(e2);
//		}
//	}
//
//	return 1;
//
//}

//bool Tissue::UnpopVertices()
//{
//	std::list<Vertex>::iterator itv = ListVertex.begin();
//	
//	bool unpopped = false;
//	
//	Edge e;
//	
//	for(;itv!=ListVertex.end();itv++)
//	{
//		if((*itv).checkAndPop(0.01, &e) == 1) unpopped = true;
//	}
//    
//	return unpopped;
//}

//void Tissue::QPushT1Transition()
//{
//
//	T1Transition((*CurrentEdge.begin()));
//}
//
//void Tissue::QPushApoptosis()
//{
//	if(!CurrentCell.empty())
//	{
//		(*CurrentCell.begin())->Apoptosis();
//		RemoveCell((*CurrentCell.begin()));
//		CurrentCell.pop_front();
//	}
//	update();
//}
//
//void Tissue::T1Transition(Edge* e)
//{
//	SetTissueStructure();
//	// how many edges does every vertex have?
//	std::list<Vertex>::iterator it = ListVertex.begin();
//	for(;it!=ListVertex.end();it++)
//	{
//		std::cout << (*it).getListEdge()->size() << "  ";
//	}
//	e->RemoveFromTissue();
//	update();
//	repaint();
//	QApplication::processEvents();
//}

//void Tissue::updateProfile()
//{
//
//}

//
//bool PopEdge(Vertex* v, Cell* c1, Cell* c2)
//{
//	// cells c1 and c2 will become neighbours with common edge e, which was 'created' from vertex v
//
//	// check if the vertex is four fold
//	if(v->ListVertexEdge.size()==4)
//	{
//
//
//	}
//
//}

//Edge* Tissue::Pop(Vertex* v, Edge* e1, Edge* e3, Cell* c1, Cell* c2, Cell* c3, Cell* c4, double epsilon)
//{
//	//epsilon is the displacement imposed when forces are tested
//	Point e1Director;
//	e1Director=e1->MidPoint(v);
//	e1Director=e1Director*(1.0/e1Director.Norme());
//	Point e3Director;
//	e3Director=e3->MidPoint(v);
//	e3Director=e3Director*(1.0/e3Director.Norme());
//
//		std::list<Edge*>::iterator it2;
//
//		Point ForceDirection1=e1->ForceDirection(v);//direction of the first vertex
//		Point ForceDirection2=e3->ForceDirection(v);//direction of the second vertex
//		Point SommeForces=ForceDirection1+ForceDirection2;
//		double NormeForce=SommeForces.Norme();
//		Point deplacement=SommeForces*(epsilon/NormeForce/2.0);
//		Vertex newvtemp=Vertex(v->coord+deplacement);//creates the new vertex
//		AddVertex(&newvtemp);
//		Vertex * newv=&(ListVertex.front());
//		v->coord=v->coord-deplacement;//moves the old vertex
//
//		//connect new vertex to the relevant edges
//		(newv->ConnectedEdges).push_front(e1);
//		(newv->ConnectedEdges).push_front(e3);
//		if(*(e1->v1)==*v) {e1->v1=newv;}
//		if(*(e1->v2)==*v) {e1->v2=newv;}
//		if(*(e3->v1)==*v) {e3->v1=newv;}
//		if(*(e3->v2)==*v) {e3->v2=newv;}
//
//		for(it2=v->ConnectedEdges.begin();it2!=v->ConnectedEdges.end();)
//		{
//			if(**it2==*e1||**it2==*e3) it2=v->ConnectedEdges.erase(it2); else ++it2;
//		}
//
//	//create new edge
//	Edge newetemp=Edge(newv,v);
//	Point newEdgeDirector=newetemp.v1v2();
//	newEdgeDirector=newEdgeDirector*(1.0/newEdgeDirector.Norme());
//	double pv=(e3Director+e1Director*(-1)).ProduitVectoriel(newEdgeDirector+e3Director*(-1));
//
//			if(pv>0)
//				{
//				newetemp.c1=c2;
//				newetemp.c2=c4;
//				}
//			else
//				{
//				newetemp.c1=c4;
//				newetemp.c2=c2;
//				}
//		newetemp.Tension=Tensions[newetemp.c1->Type][newetemp.c2->Type];
//		AddEdge(&newetemp);
//		Edge* newe=&(ListEdge.front());
//		if(c2->Number!=0) c2->ListEdges.push_front(newe);
//		if(c4->Number!=0) c4->ListEdges.push_front(newe);
//		//actualiseTensions();
//		ActualiseCrossingVertex(v);
//		ActualiseCrossingVertex(newv);
//
//
//
//
//		return(newe);
//}



//void Tissue::UnPop(Edge* e)
//{
//
//	if((!e->v1->xfixed&&!e->v2->xfixed)&&(!e->v1->yfixed&&!e->v2->yfixed))
//		{
//	std::list<Edge>::iterator finde;
//	finde = find( ListEdge.begin(), ListEdge.end(), *e);//pour s'assurer qu'on prend bien celui qui est dans la liste
//
//	//move kept vertex (v2) toward the middle of the old edge
//	Point v1v2=finde->v1v2();
//	finde->v2->coord=finde->v2->coord-v1v2*(0.5);
//	if(xPeriodic||yPeriodic)
//	{
//		Point v1coordsouvenir=finde->v1->coord;
//		finde->v1->coord=finde->v1->coord+(v1v2*(0.5));
//		if(finde->yCrossBoundary12&
//		   finde->v2->coord.y>0
//		   &finde->v2->coord.y<10e-4
//		   &finde->v1->coord.y<Ly
//		   &finde->v1->coord.y>Ly-10e-4)
//		{
//			finde->v1->coord.y+=10e-3;
//		}
//
//		if(finde->yCrossBoundary12&
//		   finde->v1->coord.y>0
//		   &finde->v1->coord.y<10e-4
//		   &finde->v2->coord.y<Ly
//		   &finde->v2->coord.y>Ly-10e-4)
//		{
//			finde->v2->coord.y+=10e-3;
//		}
//		if(finde->xCrossBoundary12&
//		   finde->v2->coord.x>0
//		   &finde->v2->coord.x<10e-4
//		   &finde->v1->coord.x<Lx
//		   &finde->v1->coord.x>Lx-10e-4)
//		{
//			finde->v1->coord.x+=10e-3;
//		}
//		if(finde->xCrossBoundary21&
//		   finde->v1->coord.x>0
//		   &finde->v1->coord.x<10e-4
//		   &finde->v2->coord.x<Lx
//		   &finde->v2->coord.x>Lx-10e-4)
//		{
//			finde->v2->coord.x+=10e-3;
//		}
//
//		ActualiseCrossingVertex(finde->v1);
//		ActualiseCrossingVertex(finde->v2);
//		finde->v1->coord=v1coordsouvenir;//necessary to avoid confusion between the two vertices
//	}
//			//connect edges touching vertex v1 to vertex v2
//
//	std::list<Edge*>::iterator it2;
//	for(it2=finde->v1->ConnectedEdges.begin();it2!=finde->v1->ConnectedEdges.end();it2++)
//	{
//		if(!(**it2==*finde))
//		{
//			if(*((**it2).v1)==*(finde->v1)) (*it2)->v1=&(*(finde->v2)); else (*it2)->v2=&(*(finde->v2));//fais le remplacement du vertex v2 par v1 dans les edges touchant v2 (il faut d'abord voir si le v2 en cours est le v1 ou le v2 de l'edge avoisinant)
//			//il faut egalement completer les edges accroches a v2
//			finde->v2->ConnectedEdges.push_front(&(**it2));
//		}
//	}
//			//suppress old edge and vertex v1
//	Vertex* v1souvenir=finde->v1;
//	Vertex v2=*(finde->v2);
//	RemoveEdge(&(*finde));
//	RemoveVertex(v1souvenir);
//	}
//
//}
//


// interesting function, has to be analysed more in detail
//void Tissue::RemoveShortEdges(double minlength)
//{
//	std::list<Edge>::iterator it;
//	std::list<Edge>::iterator itFound;
//	for(it=ListEdge.begin();it!=ListEdge.end();)
//	{
//		if(it->Length()<minlength)//if the edge is too short it is erased
//		{
//			QLOG_INFO() << "Remove edge between "<<it->c1->Number<<" and "<<it->c2->Number<<" with length "<<it->Length();
//			itSouvenir=it;
//			it++;//il faut iterer avant de le supprimer
//			if(((itSouvenir->c1->A>20||itSouvenir->c1->Number==0)
//				&&(itSouvenir->c2->A>20||itSouvenir->c2->Number==0))//evite la suppression des cotes des petites cellules
//			   &&((itSouvenir->c1->NbEdges()>3||itSouvenir->c1->Number==0)
//				  &&(itSouvenir->c2->NbEdges()>3||itSouvenir->c2->Number==0))//evite qu'on se retrouve avec une cellule a deux cotes
//			   )
//			{
//				UnPop(&(*itSouvenir));
//				//UnPopAndPop(&(*itSouvenir));
//			}
//
//
//		}
//		else
//		{
//			it++;
//		}
//	}
//}

//void Tissue::updateCells()
//{
//	std::list<Cell>::iterator it;
//	for(it=ListCell.begin();it!=ListCell.end();it++)
//	{
//		//it->Print();
//		if(it->AreEdgesSorted==false) it->sortEdges();
//		it->UpdateCentroid();
//		it->A=it->Area();
//	}
//}
//
//void Tissue::Evolve(double mobility)
//{
//
//	std::list<Vertex>::iterator it;//it iterates on the list of vertices
//	//at every evolve iteration, reevaluate maximum residual
//	maxResidual=0;
//	for(it=ListVertex.begin();it!=ListVertex.end();it++)
//	{
//		it->dx=0;
//		it->dy=0;
//		std::list<Edge*>::iterator it2;//it2 iterates on the list of pointers pointing on the edges connected to a vertex
//		double AverageTension=0;
//		double residual;
//		int nbEdges=0;
//		for(it2=it->ConnectedEdges.begin();it2!=it->ConnectedEdges.end();it2++)
//		{
//			//average tension useful to evaluate normalized residuals
//			AverageTension+=(**it2).Tension;
//			nbEdges++;
//
//			//compute forces on every vertex
//			Point ForceDirection=(**it2).ForceDirection(&(*it));
//			Point ElasticEnergyForceDirection=(**it2).ElasticEnergyForceDirection(&(*it));
//			double deltaP=((*it2)->c2->ElasticEnergy()-(*it2)->c1->ElasticEnergy())/2.0;
//			if(!(it->xfixed))
//			{
//				it->dx+=mobility*((**it2).Tension*ForceDirection.x+deltaP*ElasticEnergyForceDirection.x);
//			}
//			if(!(it->yfixed))
//			{
//				it->dy+=mobility*((**it2).Tension*ForceDirection.y+deltaP*ElasticEnergyForceDirection.y);
//			}
//
//		}
//
//		//actualiseTensions();//necessary when the edges tensions depend on the position of vertices
//
//		//evaluate maximum residual
//		AverageTension/=((double)nbEdges);
//		residual=sqrt(pow(it->dx/mobility,2)+pow(it->dy/mobility,2.0))/AverageTension;
//		maxResidual=std::max(maxResidual,residual);
//	}
//
//
//	//update position of vertices
//	for(it=ListVertex.begin();it!=ListVertex.end();it++)
//	{
//		it->coord.x+=it->dx;
//		it->coord.y+=it->dy;
//
//		if(xPeriodic||yPeriodic)
//		{
//			//if periodic boundary condition, make sure that vertices are correctly placed
//			ActualiseCrossingVertex(&(*it));
//		}
//	}
//
//	//update cells areas and centroids
//	updateCells();
//	//display value of residual
//	QLOG_INFO() <<"max residual="<<maxResidual;
//
//
//}
//
//
//
//
////function which computes the forces without updating position of vertices
////required for evolution with adaptative step
//double* Tissue::EvolveNoIteration(double mobility, double CurrentStep)
//{
//
//	std::list<Vertex>::iterator it;//it iterates on the list of vertices
//
//	maxResidual=0;
//	int nbVertex=ListVertex.size();
//	int nbEdges=ListEdge.size();
//	double* ListIterations=new double[2*nbVertex];
//	int k=0;
//	for(it=ListVertex.begin();it!=ListVertex.end();it++)
//	{
//		it->dx=0;
//		it->dy=0;
//		std::list<Edge*>::iterator it2;//it2 itere sur la liste de pointeurs listant les edges d'un vertex
//		double AverageTension=0;
//		int nbEdges=0;
//		for(it2=it->ConnectedEdges.begin();it2!=it->ConnectedEdges.end();it2++)
//		{
//			//average tension useful to evaluate normalized residuals
//			AverageTension+=(**it2).Tension;
//			nbEdges++;
//
//			//compute forces on every vertex
//			Point ForceDirection=(**it2).ForceDirection(&(*it));
//			Point ElasticEnergyForceDirection=(**it2).ElasticEnergyForceDirection(&(*it));
//			double deltaP=((*it2)->c2->ElasticEnergy()-(*it2)->c1->ElasticEnergy())/2.0;
//			if(!(it->xfixed))
//			{
//				it->dx+=mobility*CurrentStep*((**it2).Tension*ForceDirection.x+deltaP*ElasticEnergyForceDirection.x);
//			}
//			if(!(it->yfixed))
//			{
//				it->dy+=mobility*CurrentStep*((**it2).Tension*ForceDirection.y+deltaP*ElasticEnergyForceDirection.y);
//			}
//
//
//		}
//		//store positions of vertices
//		ListIterations[k]=it->dx;
//		ListIterations[k+1]=it->dy;
//		k=k+2;
//		//actualiseTensions();//necessary when the edges tensions depend on the position of vertices
//	}
//
//
//	//no position of updates so no need to update cells areas and centroids
//	//updateCells();
//	return(ListIterations);
//}
//
//
//void Tissue::EvolveAdaptativeStep(double mobility, double deltat)
//{
//	//total time spent
//	double Deltat=0;
//	int nbIterations=0;
//	double erreur=0;
//	int nbVertex;
//	int nbEdges;
//
//	//first condition wait for total time
//	//second condition limits the number of steps, is for safety if something goes wrong
//	while(Deltat<deltat&nbIterations<10000)
//	{
//
//		nbVertex=ListVertex.size();
//		//extract position with two different steps
//		//1- full step
//		double* dPosition1=EvolveNoIteration(mobility,CurrentStep);
//
//		//2- half step
//		//compute first half step
//		double* dPosition2=EvolveNoIteration(mobility,CurrentStep/2.0);
//		int k=0;
//		//update
//		std::list<Vertex>::iterator itVertex;
//		for(itVertex=ListVertex.begin();itVertex!=ListVertex.end();itVertex++)
//		{
//			itVertex->coord.x+=dPosition2[k];
//			itVertex->coord.y+=dPosition2[k+1];
//			k=k+2;
//			if(xPeriodic||yPeriodic)
//			{
//				ActualiseCrossingVertex(&(*itVertex));
//			}
//		}
//
//		//compute second half-step
//		double* dPosition3=EvolveNoIteration(mobility,CurrentStep/2.0);
//
//
//		//compute error
//		k=0;
//		erreur=0;
//		for(;k<2*nbVertex+nbEdges;k++)
//		{
//			double Delta=fabs(dPosition3[k]+dPosition2[k]-dPosition1[k])*2.0;//le facteur 2 est la car on a une methode d'ordre 2 et donc l'erreur est 2 fois plus grande que la difference
//			double Abs=fmax(fabs(dPosition2[k]+dPosition3[k]),fabs(dPosition1[k]));
//			double scale=absTol+Abs*relTol;
//			erreur=fmax(erreur,Delta/scale);
//			//std::cout<<"Delta="<<Delta<<", scale="<<scale<<", erreur="<<erreur<<"\n";
//
//		}
//		std::cout<<"erreur="<<erreur<<"\n";
//
//		//if error is below tolerance, one can update positions
//		if(erreur<1)
//		{
//			k=0;
//			std::list<Vertex>::iterator itVertex;
//			for(itVertex=ListVertex.begin();itVertex!=ListVertex.end();itVertex++)
//			{
//				itVertex->coord.x+=dPosition2[k];
//				itVertex->coord.y+=dPosition2[k+1];
//				k=k+2;
//				if(xPeriodic||yPeriodic)
//				{
//					ActualiseCrossingVertex(&(*itVertex));
//				}
//			}
//
//			//in that case move Deltat
//			Deltat+=CurrentStep;
//
//			//update cells areas and centroids
//			updateCells();
//
//			//update average residuals
//
//			maxResidual=getAverageResidual(true);
//
//			//add noise possibly
//			int TpsEcoule=(int) Deltat;
//			if(AddNoise&(TpsEcoule%addNoiseFrequency==0)) addNoise(AmpNoise);
//
//			if(DoStretch&(TpsEcoule%stretchFrequency==0))
//			{
//				{
//					stretch();
//				}
//			}
//
//
//			//remove short edges
//			RemoveShortEdges(LengthShortEdge);
//			//apoptosis of small cells
//			ApoptosisSmallCells(5.0);
//
//			//pop vertices
//			PopVertices();
//
//		}
//
//		//step changing
//		double newStep=0.8*CurrentStep/pow(erreur,0.5);
//		if(newStep>10*CurrentStep)
//			CurrentStep=10*CurrentStep;
//		else {
//			if(newStep<CurrentStep/5.0) CurrentStep=CurrentStep/5.0;
//			else CurrentStep=newStep;
//		}
//		//safety limit in the new value of the step
//		CurrentStep=fmin(10,CurrentStep);
//		std::cout<<"CurrentStep="<<CurrentStep<<"\n";
//
//		//count number of iterations
//		nbIterations++;
//		//erase arrays from memory
//		delete dPosition1;
//		delete dPosition2;
//		delete dPosition3;
//
//	}
//
//}
//
//
//
//
//
//
//
//void Tissue::AddVertex(Vertex *v)
//{
//	ListVertex.push_front(*v);
//}
//
//void Tissue::RemoveVertex(Vertex* v)
//{
//	ListVertex.remove(*v);
//}
//
//void Tissue::AddEdge(Edge *e)
//{
//	//add edge in the list of edges
//	ListEdge.push_front(*e);
//
//	//find whether first vertex is already in the list
//	std::list<Vertex>::iterator findv1;
//	findv1 = find( ListVertex.begin(), ListVertex.end(), *(e->v1) );
//
//	//if not, copy it in the list of vertices
//	if( findv1 == ListVertex.end() )
//	{
//		ListVertex.push_front(*(e->v1));
//		ListVertex.front().ConnectedEdges.push_front(&ListEdge.front());
//		//update pointers
//		ListEdge.front().v1=&(ListVertex.front());
//	}
//	//if yes, connect the edge to the existing vertex
//	else {
//		//find whether e is already included in the list of edges to which v1 connects
//		std::list<Edge*>::iterator findeinv1;
//		findeinv1 = find((*findv1).ConnectedEdges.begin(), (*findv1).ConnectedEdges.end(), e );
//		//if not, copy it in the list edges
//		if( findeinv1 == (*findv1).ConnectedEdges.end() )
//		{
//		(*findv1).ConnectedEdges.push_front(&(ListEdge.front()));
//		}
//		ListEdge.front().v1=&(*findv1);
//	}
//
//	//find whether second vertex is already in the list
//	std::list<Vertex>::iterator findv2;
//	findv2 = find( ListVertex.begin(), ListVertex.end(), *(e->v2) );
//
//	//if not, copy it in the list of vertices
//	if( findv2 == ListVertex.end() ) {
//		ListVertex.push_front(*(e->v2));
//		ListVertex.front().ConnectedEdges.push_front(&ListEdge.front());
//		//update pointers
//		ListEdge.front().v2=&(ListVertex.front());
//	}
//	//if yes, connect the edge to the existing vertex
//	else {
//		//find whether e is already included in the list of edges to which v1 connects
//		std::list<Edge*>::iterator findeinv2;
//		findeinv2 = find((*findv2).ConnectedEdges.begin(), (*findv2).ConnectedEdges.end(), e );
//		//if not, copy it in the list edges
//		if( findeinv2 ==(*findv2).ConnectedEdges.end() )
//			{
//			(*findv2).ConnectedEdges.push_front(&(ListEdge.front()));
//			}
//		ListEdge.front().v2=&(*findv2);
//		}
//
//
//
//}
//
//void Tissue::RemoveEdge(Edge *e)
//{
//	std::list<Edge*>::iterator it;
//	for(it=e->v1->ConnectedEdges.begin();it!=e->v1->ConnectedEdges.end();)
//	{
//		if((**it)==*e) it=e->v1->ConnectedEdges.erase(it); else ++it;
//	}
//	for(it=e->v2->ConnectedEdges.begin();it!=e->v2->ConnectedEdges.end();)
//	{
//		if((**it)==*e) it=e->v2->ConnectedEdges.erase(it); else ++it;
//	}
//	//
//
//	//et il faut supprimer les references a l'edge supprime dans les deux cellules qu'il avoisine
//	Vertex v1=*(e->v1);
//	Vertex v2=*(e->v2);
//
//	e->c1->RemoveEdge(e);
//	e->c2->RemoveEdge(e);
//
//	Edge ecopy=*e;
//	ListEdge.remove(ecopy);//enleve l'edge de la liste
//
//}
//void Tissue::AddCell(Cell *c)
//{
//	Point center;
//	center=c->Centroid;
//	ListCell.push_front(*c);//Rajoute a la liste dans le tissu
//	ListCell.front().ListEdges.clear();//on enleve tout car on va tout remplacer par des pointeurs qui vont vers des objets du tissu
//	//il faut egalement actualiser les pointeurs, pour qu'ils pointent vers les objets qui sont dans la liste du tissu
//	std::list<Edge*>::iterator Celledge;
//	for(Celledge=((*c).ListEdges).begin(); Celledge!=((*c).ListEdges).end();Celledge++)
//	{//pour chaque cote on regarde s'il est dans la liste du tissu
//		//std::cout<<"edge recherche : \n";
//		//(**Celledge).Print();
//		//std::cout<<"etat de la liste des edges :\n";
//		//std::list<Edge>::iterator it;
//		//for(it=ListEdge.begin(); it!=ListEdge.end(); it++)
//		//{(*it).Print();
//		//}
//		std::list<Edge>::iterator findCelledge;
//		findCelledge = std::find(ListEdge.begin(), ListEdge.end(), **Celledge);
//
//		if( findCelledge == ListEdge.end() )
//		{//dans ce cas l'edge n'est pas deja dans la liste
//			//std::cout<<"en ajoutant une cellule, un nouveau cote a ete ajoute\n";
//			AddEdge(*Celledge);//copie l'edge dans le tissu s'il n'y etait pas
//			((ListCell.front())).ListEdges.push_front(&((ListEdge.front())));//ajoute le cote a la liste de la cellule qui est dans le tissu;
//			//ici je suppose que les cotes sont deja au courant de quel cellule ils cotoient, car c'est ce que j'ai mis dans le constructeur de cellule, mais peut etre il faudrait rajouter qqch ici
//			//en fait je m'assure quand meme que le cote qu'on vient d'ajouter sait bien ou est sa cellule
//			//en fait c'est probablement absolument necessaire pour etre sur, encore une fois, qu'on pointe bien dans la cellule qui est rangee dans le tissu.
//
//			//if(ListEdge.front().Director().ProduitVectoriel(ListEdge.front().v1->coord.Vector(center))>0)//teste a l'aide d'un produit vectoriel dans quel sens est oriente l'edge
//			if(!(c->Isc2(&ListEdge.front())))
//			{(ListEdge.front()).c1=&(ListCell.front());}
//			else {(ListEdge.front()).c2=&(ListCell.front());}
//		}
//		else
//		{//si on a trouve un edge qui est deja dans la liste, on l'ajoute a la cellule (nouvelle version dans la liste ListCell)
//
//			ListCell.front().ListEdges.push_front(&(*findCelledge));//ajoute le cote a la liste de la cellule qui est dans le tissu;
//			//et j'ajoute la cellule a l'edge en question
//			//std::cout<<"coordonnee du centroide : ("<<center.x<<","<<center.y<<")\n";
//			//
//			   //if((*findCelledge).Director().ProduitVectoriel((*findCelledge).v1->coord.Vector(center))>0 )//teste a l'aide d'un produit vectoriel dans quel sens est oriente l'edge
//				if(!(c->Isc2(&(*findCelledge))))
//			   {(*findCelledge).c1=&((ListCell.front()));
//				//std::cout<<"la cellule a ete ajoute en position 1"<<"\n";
//				}
//				else
//				{
//					(*findCelledge).c2=&((ListCell.front()));
//					//std::cout<<"la cellule a ete ajoute en position 2"<<"\n";
//				}
//			//std::cout<<"edge retrouve appartenant a la cellule "<<ListCell.front().Number<<":\n";
//			//(*findCelledge).Print();
//		}
//
//	}
//
//
//}
//
//void Tissue::RemoveCell(Cell *c)
//{
//	ListCell.remove(*c);//enleve de la liste dans le tissu
//
//	// si ca se revele utile il faudra ajouter qqch pour effacer les references a cette cellules dans les cotes qui l'entourent
//}
//
//
//void Tissue::clear()
//{
//	ListVertex.clear();
//	ListEdge.clear();
//	ListCell.clear();
//	xPeriodic=false;
//	yPeriodic=false;
//	Lx=0;
//	Ly=0;
//}
//
//
//
//
//void Tissue::RemoveShortEdges(double minlength)
//{
//	std::list<Edge>::iterator it;
//	std::list<Edge>::iterator itSouvenir;
//	for(it=ListEdge.begin();it!=ListEdge.end();)
//	{
//		if(it->Length()<minlength)//si le cote est trop petit il faut le supprimer :
//		{
//			QLOG_INFO() << "Remove edge between "<<it->c1->Number<<" and "<<it->c2->Number<<" with length "<<it->Length();
//			itSouvenir=it;
//			it++;//il faut iterer avant de le supprimer
//			if(((itSouvenir->c1->A>20||itSouvenir->c1->Number==0)
//				 &&(itSouvenir->c2->A>20||itSouvenir->c2->Number==0))//evite la suppression des cotes des petites cellules
//			&&((itSouvenir->c1->NbEdges()>3||itSouvenir->c1->Number==0)
//			   &&(itSouvenir->c2->NbEdges()>3||itSouvenir->c2->Number==0))//evite qu'on se retrouve avec une cellule a deux cotes
//				)
//			{
//			UnPop(&(*itSouvenir));
//			//UnPopAndPop(&(*itSouvenir));
//			}
//
//
//		}
//		else
//		{
//			it++;
//		}
//	}
//}
//
////T1 transition where opposite configuration is forced
//void Tissue::UnPopAndPop(Edge* e)
//{
//
//
//	if((!e->v1->xfixed&&!e->v2->xfixed)&&(!e->v1->yfixed&&!e->v2->yfixed))
//	{
//		std::list<Edge>::iterator finde;
//		std::list<Edge*>::iterator it2;
//
//		//make sure the edge is the one belonging to ListEdge
//		finde = find(ListEdge.begin(), ListEdge.end(), *e);
//
//		finde->c1->AreEdgesSorted=false;
//        finde->c2->AreEdgesSorted=false;
//		//find cells involved in the T1
//		Cell* c1old=finde->c1;
//		Cell* c2old=finde->c2;
//		Cell* c1new;
//		Cell* c2new;
//		Edge* e1=0;
//		Edge* e2=0;
//		Edge* e3=0;
//		Edge* e4=0;
//
//		for(it2=finde->v1->ConnectedEdges.begin();it2!=finde->v1->ConnectedEdges.end();it2++)
//		{
//			if(!(*(*it2)==(*finde)))
//			{
//
//				if((*((*it2)->v1)==*(finde->v1))&(*((*it2)->c2)==*c1old)) {c2new=(*it2)->c1;e1=*it2;}
//				if((*((*it2)->v2)==*(finde->v1))&(*((*it2)->c1)==*c1old)) {c2new=(*it2)->c2;e1=*it2;}
//
//				if((*((*it2)->v1)==*(finde->v1))&(*((*it2)->c1)==*c2old)) {e2=*it2;}
//				if((*((*it2)->v2)==*(finde->v1))&(*((*it2)->c2)==*c2old)) {e2=*it2;}
//
//
//			}
//		}
//
//		for(it2=finde->v2->ConnectedEdges.begin();it2!=finde->v2->ConnectedEdges.end();it2++)
//		{
//			if(!(*(*it2)==(*finde)))
//			{
//
//				if((*((*it2)->v1)==*(finde->v2))&(*((*it2)->c1)==*c1old)) {c1new=(*it2)->c2;e3=*it2;}
//				if((*((*it2)->v1)==*(finde->v2))&(*((*it2)->c2)==*c2old)) {e4=*it2;}
//				if((*((*it2)->v2)==*(finde->v2))&(*((*it2)->c2)==*c1old)) {c1new=(*it2)->c1;e3=*it2;}
//				if((*((*it2)->v2)==*(finde->v2))&(*((*it2)->c1)==*c2old)) {e4=*it2;}
//			}
//		}
//
//		//move kept vertex (v2) toward the middle of the old edge
//		Point v1v2=finde->v1v2();
//		finde->v2->coord=finde->v2->coord-v1v2*(0.5);
//		if(xPeriodic||yPeriodic)
//		{
//			Point v1coordsouvenir=finde->v1->coord;
//			finde->v1->coord=finde->v1->coord+(v1v2*(0.5));
//				if(fabs(fmod(finde->v2->coord.x,Lx))<0.0001)
//				{
//					finde->v2->coord.x-=fmod(finde->v2->coord.x,Lx);
//				    finde->v1->coord.x-=fmod(finde->v1->coord.x,Lx);
//				}
//				if(fabs(fmod(finde->v2->coord.y,Ly))<0.0001)
//				{
//					finde->v2->coord.y-=fmod(finde->v2->coord.y,Ly);
//					finde->v1->coord.y-=fmod(finde->v1->coord.y,Ly);
//				}
//
//			ActualiseCrossingVertex(finde->v1);
//			ActualiseCrossingVertex(finde->v2);
//			finde->v1->coord=v1coordsouvenir;//necessary to avoid confusion between the two vertices
//
//		}
//		//connect edges touching vertex v1 to vertex v2
//
//		for(it2=finde->v1->ConnectedEdges.begin();it2!=finde->v1->ConnectedEdges.end();it2++)
//		{
//			if(!(**it2==*finde))
//			{
//				if(*((**it2).v1)==*(finde->v1)) (*it2)->v1=&(*(finde->v2)); else (*it2)->v2=&(*(finde->v2));//fais le remplacement du vertex v2 par v1 dans les edges touchant v2 (il faut d'abord voir si le v2 en cours est le v1 ou le v2 de l'edge avoisinant)
//				//il faut egalement completer les edges accroches a v2
//				finde->v2->ConnectedEdges.push_front(&(**it2));
//			}
//		}
//		//suppress edge and vertex v1
//		Vertex* v1souvenir=finde->v1;
//		Vertex* v2=(finde->v2);
//		RemoveEdge(&(*finde));
//		RemoveVertex(v1souvenir);
//
//		//partie pop
//
//
//		Point ForceDirection1=e1->ForceDirection(v2);//direction du premier vertex
//		Point ForceDirection2=e3->ForceDirection(v2);//direction du second vertex : la somme des deux fixe le deplacement
//		Point SommeForces=ForceDirection1+ForceDirection2;
//		double epsilon=3.5;//norme du deplacement elementaire du vertex
//		//il faut s'assurer que ce deplacement est toujours superieur a la limite de suppression des vertex sinon ca bloque...
//		double NormeForce=SommeForces.Norme();
//		Point deplacement=SommeForces*(epsilon/NormeForce/2.0);//je divise par 2 pour chaque moitie
//		Vertex newvtemp=Vertex(v2->coord+deplacement);//cree un vertex translate dans une direction
//		AddVertex(&newvtemp);
//		Vertex * newv=&(ListVertex.front());//cree le vertex, le range dans la liste, cree le pointeur qui va bien vers la liste
//		v2->coord=v2->coord-deplacement;//et translate le vertex en question dans une autre direction
//		//si c'est periodique il faut voir si l'un des deux nouveaux vertex sort du cadre (fait plus bas)
//
//
//		//actualise les edges qui recoivent le nouveau vertex
//		(newv->ConnectedEdges).push_front(e1);
//		(newv->ConnectedEdges).push_front(e3);
//		if((e1->c1)==c1old) {e1->v2=newv;}
//		if((e1->c2)==c1old) {e1->v1=newv;}
//		if((e3->c1)==c1old) {e3->v1=newv;}
//		if((e3->c2)==c1old) {e3->v2=newv;}
//
//		for(it2=v2->ConnectedEdges.begin();it2!=v2->ConnectedEdges.end();)
//		{
//			if(**it2==*e1||**it2==*e3) it2=v2->ConnectedEdges.erase(it2); else ++it2;
//		}
//
//		Edge newetemp=Edge(newv,v2);
//		newetemp.c1=c1new;
//		newetemp.c2=c2new;
//		newetemp.Tension=Tensions[newetemp.c1->Type][newetemp.c2->Type];
//		AddEdge(&newetemp);
//		Edge * newe=&(ListEdge.front());//meme principe que pour le vertex : d'abord ajouter, ensuite prendre le pointeur vers l'element ajoute
//		if(c1new->Number!=0) c1new->ListEdges.push_front(newe);
//		if(c2new->Number!=0) c2new->ListEdges.push_front(newe);
//
//		//actualiseTensions();
//		//il me semble que partant du principe que le nouveau edge a une taille presque 0 au depart, on peut faire comme si il n'avait pas bouge et
//		//actualiser ses deux vertex pour qu'il croise correctement la ligne
//		ActualiseCrossingVertex(v2);
//		ActualiseCrossingVertex(newv);
//
//		//il faut annuler si la force n'est pas dans le bon sens
//		//temporairement je l'enleve
//		double Projection=newe->ForceExtension();
//		if(Projection<0) UnPop(newe);
//
//	}
//
//}
//
//
//Edge* Tissue::Pop(Vertex* v, Edge* e1, Edge* e3, Cell* c1, Cell* c2, Cell* c3, Cell* c4, double epsilon)
//{
//	//epsilon is the displacement imposed when forces are tested
//	Point e1Director;
//	e1Director=e1->MidPoint(v);
//	e1Director=e1Director*(1.0/e1Director.Norme());
//	Point e3Director;
//	e3Director=e3->MidPoint(v);
//	e3Director=e3Director*(1.0/e3Director.Norme());
//
//		std::list<Edge*>::iterator it2;
//
//		Point ForceDirection1=e1->ForceDirection(v);//direction of the first vertex
//		Point ForceDirection2=e3->ForceDirection(v);//direction of the second vertex
//		Point SommeForces=ForceDirection1+ForceDirection2;
//		double NormeForce=SommeForces.Norme();
//		Point deplacement=SommeForces*(epsilon/NormeForce/2.0);
//		Vertex newvtemp=Vertex(v->coord+deplacement);//creates the new vertex
//		AddVertex(&newvtemp);
//		Vertex * newv=&(ListVertex.front());
//		v->coord=v->coord-deplacement;//moves the old vertex
//
//		//connect new vertex to the relevant edges
//		(newv->ConnectedEdges).push_front(e1);
//		(newv->ConnectedEdges).push_front(e3);
//		if(*(e1->v1)==*v) {e1->v1=newv;}
//		if(*(e1->v2)==*v) {e1->v2=newv;}
//		if(*(e3->v1)==*v) {e3->v1=newv;}
//		if(*(e3->v2)==*v) {e3->v2=newv;}
//
//		for(it2=v->ConnectedEdges.begin();it2!=v->ConnectedEdges.end();)
//		{
//			if(**it2==*e1||**it2==*e3) it2=v->ConnectedEdges.erase(it2); else ++it2;
//		}
//
//	//create new edge
//	Edge newetemp=Edge(newv,v);
//	Point newEdgeDirector=newetemp.v1v2();
//	newEdgeDirector=newEdgeDirector*(1.0/newEdgeDirector.Norme());
//	double pv=(e3Director+e1Director*(-1)).ProduitVectoriel(newEdgeDirector+e3Director*(-1));
//
//			if(pv>0)
//				{
//				newetemp.c1=c2;
//				newetemp.c2=c4;
//				}
//			else
//				{
//				newetemp.c1=c4;
//				newetemp.c2=c2;
//				}
//		newetemp.Tension=Tensions[newetemp.c1->Type][newetemp.c2->Type];
//		AddEdge(&newetemp);
//		Edge* newe=&(ListEdge.front());
//		if(c2->Number!=0) c2->ListEdges.push_front(newe);
//		if(c4->Number!=0) c4->ListEdges.push_front(newe);
//		//actualiseTensions();
//		ActualiseCrossingVertex(v);
//		ActualiseCrossingVertex(newv);
//
//
//
//
//		return(newe);
//}
//
//void Tissue::PopVertices()
//{
//	std::list<Vertex>::iterator it;
//	std::list<Edge*>::iterator it2;
//	bool OnePopDone=false;
//	//test for all vertices
//	for(it=ListVertex.begin();it!=ListVertex.end();++it)
//	{
//		if(!(it->xfixed)&&!(it->yfixed))//fixed vertices are not considered
//		{
//			if(it->ConnectedEdges.size()>=4)
//			{
//				//4 edges involved in T1 transition
//				Edge* e1;
//				Edge* e2;
//				Edge* e3;
//				Edge* e4;
//				//4 cells involved in T1 transition
//				Cell *c1;
//				Cell *c2;
//				Cell *c3;
//				Cell *c4;
//				it2=it->ConnectedEdges.begin();
//				e1=*(it2);
//				c1=e1->c1;
//				c2=e1->c2;
//				for(;it2!=it->ConnectedEdges.end();it2++)
//				{
//					if(*((*it2)->c1)==*(c1)) {e3=(*it2);c4=(*it2)->c2;}
//					if(*((*it2)->c2)==*(c1)) {e3=(*it2);c4=(*it2)->c1;}
//					if(*((*it2)->c1)==*(c2)) {e2=(*it2);c3=(*it2)->c2;}
//					if(*((*it2)->c2)==*(c2)) {e2=(*it2);c3=(*it2)->c1;}
//				}
//				//find e4 which is the neighbor of e2
//				for(it2=it->ConnectedEdges.begin();it2!=it->ConnectedEdges.end();it2++)
//				{
//					if((*((*it2)->c1)==*c3&!(*((*it2)->c2)==*c2))||(*((*it2)->c2)==*c3&!(*((*it2)->c1)==*c2)))
//						e4=*it2;
//				}
//
//				//everything is ready for test of different pops
//
//				//first test
//				Edge * newe=Pop(&(*it),e1,e3,c1,c2,c3,c4,LengthAfterT1);
//				//compute whether the forces exerted on the new vertex tends to open or close the newly formed edge
//				double Projection1=newe->ForceExtension();
//				UnPop(newe);
//
//				//second test
//				newe=Pop(&(*it),e2,e1,c2,c3,c4,c1,LengthAfterT1);
//				double Projection2=newe->ForceExtension();
//				UnPop(newe);
//
//				//decides which T1 to do according to different forces
//
//				if(Projection1<0&Projection2>0)
//				{//if second T1 is favorable, not the first one, one does the second one
//					newe=Pop(&(*it),e2,e1,c2,c3,c4,c1,LengthAfterT1);
//					c1->AreEdgesSorted=false;
//					c2->AreEdgesSorted=false;
//					c3->AreEdgesSorted=false;
//					c4->AreEdgesSorted=false;
//					OnePopDone=true;
//				}
//
//				if((Projection1>0&Projection2<0))
//				{
//					//if first T1 is favorable, not the second one, one does the first one
//					newe=Pop(&(*it),e1,e3,c1,c2,c3,c4,LengthAfterT1);
//					c1->AreEdgesSorted=false;
//					c2->AreEdgesSorted=false;
//					c3->AreEdgesSorted=false;
//					c4->AreEdgesSorted=false;
//					OnePopDone=true;
//
//				}
//				if(Projection1<0&Projection2<0)
//				{
//					//neither of them is favorable
//					//std::cout<<"aucun favorable\n";
//					//QLOG_INFO() <<"Not Popped";
//					//UnPop(newe);//cas defavorable pour tous
//				}
//				if(Projection1>0&Projection2>0)
//				{
//					//all T1s are favorable
//					if(ProbabilisticPop)//random decision
//					{
//						double randomnumber=(rand()%1000+1);
//						if(randomnumber>500)
//						{
//							newe=Pop(&(*it),e1,e3,c1,c2,c3,c4,LengthAfterT1);
//							c1->AreEdgesSorted=false;
//							c2->AreEdgesSorted=false;
//							c3->AreEdgesSorted=false;
//							c4->AreEdgesSorted=false;
//							OnePopDone=true;
//						}
//					}
//					else
//					{
//						//decision based on comparison of magnitude of forces
//						if(Projection1>Projection2)
//						{
//							newe=Pop(&(*it),e1,e3,c1,c2,c3,c4,LengthAfterT1);
//							c1->AreEdgesSorted=false;
//							c2->AreEdgesSorted=false;
//							c3->AreEdgesSorted=false;
//							c4->AreEdgesSorted=false;
//							OnePopDone=true;
//
//						}
//						else
//						{
//						newe=Pop(&(*it),e2,e1,c2,c3,c4,c1,LengthAfterT1);
//						c1->AreEdgesSorted=false;
//						c2->AreEdgesSorted=false;
//						c3->AreEdgesSorted=false;
//						c4->AreEdgesSorted=false;
//						OnePopDone=true;
//						}
//
//
//					}
//				}
//
//
//
//			}
//
//
//
//			}
//
//	}
//
//	//updates
//	if(OnePopDone)
//	{
//
//		QLOG_INFO() <<"Popped";
//		std::list<Cell>::iterator itCell;
//		for(itCell=ListCell.begin();itCell!=ListCell.end();++itCell)
//		{
//			itCell->AreEdgesSorted=false;
//		}
//		updateCells();
//
//	}
//	else
//	{
//		//QLOG_INFO() <<"No Pop Done";
//	}
//
//}
//
//void Tissue::Mitosis(Cell* c)//mitosis with random angle
//{
//	QLOG_INFO() << "Mitosis of cell "<<c->Number;
//
//	srand ( time(NULL) );//initialisation de random
//	//genere un angle au hasard ou separer la cellule
//	double angle = (rand() % 180+ 1);//necessaire pour generer un chiffre entre 0 et PI, avec 1/180 de possibilite
//	Mitosis(c,angle);
//}
//
//int Tissue::Mitosis(Cell* c, double angle)
//{
//
//	QLOG_INFO() << "Mitosis of cell "<<c->Number;
//
//	//transforme l'angle de degre en radian
//	angle=angle*PI/180;
//
//	bool CrossBoundary=c->CrossBoundary();
//	//il faut enregistrer ces valeurs avant qu'on les modifie
//	Point center;
//	Cell csouvenir=*c;
//	std::list<Edge*>::iterator it;
//	Point p;
//	Point p1=Point(0,0);
//	Point p2=Point(0,0);
//	Edge * e1;
//	Edge * e2;
//	int flag =0;
//
//
//
//	if(CrossBoundary) c->cellSameScreen();
//
//	c->UpdateCentroid();
//	center=c->Centroid;
//	for(it=c->ListEdges.begin();it!=c->ListEdges.end();it++)
//	{
//
//		p=(**it).Intersection(center,angle);
//		if(!(p.x==0&&p.y==0))
//		{
//
//			flag++;
//			if (p1.x==0&&p1.y==0)
//			{
//				p1.x=p.x;
//				p1.y=p.y;
//				e1=*it;
//			}
//			else
//			{
//				p2.x=p.x;
//				p2.y=p.y;
//				e2=*it;
//			}
//		}
//	}
//
//
//
//	//il faut enregistrer ces valeurs avant qu'on les modifie
//
//
//	if (flag==2){
//        if(!((p1.x==0&&p1.y==0)||(p2.x==0&&p2.y==0)))//on verifie qu'on a bien trouve deux solutions
//		{
//			Vertex e1v1souvenir=*e1->v1;
//			Vertex e1v2souvenir=*e1->v2;
//			Vertex e2v1souvenir=*e2->v1;
//			Vertex e2v2souvenir=*e2->v2;
//			Vertex newv1inter=Vertex(p1);
//			Vertex newv2inter=Vertex(p2);
//			AddVertex(&newv1inter);
//			Vertex * newv1=&(ListVertex.front());
//			AddVertex(&newv2inter);
//			Vertex * newv2=&(ListVertex.front());
//			Edge neweinter=Edge(newv1,newv2);//relie les deux nouveaux vertex bien sur
//			ListEdge.push_front(neweinter);//on l'ajoute au tissu
//			Edge * newe=&(ListEdge.front());//comme d'hab il faut prendre l'adresse de ce qu'on ajoute
//			newe->Tension=Tensions[c->Type][c->Type];
//			//newe->Lx=e1->Lx;
//			//newe->Ly=e1->Ly;
//			newv1->ConnectedEdges.push_front(newe);
//			newv2->ConnectedEdges.push_front(newe);
//			std::cout<<"newv1=";
//			newv1->coord.Print();
//			std::cout<<"newv2=";
//			newv2->coord.Print();
//
//			Cell newcinter=Cell();
//			ListCell.push_front(newcinter);
//			Cell* newc=&(ListCell.front());
//			double pv1;
//			double pv2;//important pour savoir ou sont les vertex v1 de e1 et e2
//			Point e1v1local;
//			Point e1v2local;
//			Point e2v1local;
//			Point e2v2local;
//
//			double outsider=0;
//			Point cnv1= newv1->coord-center;
//
//			Point cnv2= newv2->coord-center;
//			if (cnv1.x*cnv2.x+cnv1.y*cnv2.y > 0)
//			{
//				if (cnv1.Norme()>cnv2.Norme())
//					outsider = 1;
//				else outsider = -1;
//
//			}
//
//
//			for(it=c->ListEdges.begin();it!=c->ListEdges.end();)
//			{
//				if(!(**it==*e1||**it==*e2))
//				{
//					Point v1coord=(*it)->v1->coord;
//
//					double pv=((newv1->coord-center).ProduitVectoriel(v1coord-center));
//					if(outsider!=0) pv*=outsider;
//					if(pv<0)
//					{
//						//rajouter les nouveaux
//						newc->ListEdges.push_front(*it);
//						if(((**it).c2)->Number==csouvenir.Number) (**it).c2=newc; else (**it).c1=newc;
//						//enlever les edges de l'ancienne cellule
//						it=c->ListEdges.erase(it);
//					}
//					else it++;
//				}
//				else it++;
//			}
//
//			pv1=((newv1->coord-center).ProduitVectoriel(e1->v1->coord-center));
//			pv2=((newv2->coord-center).ProduitVectoriel(e2->v1->coord-center));
//			if(outsider)
//			{ pv1*=outsider;
//				pv2*=(0-outsider);
//			}
//
//
//
//			Vertex** voldcelle1;
//			Vertex** vnewcelle1;
//			Vertex** voldcelle2;
//			Vertex** vnewcelle2;
//			//std::cout<<"e1v1=\n";
//			//	e1->v1->Print();
//			std::cout<<"pv1="<<pv1<<", pv2="<<pv2<<"\n";
//			if(pv1>=0) {voldcelle1=&(e1->v1);vnewcelle1=&(e1->v2);} else {voldcelle1=&(e1->v2);vnewcelle1=&(e1->v1);}
//			if(pv2>=0) {voldcelle2=&(e2->v2);vnewcelle2=&(e2->v1);} else {voldcelle2=&(e2->v1);vnewcelle2=&(e2->v2);}
//			//std::cout<<"vnewcelle1:\n";
//			//(**vnewcelle1).Print();
//			Edge newe1inter=Edge(newv1,*vnewcelle1);
//			ListEdge.push_front(newe1inter);
//			Edge * newe1=&(ListEdge.front());
//			newe1->Tension=e1->Tension;
//			//newe1->Lx=e1->Lx;
//			//newe1->Ly=e1->Ly;
//			for(it=(*vnewcelle1)->ConnectedEdges.begin();it!=(*vnewcelle1)->ConnectedEdges.end();)
//			{
//				if((**it)==*e1) it=(*vnewcelle1)->ConnectedEdges.erase(it); else ++it;
//			}
//			(*vnewcelle1)->ConnectedEdges.push_front(newe1);
//			*vnewcelle1=newv1; //this requires understanding.....
//
//			newv1->ConnectedEdges.push_front(newe1);
//			newv1->ConnectedEdges.push_front(e1);
//			Edge newe2inter=Edge(newv2,*vnewcelle2);
//			ListEdge.push_front(newe2inter);
//			Edge * newe2=&(ListEdge.front());
//			newe2->Tension=e2->Tension;
//			std::list<Edge*>::iterator it;
//			for(it=(*vnewcelle2)->ConnectedEdges.begin();it!=(*vnewcelle2)->ConnectedEdges.end();)
//			{
//				if((**it)==*e2) it=(*vnewcelle2)->ConnectedEdges.erase(it); else ++it;
//			}
//			(*vnewcelle2)->ConnectedEdges.push_front(newe2);
//			*vnewcelle2=newv2;
//			newv2->ConnectedEdges.push_front(newe2);
//			newv2->ConnectedEdges.push_front(e2);
//
//
//			//ensuite il faut creer une nouvelle cellule
//
//			newc->ListEdges.push_front(newe);
//			newc->ListEdges.push_front(newe1);
//			newc->ListEdges.push_front(newe2);
//			newc->K=c->K;
//			newc->A0=c->A0;
//			//c->A0=gloA0;
//			newc->Type=c->Type;//la nouvelle cellule recupere les attributs de l'ancienne
//
//			//il vaut mieux donner un numero superieur au maximum
//			int maxNumber=0;
//			std::list<Cell>::iterator itCell;
//			for(itCell=ListCell.begin();itCell!=ListCell.end();itCell++)
//			{
//				if(maxNumber<itCell->Number) maxNumber=itCell->Number;
//			}
//			newc->Number=maxNumber+1;
//			//ajoute la cellule et garde la reference de la copie
//			c->ListEdges.push_front(newe);
//
//			//ensuite il faut reaffecter les cellules correctement
//			newe->c2=c;
//			newe->c1=newc;
//
//			if((e1->c1->Number)==c->Number)
//			{
//				if(e1->c2->Number!=0)
//				{
//					e1->c2->ListEdges.push_front(newe1);
//					e1->c2->AreEdgesSorted=false;
//				}
//				if(*(e1->v1)==*newv1) {newe1->c1=e1->c2;newe1->c2=newc;} else {newe1->c2=e1->c2;newe1->c1=newc;}
//
//			}
//			else
//			{
//				if(e1->c1->Number!=0)
//				{
//					e1->c1->ListEdges.push_front(newe1);
//					e1->c1->AreEdgesSorted=false;
//				}
//				if(*(e1->v1)==*newv1) {newe1->c2=e1->c1;newe1->c1=newc;} else {newe1->c1=e1->c1;newe1->c2=newc;}
//			}
//			if((e2->c1->Number)==c->Number)
//			{
//				if(e2->c2->Number!=0)
//				{ e2->c2->ListEdges.push_front(newe2);
//					e2->c2->AreEdgesSorted=false;
//				}
//				if(*(e2->v1)==*newv2) {newe2->c1=e2->c2;newe2->c2=newc;} else {newe2->c2=e2->c2;newe2->c1=newc;}
//			}
//			else
//			{
//				if(e2->c1->Number!=0)
//				{
//					e2->c1->ListEdges.push_front(newe2);
//					e2->c1->AreEdgesSorted=false;;
//
//				}
//				if(*(e2->v1)==*newv2) {newe2->c2=e2->c1;newe2->c1=newc;} else {newe2->c1=e2->c1;newe2->c2=newc;}
//			}
//			c->sortEdges();
//			newc->sortEdges();
//
//			if(CrossBoundary)
//			{
//
//
//				Point v1memor;
//				v1memor.x=newv1->coord.x;
//				v1memor.y=newv1->coord.y;
//
//				Point v2memor=newv2->coord;
//				v2memor.x=newv2->coord.x;
//				v2memor.y=newv2->coord.y;
//				c->cellCrossScreen();
//
//				newv1->coord.x=v1memor.x;
//				newv1->coord.y=v1memor.y;
//				newv2->coord.x=v2memor.x;
//				newv2->coord.y=v2memor.y;
//				newc->cellCrossScreen();
//
//			}
//			c->A=c->Area();
//			newc->A=newc->Area();
//
//
//			std::list<Edge*>::iterator iteP;
//        }
//        else
//		{
//			//std::cout<<"pas de croisement trouve\n";
//		}
//
//        std::cout<< ("End of Mitosis/n");
//        return 1;
//    }
//	else
//	{
//		if(CrossBoundary) c->cellCrossScreen();
//		return 0;
//	}
//
//}
//
//void Tissue::Apoptosis(Cell* c)
//{
//
//	QLOG_INFO() << "Apoptosis of cell "<<c->Number;
//
//	std::list<Vertex>::iterator itVertex;
//	if (c->AreEdgesSorted==false)
//    c->sortEdges();
//	std::list<Edge*>::iterator itEdge;
//
//	itEdge=c->ListEdges.begin();
//	Vertex* v=(*itEdge)->v1;//ce vertex sera le seul a etre garde
//	//je change ses coordonnees pour qu'il soit au milieu
//	v->coord=c->Centroid;
//	//il faut lui enlever les deux voisins qu'il ne gardera pas
//	std::list<Edge*>::iterator itEdge2;
//	std::list<Edge>::iterator itEdge3;
//for(itEdge2=v->ConnectedEdges.begin();
//	itEdge2!=v->ConnectedEdges.end();)
//{
//	if((*((*itEdge2)->c1)==*c)||(*((*itEdge2)->c2)==*c))
//	{
//		itEdge2=v->ConnectedEdges.erase(itEdge2);
//	}
//	else itEdge2++;
//}
//	Vertex* Vdejatraite=(*itEdge)->v1;
//	Vertex* Vatraiter;
//	Cell* cvoisine;
//	for(;itEdge!=c->ListEdges.end();itEdge++)
//	{
//		//d'abord trouver le bon vertex
//		if(*((*itEdge)->v1)==*Vdejatraite) Vatraiter=(*itEdge)->v2; else Vatraiter=(*itEdge)->v1;
//		//ensuite faire toutes les modifications
//			//d'abord a ce stade on peut supprimer le vertex precedent
//			if(!(*Vdejatraite==*v))
//			{
//				for(itVertex=ListVertex.begin();
//					itVertex!=ListVertex.end();)
//				{
//					if(*itVertex==*Vdejatraite) itVertex=ListVertex.erase(itVertex); else itVertex++;
//				}
//			}
//			//d'abord il faut reconnecter le vertex
//		if(!(*Vatraiter==*v))//pour etre sur qu'apres avoir fait un tour on ne reconnecte pas le vertex de base
//		{
//			for(itEdge2=Vatraiter->ConnectedEdges.begin();
//				itEdge2!=Vatraiter->ConnectedEdges.end();
//				itEdge2++)
//			{
//				if(!(*((*itEdge2)->c1)==*c)&&!(*((*itEdge2)->c2)==*c))//on cherche les cotes qui seront gardes
//				{
//					if(*((*itEdge2)->v1)==*Vatraiter)
//					{(*itEdge2)->v1=v;} else (*itEdge2)->v2=v;
//					v->ConnectedEdges.push_front(*itEdge2);
//				}
//
//			}
//		}
//
//			//ensuite il faut supprimer la mention du cote dans la cellule voisine
//		if((*(*itEdge)->c1)==*c) cvoisine=(*itEdge)->c2; else cvoisine=(*itEdge)->c1;
//			cvoisine->AreEdgesSorted=false;
//			for(itEdge2=cvoisine->ListEdges.begin();
//				itEdge2!=cvoisine->ListEdges.end();)
//			{
//				if((**itEdge2)==(**itEdge)) itEdge2=cvoisine->ListEdges.erase(itEdge2); else itEdge2++;
//			}
//
//		//il faut aussi supprimer la mention de l'edge dans la liste
//		for(itEdge3=ListEdge.begin();
//			itEdge3!=ListEdge.end();)
//		{
//			if(*itEdge3==**itEdge) itEdge3=ListEdge.erase(itEdge3); else itEdge3++;
//		}
//		//avant de continuer il faut dire ou on en est
//		Vdejatraite=Vatraiter;
//	}
//
////enfin il faut supprimer la cellule elle-meme
//	//c->ListEdges=std::list<Edge*>();//je fais ca parce qu'on dirait que ca lui pose un probleme
//	std::list<Cell>::iterator itCell;
//	for(itCell=ListCell.begin();itCell!=ListCell.end();)
//	{
//		if(*itCell==*c) {ListCell.erase(itCell); return;} else itCell++;
//	}
//}
//
//void Tissue::ApoptosisSmallCells(double minArea)
//{
//	std::list<Cell>::iterator itCell;
//	bool ApoptosisDone=false;
//	for(itCell=ListCell.begin();itCell!=ListCell.end();)
//	{
//		if(itCell->A<minArea)
//		{
//			if(!(itCell->CrossBoundary()))
//			{
//				QLOG_INFO() << "Small cell number"<<itCell->Number<<" found with area "<<itCell->A;
//				Cell* c=&(*itCell);itCell++;Apoptosis(c); ApoptosisDone=true;
//			}
//			else itCell++;
//		}
//			else itCell++;
//
//	}
//
//	if(ApoptosisDone) updateCells();
//}
//void Tissue::xStretch(double Tension, double Text,double mobility)
//{
//	if(StretchUniformly)
//	{
//
//		std::list<Vertex>::iterator it;
//		for(it=ListVertex.begin();it!=ListVertex.end();it++)
//		{
//		it->coord.x=it->coord.x+mobility*(Text-Tension)*(it->coord.x);
//		}
//	}
//
//	if(xPeriodic)
//		{
//			Lx=Lx*(1+mobility*(Text-Tension));
//			actualiseLx(Lx);
//		}
//	//QLOG_INFO() << "Stretch in x by a factor "<<mobility*(Text-Tension)<<", ext. stress "<<Text<<", int. stress "<<Tension;
//
//}
//
//void Tissue::yStretch(double Tension, double Text,double mobility)
//{
//	if(StretchUniformly)
//	{
//	std::list<Vertex>::iterator it;
//	for(it=ListVertex.begin();it!=ListVertex.end();it++)
//	{
//		it->coord.y=it->coord.y+mobility*(Text-Tension)*(it->coord.y);
//	}
//	}
//	if(yPeriodic)
//	{Ly=Ly*(1+mobility*(Text-Tension));
//	actualiseLy(Ly);
//	}
//	QLOG_INFO() << "Stretch in y by a factor "<<mobility*(Text-Tension);
//
//	//Validite();
//}
//
//
//double Tissue::xm()//renvoie l'abscisse moyenne du tissu
//{
//	if(xPeriodic) return(Lx/2.0);
//
//	double xm=0;
//	double compteur=0;
//	std::list<Vertex>::iterator it;
//	for(it=ListVertex.begin();it!=ListVertex.end();it++)
//	{
//		xm+=it->coord.x;
//		compteur++;
//	}
//	return(xm/compteur);
//}
//
//double Tissue::ym()//renvoie l'abscisse moyenne du tissu
//{
//	if(yPeriodic) return(Ly/2.0);
//	double ym=0;
//	double compteur=0;
//	std::list<Vertex>::iterator it;
//	for(it=ListVertex.begin();it!=ListVertex.end();it++)
//	{
//		ym+=it->coord.y;
//		compteur++;
//	}
//	return(ym/compteur);
//}
//
//double Tissue::longueurx()//renvoie la longueur en x du tissu (approximativement)
//{
//	if(xPeriodic) return(Lx);
//
//	else
//	{
//		double xmax=-10e10;
//		double xmin=10e10;
//		std::list<Vertex>::iterator it;
//		for(it=ListVertex.begin();it!=ListVertex.end();it++)
//		{
//			if(it->coord.x>xmax) xmax=it->coord.x;
//			if(it->coord.x<xmin) xmin=it->coord.x;
//		}
//		return(xmax-xmin);
//	}
//}
//
//double Tissue::longueury()//renvoie la longueur en y du tissu (approximativement)
//{
//	if(yPeriodic) return(Ly);
//	else
//	{
//		double ymax=-10e10;
//		double ymin=10e10;
//		std::list<Vertex>::iterator it;
//		for(it=ListVertex.begin();it!=ListVertex.end();it++)
//		{
//			if(it->coord.y>ymax) ymax=it->coord.y;
//			if(it->coord.y<ymin) ymin=it->coord.y;
//		}
//		return(ymax-ymin);
//	}
//}
//
//
//void Tissue::ActualiseCrossingVertex(Vertex *v)
//{
//	if(xPeriodic&&v->coord.x>=Lx)
//	{
//		std::list<Edge*>::iterator itEdge;
//		for(itEdge=v->ConnectedEdges.begin();itEdge!=v->ConnectedEdges.end();itEdge++)
//		{
//			if((**itEdge).Isv1(v)&&!((**itEdge).xCrossBoundary12||(**itEdge).xCrossBoundary21))
//
//				(*itEdge)->xCrossBoundary21=true;
//				else{
//			if((**itEdge).Isv2(v)&&!((**itEdge).xCrossBoundary12||(**itEdge).xCrossBoundary21))
//				(*itEdge)->xCrossBoundary12=true;
//							else{
//			if((**itEdge).Isv1(v)&&((**itEdge).xCrossBoundary12))
//				(*itEdge)->xCrossBoundary12=false;
//								else{
//			if((**itEdge).Isv2(v)&&((**itEdge).xCrossBoundary21))
//				(*itEdge)->xCrossBoundary21=false;
//								}
//							}
//						}
//		}
//		v->coord.x-=Lx;
//	}
//
//	if(xPeriodic&&(v->coord.x<0))
//	{
//		std::list<Edge*>::iterator itEdge;
//		for(itEdge=v->ConnectedEdges.begin();itEdge!=v->ConnectedEdges.end();itEdge++)
//		{
//			if((**itEdge).Isv1(v)&&!((**itEdge).xCrossBoundary12||(**itEdge).xCrossBoundary21))
//				(*itEdge)->xCrossBoundary12=true;
//			else{
//				if((**itEdge).Isv2(v)&&!((**itEdge).xCrossBoundary12||(**itEdge).xCrossBoundary21))
//					(*itEdge)->xCrossBoundary21=true;
//				else{
//					if((**itEdge).Isv1(v)&&((**itEdge).xCrossBoundary21))
//						(*itEdge)->xCrossBoundary21=false;
//					else{
//						if((**itEdge).Isv2(v)&&((**itEdge).xCrossBoundary12))
//							(*itEdge)->xCrossBoundary12=false;
//					}
//				}
//			}
//			//**itEdge).Print();
//		}
//
//		//it->coord.x=fmod(it->coord.x+Lx,Lx);
//		v->coord.x+=Lx;
//	}
//
//	if(yPeriodic&&v->coord.y>=Ly)
//	{
//		std::list<Edge*>::iterator itEdge;
//		for(itEdge=v->ConnectedEdges.begin();itEdge!=v->ConnectedEdges.end();itEdge++)
//		{
//			if((**itEdge).Isv1(v)&&!((**itEdge).yCrossBoundary12||(**itEdge).yCrossBoundary21))
//				(*itEdge)->yCrossBoundary21=true;
//				else{
//					if((**itEdge).Isv2(v)&&!((**itEdge).yCrossBoundary12||(**itEdge).yCrossBoundary21))
//						(*itEdge)->yCrossBoundary12=true;
//					else{
//						if((**itEdge).Isv1(v)&&((**itEdge).yCrossBoundary12))
//						(*itEdge)->yCrossBoundary12=false;
//							else{
//								if((**itEdge).Isv2(v)&&((**itEdge).yCrossBoundary21))
//									(*itEdge)->yCrossBoundary21=false;
//								}
//						}
//					}
//		}
//		v->coord.y-=Ly;
//	}
//
//	if(yPeriodic&&(v->coord.y<0))
//	{
//		std::list<Edge*>::iterator itEdge;
//		for(itEdge=v->ConnectedEdges.begin();itEdge!=v->ConnectedEdges.end();itEdge++)
//		{
//			if((**itEdge).Isv1(v)&&!((**itEdge).yCrossBoundary12||(**itEdge).yCrossBoundary21))
//				(*itEdge)->yCrossBoundary12=true;
//			else{
//				if((**itEdge).Isv2(v)&&!((**itEdge).yCrossBoundary12||(**itEdge).yCrossBoundary21))
//					(*itEdge)->yCrossBoundary21=true;
//				else{
//					if((**itEdge).Isv1(v)&&((**itEdge).yCrossBoundary21))
//						(*itEdge)->yCrossBoundary21=false;
//					else{
//						if((**itEdge).Isv2(v)&&((**itEdge).yCrossBoundary12))
//							(*itEdge)->yCrossBoundary12=false;
//						}
//					}
//				}
//			//**itEdge).Print();
//		}
//
//		//it->coord.x=fmod(it->coord.x+Lx,Lx);
//		v->coord.y+=Ly;
//	}
//
//}
//
//void Tissue::PrintEdges()
//{
//
//	std::cout<<"Edges\n";
//	std::list<Edge>::iterator it;
//	for(it=ListEdge.begin();it!=ListEdge.end();it++)
//	{
//	(*it).Print();}
//}
//
//double Tissue::TensionAverage()
//{
//	double sumTensions=0;
//	int nbEdges=0;
//	std::list<Edge>::iterator it;
//	for(it=ListEdge.begin();it!=ListEdge.end();it++)
//	{
//	sumTensions+=it->Tension;
//	nbEdges++;
//	}
//	return(sumTensions/nbEdges);
//}
//
//double Tissue::TensionMax()
//{
//	double TensionMax=-1000;
//	std::list<Edge>::iterator it;
//	for(it=ListEdge.begin();it!=ListEdge.end();it++)
//	{
//		if((it->Tension)>TensionMax) TensionMax=it->Tension;
//	}
//	return(TensionMax);
//}
//
//void Tissue::PrintVertices()
//{
//	std::cout<<"Vertices\n";
//	std::list<Vertex>::iterator it;
//	for(it=ListVertex.begin();it!=ListVertex.end();it++)
//	{
//		it->Print();
//		it->PrintConnectedEdges();
//	}
//}
//
//void Tissue::PrintCells()
//{
//	std::cout<<"Cells\n";
//	std::list<Cell>::iterator it;
//	for(it=ListCell.begin();it!=ListCell.end();it++)
//	{
//		it->Print();
//	}
//}
//
//bool Tissue::Validite()
//{
//	std::cout<<"Test de la validite des liens :\n";
//	std::list<Cell>::iterator it;
//	bool result=true;
//	for(it=ListCell.begin();it!=ListCell.end();it++)
//	{
//		//if(it->Number==93) it->Print();
//		//if(it->Number==106) it->Print();
//		if(it->A<0) {std::cout<<"erreur : aire negative pour la cellule\n";
//			it->Print();
//		}
//		std::list<Edge*>::iterator itEdge;
//		for(itEdge=it->ListEdges.begin();itEdge!=it->ListEdges.end();itEdge++)
//		{
//			if(!((*((*itEdge)->c1)==*it)||(*((*itEdge)->c2)==*it)))
//			{result=false;
//				std::cout<<"dans la cellule "<<it->Number<<"\n l'edge suivant est mal connecte :\n";
//				(**itEdge).Print();
//			}
//			if((it->Isc2(*itEdge))&&(((*((*itEdge)->c1)==*it))))
//			{result=false;
//			std::cout<<"dans la cellule "<<it->Number<<"\n l'edge suivant est mal oriente :\n";
//			(**itEdge).Print();
//			}
//			if(!(it->Isc2(*itEdge))&&(((*((*itEdge)->c2)==*it))))
//			{result=false;
//			std::cout<<"dans la cellule "<<it->Number<<"\n l'edge suivant est mal oriente :\n";
//			(**itEdge).Print();
//			}
//		}
//	}
//
//	std::list<Edge>::iterator it2;
//	for(it2=ListEdge.begin();it2!=ListEdge.end();it2++)
//	{
//		bool c1trouve=false;
//		bool c2trouve=false;
//		std::list<Edge*>::iterator itEdge;
//		if(it2->c1->Number==0) c1trouve=true; else{
//		for(itEdge=it2->c1->ListEdges.begin();itEdge!=it2->c1->ListEdges.end();itEdge++)
//		{
//			if(!((**itEdge)==(*it2))) c1trouve=true;
//		}
//		}
//		if(!c1trouve)
//		{
//			std::cout<<"dans l'edge ";
//			it2->Print();
//			std::cout<<"la cellule c1 n'a pas l'edge correspondant\n";
//		}
//			if(it2->c2->Number==0) c2trouve=true; else{
//		for(itEdge=it2->c2->ListEdges.begin();itEdge!=it2->c2->ListEdges.end();itEdge++)
//		{
//			if(!((**itEdge)==(*it2))) c2trouve=true;
//		}
//			}
//		if(!c2trouve)
//		{
//			std::cout<<"dans l'edge ";
//			it2->Print();
//			std::cout<<"la cellule c2 n'a pas l'edge correspondant\n";
//		}
//		result=result*c1trouve*c2trouve;
//	}
//
//
//	if(result) std::cout<<"Validite Cell-Edge satisfaite\n";
//	bool result2=false;
//
//	for(it2=ListEdge.begin();it2!=ListEdge.end();it2++)
//	{
//		bool v1trouve=false;
//		bool v2trouve=false;
//		std::list<Edge*>::iterator itEdge;
//		for(itEdge=it2->v1->ConnectedEdges.begin();itEdge!=it2->v1->ConnectedEdges.end();itEdge++)
//		{
//			if(!((**itEdge)==(*it2))) v1trouve=true;
//		}
//		if(!v1trouve)
//		{
//		std::cout<<"dans l'edge ";
//		it2->Print();
//		std::cout<<"le vertex v1 est mal connecte\n";
//		}
//
//		for(itEdge=it2->v2->ConnectedEdges.begin();itEdge!=it2->v2->ConnectedEdges.end();itEdge++)
//		{
//			if(!((**itEdge)==(*it2))) v2trouve=true;
//		}
//		if(!v2trouve)
//		{
//			std::cout<<"dans l'edge ";
//			it2->Print();
//			std::cout<<"le vertex v2 est mal connecte\n";
//		}
//		result2=v1trouve*v2trouve;
//	}
//	if(result2) std::cout<<"Validite Edge-Vertex satisfaite\n";
//
//
//	return (result*result2);
//}
//
////ce type doit servir a classer les edges rencontres
//typedef struct yEdge {
//	double yinter;
//	Edge * e;
//}yEdge;
//
//
//double Tissue::VerticalTension(double x)
//{
//	std::list<yEdge> EdgeCrossed;
//	std::list<Edge>::iterator itEdge;
//	double x1;
//	double x2;
//	double y1;
//	double y2;
//	double yinter;
//	Point Center;
//	//fabrication de la liste
//
//	for(itEdge=ListEdge.begin();itEdge!=ListEdge.end();itEdge++)
//	{
//		x1=itEdge->v1->coord.x;
//		x2=itEdge->v2->coord.x;
//		y1=itEdge->v1->coord.y;
//		y2=itEdge->v2->coord.y;
//		if(itEdge->xCrossBoundary12) x2+=Lx;
//		if(itEdge->xCrossBoundary21) x1+=Lx;
//		if(itEdge->yCrossBoundary12) y2+=Ly;
//		if(itEdge->yCrossBoundary21) y1+=Ly;
//
//		if(((x1<x)&&(x2>x))||((x1>x)&&(x2<x)))
//		{
//
//			yinter=(y2-y1)/(x2-x1)*(x-x1)+y1;
//
//
//			yEdge eyEdge;
//			eyEdge.yinter=yinter;
//			eyEdge.e=(&(*itEdge));
//			EdgeCrossed.push_front(eyEdge);
//		}
//	}
//		std::list<yEdge>::iterator it;
//		//std::cout<<"Liste avant classement\n";
//	/*for(it=EdgeCrossed.begin();it!=EdgeCrossed.end();it++)
//	{
//		//std::cout<<"yinter="<<it->yinter<<"\n";
//		//it->e->Print();
//	}*/
//
//	//il faut ensuite classer cette liste (je fais un bubble sort comme ca)
//	bool swapped=true;
//
//	yEdge esouvenir;
//	std::list<yEdge>::iterator itsuivant;
//	while(swapped) {
//		swapped = false;
//		for(it=EdgeCrossed.begin();it!=EdgeCrossed.end();)
//		{
//			itsuivant=it;
//			itsuivant++;
//			if ((itsuivant!=EdgeCrossed.end())&&(it->yinter>itsuivant->yinter))
//			{
//				esouvenir=*itsuivant;
//				(*itsuivant)=*it;
//				(*it)=esouvenir;
//				swapped=true;
//			}
//			else it++;
//		}
//	}
//
//	//Normalement la la liste est classee, il faut calculer les forces
//	//std::cout<<"Liste apres classement\n";
//
//	double Tension=0;
//	for(it=EdgeCrossed.begin();it!=EdgeCrossed.end();it++)
//	{
//
//		//std::cout<<"xinter="<<it->yinter<<"\n";
//		//it->e->Print();
//		itsuivant=it;
//		itsuivant++;
//		x1=it->e->v1->coord.x;
//		x2=it->e->v2->coord.x;
//		if(it->e->xCrossBoundary12) x2+=Lx;
//		if(it->e->xCrossBoundary21) x1+=Lx;
//		//partie qui calcule la contribution de l'edge
//		Point LocalTension;
//
//			LocalTension=(it->e->Director())*(it->e->Tension);
//		if (x2<x) LocalTension=LocalTension*(-1);
//		Tension+=(LocalTension.x);
//		//std::cout<<"on vient d'ajouter pour un edge "<<LocalTension.x<<"\n";
//
//		//et maintenant il faut calculer la pression de la cellule
//		double LocalForce;
//		if(itsuivant!=EdgeCrossed.end()) {
//		double L=itsuivant->yinter-it->yinter;
//		if(x2>x) LocalForce=L*(it->e->c1->ElasticEnergy());
//		else LocalForce=L*(it->e->c2->ElasticEnergy());
//		Tension+=-LocalForce;
//		//std::cout<<"on vient d'ajouter pour une cellule la pression "<<LocalForce/L<<"\n";
//		//	std::cout<<"integree sur une longueur "<<L<<"\n";
//		//	std::cout<<"soit un total "<<-LocalForce<<"\n";
//		}
//	}
//
//	//si c'est periodique il doit rester une cellule a compter
//	if(yPeriodic)
//	{
//		if(EdgeCrossed.size()>0)
//		{
//			double LocalForce;
//			it=EdgeCrossed.begin();
//			double L=it->yinter+(Ly-EdgeCrossed.rbegin()->yinter);
//			double x2=it->e->v2->coord.x;
//			if(it->e->xCrossBoundary12) x2+=Lx;
//			if(x2>x)  LocalForce=L*(it->e->c2->ElasticEnergy());
//			else LocalForce=L*(it->e->c1->ElasticEnergy());
//			Tension+=-LocalForce;
//		}
//	}
//
//
//
//	return(Tension);//car c'est une grande valeur souvent
//
//}
//
//double Tissue::HorizontalTension(double y)
//{
//	std::list<yEdge> EdgeCrossed;
//	std::list<Edge>::iterator itEdge;
//	double x1;
//	double x2;
//	double y1;
//	double y2;
//	double xinter;
//	Point Center;
//	//fabrication de la liste
//
//	for(itEdge=ListEdge.begin();itEdge!=ListEdge.end();itEdge++)
//	{
//		x1=itEdge->v1->coord.x;
//		x2=itEdge->v2->coord.x;
//		y1=itEdge->v1->coord.y;
//		y2=itEdge->v2->coord.y;
//		if(itEdge->xCrossBoundary12) x2+=Lx;
//		if(itEdge->xCrossBoundary21) x1+=Lx;
//		if(itEdge->yCrossBoundary12) y2+=Ly;
//		if(itEdge->yCrossBoundary21) y1+=Ly;
//		if(((y1<y)&&(y2>y))||((y1>y)&&(y2<y)))
//		{
//			//comme d'hab c'est mieux de distinguer si c'est tres plat ou non
//				xinter=(x2-x1)/(y2-y1)*(y-y1)+x1;
//
//
//			yEdge eyEdge;
//			eyEdge.yinter=xinter;
//			eyEdge.e=(&(*itEdge));
//			EdgeCrossed.push_front(eyEdge);
//		}
//	}
//	std::list<yEdge>::iterator it;
//	//std::cout<<"Liste avant classement\n";
//	/*for(it=EdgeCrossed.begin();it!=EdgeCrossed.end();it++)
//	 {
//	 //std::cout<<"yinter="<<it->yinter<<"\n";
//	 //it->e->Print();
//	 }*/
//
//	//il faut ensuite classer cette liste (je fais un bubble sort comme ca)
//	bool swapped=true;
//
//	yEdge esouvenir;
//	std::list<yEdge>::iterator itsuivant;
//	while(swapped) {
//		swapped = false;
//		for(it=EdgeCrossed.begin();it!=EdgeCrossed.end();)
//		{
//			itsuivant=it;
//			itsuivant++;
//			if ((itsuivant!=EdgeCrossed.end())&&(it->yinter>itsuivant->yinter))
//			{
//				esouvenir=*itsuivant;
//				(*itsuivant)=*it;
//				(*it)=esouvenir;
//				swapped=true;
//			}
//			else it++;
//		}
//	}
//
//	//Normalement la la liste est classee, il faut calculer les forces
//
//	//std::cout<<"Liste apres classement\n";
//
//	double Tension=0;
//	for(it=EdgeCrossed.begin();it!=EdgeCrossed.end();it++)
//	{
//
//		//std::cout<<"yinter="<<it->yinter<<"\n";
//		//it->e->Print();
//		itsuivant=it;
//		itsuivant++;
//		y1=it->e->v1->coord.y;
//		y2=it->e->v2->coord.y;
//		if(it->e->yCrossBoundary12) y2+=Ly;
//		if(it->e->yCrossBoundary21) y1+=Ly;
//		//partie qui calcule la contribution de l'edge
//		Point LocalTension;
//
//			LocalTension=(it->e->Director())*(it->e->Tension);
//
//		if (y2<y) LocalTension=LocalTension*(-1);
//		Tension+=(LocalTension.y);
//		//std::cout<<"y="<<y<<", y2="<<y2<<"\n";
//		//std::cout<<"on vient d'ajouter pour un edge "<<LocalTension.y<<"\n";
//
//		//et maintenant il faut calculer la pression de la cellule
//		double LocalForce;
//		if(itsuivant!=EdgeCrossed.end()) {
//			double L=itsuivant->yinter-it->yinter;
//			if(y2>y) LocalForce=L*(it->e->c2->ElasticEnergy());
//			else LocalForce=L*(it->e->c1->ElasticEnergy());
//			Tension+=-LocalForce;
//			//std::cout<<"on vient d'ajouter pour une cellule la pression "<<LocalForce/L<<"\n";
//			//	std::cout<<"integree sur une longueur "<<L<<"\n";
//			//	std::cout<<"soit un total "<<-LocalForce<<"\n";
//		}
//	}
//	if(xPeriodic)
//	{
//		if(EdgeCrossed.size()>0)
//		{
//
//			double LocalForce;
//			it=EdgeCrossed.begin();
//			std::list<yEdge>::reverse_iterator itfin=EdgeCrossed.rbegin();
//			double xdebut=it->yinter;
//			double xfin=itfin->yinter;
//			double L;
//			if(xfin>Lx) L=xdebut-(xfin-Lx);
//			if(xfin<Lx) L=xdebut+Lx-xfin;
//			double y2=itfin->e->v2->coord.y;
//			if(itfin->e->yCrossBoundary12) y2+=Ly;
//			if(y2>y)  LocalForce=L*(itfin->e->c2->ElasticEnergy());
//			else LocalForce=L*(itfin->e->c1->ElasticEnergy());
//			Tension+=-LocalForce;
//			//std::cout<<"on vient d'ajouter pour une cellule la pression "<<LocalForce/L<<"\n";
//			//std::cout<<"integree sur une longueur "<<L<<"\n";
//			//std::cout<<"soit un total "<<-LocalForce<<"\n";
//		}
//	}
//
//	return(Tension);
//
//}
//void Tissue::actualiseTensions()
//{
//	std::list<Edge>::iterator it;
//	for(it=ListEdge.begin();it!=ListEdge.end();it++)
//	{
//		{
//			//it->Tension=Tensions[it->c1->Type][it->c2->Type]
//			//+(QuadraticTensions[it->c1->Type]
//			//+QuadraticTensions[it->c2->Type])*it->Length();
//			//autre essai avec une organisation differente des dependences
//			//it->Tension=Tensions[it->c1->Type][it->c2->Type]*it->Length()
//			//+(QuadraticTensions[it->c1->Type]
//			//+QuadraticTensions[it->c2->Type]);
//			//le terme quadratique est dans la longueur totale de la cellule
//			it->Tension=Tensions[it->c1->Type][it->c2->Type]
//			+QuadraticTensions[it->c1->Type]*it->c1->Perimeter()
//			+QuadraticTensions[it->c2->Type]*it->c2->Perimeter();
//		}
//	}
//}
//
//void Tissue::actualiseK(double K)
//{
//	std::list<Cell>::iterator it;
//	for(it=ListCell.begin();it!=ListCell.end();it++)
//	{
//		{it->K=K;
//		}
//	}
//}
//
//void Tissue::actualiseLx(double Lx)
//{
//	Edge::Lx=Lx;
//}
//
//void Tissue::actualiseLy(double Ly)
//{
//	Edge::Ly=Ly;
//}
//


void Tissue::UpdateExtPositions()
{
	// iterate through all the vertices and find the min and max values
    if(ListVertex.empty()) return;
 
    if(1)
    {
    
	std::list<Vertex>::iterator it=ListVertex.begin();
	
	Vertex* v=&(*it);
	
	extPosition_bottomLeft=v->getCoord();
	extPosition_topRight=v->getCoord();
	
	
	// check if x values of apical and basal vertex are bigger than max or smaller than min
	if(v->getBasalCoord().x<extPosition_bottomLeft.x) extPosition_bottomLeft.x=v->getBasalCoord().x;
	if(v->getBasalCoord().x>extPosition_topRight.x) extPosition_topRight.x=v->getBasalCoord().x;
	
	// check if y values are bigger than max or smaller than min
	if(v->getBasalCoord().y<extPosition_bottomLeft.y) extPosition_bottomLeft.y=v->getBasalCoord().y;
	if(v->getBasalCoord().y>extPosition_topRight.y) extPosition_topRight.y=v->getBasalCoord().y;
	
	// check if z values are bigger than max or smaller than min
	if(v->getBasalCoord().z<extPosition_bottomLeft.z) extPosition_bottomLeft.z=v->getBasalCoord().z;
	if(v->getBasalCoord().z>extPosition_topRight.z) extPosition_topRight.z=v->getBasalCoord().z;
	
	
	it++;
	
	for(;it!=ListVertex.end();it++)
	{
		v=&(*it);
		
		// check if x values of apical and basal vertex are bigger than max or smaller than min
		if(v->getCoord().x<extPosition_bottomLeft.x) extPosition_bottomLeft.x=v->getCoord().x;
		if(v->getBasalCoord().x<extPosition_bottomLeft.x) extPosition_bottomLeft.x=v->getBasalCoord().x;
		if(v->getCoord().x>extPosition_topRight.x) extPosition_topRight.x=v->getCoord().x;
		if(v->getBasalCoord().x>extPosition_topRight.x) extPosition_topRight.x=v->getBasalCoord().x;
		
		// check if y values are bigger than max or smaller than min
		if(v->getCoord().y<extPosition_bottomLeft.y) extPosition_bottomLeft.y=v->getCoord().y;
		if(v->getBasalCoord().y<extPosition_bottomLeft.y) extPosition_bottomLeft.y=v->getBasalCoord().y;
		if(v->getCoord().y>extPosition_topRight.y) extPosition_topRight.y=v->getCoord().y;
		if(v->getBasalCoord().y>extPosition_topRight.y) extPosition_topRight.y=v->getBasalCoord().y;
		
		// check if z values are bigger than max or smaller than min
		if(v->getCoord().z<extPosition_bottomLeft.z) extPosition_bottomLeft.z=v->getCoord().z;
		if(v->getBasalCoord().z<extPosition_bottomLeft.z) extPosition_bottomLeft.z=v->getBasalCoord().z;
		if(v->getCoord().z>extPosition_topRight.z) extPosition_topRight.z=v->getCoord().z;
		if(v->getBasalCoord().z>extPosition_topRight.z) extPosition_topRight.z=v->getBasalCoord().z;
		
	}
    }
    
    else
    {
        // iterate through all the vertices and find the min and max values
        if(ListVertex.empty()) return;
        
        std::list<Vertex>::iterator it=ListVertex.begin();
        
        Vertex* v=&(*it);
        
        double minz=v->coord.z;
        double maxz=v->coord.z;
        
        // check if z values are bigger than max or smaller than min
        if(v->coord_b.z<minz) minz=v->coord_b.z;
        if(v->coord_b.z>maxz) maxz=v->coord_b.z;
        
        
        it++;
        
        for(;it!=ListVertex.end();it++)
        {
            v=&(*it);
            if(v->coord_b.z<minz) minz=v->coord_b.z;
            if(v->coord_b.z>maxz) maxz=v->coord_b.z;
            if(v->coord.z<minz) minz=v->coord.z;
            if(v->coord.z>maxz) maxz=v->coord.z;
        }
        
        extPosition_topRight =  Point(SystemSize.x,SystemSize.y,maxz);
        extPosition_bottomLeft = Point(0,0,minz);
    }
	
}

// positional shot noise on all vertices, either N(0,stdDev) for distribution==1 or U(-stdDev,stdDev) for distribution==2
void Tissue::applyPositionalNoise(double stdDev, unsigned distribution)
{
    //std::cout << "Noise" << std::endl;
    
    // set seed for random number generator
    srand(time(NULL));
    
    if(distribution==1){ // normal distribution
        
        // add a random contribution to all the vertices
        for(std::list<Vertex>::iterator it_vertex = ListVertex.begin(); it_vertex!=ListVertex.end(); it_vertex++)
        {
            Vertex* v = &(*it_vertex);
            Point P1 = Point(    random_normal()  ,  random_normal() , random_normal())*stdDev;
            Point P2;
            if(!STIFF_BASAL_MEMBRANE){ // movement in z direction for the basal vertices is allowed
                P2 = Point(    random_normal()  ,  random_normal() , random_normal())*stdDev;
            }  else  { // no movement in z direction for the basal vertices due to the basal membrane
                P2 = Point(    random_normal()  ,  random_normal() , 0)*stdDev;
            }
            
            v->moveVertex(P1,P2);
        }
        
        // initialize random number generator
        //    std::default_random_engine generator;
        //    std::normal_distribution<double> distribution(0.0,stdDev);
        //
        //    // add a random contribution to all the vertices
        //    for(std::list<Vertex>::iterator it_vertex = ListVertex.begin(); it_vertex!=ListVertex.end(); it_vertex++)
        //    {
        //        Vertex* v = &(*it_vertex);
        //
        //        Point Pa = Point(distribution(generator),distribution(generator),distribution(generator));
        //        Point Pb;
        //        if(!STIFF_BASAL_MEMBRANE) { // movement in z direction of basal vertex
        //            Pb = Point(distribution(generator),distribution(generator),distribution(generator));
        //        } else { // z coordinate of basal vertex fixed
        //            Pb = Point(distribution(generator),distribution(generator),0);
        //        }
        //        
        //        v->moveVertex(P1,Pb);
        //    }
        //    
        //    Changed = true;
    } else { // uniform distribution
        // initialize random number generator with seed being the current time
        srand(time(NULL));
        
        // add a random contribution to all the vertices
        for(std::list<Vertex>::iterator it_vertex = ListVertex.begin(); it_vertex!=ListVertex.end(); it_vertex++)
        {
            Vertex* v = &(*it_vertex);
            Point P1 = Point(    ((rand()%200) - 100) / 100.  ,  ((rand()%200) - 100) / 100.,((rand()%200) - 100) / 100.)*stdDev;
            Point P2;
            if(!STIFF_BASAL_MEMBRANE){ // movement in z direction for the basal vertices is allowed
                P2 = Point(    ((rand()%200) - 100) / 100.  ,  ((rand()%200) - 100) / 100.,((rand()%200) - 100) / 100.)*stdDev;
            }  else  { // no movement in z direction for the basal vertices due to the basal membrane
                P2 = Point(    ((rand()%200) - 100) / 100.  ,  ((rand()%200) - 100) / 100.,0)*stdDev;
            }
            v->moveVertex(P1,P2);
        }
    }
    
    Changed = true;
}

//// add positional noise on all vertices, which is distributed normally around 0 with std dev stdDev
//void Tissue::applyPositionalNoise(double stdDev)
//{
//    // initialize random number generator
//    std::default_random_engine generator;
//    std::normal_distribution<double> distribution(0.0,stdDev);
//    
//    // add a random contribution to all the vertices
//    for(std::list<Vertex>::iterator it_vertex = ListVertex.begin(); it_vertex!=ListVertex.end(); it_vertex++)
//    {
//        Vertex* v = &(*it_vertex);
//        
//        Point Pa = Point(distribution(generator),distribution(generator),distribution(generator));
//        Point Pb;
//        if(!STIFF_BASAL_MEMBRANE) { // movement in z direction of basal vertex
//            Pb = Point(distribution(generator),distribution(generator),distribution(generator));
//        } else { // z coordinate of basal vertex fixed
//            Pb = Point(distribution(generator),distribution(generator),0);
//        }
//        
//        v->moveVertex(P1,Pb);
//    }
//    
//    Changed = true;
//}

void Tissue::ResetVolumeDerivatives()
{
//    for(unsigned int c=0;c<CELL_NUMBER;c++)
//    {
//        for(unsigned int v=0;v<VERTEX_NUMBER;v++)
//        {
//            VolumeDerivativesApical_analytically[c*VERTEX_NUMBER + v] = Point();
//            VolumeDerivativesBasal_analytically[c*VERTEX_NUMBER + v] = Point();
//        }
//    }
};

void Tissue::UpdateVolumeDerivatives_numerical()
{
//    double d(10e-5), d_2(d/2);
//    Point P_x_d_2(d_2,0,0),P_x_d(-d,0,0),P_y_d_2(0,d_2,0),P_y_d(0,-d,0),P_z_d_2(0,0,d_2),P_z_d(0,0,-d),P_null(0,0,0);
//    double px,mx,py,my,pz,mz;
//    Point plus, minus;
//
//    // loop through vertices
//    for(std::list<Vertex>::iterator it_v = ListVertex.begin(); it_v != ListVertex.end(); it_v++)
//    {
//        Vertex* v = &(*it_v);
//        
//        std::list<Cell*> neighbouringCells = v->getNeighbouringCells();
//        for(std::list<Cell*>::iterator it_c = neighbouringCells.begin(); it_c != neighbouringCells.end(); it_c++)
//        {
//            Cell* c = (*it_c);
//            
//            v->moveVertex(P_x_d_2,P_null);
//            UpdateBarycenters();
//            c->UpdateVolume();
//            px = c->Volume;
//            v->moveVertex(P_x_d,P_null);
//            UpdateBarycenters();
//            c->UpdateVolume();
//            mx = c->Volume;
//            v->moveVertex(P_x_d_2,P_null);
//
//            v->moveVertex(P_y_d_2,P_null);
//            UpdateBarycenters();
//            c->UpdateVolume();
//            py = c->Volume;
//            v->moveVertex(P_y_d,P_null);
//            UpdateBarycenters();
//            c->UpdateVolume();
//            my = c->Volume;
//            v->moveVertex(P_y_d_2,P_null);
//
//            v->moveVertex(P_z_d_2,P_null);
//            UpdateBarycenters();
//            c->UpdateVolume();
//            pz = c->Volume;
//            v->moveVertex(P_z_d,P_null);
//            UpdateBarycenters();
//            c->UpdateVolume();
//            mz = c->Volume;
//            v->moveVertex(P_z_d_2,P_null);
//            
//            VolumeDerivativesApical_numerically[c->Number*VERTEX_NUMBER + v->Number]= Point((px-mx)/d,(py-my)/d,(pz-mz)/d);
//            
//            v->moveVertex(P_null,P_x_d_2);
//            UpdateBarycenters();
//            c->UpdateVolume();
//            px = c->Volume;
//            v->moveVertex(P_null,P_x_d);
//            UpdateBarycenters();
//            c->UpdateVolume();
//            mx = c->Volume;
//            v->moveVertex(P_null,P_x_d_2);
//            
//            v->moveVertex(P_null,P_y_d_2);
//            UpdateBarycenters();
//            c->UpdateVolume();
//            py = c->Volume;
//            v->moveVertex(P_null,P_y_d);
//            UpdateBarycenters();
//            c->UpdateVolume();
//            my = c->Volume;
//            v->moveVertex(P_null,P_y_d_2);
//            
//            v->moveVertex(P_null,P_z_d_2);
//            UpdateBarycenters();
//            c->UpdateVolume();
//            pz = c->Volume;
//            v->moveVertex(P_null,P_z_d);
//            UpdateBarycenters();
//            c->UpdateVolume();
//            mz = c->Volume;
//            v->moveVertex(P_null,P_z_d_2);
//            UpdateBarycenters();
//            
//            VolumeDerivativesBasal_numerically[c->Number*VERTEX_NUMBER + v->Number]= Point((px-mx)/d,(py-my)/d,(pz-mz)/d);
//        }
//    }
    
//    for(unsigned cell = 0; cell<CELL_NUMBER; cell++)
//    {
//        std::cout << "cell: " << cell << std::endl;
//        for(unsigned vertex = 0; vertex<VERTEX_NUMBER; vertex++)
//        {
//            //VolumeDerivativesApical_numerically[cell*72*sizeof(Point()) + vertex*sizeof(Point())].Print();
//            if(Point::Norm(VolumeDerivativesApical_numerically[cell*VERTEX_NUMBER + vertex])>0)
//            {
//                std::cout << "      vertex: " << vertex << std::endl << "            ";
//                VolumeDerivativesApical_numerically[cell*VERTEX_NUMBER + vertex].Print();
//            }
//        }
//    }
    
};


Cell* Tissue::CellFromNumber(int no)
{
    for(std::list<Cell>::iterator it_c = ListCell.begin(); it_c != ListCell.end(); it_c++)
    {
        if(it_c->Number == no) return &(*it_c);
    }
    
    return NULL;
}

void Tissue::slot_addCurrentPoint(double tx, double ty)
{
    //if(mainWindow_tab==1)
    {
        addCurrentCell(tx, ty);
    }
    //else if(mainWindow_tab==2)
    {
        addCurrentEdge(tx, ty);
    }
    //else if(mainWindow_tab==3)
    {
        addCurrentVertex(tx, ty);
    }
    
    updateMainWindow();
    
}

void Tissue::slot_currentPoint(double tx, double ty)
{
    //if(mainWindow_tab==2)
    {
        currentCell(tx, ty);
    }
    //else if(mainWindow_tab==3)
    {
        currentEdge(tx, ty);
    }
    //else if(mainWindow_tab==4)
    {
        currentVertex(tx, ty);
    }
    
    updateMainWindow();
    //
    //    switch(mainWindow_tab)
    //    {
    //        case 1:
    //            currentCell(tx, ty);
    //        case 2:
    //            currentEdge(tx, ty);
    //        case 3:
    //            currentVertex(tx, ty);
    //    }
    
}

std::list<Cell*> Tissue::pointerCellList()
{
    std::list<Cell*> p_cellList;
    for(std::list<Cell>::iterator it_c = ListCell.begin(); it_c != ListCell.end(); it_c++) p_cellList.push_back(&(*it_c));
    
    return p_cellList;
}

void Tissue::addCurrentCell(double tx, double ty)
{
    Cell *c = closestCell(Point(tx,ty,0.0), pointerCellList());
    currentCells.push_front(c);
    // calculate aspect ratio
    double Aa = c->CalculateApicalSurfaceArea();
    double h = Point::Norm(c->BC_a-c->BC_b);
    double AR = h/std::sqrt(Aa);
    
    signal_currentCell(c->Number,c->V0,c->Volume,c->Type,c->K,c->Pressure,c->T_a(),c->T_b(), Aa, AR);
}

void Tissue::addCurrentEdge(double tx, double ty)
{
    Edge *e = closestEdge(Point(tx,ty,0.0));
    currentEdges.push_front(e);
    if(e->c1 && e->c2)
    {
        signal_currentEdge(e->Number,e->c1->Number,e->c2->Number,e->T_l,e->G_a,e->G_b,e->l_a(), e->l_b());
    }
    else if(!e->c1)
    {
        signal_currentEdge(e->Number, 0, e->c2->Number, e->T_l, e->G_a, e->G_b, e->l_a(), e->l_b());
    }
    else
    {
        signal_currentEdge(e->Number, e->c1->Number, 0, e->T_l, e->G_a, e->G_b, e->l_a(), e->l_b());
    }
}

void Tissue::addCurrentVertex(double tx, double ty)
{
    Vertex *v = closestVertex(Point(tx,ty,0.0));
    currentVertices.push_front(v);
    signal_currentVertex(v->Number, v->coord, v->coord_b,v->dW_dva,v->dW_dvb,v->G_l);

    //v->updateNeighbouringEdges();
    //std::cout << "Surrounding Energy = " << v->SurroundingEnergy() << std::endl;

}

void Tissue::currentCell(double tx, double ty)
{
    Cell *c = closestCell(Point(tx,ty,0.0),pointerCellList());
    currentCells = std::list<Cell*>();
    currentCells.push_front(c);
    //std::cout << std::setprecision(20) << c->Volume << " " << c->Pressure << "  " << c->V0 << "  " << c->K*(c->Volume-c->V0) << std::endl;
    // calculate aspect ratio
    double Aa = c->CalculateApicalSurfaceArea();
    double h = Point::Norm(c->BC_a-c->BC_b);
    double AR = h/std::sqrt(Aa);
    
    signal_currentCell(c->Number,c->V0,c->Volume,c->Type,c->K,c->Pressure,c->T_a(),c->T_b(), Aa, AR);
}

void Tissue::currentEdge(double tx, double ty)
{
    Edge *e = closestEdge(Point(tx,ty,0.0));
    currentEdges = std::list<Edge*>();
    currentEdges.push_front(e);

    if(e->c1 && e->c2)
    {
        signal_currentEdge(e->Number, e->c1->Number, e->c2->Number, e->T_l, e->G_a, e->G_b, e->la, e->lb);
    }
    else if(!e->c1 && e->c2)
    {
        signal_currentEdge(e->Number, 0, e->c2->Number, e->T_l, e->G_a, e->G_b, e->la, e->lb);
    }
    else if(e->c1 && !e->c2)
    {
        signal_currentEdge(e->Number, e->c1->Number, 0, e->T_l, e->G_a, e->G_b, e->la, e->lb);
    }
}

void Tissue::currentVertex(double tx, double ty)
{
    Vertex *v = closestVertex(Point(tx,ty,0.0));
    currentVertices = std::list<Vertex*>();
    currentVertices.push_front(v);
    
    // v->neighbouringEdges();
    
    signal_currentVertex(v->Number, v->coord, v->coord_b,v->dW_dva,v->dW_dva,v->G_l);
    
    //std::cout << "Number: " << v->Number << std::endl;
    //v->dW_dva_numerical.Print();
    //v->dW_dva.Print();
    
    //v->updateNeighbouringEdges();
    //std::cout << "Vertex Surrounding Energy = " << v->SurroundingEnergy() << std::endl;

}

void Tissue::slot_setCellsType(int c_type)
{
    for(std::list<Cell*>::iterator it_c = currentCells.begin(); it_c!=currentCells.end(); ++it_c)
    {
        (*it_c)->Type=c_type;
    }
    
    setMechanicalProperties();
}


void Tissue::updateMainWindow()
{
    if(!currentCells.empty())
    {
        Cell *c = currentCells.front();
        
        // calculate aspect ratio
        double Aa = c->CalculateApicalSurfaceArea();
        double h = Point::Norm(c->BC_a-c->BC_b);
        double AR = h/std::sqrt(Aa);
        
        signal_currentCell(c->Number,c->V0,c->Volume,c->Type,c->K,c->Pressure,c->T_a(),c->T_b(), Aa, AR);

    }
    else
    {
        //signal_currentCell(0,0,0,0,0,0,0,0);
    }
    
    if(!currentEdges.empty())
    {
        Edge *e = currentEdges.front();
        if(e->c1 && e->c2)
        {
            signal_currentEdge(e->Number,e->c1->Number,e->c2->Number,e->T_l,e->G_a,e->G_b,e->l_a(),e->l_b());
        }
        else if(e->c1 && !e->c2)
        {
            signal_currentEdge(e->Number,e->c1->Number,0,e->T_l,e->G_a,e->G_b,e->l_a(),e->l_b());
        }
        if(!e->c1 && e->c2)
        {
            signal_currentEdge(e->Number,0,e->c2->Number,e->T_l,e->G_a,e->G_b,e->l_a(),e->l_b());
        }
    }
    else
    {
        //signal_currentEdge(0,0,0,0,0,0,0,0);
    }
    
    if(!currentVertices.empty())
    {
        Vertex *v = currentVertices.front();
        signal_currentVertex(v->Number, v->coord, v->coord_b,v->dW_dva,v->dW_dvb,v->G_l);
    }
    else
    {
        Point Zero = Point();
        signal_currentVertex(0,Zero,Zero,Zero,Zero,0);
    }
    
    signal_currentTissue(0, SystemSize.x, SystemSize.y, dW_dL.x, dW_dL.y, T_ext, Lx_fix, Ly_fix,  Energy);
    signal_currentLumen_a(V_Lumen_a);
    signal_currentLumen_b(V_Lumen_b);
}

void Tissue::update_pos(bool including_L)
{
    if(!including_L || !isPeriodic || (Lx_fix && Ly_fix)) {
        // initialize the vector
        unsigned N = ListVertex.size();
        unsigned N_times_3 = 3*N;
        pos = std::vector<double>(6*N);
        unsigned i=0;
        
        // fill in the coordinates for every vertex
        for(std::list<Vertex>::iterator it_v = ListVertex.begin(); it_v!=ListVertex.end(); it_v++) {
            pos[i]=it_v->coord.x;
            pos[i+N_times_3]=it_v->coord_b.x;
            pos[i+1]=it_v->coord.y;
            pos[i+N_times_3+1]=it_v->coord_b.y;
            pos[i+2]=it_v->coord.z;
            pos[i+N_times_3+2]=it_v->coord_b.z;
            i+=3;
        }
    } else {
        // initialize the vector
        unsigned N = ListVertex.size();
        unsigned N_times_3 = 3*N;
        pos = std::vector<double>(6*N+2);
        unsigned i=0;
        
        // fill in the coordinates for every vertex
        for(std::list<Vertex>::iterator it_v = ListVertex.begin(); it_v!=ListVertex.end(); it_v++) {
            pos[i]=it_v->coord.x;
            pos[i+N_times_3]=it_v->coord_b.x;
            pos[i+1]=it_v->coord.y;
            pos[i+N_times_3+1]=it_v->coord_b.y;
            pos[i+2]=it_v->coord.z;
            pos[i+N_times_3+2]=it_v->coord_b.z;
            i+=3;
        }
        
        // set derivatives wrt system size
        pos[6*N] = SystemSize.x;
        pos[6*N+1] = SystemSize.y;
    }

};

void Tissue::update_der(bool including_L)
{
    if(!including_L || isPeriodic || (Lx_fix && Ly_fix)) {
        // initialize the vector
        unsigned N = ListVertex.size();
        unsigned N_times_3 = 3*N;
        der = std::vector<double>(6*N);
        unsigned i=0;
        
        // fill in the coordinates for every vertex
        for(std::list<Vertex>::iterator it_v = ListVertex.begin(); it_v!=ListVertex.end(); it_v++) {
            der[i]=it_v->dW_dva.x;
            der[i+N_times_3]=it_v->dW_dvb.x;
            der[i+1]=it_v->dW_dva.y;
            der[i+N_times_3+1]=it_v->dW_dvb.y;
            der[i+2]=it_v->dW_dva.z;
            der[i+N_times_3+2]=it_v->dW_dvb.z;
            
            i+=3;
        }
    } else {
        // initialize the vector
        unsigned N = ListVertex.size();
        unsigned N_times_3 = 3*N+2;
        der = std::vector<double>(6*N+2);
        unsigned i=0;
        
        // fill in the coordinates for every vertex
        for(std::list<Vertex>::iterator it_v = ListVertex.begin(); it_v!=ListVertex.end(); it_v++) {
            der[i]=it_v->dW_dva.x;
            der[i+N_times_3]=it_v->dW_dvb.x;
            der[i+1]=it_v->dW_dva.y;
            der[i+N_times_3+1]=it_v->dW_dvb.y;
            der[i+2]=it_v->dW_dva.z;
            der[i+N_times_3+2]=it_v->dW_dvb.z;
            
            i+=3;
        }
        
        // set derivatives wrt system size
        der[6*N] = dW_dL.x;
        der[6*N+1] = dW_dL.y;
    }
};

void Tissue::set_pos(bool including_L)
{
    if(!including_L || !isPeriodic || (Lx_fix && Ly_fix)) {
        // check if the dimensions are compatible
        if(pos.size()!=6*ListVertex.size()) throw;
        const unsigned N_times_3 = 3*ListVertex.size();
        
        unsigned i=0;
        
        // change the coordinates for every vertex
        for(std::list<Vertex>::iterator it_v = ListVertex.begin(); it_v!=ListVertex.end(); it_v++)
        {
            //if(isPeriodic && (pos[i]>SystemSize.x || pos[i+1]>SystemSize.y || pos[i+N_times_3]>SystemSize.x || pos[i+N_times_3+1]>SystemSize.y || pos[i]<0 || pos[i+1]<0 || pos[i+N_times_3]<0 || pos[i+N_times_3+1]<0) ){ // if the vertex crosses the boundary
                Point apicalMovement(  pos[i]-it_v->coord.x  ,  pos[i+1]-it_v->coord.y  ,  pos[i+2]-it_v->coord.z  );
                Point basalMovement(  pos[i+N_times_3]-it_v->coord_b.x  ,  pos[i+1+N_times_3]-it_v->coord_b.y  ,  pos[i+2+N_times_3]-it_v->coord_b.z  );
            
                it_v->moveVertex(apicalMovement, basalMovement);
                
                pos[i]=it_v->coord.x;
                pos[i+1]=it_v->coord.y;
                pos[i+2]=it_v->coord.z;
                
                pos[i+N_times_3]=it_v->coord_b.x;
                pos[i+1+N_times_3]=it_v->coord_b.y;
                pos[i+2+N_times_3]=it_v->coord_b.z;
                
                i+=3;
//
//            } else {
//            
////                it_v->coord.x   = pos[i];
////                it_v->coord_b.x = pos[i+N_times_3];
////                it_v->coord.y   = pos[i+1];
////                it_v->coord_b.y = pos[i+N_times_3+1];
////                
////                if(apical_wall && pos[i+2]>apical_wall_position){ // if we want to exclude the ECM as a basal wall at apical_wall_position
////                    pos[i+2] = apical_wall_position;
////                    der[i+2] = 0;
////                }
////                
////                if(basal_wall && pos[i+N_times_3+2]<basal_wall_position){ // if we want to exclude the ECM as a basal wall at basal_wall_position
////                    pos[i+N_times_3+2] = basal_wall_position;
////                    der[i+N_times_3+2] = 0;
////                }
////                
////                it_v->coord.z   = pos[i+2];
////                it_v->coord_b.z = pos[i+N_times_3+2];
//                
//                i+=3;
//            }
        }
    } else {
        // check if the dimensions are compatible
        if(pos.size()!=6*ListVertex.size()+2) throw;
        const unsigned N_times_3 = 3*ListVertex.size();
        
        unsigned i=0;
        
        // prefactors of the new coordinates due to the change in system size
        double scaling_x = pos[6*ListVertex.size()] / SystemSize.x;
        double scaling_y = pos[6*ListVertex.size()+1] / SystemSize.y;
        
        
        // change the coordinates for every vertex
        for(std::list<Vertex>::iterator it_v = ListVertex.begin(); it_v!=ListVertex.end(); it_v++)
        {
            if(isPeriodic && (pos[i]>SystemSize.x || pos[i+1]>SystemSize.y || pos[i+N_times_3]>SystemSize.x || pos[i+N_times_3+1]>SystemSize.y || pos[i]<0 || pos[i+1]<0 || pos[i+N_times_3]<0 || pos[i+N_times_3+1]<0) ){ // if the vertex crosses the boundary
                
                // calculate the new positions of the vertices and take the periodicity into account
                Point apicalMovement(  pos[i]-it_v->coord.x  ,  pos[i+1]-it_v->coord.y  ,  pos[i+2]-it_v->coord.z  );
                Point basalMovement(  pos[i+N_times_3]-it_v->coord_b.x  ,  pos[i+1+N_times_3]-it_v->coord_b.y  ,  pos[i+2+N_times_3]-it_v->coord_b.z  );
                
                it_v->moveVertex(apicalMovement, basalMovement);

                // rescale the coordinates
                it_v->coord.x   *= scaling_x;
                it_v->coord_b.x *= scaling_x;
                it_v->coord.y   *= scaling_y;
                it_v->coord_b.y *= scaling_y;
                
                pos[i]   = it_v->coord.x;
                pos[i+1] = it_v->coord.y;
                pos[i+2] = it_v->coord.z;
                
                pos[i+N_times_3]   = it_v->coord_b.x;
                pos[i+1+N_times_3] = it_v->coord_b.y;
                pos[i+2+N_times_3] = it_v->coord_b.z;
                
                i+=3;
            }
            else {
                
                pos[i+N_times_3]   *= scaling_x;
                pos[i]             *= scaling_x;
                pos[i+N_times_3+1] *= scaling_y;
                pos[i+1]           *= scaling_y;
                
                it_v->coord.x   = pos[i];
                it_v->coord_b.x = pos[i+N_times_3];
                it_v->coord.y   = pos[i+1];
                it_v->coord_b.y = pos[i+N_times_3+1];

                if(apical_wall && pos[i+2]>apical_wall_position){ // if we want to exclude the ECM as a basal wall at basal_wall_position
                    pos[i+2] = apical_wall_position;
                    der[i+2] = 0;
                }
                
                if(basal_wall && pos[i+N_times_3+2]<basal_wall_position){ // if we want to include the ECM as a basal wall at basal_wall_position
                    pos[i+N_times_3+2] = basal_wall_position;
                    der[i+N_times_3+2] = 0;
                }
                
                it_v->coord.z   = pos[i+2];
                it_v->coord_b.z = pos[i+N_times_3+2];

                i+=3;
            }
        }
        
        // set system size
        SystemSize.x = pos[6*ListVertex.size()];
        SystemSize.y = pos[6*ListVertex.size()+1];
    }
};

void Tissue::update_pos_and_der(bool including_L)
{
    if(!including_L || !isPeriodic || (Lx_fix && Ly_fix)) {
        // initialize the vector
        unsigned N = ListVertex.size();
        unsigned N_times_3 = 3*N;
        pos = std::vector<double>(6*N);
        der = std::vector<double>(6*N);

        unsigned i=0;
        
        // fill in the coordinates and derivatives for every vertex
        for(std::list<Vertex>::iterator it_v = ListVertex.begin(); it_v!=ListVertex.end(); it_v++)
        {
            der[i]=it_v->dW_dva.x;
            pos[i]=it_v->coord.x;
            der[i+N_times_3]=it_v->dW_dvb.x;
            pos[i+N_times_3]=it_v->coord_b.x;
            
            der[i+1]=it_v->dW_dva.y;
            pos[i+1]=it_v->coord.y;
            der[i+N_times_3+1]=it_v->dW_dvb.y;
            pos[i+N_times_3+1]=it_v->coord_b.y;
            
            der[i+2]=it_v->dW_dva.z;
            pos[i+2]=it_v->coord.z;
            der[i+N_times_3+2]=it_v->dW_dvb.z;
            pos[i+N_times_3+2]=it_v->coord_b.z;
            
            i+=3;
        }
    } else {
        // initialize the vector
        unsigned N = ListVertex.size();
        unsigned N_times_3 = 3*N;
        pos = std::vector<double>(6*N+2);
        der = std::vector<double>(6*N+2);
        
        unsigned i=0;
        
        // fill in the coordinates and derivatives for every vertex
        for(std::list<Vertex>::iterator it_v = ListVertex.begin(); it_v!=ListVertex.end(); it_v++) {
            der[i]=it_v->dW_dva.x;
            pos[i]=it_v->coord.x;
            der[i+N_times_3]=it_v->dW_dvb.x;
            pos[i+N_times_3]=it_v->coord_b.x;
            
            der[i+1]=it_v->dW_dva.y;
            pos[i+1]=it_v->coord.y;
            der[i+N_times_3+1]=it_v->dW_dvb.y;
            pos[i+N_times_3+1]=it_v->coord_b.y;
            
            der[i+2]=it_v->dW_dva.z;
            pos[i+2]=it_v->coord.z;
            der[i+N_times_3+2]=it_v->dW_dvb.z;
            pos[i+N_times_3+2]=it_v->coord_b.z;
            
            i+=3;
        }
        
        // set system size and derivatives wrt system size
        pos[6*N]   = SystemSize.x;
        pos[6*N+1] = SystemSize.y;
        der[6*N]   = dW_dL.x;
        der[6*N+1] = dW_dL.y;
    }
};

void Tissue::change_size(double new_Lx, double new_Ly)
{
    if(new_Lx==SystemSize.x && new_Ly==SystemSize.y) return; // nothing to be done
    
    double factor_x = new_Lx/SystemSize.x;
    double factor_y = new_Ly/SystemSize.y;
    
    for(std::list<Vertex>::iterator it_v=ListVertex.begin(); it_v!=ListVertex.end(); it_v++)
    {
        it_v->coord.x *= factor_x;
        it_v->coord.y *= factor_y;
    }
    
    SystemSize.x=new_Lx;
    SystemSize.y=new_Ly;
}

Cell* Tissue::get_random_cell()
{
    unsigned numberOfCells = ListCell.size();
    
    // draw random numbers to determine two edges of the cell that will be splitted
    int rand_cell = rand() % numberOfCells;
    
    int counter = 0;
    // find the edges that correspond to those numbers
    for(std::list<Cell>::iterator iter_c = ListCell.begin(); iter_c!=ListCell.end(); ++iter_c)
    {
        if(counter==rand_cell) return &(*iter_c);
        counter++;
    }
    
    std::cout << "Tissue::get_random_cell() failed" << std::endl;
    throw;
}

// returns the biggest cell, that does not cross the boundary
Cell* Tissue::get_biggest_cell()
{
    //update();
    
    double maxVolume = 0;
    
    Cell* c=NULL;

    // find the edges that correspond to those numbers
    for(std::list<Cell>::iterator iter_c = ListCell.begin(); iter_c!=ListCell.end(); ++iter_c)
    {
        iter_c->UpdateVolume();
        std::cout << "V = " << iter_c->Volume << std::endl;
        if(iter_c->Volume>maxVolume) {c = &(*iter_c); maxVolume=c->Volume;std::cout<<"V_max= " << c->Volume << ", ";};
    }
    
    std::cout << std::endl;
    
    return c;
}

void Tissue::mutate_cells_in_middle(int number_of_closest_Cells)
{
    if(ListCell.size()==0) return; // nothing to be done if the tissue is not loaded!
    
    if(number_of_closest_Cells>ListCell.size()) number_of_closest_Cells = ListCell.size();
    
    std::list<Cell*> listCell = pointerCellList();
    
    // Midpoint of the tissue
    Point M = Point(SystemSize.x/2, SystemSize.y/2, 0);

    // closest vertex to midpoint is set to be reference and save the number to simplify calculation of clone properties
    Vertex* v1 = closestVertex(M);
    midVertex_number = v1->Number;
    M = v1->coord;
    
    // the closest number_of_closest_Cells of cells to v1 is set to have type 2
    for(int i=0; i < number_of_closest_Cells; i++) {
        Cell* c = closestCell(M, listCell);
        c->Type=2;
        listCell.remove(c);
    }
    
    // all other cells are type 1
    for(std::list<Cell*>::iterator it_c = listCell.begin(); it_c != listCell.end(); it_c++)  (*it_c)->Type=1;
    
    // now set the mechanical properties
    setMechanicalProperties();
    
    Cell* midCell = cell_in_middle(2);
    std::cout << "The cell in the middle of the clone has number " << midCell->Number << std::endl;
}

void Tissue::mutate_cells_in_circle(int number_of_closest_Cells)
{
    if(ListCell.size()==0) return; // nothing to be done if the tissue is not loaded!
    
    if(number_of_closest_Cells>ListCell.size()) number_of_closest_Cells = ListCell.size();
    
    std::list<Cell*> listCell = pointerCellList();
    
    // Midpoint of the tissue
    // choose midpoint of first cell as reference, if not choose closest cell to midpoint of tissue
    Point M;
    if(!currentCells.empty()){
        M = (*currentCells.begin())->BC_a;
    } else {
        M = closestCell(Point(SystemSize.x/2, SystemSize.y/2, 0),listCell)->BC_a;
    }
    
    // the number_of_closest_Cells of cells is set to have type 2
    for(int i=0; i < number_of_closest_Cells; i++) {
        Cell* c = closestCell(M, listCell);
        c->Type=2;
        listCell.remove(c);
    }
    
    // all other cells are type 1
    for(std::list<Cell*>::iterator it_c = listCell.begin(); it_c != listCell.end(); it_c++)  (*it_c)->Type=1;
    
    // now set the mechanical properties
    setMechanicalProperties();

    if(number_of_closest_Cells>0){
        Cell* midCell = cell_in_middle(2);
        std::cout << "The cell in the middle of the clone has number " << midCell->Number << std::endl;
    }
}


void Tissue::mutate_cells_in_radius(double radius)
{
    if(ListCell.size()==0) return; // nothing to be done if the tissue is not loaded!

    std::list<Cell*> listCell = pointerCellList();

    // Midpoint
    Point M = Point(SystemSize.x/2,SystemSize.y/2,0);
    Cell* c = closestCell(M, listCell);
    
    M = Point(c->BC_a.x,c->BC_a.y,0);
    
    for(std::list<Cell>::iterator it_c = ListCell.begin(); it_c != ListCell.end(); it_c++) { // loop through all cells
        Point M_c = Point(it_c->BC_a.x,it_c->BC_a.y,0); // project cell's apical BC on xy plane
        if(Point::Norm(M-M_c)<radius) it_c->Type=2;
        else it_c->Type=1;
    }
    
    // now set the mechanical properties
    setMechanicalProperties();
}

void Tissue::mutate_cells_in_annulus(double r1, double r2)
{
    if(ListCell.size()==0) return; // nothing to be done if the tissue is not loaded!
    
    std::list<Cell*> listCell = pointerCellList();
    
    // Midpoint
    Point M = Point(SystemSize.x/2,SystemSize.y/2,0);
    Cell* c = closestCell(M, listCell);
    
    M = Point(c->BC_a.x,c->BC_a.y,0);
    
    for(std::list<Cell>::iterator it_c = ListCell.begin(); it_c != ListCell.end(); it_c++) { // loop through all cells
        Point M_c = Point(it_c->BC_a.x,it_c->BC_a.y,0); // project cell's apical BC on xy plane
        double d = Point::Norm(M-M_c);
        if(d>r1 && d<r2) it_c->Type=2;
        else if(d<r1) it_c->Type=3;
        else it_c->Type=1;
    }
    
    // update the mechanical properties
    setMechanicalProperties();
}


void Tissue::mutate_cells_in_stripe(double distance)
{
    if(ListCell.size()==0) return; // nothing to be done if the tissue is not loaded!

    for(std::list<Cell>::iterator it_c = ListCell.begin(); it_c != ListCell.end(); it_c++) { // loop through all cells
        if(std::abs(SystemSize.x/2-it_c->BC_a.x)<distance) it_c->Type=2;
        else it_c->Type=1;
    }

    // set mechanical properties
    setMechanicalProperties();
}

void Tissue::mutate_cells_in_two_stripes(double width)
{
    if(ListCell.size()==0) return; // nothing to be done if the tissue is not loaded!

    for(std::list<Cell>::iterator it_c = ListCell.begin(); it_c != ListCell.end(); it_c++) { // loop through all cells
        if(std::abs(SystemSize.x/4.-it_c->BC_a.x)<width) it_c->Type=2;
        else if(std::abs(3*SystemSize.x/4.-it_c->BC_a.x)<width) it_c->Type=3;
        else it_c->Type=1;
    }
    
    // set mechanical properties
    setMechanicalProperties();
}

// routine to rectify possibly false orientations of the edges towards the cells in order to find a mistake. dont use in regular code!!
void Tissue::correct_orientations()
{
    for(std::list<Edge>::iterator it_e = ListEdge.begin(); it_e != ListEdge.end(); it_e++) {
        
        Point BC_a = it_e->c1->BC_a;
        
        bool test = BC_a.dot((BC_a-it_e->v1->coord) * (BC_a-it_e->v2->coord))>0;
        
        if(test) {
            Cell *c2 = it_e->c1;
            it_e->c1 = it_e->c2;
            it_e->c2 = c2;
        }
    }
}

// function to create a normally distributed number from a uniform distribution
double Tissue::random_normal()
{
	double U1,U2,V1,V2;
	double S=2;
	while(S>=1)
	{
		U1=((double)rand()/(double)RAND_MAX);
		U2=((double)rand()/(double)RAND_MAX);
		V1=2.0*U1-1.0;
		V2=2.0*U2-1.0;
		S=pow(V1,2)+pow(V2,2);
	}
	double X1=V1*sqrt(-2.0*log(S)/S);
    //std::cout << X1 << "  ";
	return X1;
}

int Tissue::maxVertexNumber()
{
    if(ListVertex.empty()) return 1;
    
    int maxNumber = ListVertex.front().Number;
    for(std::list<Vertex>::iterator itv = ListVertex.begin(); itv != ListVertex.end(); itv++)
    {
        if(itv->Number>maxNumber) maxNumber = itv->Number;
    }
    
    return maxNumber;
}

int Tissue::maxEdgeNumber()
{
    int maxNumber = ListEdge.front().Number;
    for(std::list<Edge>::iterator itv = ListEdge.begin(); itv != ListEdge.end(); itv++)
    {
        if(itv->Number>maxNumber) maxNumber = itv->Number;
    }
    
    return maxNumber;
}

int Tissue::maxCellNumber()
{
    int maxNumber = ListCell.front().Number;
    for(std::list<Cell>::iterator itv = ListCell.begin(); itv != ListCell.end(); itv++)
    {
        if(itv->Number>maxNumber) maxNumber = itv->Number;
    }
    
    return maxNumber;
}


void Tissue::slot_setFixLx(int x)
{
    Lx_fix = (x==2);
    //std::cout<<"Fix Lx = " << x << std::endl;
}
void Tissue::slot_setFixLy(int y)
{
    Ly_fix = (y==2);
}

bool Tissue::isFlat()
{
    std::vector<double> z_a, z_b;
    
    for(std::list<Cell>::iterator it_cell = ListCell.begin();it_cell != ListCell.end();it_cell++) {
        if(it_cell->Type==1) {
            z_a.push_back(it_cell->BC_a.z);
            z_b.push_back(it_cell->BC_b.z);
        }
    }
    
    double max_z_a, max_z_b, min_z_a, min_z_b, mean_z_a, mean_z_b, var_z_a, var_z_b;
    
    vector_statistics(z_a, mean_z_a, var_z_a, max_z_a, min_z_a);
    vector_statistics(z_b, mean_z_b, var_z_b, max_z_b, min_z_b);
    
    double maxRelDeviation = 0.01;
    
    double relDeviation = std::max((max_z_a-min_z_a), (max_z_b-min_z_b));
    
    //std::cout << "error = " << relDeviation << std::endl;
    
    return relDeviation < maxRelDeviation;
}

// calculate the volume of the volume enclosed by all the cells' apical sides in the tissue
void Tissue::update_V_Lumen_a()
{
    
    // strategy: iterate through all edges and add volume contributions from the two apical triangles
    V_Lumen_a = 0;
    
    for(std::list<Edge>::iterator ite = ListEdge.begin(); ite != ListEdge.end(); ite++)
    {
        if(ite->c1) V_Lumen_a += Point::TetrahedronSignedVolume(ite->v1->coord, ite->v2->coord, ite->c1->BC_a);
        if(ite->c2) V_Lumen_a += Point::TetrahedronSignedVolume(ite->v2->coord, ite->v1->coord, ite->c2->BC_a);
    }
    
    return;
    
    // strategy: calculate the volume enclosed by the apical side and add the volume of all the cells
    update_V_Lumen_b();
    
    V_Lumen_a = V_Lumen_b;
    
    for(std::list<Cell>::iterator it_cell = ListCell.begin(); it_cell!=ListCell.end(); it_cell++)
    {
        it_cell->UpdateVolume();
        V_Lumen_a += it_cell->Volume;
    }
}


// calculate the volume of the volume enclosed by all the cells basal sides in the tissue
void Tissue::update_V_Lumen_b()
{
    // strategy: iterate through all edges and add volume contributions from the two apical triangles
    V_Lumen_a = 0;
    
    for(std::list<Edge>::iterator ite = ListEdge.begin(); ite != ListEdge.end(); ite++)
    {
        if(ite->c1) V_Lumen_a += Point::TetrahedronSignedVolume(ite->v2->coord_b, ite->v1->coord_b, ite->c1->BC_b);
        if(ite->c2) V_Lumen_a += Point::TetrahedronSignedVolume(ite->v1->coord_b, ite->v2->coord_b, ite->c2->BC_b);
    }
    
    return;

    // iterate through all cells and add the contribution from the basal tetrahedrons
    for(std::list<Cell>::iterator it_cell = ListCell.begin(); it_cell != ListCell.end(); it_cell++)
    {
        for(std::list<Edge*>::iterator it=it_cell->ListEdges.begin();it!=it_cell->ListEdges.end();it++)	{	// iterate through all the edges of the cell
            Edge *e = (*it);
			
			// every edge yields for a basal contribution
			V_Lumen_b-=Point::TetrahedronSignedVolume(e->v2->coord_b,e->v1->coord_b,it_cell->BC_b)*it_cell->DirectionOfEdge(*it);
        }
        
    }

}

void Tissue::update_P_Lumen_a()
{
    update_V_Lumen_a();
    
    P_Lumen_a = K_Lumen_a*(V0_Lumen_a-V_Lumen_a) + P0_Lumen_a;
}

void Tissue::update_P_Lumen_b()
{
    update_V_Lumen_b();
    
    P_Lumen_b = K_Lumen_b*(V0_Lumen_b-V_Lumen_b) + P0_Lumen_b;
}

void Tissue::setBasalTensionGradient(double basalTensionLeft, double basalTensionRight)
{
    // update the external positions in order to obtain outmost points in x
    UpdateExtPositions();
    double x_min = extPosition_bottomLeft.x;
    double x_max = extPosition_topRight.x;
    
    // the lateral tension in the cells will be given by a linear function inbetween the two outmost points
    for(std::list<Cell>::iterator it_cell = ListCell.begin(); it_cell != ListCell.end(); it_cell++)
    {
        double x_curr = it_cell->BC_a.x;
        
        it_cell->Tb = basalTensionLeft + (x_curr - x_min)/(x_max-x_min)*(basalTensionRight-basalTensionLeft);
        
    }
    
}

void Tissue::randomCellDivisions(unsigned N)
{
    for(unsigned i=0;i<N;i++)  {
        get_random_cell()->RandomDivision();
    }
}

void Tissue::allCellsDivide()
{
    // obtain a list of pointers on all cells in the tissue to iterate through all of them without
    //std::list<Cell*> listCell_old

    // current number of cells
    unsigned N = ListCell.size();
    unsigned i=0;

    for(std::list<Cell>::iterator it_c = ListCell.begin(); it_c != ListCell.end(); it_c++)  {
        it_c->RandomDivision();
        i++;
        if(i==N) break; // make sure that no cell is divided twice
    }
}

void Tissue::divideCells(int cellType, int maxDivisions)
{
    // cell list is varied in the course of the divisions, so save a copy in order to not divide cells twice
    std::list<Cell*> listCell_old;
    
    int divisionCounter = 0;
    
    for(std::list<Cell>::iterator itc = ListCell.begin(); itc!=ListCell.end();itc++)  listCell_old.push_back(&(*itc));
    
    for(std::list<Cell*>::iterator it_c = listCell_old.begin(); it_c != listCell_old.end(); it_c++)  {
        if((*it_c)->Type != cellType) continue;
        (*it_c)->RandomDivision();
        
        // check if the maximum number of allowed cell divisions has been reached
        if(maxDivisions>0 && ++divisionCounter>=maxDivisions) return;
    }
}

void Tissue::mutate_cells_in_box(int cell_type, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max)
{
    // works only for non periodic tissues
    if(isPeriodic) return;
    
    // non valid cell type
    if(cell_type<0 || cell_type> NO_CELL_TYPES-1) return;
    
    for(std::list<Cell>::iterator it_c = ListCell.begin(); it_c != ListCell.end(); it_c++)  {
        // calculate the midpoint between the apical and the basal sides
        Point M = it_c->BC_a;
        
        // if the midpoint is in the box mutate the cell
        if(M.x>=x_min && M.x<=x_max && M.y>=y_min && M.y<=y_max && M.z>=z_min && M.z<=z_max) {
            it_c->Type = cell_type;
        }
    }
}

Cell* Tissue::cell_in_middle(int type)
{
    double minDistanceSum = 999999999999;
    Cell *bestCell = NULL;
    
    // find first cell of given type in the list
    for(std::list<Cell>::iterator itc = ListCell.begin(); itc != ListCell.end(); itc++)
    {
        Cell* c = &(*itc);
        
        if(c->Type == type) {
            double distanceSum = 0;
            
            for(std::list<Cell>::iterator itc2 = ListCell.begin(); itc2 != ListCell.end(); itc2++)
            {
                if(itc2->Type==type) distanceSum += Point::Norm(Point(c->BC_a.x,c->BC_a.y,0) - Point(itc2->BC_a.x,itc2->BC_a.y,0));
            }
            
            if(distanceSum<minDistanceSum) {
                minDistanceSum = distanceSum;
                bestCell = c;
            }
        }
    }
    
    return bestCell;
}