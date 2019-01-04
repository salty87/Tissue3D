#include "Cell.h"
#include <algorithm>
#include <time.h>
#include <iomanip> 
#include <math.h>

#define PI 3.1415926535897932384
#define EPSILON 10e-8

Cell::Cell()
{
	static int NumberOfCreatedCells(1);
	Number=NumberOfCreatedCells++;
	crossBoundary=false;
    q_BC_b = Quadrant();
	Type=1;
    positionFixed = false;
}

Cell::Cell(std::list<Edge*> ListEdges,double K, double Ta, double Tb, double V0)
:  ListEdges(ListEdges), K(K), Ta(Ta), Tb(Tb), V0(V0)
{
	static int NumberOfCreatedCells(1);
	Number=NumberOfCreatedCells++;
	
    crossBoundary=false;
    q_BC_b = Quadrant();
	Type=1;
    positionFixed = false;
}

bool Cell::operator==(Cell c)
{
	return(Number==c.Number);
}

int Cell::DirectionOfEdge(Edge* e)
{
	if((e->getCell1())==this){
		// edge	e is directed clockwise towards this cell
		return 1;
	}
	else if((e->getCell2())==this){
		// edge e is counterclockwise towards this cell
		return -1;
	}
	else{
		// the cell should not include the edge e
		return 0;
	}
}


void Cell::updatePerimeter()
{
	Perimeter_a = 0;
	Perimeter_b = 0;
	
	std::list<Edge*>::iterator ite=ListEdges.begin();
	for(;ite!=ListEdges.end();ite++)
	{
		Perimeter_a += (*ite)->l_a();
		Perimeter_b += (*ite)->l_b();
	}
}

// BC_as here are given by the weighted middle of the edges
void Cell::UpdateCenter()
{
    if(T->onlyApicalContributions)
    {
        BC_a = Point();
        
        for(std::list<Edge*>::iterator it_e = ListEdges.begin();it_e!=ListEdges.end();it_e++)
        {
            Edge *e = (*it_e);
            if(e->c1 == this)
            {
                BC_a += e->v1->coord;
            }
            else
            {
                BC_a += e->v2->coord;
            }
        }
        
        BC_a = BC_a / ListEdges.size();
        
        return;
    }
    

    //reorientAllEdges();
    
    if(T->BaryType == PointCenter)
    {
        if(!T->isPeriodic)
        {
            BC_a = Point();
            BC_b = Point();
            
            for(std::list<Edge*>::iterator it_e = ListEdges.begin();it_e!=ListEdges.end();it_e++)
            {
                Edge *e = (*it_e);
                if(e->c1 == this)
                {
                    BC_a += e->v1->coord;
                    BC_b += e->v1->coord_b;
                }
                else
                {
                    BC_a += e->v2->coord;
                    BC_b += e->v2->coord_b;
                }
            }
        
            BC_a = BC_a / ListEdges.size();
            BC_b = BC_b / ListEdges.size();
            
        }
        
        else
        {
            BC_a = Point();
            BC_b = Point();
            
            // the current quadrant is needed to keep in mind in which quadrant the current vertex lies
            Quadrant currentQuadrant=Quadrant();
            
            for(std::list<Edge*>::iterator ite = ListEdges.begin();ite!=ListEdges.end();ite++)  {
                Edge *e = (*ite);
                
                if(e->c1 == this )  { // get the vertex that is located in the counter clockwise direction of the edge towards the cell
                    BC_a += shift(e->v1->coord, currentQuadrant*T->SystemSize);
                    BC_b += shift(e->v1->coord_b, (e->v1->q_b+currentQuadrant)*T->SystemSize);
                                  
                    // save the relative position to the first vertex to be able to later locate the position towards the BC_a for the edge
                    e->q_c1 = currentQuadrant;
                    
                    // update the current quadrant = quadrant of the second vertex of this edge (towards the cell) with respect to the first vertex of the first edge of the cell
                    currentQuadrant += e->q_v2;
                } else {
                    BC_a += shift(e->v2->coord, currentQuadrant*T->SystemSize);
                    BC_b += shift(e->v2->coord_b, (e->v2->q_b+currentQuadrant)*T->SystemSize);
                    
                    // save the relative position to the secpmd vertex to be able to later locate the position towards the BC_a for the edge
                    e->q_c2 = currentQuadrant;

                    // update the current quadrant = quadrant of the second vertex of this edge (towards the cell) with respect to the first vertex of the first edge of the cell
                    currentQuadrant -= e->q_v2;                    
                }
            }
            
            // divide by the number of connected vertices
            BC_a = BC_a/ListEdges.size();
            BC_b = BC_b/ListEdges.size();
            
            
            // get the quadrant of the new barycenter wrt the previous position
            Quadrant q_BC_a_new = Quadrant(floor(BC_a.x/T->SystemSize.x),floor(BC_a.y/T->SystemSize.y));
            BC_a = shift(BC_a, T->SystemSize*q_BC_a_new*(-1));
            Quadrant q_BC_b_new = Quadrant(floor(BC_b.x/T->SystemSize.x),floor(BC_b.y/T->SystemSize.y));
            BC_b = shift(BC_b, T->SystemSize*q_BC_b_new*(-1));
            
            q_BC_b = q_BC_b_new - q_BC_a_new;
            
            for(std::list<Edge*>::iterator ite = ListEdges.begin();ite!=ListEdges.end();ite++) { // set the quadrant of BC_a relative to the edges
                if(this == (*ite)->getCell1())  {(*ite)->q_c1 = q_BC_a_new - (*ite)->q_c1;}
                else {(*ite)->q_c2 = q_BC_a_new - (*ite)->q_c2;}
            }
        }
        
    }
    else if(T->BaryType == ContourCenter)
    {
        
        double totalLength_a=0, totalLength_b=0;

        if(!T->isPeriodic || !CrossBoundary())
        {
            BC_a = Point();
            BC_b = Point();
            
            std::list<Edge*>::iterator it_e = ListEdges.begin();
            
            for(;it_e!=ListEdges.end();it_e++)
            {
                Edge *e = (*it_e);

                totalLength_a += e->l_a();
                totalLength_b += e->l_b();
                
                BC_a += e->MP_a()*e->l_a();
                BC_b += e->MP_b()*e->l_b();
            }
            
            // divide by the length of the perimeter
            BC_a = BC_a/(2.*totalLength_a);
            BC_b = BC_b/(2.*totalLength_b);
            
            
        }
        
        else
        {
            BC_a = Point();
            BC_b = Point();
            
            
            // the current quadrant is needed to keep in mind in which quadrant the current vertex lies
            Quadrant currentQuadrant=Quadrant();
            
            for(std::list<Edge*>::iterator ite = ListEdges.begin();ite!=ListEdges.end();ite++)
            {
                
                Edge *e = (*ite);
                
                totalLength_a += e->la;
                totalLength_b += e->lb;
                
                BC_a += e->MP_a(this) * e->la;
                BC_b += e->MP_b(this) * e->lb;
                
                // save the relative position to the first vertex to be able to later locate the position towards the BC_a for the edge
                e->q_c1 = currentQuadrant;
                
                if(this == e->getCell1())
                {
                                    
                    // update the current quadrant = quadrant of the second vertex of this edge (towards the cell) with respect to the first vertex of the first edge of the cell
                    currentQuadrant += e->q_v2;
                }
                else
                {
                    // update the current quadrant = quadrant of the second vertex of this edge (towards the cell) with respect to the first vertex of the first edge of the cell
                    currentQuadrant -= e->q_v2;
                }
            }
            
            // divide by the length of the perimeter
            BC_a = BC_a/(2*totalLength_a);
            BC_b = BC_b/(2*totalLength_b);
            
            
            // get the quadrant of the new barycenter wrt the previous position
            Quadrant q_BC_a = Quadrant(floor(BC_a.x/T->SystemSize.x),floor(BC_a.y/T->SystemSize.y));
            BC_a = shift(BC_a, T->SystemSize*q_BC_a*(-1));
            //q_BC_b = Quadrant(floor(BC_b.x/T->SystemSize.x),floor(BC_b.y/T->SystemSize.y));
            BC_b = shift(BC_b, T->SystemSize*q_BC_b*(-1));
            
            // set the quadrant of the BC_a relative to the edges of the cell
            for(std::list<Edge*>::iterator ite = ListEdges.begin();ite!=ListEdges.end();ite++)
            {
                if(this == (*ite)->getCell1())  {(*ite)->q_c1 = q_BC_a - (*ite)->q_c1;}
                else {(*ite)->q_c2 = q_BC_a - (*ite)->q_c2;}
            }

            
        }
        
    }
    
    
}

/// Ta as a function of Aa -> Ta = dW_dA!
double Cell::T_a()
{
    // elasticity is turned off
    if(Ta==T0_a || A0_a<=0) return Ta;
    else {
        double Aa = CalculateApicalSurfaceArea();
        if(Aa>=A0_a) return Ta;
        else return T0_a + Aa/A0_a*(Ta-T0_a);
    }
}

// Tb as a function of Ab
double Cell::T_b()
{
    // elasticity is turned off
    if(Tb==T0_b || A0_b<=0) return Tb;

    else {
      
        double Ab = CalculateBasalSurfaceArea();
        if(Ab>=A0_b) return Tb;
        else return T0_b + Ab/A0_b*(Tb-T0_b);

    }
}


double Cell::height()
{
    if(!T->isPeriodic || q_BC_b.isZero()) return BC_a.Distance(BC_b);
    else return BC_a.Distance(shift(BC_b,q_BC_b*T->SystemSize));
}


/// angle between the cell and (0,0,1)
double Cell::angleWithZPlane()
{
    Point direction;
    
    if(!T->isPeriodic || q_BC_b.isZero()) direction = BC_a-BC_b;
    else direction = BC_a-shift(BC_b,q_BC_b*T->SystemSize);
    
    return std::acos(Point(0,0,1).dot(direction/(Point::Norm(direction))));
}

double Cell::UpdateVolume()
{
	// first calculate current BC_as
	UpdateCenter();
    
	Volume=0;
    
	if (!T->isPeriodic || !crossBoundary)
    {
		std::list<Edge*>::iterator it;
		int direction;

		for(it=ListEdges.begin();it!=ListEdges.end();it++)	{	// iterate through all the edges of the cell
            Edge *e = (*it);

			// contributions from any edge to the cell:
			
			// 1st get the direction of the edge towards the cell, 1 = clockwise, -1 = counterclockwise
			direction=DirectionOfEdge(*it);
			
			// every edge yields for 6 contributions: four coming from the inside, and one each from apical and basal side
			// apical:
			double VolumeEdge = Point::TetrahedronSignedVolume(e->v1->coord,e->v2->coord,BC_a);
			//basal
			VolumeEdge+=Point::TetrahedronSignedVolume(e->v2->coord_b,e->v1->coord_b,BC_b);
			            
			// now calculate the contributions from the interior parts
			// top
			VolumeEdge+=Point::TetrahedronSignedVolume(e->v2->coord,e->v1->coord,e->BC_l);
			// bottom
			VolumeEdge+=Point::TetrahedronSignedVolume(e->v1->coord_b,e->v2->coord_b,e->BC_l);
			// towards v1
			VolumeEdge+=Point::TetrahedronSignedVolume(e->v1->coord,e->v1->coord_b,e->BC_l);
			// towards v2
			VolumeEdge+=Point::TetrahedronSignedVolume(e->v2->coord_b,e->v2->coord,e->BC_l);
        
			Volume+=VolumeEdge*direction;
		}
	} else {
		
		std::list<Edge*>::iterator it = ListEdges.begin();
		        
		// iterate through all the edges of the cell
		for(;it!=ListEdges.end();++it) {
			
            Edge *e = (*it);
						
			if((*it)->getCell1()==this)	{
				// position the vertices such that they are all in the same quadrant as the BC_a of the cell
				Point v1a_rel = shift( e->v1->coord, T->SystemSize*e->q_c1*(-1));
				Point v1b_rel = shift( e->v1->coord_b , T->SystemSize*(e->v1->q_b - e->q_c1));
				
				Point v2a_rel = shift( e->v2->coord, T->SystemSize*(e->q_v2 - e->q_c1));
				Point v2b_rel = shift( e->v2->coord_b, T->SystemSize*(e->q_v2 + e->v2->q_b - e->q_c1));
				
				Point BC_b_rel = shift( BC_b, T->SystemSize*q_BC_b);
                Point BC_l_rel = shift( e->BC_l, T->SystemSize*(e->q_c1 - e->q_BC_l)*(-1));
                
				
				// every edge yields for 6 contributions: four coming from the inside, and one each from apical and basal side
				// apical:
				double VolumeEdge = Point::TetrahedronSignedVolume(v1a_rel,v2a_rel,BC_a);
				//basal
				VolumeEdge+=Point::TetrahedronSignedVolume(v2b_rel,v1b_rel,BC_b_rel);

				// now calculate the contributions from the interior parts
				// top
				VolumeEdge+=Point::TetrahedronSignedVolume(v2a_rel,v1a_rel,BC_l_rel);
				// bottom
				VolumeEdge+=Point::TetrahedronSignedVolume(v1b_rel,v2b_rel,BC_l_rel);
				// towards v1
				VolumeEdge+=Point::TetrahedronSignedVolume(v1a_rel,v1b_rel,BC_l_rel);
				// towards v2
				VolumeEdge+=Point::TetrahedronSignedVolume(v2b_rel,v2a_rel,BC_l_rel);
                
                Volume+=VolumeEdge;
			
            } else {
                
				Point v1a_rel = shift( e->v1->coord,  T->SystemSize*( e->q_c2 + e->q_v2)*(-1));
				Point v1b_rel = shift( e->v1->coord_b,T->SystemSize*(e->v1->q_b - e->q_c2 - e->q_v2));
				
				Point v2a_rel = shift( e->v2->coord,   T->SystemSize*e->q_c2*(-1));
				Point v2b_rel = shift( e->v2->coord_b, T->SystemSize*(e->q_c2 - e->v2->q_b)*(-1));
				
				Point BC_b_rel = shift( BC_b, T->SystemSize*q_BC_b);
				Point BC_l_rel = shift( e->BC_l, T->SystemSize* (e->q_c2 + e->q_v2 - e->q_BC_l)*(-1));
                
				// every edge yields 6 contributions: four coming from the inside, and one each from apical and basal side
				// apical:
				double VolumeEdge = Point::TetrahedronSignedVolume(v1a_rel,v2a_rel,BC_a);
				//basal
				VolumeEdge+=Point::TetrahedronSignedVolume(v2b_rel,v1b_rel,BC_b_rel);
				
				// now calculate the contributions from the interior parts
				// top
				VolumeEdge+=Point::TetrahedronSignedVolume(v2a_rel,v1a_rel,BC_l_rel);
				// bottom
				VolumeEdge+=Point::TetrahedronSignedVolume(v1b_rel,v2b_rel,BC_l_rel);
				// towards v1
				VolumeEdge+=Point::TetrahedronSignedVolume(v1a_rel,v1b_rel,BC_l_rel);
				// towards v2
				VolumeEdge+=Point::TetrahedronSignedVolume(v2b_rel,v2a_rel,BC_l_rel);
				
				
				Volume-=VolumeEdge;
			}
			
		}
	}
    
    return Volume;
}


void Cell::UpdateElasticEnergy()
{
	UpdateVolume();
	   
    Pressure = K*(V0-Volume);
    ElasticEnergy = K/2.0*(V0-Volume)*(V0-Volume);
}

//
void Cell::assignCellType(int typeNo)
{
//    if(typeNo>=T->MechProp.CellPropVector.size()) {
//        std::cout << "No such cell type number (in Cell::assignCellType(int typeNo)) !" << std::endl;
//        return;
//    }
    
    Ta = T->MechProp.CellPropVector[typeNo].T_a;
    Tb = T->MechProp.CellPropVector[typeNo].T_b;
    V0 = T->MechProp.CellPropVector[typeNo].V0;
    K = T->MechProp.CellPropVector[typeNo].K;
    A0_a = T->MechProp.CellPropVector[typeNo].A0_a;
    A0_b = T->MechProp.CellPropVector[typeNo].A0_b;
    T0_a = T->MechProp.CellPropVector[typeNo].T0_a;
    T0_b = T->MechProp.CellPropVector[typeNo].T0_b;
}


bool Cell::CrossBoundary()
{
	std::list<Edge*>::iterator it;
	
	for(it=ListEdges.begin();it!=ListEdges.end();it++)
	{
		if((*it)->CrossBoundary())
		{
			crossBoundary=true;
			return true;
		}
	}
	
    crossBoundary = false;
	return false;
	
}

void Cell::Print()
{
	//Point Centre=Centroid;
	//std::cout<<"Edges of the cell\n";
    //	std::list<Edge*>::iterator it;
    //
    //	std::cout << ListEdges.size();
    //
    //	if(!ListEdges.empty())
    //	{
    //	for(it=ListEdges.begin();it!=ListEdges.end();it++)
    //	{
    //		(**it).Print();
    //	}
    //	}
    //
    //	std::cout<<"Number="<<Number<<", V0="<<V0<<", V="<<
    //	Volume<<", ElasticEnergy="<<ElasticEnergy<<", K="<<K;
    //	//", Type="<<Type<<", Centroid=("<<Centre.x<<","<<Centre.y<<")";
    //	std::cout<<"\n";
    //
	std::cout << "Cell:" << Number << ", SurfaceTensions:" << T_a() << ", " << T_b() << ", K:" << K << std::endl;
    
}

void Cell::Division(Edge* e1, Edge* e2, Edge* e_new, double ratio_e1, double ratio_e2, Cell *c_new)
{
    // if the cell crosses the boundary return -3, as this shall not be allowed for the moment
    if(T->isPeriodic && CrossBoundary()) return;
    
    // turn around all edges of the cell such that they are oriented all in the same way towards the cell
    orderEdges();
    
    // check if the two ratios are inbetween 0 and 1 and else return -1
    if(ratio_e1<0 || ratio_e1>1 || ratio_e2<0 || ratio_e2>1) return;
    
    // check if the two edges are parts of the cell and if not return -2
    std::list<Edge*>::iterator iter_e1 =  std::find(ListEdges.begin(),ListEdges.end(),e1);
    std::list<Edge*>::iterator iter_e2 =  std::find(ListEdges.begin(),ListEdges.end(),e2);
    
    if(iter_e1==ListEdges.end() || iter_e2==ListEdges.end()) return;
    
    // we obtain two new vertives, three new edges and one new cell
    
    // add all the newly created structures to the tissue lists
    T->ListVertex.push_back(*( new Vertex(T)));
    Vertex *v_e1 = &T->ListVertex.back();
    T->ListVertex.push_back(*(new Vertex(T)));
    Vertex *v_e2 = &T->ListVertex.back();
    
    T->ListEdge.push_back(*(new Edge(T,v_e1,e1->v2)));
    Edge *e_new1 = &T->ListEdge.back();
    T->ListEdge.push_back(*(new Edge(T,v_e2,e2->v2)));
    Edge *e_new2 = &T->ListEdge.back();
    T->ListEdge.push_back(*(new Edge(T,v_e1,v_e2)));
    e_new = &T->ListEdge.back();
    
    e_new->T = T;
    e_new->Number = T->maxEdgeNumber()+1;
    
    T->ListCell.push_back(Cell());
    c_new = &T->ListCell.back();
    
    // assign the newly created daughter cell the same mechanical properties as the old one
//    c_new->Ta = Ta;
//    c_new->Tb = Tb;
//    c_new->V0 = V0;
//    c_new->K = K;

    c_new->T = T;
    c_new->Type = Type;
    c_new->Number = T->maxCellNumber()+1;
    
    // prepare the new vertices by calculating the apical and basal coordinates and make them have the same line tension as the v1's of the edges
    v_e1->coord = e1->v1->coord*ratio_e1 + e1->v2->coord*(1-ratio_e1);
    v_e1->coord_b = e1->v1->coord_b*ratio_e1 + e1->v2->coord_b*(1-ratio_e1);
    v_e1->T = T;
    v_e1->Number = T->maxVertexNumber()+1;
    v_e2->coord = e2->v1->coord*ratio_e2 + e2->v2->coord*(1-ratio_e2);
    v_e2->coord_b = e2->v1->coord_b*ratio_e2 + e2->v2->coord_b*(1-ratio_e2);
    v_e1->G_l = e1->v1->G_l;
    v_e2->G_l = e2->v1->G_l;
    v_e2->T = T;
    v_e2->Number = v_e1->Number+1;
    // make the edges have the same mechanical properties as the whole edge had before
//    e_new1->T_l = e_new2->T_l = e1->T_l;
//    e_new1->G_a = e_new2->G_a = e1->G_a;
//    e_new1->G_b = e_new2->G_b = e1->G_b;

    e_new1->T = e_new2->T = T;

    e_new1->c1 = c_new;
    e_new1->c2 = e1->c2;

    e_new2->c1 = this;
    std::cout << e2->c2->Number << std::endl;
    e_new2->c2 = e2->c2;
    
    // assign the new vertices to the existing edges to make them the first half of their previous form
    e1->v2 = v_e1;
    e2->v2 = v_e2;
    e2->c1 = c_new;
    
    // now create the new edge that divides the cell and give it the same properties as edges inbetween two cells of this cell type have
//    e_new->T_l = T->MechProp.Intercell_SurfaceTension[Type][Type];
//    e_new->G_a = T->MechProp.Intercell_ApicalLineTension[Type][Type];
//    e_new->G_b = T->MechProp.Intercell_BasalLineTension[Type][Type];
    e_new->c1 = this;
    e_new->c2 = c_new;

    std::list<Edge*>::iterator iter_e;
    
    for(iter_e = ListEdges.begin(); iter_e != ListEdges.end(); ++iter_e)
    {
        std::cout<<(*iter_e)->Number << std::endl;
    }
    
    // now arrange the edges of the two cells after the division
    std::list<Edge*> ListEdges_new;
    ListEdges_new.push_back(e1); // first part of divided edge1
    ListEdges_new.push_back(e_new); // dividing edge
    ListEdges_new.push_back(e_new2); // second part of divided edge2

    std::cout<<e1->Number<<" "<< e_new->Number<<" "<< e_new2->Number<< std::endl;
    
    iter_e = iter_e2; iter_e++;
    
    for(; iter_e != iter_e1 && iter_e != ListEdges.end(); ++iter_e)
    {
        std::cout<<(*iter_e)->Number << std::endl;
        ListEdges_new.push_back(*iter_e);
    }
    
    if(iter_e==ListEdges.end())
    {
        for(iter_e = ListEdges.begin(); iter_e != iter_e1 ; ++iter_e)
        {
            ListEdges_new.push_back(*iter_e);
        }
    }
    
    c_new->ListEdges.push_back(e2);
    c_new->ListEdges.push_back(e_new);
    c_new->ListEdges.push_back(e_new1);
    
    iter_e = iter_e1; iter_e++;
    for(; iter_e != iter_e2 && iter_e != ListEdges.end(); ++iter_e)
    {
        std::cout<<(*iter_e)->Number<<std::endl;
        c_new->ListEdges.push_back(*iter_e);
        (*iter_e)->c1 = c_new;
    }
    
    if(iter_e==ListEdges.end())
    {
        for(iter_e = ListEdges.begin(); iter_e != iter_e2 ; ++iter_e)
        {
            c_new->ListEdges.push_back(*iter_e);
            (*iter_e)->c1 = c_new;
        }
    }

    ListEdges = ListEdges_new;
    
    // identify the two cells that are adjacent to the edges 1 and 2 (c_e1 / c_e2) who are also changing their lists
    Cell *c_e1 = e1->c2;
    Cell *c_e2 = e2->c2;
    
    // insert the edge e_new1 in c_e1 before e1 and e_new2 in c_e2 after e2
    c_e1->ListEdges.insert(std::find(c_e1->ListEdges.begin(),c_e1->ListEdges.end(),e1),e_new1);
    c_e2->ListEdges.insert(std::find(c_e2->ListEdges.begin(),c_e2->ListEdges.end(),e2),e_new2);
    
    // set the mechanical properties of the newly created structures (known from type / type of constituting cells)
    c_new->setMechanicalProperties();
    e_new->setMechanicalProperties();
    e_new1->setMechanicalProperties();
    e_new2->setMechanicalProperties();
    
    T->update();
    
    //std::cout << e_new1->Number;
}

void Cell::RandomDivision()
{
    unsigned numberOfEdges = ListEdges.size();
    
    // draw random numbers to determine two edges of the cell that will be splitted
    int rand1 = rand() % numberOfEdges;
    int rand2 = (rand1 + numberOfEdges/2) % numberOfEdges;
    
    
    Edge *e1, *e2;
    
    int counter = 0;
    // find the edges that correspond to those numbers
    for(std::list<Edge*>::iterator iter_e = ListEdges.begin(); iter_e!=ListEdges.end(); ++iter_e)
    {
        if(counter==rand1) e1 = (*iter_e);
        else if(counter==rand2) e2 = (*iter_e);
        counter++;
    }
    
    // call the cell division routine with 50-50 splits
    Division(e1,e2,NULL, 0.5,0.5, NULL);
        
}

void Cell::reorientAllEdges()
{
    // get all edges that are neighbours to the cell
    ListEdges.clear();

    for(std::list<Edge>::iterator iter_e = T->ListEdge.begin(); iter_e != T->ListEdge.end(); ++iter_e)
    {
        Edge *e = &(*iter_e);
        if(e->c1 == this)
        {
            ListEdges.push_back(e);
        } else if(e->c2 == this) {
            ListEdges.push_back(e);
            
            // turn around the edge
            Vertex *v_temp = e->v1;
            e->v1 = e->v2;
            e->v2 = v_temp;
            
            // swap the cells
            e->c2 = e->c1;
            e->c1 = this;
            
            // swap the quadrants of the cell
            Quadrant q_temp = e->q_c1;
            e->q_c1 = e->q_c2;
            e->q_c2 = q_temp;
            
            // turn around the quadrant of the edge
            e->q_v2 = e->q_v2*(-1);
        }
    }
    
//    for(std::list<Edge*>::iterator iter_e = ListEdges.begin(); iter_e != ListEdges.end(); ++iter_e)
//    {
//        Edge *e = (*iter_e);
//        if(!(e->c1 == this))
//        {
//            // turn around the edge
//            Vertex *v_temp = e->v1;
//            e->v1 = e->v2;
//            e->v2 = v_temp;
//            
//            // swap the cells
//            e->c2 = e->c1;
//            e->c1 = this;
//            
//            // swap the quadrants of the cell
//            Quadrant q_temp = e->q_c1;
//            e->q_c1 = e->q_c2;
//            e->q_c2 = q_temp;
//            
//            // turn around the quadrant of the edge
//            e->q_v2 = e->q_v2*(-1);
//        }
//    }
    
//    for(std::list<Edge*>::iterator iter_e = ListEdges.begin(); iter_e != ListEdges.end(); ++iter_e)
//    {
//        std::cout << (*iter_e)->v1->Number << "  " << (*iter_e)->v2->Number << std::endl;
//    }
    
    //return;
}

void Cell::orderEdges()
{
    reorientAllEdges();
    
    // now order them correctly in ListEdges
    std::list<Edge*> new_list_edge;
    std::list<Edge*> old_list_edge = ListEdges;
    
    // add starting element
    new_list_edge.push_back(old_list_edge.front());
    old_list_edge.pop_front();
    
    while(!old_list_edge.empty())
    {
        for(std::list<Edge*>::iterator ite = old_list_edge.begin(); ite!=old_list_edge.end();ite++)
        {
            if((*ite)->v1 == new_list_edge.back()->v2)// || (*ite)->v2 == new_list_edge.back()->v2 || (*ite)->v1 == new_list_edge.back()->v1 || (*ite)->v2 == new_list_edge.back()->v1 )
            {
                new_list_edge.push_back(*ite);
                old_list_edge.erase(ite);
                break;
            }
        }
    }
    
    ListEdges = new_list_edge;
}


Vertex* Cell::Apoptosis()
{
	/*  three steps:  1. move all the vertices of the cell to the current (apical/basal) BC_as
     2. pop all the edges of the cell
     3. destroy the cell
	 */
	
	// save initial barycenters of cell
	UpdateCenter();
    Point barycenter_a = BC_a;
    Point barycenter_b = BC_b;
	
    //
    Vertex *v;
	// now iterate through all edges that form the cell
    while(!ListEdges.empty()) {
        v = (*ListEdges.begin())->T1_remove();
    }
    
//    for(std::list<Edge*>::iterator it_edge = ListEdges.begin(); it_edge != ListEdges.end(); it_edge++)
//    {
//        v = (*it_edge)->T1_remove();
//    }
    
    //v->coord = barycenter_a; v->coord_b = barycenter_b;
//    v->T = T;
//    
//    T->SetTissueStructure();
//    
//    v->updateNeighbouringEdges();
//    
    return v;
//	
//	// iterate through all the vertices of the cell to get all the edges the cell is connected to
//	std::list<Vertex*>::iterator it_vertex=T->ListVertex.begin();
//	std::list<Edge*>::iterator it_edge;
//	for(;it_vertex!=ListVertex.end();it_vertex++)
//	{
//		Vertex *vit = *it_vertex;
//		// iterate through all the edges that are connected to vit
//		for(it_edge=vit->getListEdgeBegin();it_edge!=vit->getListEdgeEnd();it_edge++)
//		{
//			std::list<Vertex>::iterator itv2 = T->ListVertex.begin();
//			
//			if((*it_edge)->getCell1()!=NULL && (*it_edge)->getCell2()!=NULL) // cell is not on the boundary
//			{
//				if((*it_edge)->getCell1()!=this && (*it_edge)->getCell2()!=this) // if the edge is not part of this cell
//				{
//					if((*it_edge)->VertexIsPart(vit)==1)
//						(*it_edge)->setVertex1(v);
//					else
//						(*it_edge)->setVertex2(v);
//					
//					v->addEdge(*it_edge);
//                    
//				}
//			}
//			
//			else if((*it_edge)->getCell1()==NULL) // cell is not on the boundary
//			{
//				if((*it_edge)->getCell2()!=this) // if the edge is not part of this cell
//				{
//					if((*it_edge)->VertexIsPart(vit)==1)
//						(*it_edge)->setVertex1(v);
//					else
//						(*it_edge)->setVertex2(v);
//					
//					v->addEdge(*it_edge);
//                    
//				}
//			}
//			
//			else if((*it_edge)->getCell2()==NULL) // cell is not on the boundary
//			{
//				if((*it_edge)->getCell1()!=this) // if the edge is not part of this cell
//				{
//					if((*it_edge)->VertexIsPart(vit)==1)
//						(*it_edge)->setVertex1(v);
//					else
//						(*it_edge)->setVertex2(v);
//					
//					v->addEdge(*it_edge);
//                    
//				}
//			}
//			
//			
//		}
//	}
//	
//	// at this stage all the edges that touched this cell point to the new vertex
//	
//	// now the edges of the cell have to be removed from the tissue and so have the vertices
//	while(ListEdges.size()!=0)
//	{
//		T->RemoveEdge(*ListEdges.begin());
//	}
//	
//	while(ListVertex.size()!=0)
//	{
//		//Vertex *e = ListEdges.pop_front();
//		T->RemoveVertex(*ListVertex.begin());
//		ListVertex.pop_front();
//	}
//	
//	std::list<Vertex>::iterator itv = T->ListVertex.begin();
//	
//	for(; itv != T->ListVertex.end();itv++)
//	{
//		//std::cout << (*itv).getNumber() << "  ";
//	}
//    
//	//T->PrintCells();
	
}


/// returns 0 if the two cells dont have a common edge or 1 and the first common edge that has been found
bool Cell::hasCommonEdgeWith(Cell * c, Edge ** e)
{
	std::list<Edge*>::iterator ite = ListEdges.begin();
	std::list<Edge*>::iterator found;
	
	if(c==this) return 0;
	
	for(;ite!=ListEdges.end();ite++)
	{
		found = std::find(c->ListEdges.begin(), c->ListEdges.end(), (*ite));
        
		if(found!=c->ListEdges.end())
		{
			*e = (*found);
			return 1;
		}
		
	}
	
	return 0;
}

double Cell::energy()
{
    UpdateElasticEnergy();
    double E = ElasticEnergy;
    
    for(std::list<Edge*>::iterator it_e = ListEdges.begin(); it_e!=ListEdges.end(); it_e++)
    {
        E += (*it_e)->energy(this);
    }
    
    return E;
}

double Cell::W()
{
    if(positionFixed) return 0;
    
    // in the case without cellularization there is only a contribution from apical surface tension
    if(T->onlyApicalContributions) {
        double Aa = CalculateApicalSurfaceArea();
        
        if(A0_a>0 && T0_a!=Ta) { // area dependend apical tension
            if(Aa>A0_a) return A0_a*(Ta+T0_a)/2 + (Aa-A0_a)*Ta;
            else return Aa*T0_a + pow(Aa,2)*(Ta-T0_a)/(2*A0_a);
        } else {
            return Ta*Aa;
        }
    };
    
    UpdateVolume();
	
    Pressure = K*(V0-Volume);
    
	ElasticEnergy = K/2.0*(V0-Volume)*(V0-Volume);
    
    double energy = ElasticEnergy;
    
    // contribution from apical surface
    double Aa = CalculateApicalSurfaceArea();
    double Ab = CalculateBasalSurfaceArea();
    
    if(A0_a>0 && T0_a!=Ta) { // area dependend apical tension
        if(Aa>A0_a) energy += A0_a*(Ta+T0_a)/2 + (Aa-A0_a)*Ta;
        else energy += Aa*T0_a + pow(Aa,2)*(Ta-T0_a)/(2*A0_a);
    } else {
        energy += Ta*Aa;
    }
    
    if(A0_b>0 && T0_b!=Tb) { // area dependend basal tension
        if(Ab>A0_b) energy += A0_b*(Tb+T0_b)/2 + (Ab-A0_b)*Tb;
        else energy += Ab*T0_b + pow(Ab,2)*(Tb-T0_b)/(2*A0_b);
    } else {
        energy += Tb*Ab;
    }

    
    return energy;
}

double Cell::CalculateApicalSurfaceArea()
{
    if(!CrossBoundary()) {
    
        double apicalArea = 0;
        
        for(std::list<Edge*>::iterator it_e = ListEdges.begin(); it_e!=ListEdges.end(); it_e++) {   // iterate over all connected edges
            apicalArea+=Point::TriangleArea((*it_e)->v1->coord, (*it_e)->v2->coord, BC_a);
        }
        
        return apicalArea;
    } else {
        
        double apicalArea = 0;

        for(std::list<Edge*>::iterator it_e = ListEdges.begin(); it_e!=ListEdges.end(); it_e++) {   // iterate over all connected edges
            
            Edge* e = *it_e;
            
            if(e->c1==this) {
                // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE APICAL BARYCENTER OF CELL 1
                Point v1a_rel = shift( e->v1->coord, T->SystemSize*e->q_c1*(-1));
                Point v2a_rel = shift( e->v2->coord, T->SystemSize*( e->q_v2 -  e->q_c1));

                apicalArea+=Point::TriangleArea(v1a_rel,v2a_rel, BC_a);
            } else {
                Point v1a_rel = shift( e->v1->coord, T->SystemSize*(  e->q_c2 +  e->q_v2 ) *(-1));
                Point v2a_rel = shift( e->v2->coord, T->SystemSize* e->q_c2*(-1));

                apicalArea+=Point::TriangleArea(v1a_rel,v2a_rel, BC_a);
            }
        
        }

        return apicalArea;
    }
}

double Cell::CalculateBasalSurfaceArea()
{
    if(!CrossBoundary()) {
        
        double basalArea = 0;
        
        for(std::list<Edge*>::iterator it_e = ListEdges.begin(); it_e!=ListEdges.end(); it_e++) {   // iterate over all connected edges
            basalArea+=Point::TriangleArea((*it_e)->v1->coord_b, (*it_e)->v2->coord_b, BC_b);
        }
        
        return basalArea;
    } else {
        
        double basalArea = 0;
        
        for(std::list<Edge*>::iterator it_e = ListEdges.begin(); it_e!=ListEdges.end(); it_e++) {   // iterate over all connected edges
            
            Edge* e = *it_e;
            
            //if(e->v1->q_b.x!=0 || e->v1->q_b.y!=0 || e->v2->q_b.x!=0 ||e->v2->q_b.y!=0) std::cout << "cross Boundary!" << std::endl;
            
            if(e->c1==this) {
                // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE basal BARYCENTER OF CELL 1
                Point v1b_rel = shift( e->v1->coord_b, T->SystemSize*(e->q_c1-e->v1->q_b)*(-1));
                Point v2b_rel = shift( e->v2->coord_b, T->SystemSize*( e->q_v2 -  e->q_c1 + e->v2->q_b));
                
                Point BC_b_rel = shift(BC_b, T->SystemSize*q_BC_b);
                
                basalArea+=Point::TriangleArea(v1b_rel,v2b_rel, BC_b_rel);
            } else {
                Point v1b_rel = shift( e->v1->coord_b, T->SystemSize*(  e->q_c2 +  e->q_v2 - e->v1->q_b) *(-1));
                Point v2b_rel = shift( e->v2->coord_b, T->SystemSize* (e->q_c2 - e->v2->q_b)*(-1));
                
                Point BC_b_rel = shift(BC_b, T->SystemSize*q_BC_b);
                
                basalArea+=Point::TriangleArea(v1b_rel,v2b_rel, BC_b_rel);
            }
            
        }
        
//        //std::cout << "basal area = "<< basalArea<<std::endl;
//        
//        if(basalArea>500 || basalArea<-500)
//        {
//            std::cout << "Basal area = "<< basalArea<<std::endl << std::endl << std::endl;
//
//            
//            for(std::list<Edge*>::iterator it_e = ListEdges.begin(); it_e!=ListEdges.end(); it_e++) {   // iterate over all connected edges
//                
//                Edge* e = *it_e;
//                
//                //if(e->v1->q_b.x!=0 || e->v1->q_b.y!=0 || e->v2->q_b.x!=0 ||e->v2->q_b.y!=0) std::cout << "cross Boundary!" << std::endl;
//                
////                std::cout << "BC_a = " << std::endl;
////                BC_a.Print();
//
//                
//                if(e->c1==this) {
//                    // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE basal BARYCENTER OF CELL 1
//                    Point v1b_rel = shift( e->v1->coord_b, T->SystemSize*(e->q_c1-e->v1->q_b)*(-1));
//                    Point v2b_rel = shift( e->v2->coord_b, T->SystemSize*( e->q_v2 -  e->q_c1 + e->v2->q_b));
//                    
//                    Point BC_b_rel = shift(BC_b, T->SystemSize*(q_BC_b));
//                    
//                    std::cout << "v1b_rel 1 = " << std::endl;
//                    v1b_rel.Print();
//                    std::cout << "v2b_rel 1 = " << std::endl;
//                    v2b_rel.Print();
//                    std::cout << "BC_b_rel 1 = " << std::endl;
//                    BC_b_rel.Print();
//                    std::cout << "Added Area =" << Point::TriangleArea(v1b_rel,v2b_rel, BC_b_rel) << std::endl;
//
//                } else {
//                    Point v1b_rel = shift( e->v1->coord_b, T->SystemSize*(  e->q_c2 +  e->q_v2 - e->v1->q_b) *(-1));
//                    Point v2b_rel = shift( e->v2->coord_b, T->SystemSize* (e->q_c2 - e->v2->q_b)*(-1));
//                    
//                    Point BC_b_rel = shift(BC_b, T->SystemSize*(q_BC_b));
//                    
//                    std::cout << "v1b_rel 2 = " << std::endl;
//                    v1b_rel.Print();
//                    std::cout << "v2b_rel 2 = " << std::endl;
//                    v2b_rel.Print();
//                    std::cout << "BC_b_rel 2 = " << std::endl;
//                    BC_b_rel.Print();
//                    std::cout << "Added Area =" << Point::TriangleArea(v1b_rel,v2b_rel, BC_b_rel)<< std::endl;
//                }
//                
//            }
//            
//        }
        return basalArea;
    }
}


void Cell::setMechanicalProperties()
{
    Ta =    T->MechProp.CellPropVector[Type].T_a;
    T0_a =  T->MechProp.CellPropVector[Type].T0_a;
    A0_a =  T->MechProp.CellPropVector[Type].A0_a;
    
    Tb =    T->MechProp.CellPropVector[Type].T_b;
    T0_b =  T->MechProp.CellPropVector[Type].T0_b;
    A0_b =  T->MechProp.CellPropVector[Type].T0_b;
    
    K =     T->MechProp.CellPropVector[Type].K;
    V0 =    T->MechProp.CellPropVector[Type].V0;
}

//
//double Cell::CalculateBasalSurfaceArea()
//{
//    if(!CrossBoundary()) {
//        
//        double basalArea = 0;
//        
//        for(std::list<Edge*>::iterator it_e = ListEdges.begin(); it_e!=ListEdges.end(); it_e++) {   // iterate over all connected edges
//            basalArea+=Point::TriangleArea((*it_e)->v1->coord_b, (*it_e)->v2->coord_b, BC_b);
//        }
//        
//        return basalArea;
//        
//    } else {
//        
//        double basalArea = 0;
//        
//        for(std::list<Edge*>::iterator it_e = ListEdges.begin(); it_e!=ListEdges.end(); it_e++) {   // iterate over all connected edges
//            
//            Edge* e = *it_e;
//            
//            if(e->c1==this) {
//                // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE APICAL BARYCENTER OF CELL 1
//                Point v1b_rel = shift( e->v1->coord_b, T->SystemSize*( e->q_c1+e->v1->q_b)*(-1));
//                Point v2b_rel = shift( e->v2->coord_b, T->SystemSize*( e->q_v2 - e->q_c1+e->v2->q_b));
//                
//                basalArea+=Point::TriangleArea(v1b_rel,v2b_rel, BC_b);
//            } else {
//                Point v1b_rel = shift( e->v1->coord_b, T->SystemSize* (e->q_c2 + e->q_v2 +e->v1->q_b) *(-1));
//                Point v2b_rel = shift( e->v2->coord_b, T->SystemSize* (e->q_c2 + e->v2->q_b)*(-1));
//                
//                basalArea+=Point::TriangleArea(v1b_rel,v2b_rel, BC_b);
//            }
//
//        }
//        
//        std::cout << "Basal area = "<< basalArea<<std::endl;
//        
//        return basalArea;
//    }
//}

// works so far only if the cell does not cross the boundary!
//double Cell::CalculateBasalSurfaceArea()
//{
//    if(CrossBoundary()) return -1;
//    
//    double basalArea = 0;
//    
//    for(std::list<Edge*>::iterator it_e = ListEdges.begin(); it_e!=ListEdges.end(); it_e++) {   // iterate over all connected edges
//        Edge* e = (*it_e);
//        basalArea+=Point::TriangleArea(e->v1->coord_b, e->v2->coord_b, BC_b);
//    }
//    
//    return basalArea;
//}

//// calculate the tension given for a certain apical surface area if we apply a Michaelis-Menten type of tension-surface area relationship
//double Cell::T_a_minArea()
//{
//    double apicalSurfaceArea = CalculateApicalSurfaceArea();
//    
//    if(apicalSurfaceArea>150) return T_a;
//    
//    double A_min = 60;
//    int n = 5;
//    double T_greater0 = T_a;
//    double T_smaller0 = 10;
//    
//    double ret = T_smaller0*((T_smaller0-T_greater0)/(T_greater0*pow(A_min/apicalSurfaceArea, n ) - T_smaller0)+1.0);
//    
//    //double ret = T_a*(2.0/(1+pow(A_min/apicalSurfaceArea, n )) - 1);
//    
//    //if(apicalSurfaceArea<140)
//    std::cout << "Apical SurfaceArea = " << apicalSurfaceArea << " , T_a_minArea() = " << ret <<std::endl;
//    
//    return ret;
//    
//}

// calculate the tension given for a certain apical surface area if we suppose a quadratic potential around a preferred apical area A0
//double Cell::T_a_minArea()
//{
//    double apicalSurfaceArea = CalculateApicalSurfaceArea();
//    
//    return T_a*(apicalSurfaceArea-T->cells_A0);
//
//}