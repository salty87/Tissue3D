#include "Edge.h"

#define PI 3.1415926535897932384
#define EPSILON 10e-5

Edge::Edge(Tissue *T)
: v1(NULL), v2(NULL), c1(NULL), c2(NULL), T(T)
{
	Number= T->maxEdgeNumber();
	crossBoundary=false;
    
    q_v2 = Quadrant();
    q_c1 = Quadrant();
    q_c2 = Quadrant();
    
    l0_a = 3;
    l0_b = 3;
};

Edge::Edge(Tissue *T, Vertex *v1, Vertex *v2)
: v1(v1), v2(v2), c1(NULL), c2(NULL), T(T)
{
	Number= T->maxEdgeNumber()+1;
	crossBoundary=false;
   
    q_v2 = Quadrant();
    q_c1 = Quadrant();
    q_c2 = Quadrant();
    
    l0_a = 3;
    l0_b = 3;
}


Edge::Edge(Tissue *T, Vertex *v1, Vertex *v2, double G_a, double G_b, double T_l)
: v1(v1), v2(v2), c1(NULL), c2(NULL), G_a(G_a), G_b(G_b), T_l(T_l), T(T)
{
	Number= T->maxEdgeNumber()+1;
	crossBoundary=false;
    
    q_v2 = Quadrant();
    q_c1 = Quadrant();
    q_c2 = Quadrant();
    
    l0_a = 3;
    l0_b = 3;
}

// complete constructor
Edge::Edge(Tissue *T, Vertex *v1, Vertex *v2, Cell* c1, Cell* c2, double G_a, double G_b, double T_l)
: v1(v1), v2(v2), c1(c1), c2(c2), G_a(G_a), G_b(G_b), T_l(T_l), T(T)
{
	Number= T->maxEdgeNumber()+1;
	
    crossBoundary=false;
    q_v2 = Quadrant(0,0);
	Type=0;
    q_c1 = Quadrant();
    q_c2 = Quadrant();
    
    l0_a = 3;
    l0_b = 3;
}

// operator overloadings
bool Edge::operator==(const Edge & e)
{
	return((*e.v1==*v1)&&(*e.v2==*v2));
}


void Edge::UpdateLengths()
{
    la = l_a();
    if(T->onlyApicalContributions) return; // no need to calculate basal lengths
    lb = l_b();
}

double Edge::l_a()
{
	// if the apical side doesnt cross the boundary
	if(!T->isPeriodic || q_v2.isZero())
	{
		return Point::Norm(v1->getCoord()-v2->getCoord());
	}
	else
	{
		return Point::Norm(shift( v1->getCoord()-v2->getCoord(),q_v2*T->SystemSize*(-1)));
	}
}
	
double Edge::l_b()
{
	if(!T->isPeriodic)
	{
		return Point::Norm(v1->getBasalCoord()-v2->getBasalCoord());
	}

	else 
	{
		Quadrant Quadrant_v2_basal = v1->q_b - q_v2 - v2->q_b;
		
		return Point::Norm(shift( v1->getBasalCoord()-v2->getBasalCoord(), Quadrant_v2_basal*T->SystemSize));
		
	}

}


/// length dependent tension: dW/dla
double Edge::Lambda_a()
{
    double la = l_a();
    
    if(la>l0_a) return G_a;
    else return la*G_a/l0_a;
}

/// length dependent tension: dW/dlb
double Edge::Lambda_b()
{
    double lb = l_b();
    
    if(lb>l0_b) return G_b;
    else return lb*G_b/l0_b;
}

// update barycenter of the face inbetween the apical and the basal edge
bool Edge::UpdateBarycenter()
{
    if(T->BaryType == PointCenter)
    {
        
        if(!T->isPeriodic)
        {
            //average the 4 points of the two vertices
            BC_l=(v1->coord + v2->coord + v1->coord_b + v2->coord_b)/4.;
        }
        else
        {
            //Point SystemSize=T->SystemSize;
            BC_l = shift(v1->getCoord()+v1->getBasalCoord()+v2->getCoord()+v2->getBasalCoord() , T->SystemSize*(v1->getBasalQuadrant()+q_v2*2+v2->getBasalQuadrant()))/4.;
            
            //where, i.e. in which quadrant, is the barycenter in relation to v1?
            q_BC_l = Quadrant((BC_l.x>T->SystemSize.x) - (BC_l.x<0),(BC_l.y>T->SystemSize.y) - (BC_l.y<0));
            
            BC_l=Quadrant::Modulus(BC_l,T->SystemSize);
        }
    }
    else if(T->BaryType == ContourCenter)
    {        
        if(!T->isPeriodic)
        {
            double l1 = v1->l_l;
            double l2 = v2->l_l;

            // weighted average of the 4 midpoints
            BC_l= (v1->coord*(la+l1) + v2->coord*(la+l2) + v1->coord_b*(lb+l1) + v2->coord_b*(lb+l2))/(2*(la + lb + l1 + l2));
        }
        else
        {
            double l1 = v1->l_l;
            double l2 = v2->l_l;
            
            Point v2a = shift(v2->coord,q_v2*T->SystemSize);
            Point v1b = shift(v1->coord_b,v1->q_b*T->SystemSize);
            Point v2b = shift(v2->coord_b,(q_v2+v2->q_b)*T->SystemSize);
            
            // weighted average of the 4 midpoints
            BC_l= (v1->coord*(la+l1) + v2a*(la+l2) + v1b*(lb+l1) + v2b*(lb+l2))/(2*(la + lb + l1 + l2));
            
            //where, i.e. in which quadrant, is the barycenter in relation to v1?
            q_BC_l = Quadrant((BC_l.x>T->SystemSize.x) - (BC_l.x<0),(BC_l.y>T->SystemSize.y) - (BC_l.y<0));
            
            BC_l=Quadrant::Modulus(BC_l,T->SystemSize);
        }
    }
	
	return true;
}

// what are the forces acting on v1/v1_b of the edge, that try to minimize the length of the edge
Point Edge::ApicalLineForce()
{
	// save calculation time if there is no line tension anyway
	if(Lambda_a()==0) return Point();
	
	// if the edge does not cross the boundary
	if(!T->isPeriodic || !crossBoundary)
	{
		return (v2->getCoord()-v1->getCoord())*(Lambda_a()/l_a());
	}
	// if the edge crosses the periodic boundary
	else
	{
		// vector between the two points
		Point connection = shift(v2->getCoord(),T->SystemSize*q_v2) - v1->getCoord();
		// norm of the connection
		double la = Point::Norm(connection);
		
        return connection*(Lambda_a()/la);
	}
	
}

Point Edge::BasalLineForce()
{
	// save calculation time if there is no line tension anyway
	if(Lambda_b()==0) return Point();
	
	if(!T->isPeriodic || !crossBoundary)
	{
        return (v2->getBasalCoord()-v1->getBasalCoord())*(Lambda_b()/l_b());
	}
	else
	{
		Point connection = shift( v2->coord_b, T->SystemSize*(q_v2+v2->q_b-v1->q_b)) - v1->coord_b;
		
        return connection*(Lambda_b()/Point::Norm(connection));
		
	}
}

void Edge::UpdateForcesFromBarycenters()
{
	if(!T->isPeriodic) { // non-periodic tissue
		// add the contributions to the cells 
		if(c1 && c2) // save computational time, if edge has two neighbouring cells
		{
            // apical contributions
			c1->dAa_dBCa += Point::TriangleAreaChange(c1->BC_a,v1->getCoord(),v2->getCoord());
			c2->dAa_dBCa += Point::TriangleAreaChange(c2->BC_a,v1->getCoord(),v2->getCoord());
            Point volumeForce= Point::dV_dp(v1->getCoord(),v2->getCoord());
			c1->dV_dBCa += volumeForce;
			c2->dV_dBCa -= volumeForce;
			
            if(T->onlyApicalContributions) return;
            
            // basal contributions
			c1->dAb_dBCb  += Point::TriangleAreaChange(c1->BC_b,v1->getBasalCoord(),v2->getBasalCoord());
			c2->dAb_dBCb  += Point::TriangleAreaChange(c2->BC_b,v1->getBasalCoord(),v2->getBasalCoord());
			volumeForce = Point::dV_dp(v2->getBasalCoord(),v1->getBasalCoord());
			c1->dV_dBCb += volumeForce;
			c2->dV_dBCb -= volumeForce;
		}
		else 
		{
			if(c1)
			{
				c1->dAa_dBCa += Point::TriangleAreaChange(c1->BC_a,v1->getCoord(),v2->getCoord());
				c1->dAb_dBCb += Point::TriangleAreaChange(c1->BC_b,v1->getBasalCoord(),v2->getBasalCoord());
				
				// add the volume contributions for cell 1
				c1->dV_dBCa += Point::dV_dp(v1->getCoord(),v2->getCoord());
				c1->dV_dBCb  += Point::dV_dp(v2->getBasalCoord(),v1->getBasalCoord());
			}
			if(c2)
			{
				c2->dAa_dBCa += Point::TriangleAreaChange(c2->BC_a,v1->getCoord(),v2->getCoord());
				c2->dAb_dBCb += Point::TriangleAreaChange(c2->BC_b,v1->getBasalCoord(),v2->getBasalCoord());
				
				// add the volume contributions for cell 2
				c2->dV_dBCa += Point::dV_dp(v2->getCoord(),v1->getCoord());
				c2->dV_dBCb  += Point::dV_dp(v1->getBasalCoord(),v2->getBasalCoord());
			}
		}
	} else { // periodic tissue
		// if both neighbouring cells don't cross the boundary
		if(!c1->crossBoundary && !c2->crossBoundary)
		{
            // area derivatives with respect to movement in barycenter (IS LATER IN TISSUE MULTIPLIED WITH THE SURFACE TENSIONS OF THE CELLS)
			c1->dAa_dBCa += Point::TriangleAreaChange(c1->BC_a,v1->getCoord(),v2->getCoord());
			c1->dAb_dBCb += Point::TriangleAreaChange(c1->BC_b,v1->getBasalCoord(),v2->getBasalCoord());
			c2->dAa_dBCa += Point::TriangleAreaChange(c2->BC_a,v1->getCoord(),v2->getCoord());
			c2->dAb_dBCb += Point::TriangleAreaChange(c2->BC_b,v1->getBasalCoord(),v2->getBasalCoord());
			
			// volume change from barycenter movement
			Point dummy_a = Point::dV_dp(v1->getCoord(), v2->getCoord());
			Point dummy_b = Point::dV_dp(v2->getBasalCoord(), v1->getBasalCoord());
			
			c1->dV_dBCa += dummy_a;
			c1->dV_dBCb += dummy_b;
			c2->dV_dBCa -= dummy_a;
			c2->dV_dBCb -= dummy_b;
		}
		else 
		{
			// first define all the points in the same relation to the center of cell number 1
			Point barycenter_rel = shift(BC_l,(q_BC_l-q_c1)*T->SystemSize);
			Point v1_rel		 = shift(v1->getCoord(), q_c1*T->SystemSize*(-1.0));
			Point v2_rel		 = shift( v2->getCoord(), (q_v2-q_c1)*T->SystemSize);
			Point v1b_rel		 = shift( v1->getBasalCoord(), (v1->getBasalQuadrant()-q_c1)*T->SystemSize);
			Point v2b_rel		 = shift( v2->getBasalCoord(), (q_v2+v2->getBasalQuadrant()-q_c1)*T->SystemSize);
			Point barycenter_c_rel = shift( c1->BC_b, c1->q_BC_b*T->SystemSize);
			           			
			// add the force contributions of the barycenter of cell 1
			c1->dAa_dBCa += Point::TriangleAreaChange(c1->BC_a,v1_rel,v2_rel);
			c1->dAb_dBCb  += Point::TriangleAreaChange(barycenter_c_rel,v1b_rel,v2b_rel);
			c1->dV_dBCa  += Point::dV_dp(v1_rel, v2_rel);
			c1->dV_dBCb   += Point::dV_dp(v2b_rel, v1b_rel);
			
			
			// define all the points in the same relation to the center of cell number 2
			barycenter_rel = shift( BC_l, T->SystemSize*(q_BC_l-q_c2-q_v2) );
			v1_rel		 = shift( v1->getCoord(), (q_c2+q_v2)*T->SystemSize*(-1));
			v2_rel		 = shift( v2->getCoord(), q_c2*T->SystemSize*(-1));
			v1b_rel		 = shift( v1->getBasalCoord(), (v1->getBasalQuadrant()-q_c2 - q_v2)*T->SystemSize);
			v2b_rel		 = shift( v2->getBasalCoord(), (v2->getBasalQuadrant()-q_c2)*T->SystemSize);
			barycenter_c_rel = shift( c2->BC_b, c2->q_BC_b*T->SystemSize);
			
			// add the force contributions of the barycenter of cell 1
			c2->dAa_dBCa += Point::TriangleAreaChange(c2->BC_a,v1_rel,v2_rel);
			c2->dAb_dBCb += Point::TriangleAreaChange(barycenter_c_rel,v1b_rel,v2b_rel);
			c2->dV_dBCa  -= Point::dV_dp(v1_rel, v2_rel);
			c2->dV_dBCb  -= Point::dV_dp(v2b_rel, v1b_rel);
        }
	}
}

	
int Edge::VertexIsPart(Vertex *v)
{
	if(v==v1)
	{
		return 1;
	}
	else if(v==v2)
	{
		return 2;
	}
	else 
	{
		return 0;
	}
};

void Edge::Print()
{
	//v1->Print();
	//%v2->Print();
	//std::cout << std::endl<< "Edge:" << Number << ", T_l:" << T_l << ", LineTensions:" << G_a << "," << G_b << std::endl;
	//std::cout << "v1=" << v1->Number << ", v2=" << v2->Number << ", c1=" << c1->Number<< ", c2=" << c2->Number;
};

Point Edge::MP_a()
{

    //return (v1->coord+v2->coord)/2;

    if(!T->isPeriodic || !CrossBoundary())
    {
        return (v1->coord+v2->coord)/2;
    }
    else
	{
		return Quadrant::Modulus(shift( v1->coord+v2->coord, q_v2*T->SystemSize)/2, T->SystemSize);
	}

}

Point Edge::MP_a(Cell *c)
{
    if(!T->isPeriodic)
    {
        return ((v1->coord+v2->coord)/2);
    }
    
    if(c == c1)
    {
        return shift(v1->coord+v2->coord, (q_c1*2 + q_v2)*T->SystemSize) * 0.5;
    }
    else if(c == c2)
    {
        return shift(v1->coord+v2->coord, (q_c2*2 - q_v2)*T->SystemSize) * 0.5;
    }
    else
    {
        std::cout << "Error in function Edge::MP_a" << std::endl;
        throw("Error in function Edge::MP_a");
    }    
}

Point Edge::MP_b()
{
    if(!T->isPeriodic || !crossBoundary)
    {
        return ((v1->coord_b+v2->coord_b)/2);
    }
   	else
	{
		return Quadrant::Modulus(shift( v1->getBasalCoord()+v2->getBasalCoord(), (q_v2+v2->getBasalQuadrant()* T->SystemSize))/2, T->SystemSize);
	}
}

Point Edge::MP_b(Cell *c)
{
    if(!T->isPeriodic)
    {
        return ((v1->coord_b+v2->coord_b)/2);
    }
    
    if(c == c1)
    {
        return shift(v1->coord_b+v2->coord_b, (q_c1*2+v1->q_b+q_v2+v2->q_b)*T->SystemSize) * 0.5;
    }
    else if(c == c2)
    {
        return shift(v1->coord_b+v2->coord_b, (q_c2*2+v1->q_b-q_v2+v2->q_b)*T->SystemSize) * 0.5;
    }
    else
    {
        std::cout << "Error in function Edge::MP_b" << std::endl;
        throw("Error in function Edge::MP_b");
    }
}
					
bool Edge::CrossBoundary()
{
	crossBoundary = !(q_v2.isZero() & !v1->CrossBoundary() & !v2->CrossBoundary());
    return crossBoundary;
}



// function that removes an edge from the tissue - up to now just for the nonperiodic case!!
Vertex* Edge::T1_remove()
{	
	//version if e does NOT cross a boundary
	if(!T->isPeriodic || !CrossBoundary())
	{
		// move v1 to the middle between v1 and v2
		v1->coord = (v1->coord+v2->coord)/2;
		v1->coord_b = (v1->coord_b+v2->coord_b)/2;
		v1->updateNeighbouringEdges();
        v2->updateNeighbouringEdges();
        
        
		// iterate through all neighbouring edges of v2 and tell them to point at v1
		for(std::list<Edge*>::iterator iter_e = v2->NeighbouringEdges_sorted.begin(); iter_e!=v2->NeighbouringEdges_sorted.end();++iter_e)
		{
			Edge* e = (*iter_e);
            
            if(e==this) continue;
            
            if(e->v1 == v2)
            {
                e->v1 = v1;
            }
            else
            {
                e->v2 = v1;
            }
		}
        		 
		// erase e and  v2 from the lists of the surrounding edges and vertices of  c1 and  c2				
		if(c1!=NULL)
		{
			// erase v2 from c1
			c1->ListEdges.remove(this);
		}
		if(c2!=NULL)
		{
			// erase v2 from c2
			c2->ListEdges.remove(this);
		}
		
        // remove this edge and the vertex v2 from the list in Tissue
		T->ListEdge.remove(*this);
        T->ListVertex.remove(*v2);
	}
	
	// implement the same for the periodic case!!
	else 
	{
		
	}

	T->UpdateBarycenters();
    
	return v1;
}

Vertex* Edge::T1_remove_v2()
{
	//version if e does NOT cross a boundary
	if(!T->isPeriodic || !CrossBoundary())
	{
		// move v1 to the middle between v1 and v2
		v2->coord = (v1->coord+v2->coord)/2;
		v2->coord_b = (v1->coord_b+v2->coord_b)/2;
		//v1->updateNeighbouringEdges();
        //v2->updateNeighbouringEdges();
        
        
		// iterate through all neighbouring edges of v1 and tell them to point at v2
		for(std::list<Edge*>::iterator iter_e = v1->NeighbouringEdges_sorted.begin(); iter_e!=v1->NeighbouringEdges_sorted.end();++iter_e)
		{
			Edge* e = (*iter_e);
            
            if(e==this) continue;
            
            if(e->v1 == v1)
            {
                e->v1 = v2;
            }
            else
            {
                e->v2 = v2;
            }
		}
        
		// erase e and  v2 from the lists of the surrounding edges and vertices of  c1 and  c2
		if(c1!=NULL)
		{
			// erase v2 from c1
			c1->ListEdges.remove(this);
		}
		if(c2!=NULL)
		{
			// erase v2 from c2
			c2->ListEdges.remove(this);
		}
		
        // remove this edge and the vertex v2 from the list in Tissue
		T->ListEdge.remove(*this);
        T->ListVertex.remove(*v1);
	}
	
	// implement the same for the periodic case!!
	else
	{
		
	}
    
	T->UpdateBarycenters();
    
	return v2;
}

std::list<Edge*>* Edge::neighbouringEdges()
{
    // create the list that will contain all the neighbouring edges
    std::list<Edge*> *NEdges = new std::list<Edge*>();
    
    // iterator through all the edges of the tissue and add them to the list if they have one vertex in common
    for(std::list<Edge>::iterator it_edge = T->ListEdge.begin();it_edge!=T->ListEdge.end();it_edge++)
    {
        Edge *e = &(*it_edge);
        if(e->v1 == v1 || e->v1 == v2 || e->v2 == v1 || e->v2 == v2) NEdges->push_back(e);
    }
    
    return NEdges;
    
}

/// checks if the edge should be removed from the tissue and removes it if neccessary
bool Edge::CheckAndRemoveFromTissue()
{
	// the tissue should be updated
	
	// direction of force on apical edge
	double directiondW_dva = (v2->getCoord() - v1->getCoord()).dot(v2->getdW_dva() - v1->getdW_dva());
	
    // check if my criterion makes sense!
    double TINY = 10e-10;
    double NormBefore = Point::Norm(v2->coord-v1->coord);
    double NormAfter = Point::Norm((v2->coord-v2->dW_dva*TINY)-(v1->coord-v1->dW_dva*TINY));
    bool popTest = (NormBefore-NormAfter)>0;
    bool pop = directiondW_dva<0;
    
    if(!popTest==pop) std::cout << "WRONG CALCULATION!!" << std::endl;
    
	// direction of force on basal edge
	double directiondW_dvb = (v2->getBasalCoord() - v1->getBasalCoord()).dot(v2->getdW_dvb() - v1->getdW_dvb());

	// the edge will be removed if both sides "want" to shrink
	if(directiondW_dva<0 || directiondW_dvb<0)
	{
		T1_remove();
		return true;
	}
	
	return false;
}

bool Edge::isAtBoundary()
{
    if(!T->isPeriodic)
    {
        // on one side no neighbour
        return c1==NULL || c2==NULL;
    }
    else
    {
        // the barycenters of the two neighbours are not in the same quadrant
        return !(q_c1 == q_v2+q_c2);
    }
}


//void Edge::getXcrossing(double xCross, Point* ApicalCross, bool* doesApicalCross, Point* BasalCross, bool* doesBasalCross, Point* SecondApicalCross, bool* isApicallyAligned, Point* SecondBasalCross, bool* isBasallyAligned)
//{
//	Point P1 = getVertex1()->getCoord(); Point P1b = getVertex1()->getBasalCoord();
//	Point P2 = getVertex2()->getCoord(); Point P2b = getVertex2()->getBasalCoord();
//	
//	// check if the apical side of the edge lies in the y-z-plane
//	if(abs(P1.x-xCross)<EPSILON && abs(P2.x-xCross)<EPSILON)
//	{
//		*ApicalCross = P1;
//		*SecondApicalCross = P2;
//		*doesApicalCross = true;
//		*isApicallyAligned = true;
//	}
//
//	// check if the basal side of the edge lies in the y-z-plane
//	if(abs(P1b.x-xCross)<EPSILON && abs(P2b.x-xCross)<EPSILON)
//	{
//		*BasalCross = P1b;
//		*SecondBasalCross = P2b;
//		*doesBasalCross = true;
//		*isBasallyAligned = true;
//	}
//	
//	double CrossCheckApical = (P1.x-xCross)*(P2.x-xCross);
//	double CrossCheckBasal = (P1b.x-xCross)*(P2b.x-xCross);
//	
//	if(CrossCheckApical>EPSILON)
//	{
//		*doesApicalCross = false;
//		*ApicalCross = Point();
//	}
//	else if(CrossCheckApical<-EPSILON)
//	{
//		double ratio = (xCross - P1.x)/(P2.x-P1.x);
//		*ApicalCross=Point(xCross, P1.y + ratio * (P2.y-P1.y), P1.z + ratio * (P2.z-P1.z));
//		*doesApicalCross = true;
//		*isApicallyAligned = false;
//		
//	}
//	
//
//	if(CrossCheckBasal>=EPSILON)
//	{
//		*doesBasalCross = false;
//		*BasalCross = Point();
//	}
//	
//	else if (CrossCheckBasal<=-EPSILON) 
//	{
//		double ratio = (xCross - P1b.x)/(P2b.x-P1b.x);
//		*BasalCross=Point(xCross, P1b.y + ratio * (P2b.y-P1b.y), P1b.z + ratio * (P2b.z-P1b.z));
//		*doesBasalCross = true;
//	}
//	
//}
//	
//void Edge::getYcrossing(double yCross, Point* ApicalCross, bool* doesApicalCross, Point* BasalCross, bool* doesBasalCross, Point* SecondApicalCross, bool* isApicallyAligned, Point* SecondBasalCross, bool* isBasallyAligned)
//{
//	Point P1 = getVertex1()->getCoord(); Point P1b = getVertex1()->getBasalCoord();
//	Point P2 = getVertex2()->getCoord(); Point P2b = getVertex2()->getBasalCoord();
//	
//	// check if the apical side of the edge lies in the y-z-plane
//	if(abs(P1.y-yCross)<EPSILON && abs(P2.y-yCross)<EPSILON)
//	{
//		*ApicalCross = P1;
//		*SecondApicalCross = P2;
//		*doesApicalCross = true;
//		*isApicallyAligned = true;
//	}
//	
//	// check if the basal side of the edge lies in the y-z-plane
//	if(abs(P1b.y-yCross)<EPSILON && abs(P2b.y-yCross)<EPSILON)
//	{
//		*BasalCross = P1b;
//		*SecondBasalCross = P2b;
//		*doesBasalCross = true;
//		*isBasallyAligned = true;
//	}	
//	
//	double CrossCheckApical = (P1.y-yCross)*(P2.y-yCross);
//	double CrossCheckBasal = (P1b.y-yCross)*(P2b.y-yCross);
//	
//	
//	if(CrossCheckApical>EPSILON)
//	{
//		*doesApicalCross = false;
//		*ApicalCross = Point();
//		*isApicallyAligned = false;
//	}
//	else if(CrossCheckApical<-EPSILON)
//	{
//		double ratio = (yCross - P1.y)/(P2.y-P1.y);
//		*ApicalCross=Point( P1.x + ratio * (P2.x-P1.x),yCross, P1.z + ratio * (P2.z-P1.z));
//		*doesApicalCross = true;
//		*isApicallyAligned = false;
//	}
//	
//	if(CrossCheckBasal>EPSILON)
//	{
//		*doesBasalCross = false;
//		*BasalCross = Point();
//	}
//	else if(CrossCheckBasal<-EPSILON)
//	{
//		double ratio = (yCross - P1b.y)/(P2b.y-P1b.y);
//		*BasalCross=Point( P1b.x + ratio * (P2b.x-P1b.x),yCross, P1b.z + ratio * (P2b.z-P1b.z));
//		*doesBasalCross = true;
//		*isApicallyAligned = false;
//	}
//	
//}
//
//void Edge::updatePlaneSection() //Point positionVectorPlane, Point normalVectorPlane);
//{
//	// get the position of the edge after shifting it by the position vector of the plane (needs to be shifted back!!)
//	Point P1 = getVertex1()->getCoord() - T->SectionPlane_PositionVector;
//	Point P1b = getVertex1()->getBasalCoord() - T->SectionPlane_PositionVector;
//	Point P2 = getVertex2()->getCoord() - T->SectionPlane_PositionVector; 
//	Point P2b = getVertex2()->getBasalCoord() - T->SectionPlane_PositionVector;
//
//	/* ********* start with the apical side (use notation from http://en.wikipedia.org/wiki/Line-plane_intersection )  ******* */
//	
//	// l = direction vector of the apical edge
//	Point l = P2 - P1;
//	
//	double enumerator = P1.dot(T->SectionPlane_NormalVector);
//	double denominator = l.dot(T->SectionPlane_NormalVector);
//	
//	if(denominator==0 && enumerator!=0)
//	{
//		// line is parallel to the plane but not inside the plane
//		planeCrosses_apical = false;
//	}
//	else if(denominator==0 && enumerator==0)
//	{
//		// line lies inside the plane -> dont contribute any points since they will be given by the two vertices v1 and v2 whose apical sides both lie in the plane
//		planeCrosses_apical = false;		
//	}
//	else 
//	{
//		double d = -enumerator/denominator;
//		if(d>=0 && d<=1)
//		{
//			planeSection_apical = P1 + l*d + T->SectionPlane_PositionVector;
//			planeCrosses_apical = true;
//		}
//		else
//		{
//			planeCrosses_apical = false;
//		}
//	}
//	
//	/* ********* now do section with the basal side ******* */
//
//	l = P2b - P1b;
//	
//	enumerator = P1b.dot(T->SectionPlane_NormalVector);
//	denominator = l.dot(T->SectionPlane_NormalVector);
//	
//	if(denominator==0 && enumerator!=0)
//	{
//		// line is parallel to the plane but not inside the plane
//		planeCrosses_basal = false;
//	}
//	else if(denominator==0 && enumerator==0)
//	{
//		// line lies inside the plane -> dont contribute any points since they will be given by the two vertices v1 and v2 whose apical sides both lie in the plane
//		planeCrosses_basal = false;		
//	}
//	else 
//	{
//		double d = -enumerator/denominator;
//		if(d>=0 && d<=1)
//		{
//			planeSection_basal = P1b + l*d + T->SectionPlane_PositionVector;
//			planeCrosses_basal = true;
//		}
//		else
//		{
//			planeCrosses_basal = false;
//		}
//	}
//			
//}


// calculate the force coming from Cell c that acts towards a T2 transition of Vertex v 
Point Edge::Force_PressureT2Transition(Cell *c, Vertex *v)
{
	int directionEdge;
	
	if(getCell1() == c )
	{
		directionEdge = 1;
	}
	else 
	{
		directionEdge = -1;
	}
	
	int noVertex = VertexIsPart(v);
	
	
	Point A,B;
	Point C = v->getCoord();
	Point D = v->getBasalCoord();
	
	if(noVertex==1)
	{
		A = v2->getCoord();
		B = v2->getBasalCoord();
	}
	else
	{
		A = v1->getCoord();
		B = v1->getBasalCoord();
	}	
	
	// the volume increases by ((A-D)*(C-B) + (B+D)*(A+C)*4 )/12
	return  ( (A-D)*(C-B) + (B+D)*(A+C)*4 )*(directionEdge*c->Pressure /12) ; 
	

}

// calculate the force coming from Cell c that acts toward a T2 transition of Vertex v 
Point Edge::Force_SurfaceTensionT2Transition(Cell *c, Vertex *v, Point shift)
{
	int vertexNo = VertexIsPart(v);
	
	Point v_edge_a,v_edge_b;
	
	if(vertexNo==1)
	{
		v_edge_a = v2->getCoord() - v1->getCoord();
		v_edge_b = v2->getBasalCoord() - v1->getBasalCoord();
	}
	else if(vertexNo==2)
	{
		v_edge_a = v1->getCoord() - v2->getCoord();
		v_edge_b = v1->getBasalCoord() - v2->getBasalCoord();
	}
	else
	{
		throw;
	}
	
	// normalize:
	v_edge_a = v_edge_a/Point::Norm(v_edge_a);
	v_edge_b = v_edge_b/Point::Norm(v_edge_b);
	
	// project the shift on the orthogonal space of the edge as only that contribution changes the area of the triangle
	Point NormalToEdge_a = shift - v_edge_a*shift.dot(v_edge_a);
	Point NormalToEdge_b = shift - v_edge_b*shift.dot(v_edge_b);
	
	return NormalToEdge_a*(c->T_a()/2) + NormalToEdge_b*(c->T_b()/2);
	
}

//void Edge::UpdateForcesAndEnergy()
//{
//    if(!T->isPeriodic) {
//        // SET THE FORCE CONTRIBUTIONS FROM THE EDGE ON EVERY OF ITS VERTICES
//        Point Fa_v1;
//        Point Fb_v1;
//        Point Fa_v2;
//        Point Fb_v2;
//        
//        // positions
//        Point v1a = v1->coord;
//        Point v1b = v1->coord_b;
//        Point v2a = v2->coord;
//        Point v2b = v2->coord_b;
//        Point b = BC_l;
//        
//        // APICAL AND BASAL LINE TENSION
//        // force on vertices
//        Fa_v1 = (v2a-v1a)*(Lambda_a()/l_a());
//        Fa_v2 = Fa_v1*(-1);
//
//        Fb_v1 = (v2b-v1b)*(Lambda_b()/lb);
//        Fb_v2 = Fb_v1*(-1);
//        
//        // energy contribution
//        T->Energy_Line+= Lambda_a()*la + Lambda_b()*lb;
//        
//        
//        // LATEREAL SURFACE TENSION
//        // force on vertices
//        Fa_v1 += ForceFromBarycenter_Area - (Point::TriangleAreaChange(v1a,v2a,b)+Point::TriangleAreaChange(v1a,v1b,b))*T_l;
//        Fa_v2 += ForceFromBarycenter_Area - (Point::TriangleAreaChange(v2a,v1a,b)+Point::TriangleAreaChange(v2a,v2b,b))*T_l;
//        Fb_v1 += ForceFromBarycenter_Area - (Point::TriangleAreaChange(v1b,v2b,b)+Point::TriangleAreaChange(v1b,v1a,b))*T_l;
//        Fb_v2 += ForceFromBarycenter_Area - (Point::TriangleAreaChange(v2b,v1b,b)+Point::TriangleAreaChange(v2b,v2a,b))*T_l;
//        
//        // energy contribution
//        T->Energy_Surface +=  Point::TriangleArea(v1a,v2a,c1->BC_a)*c1->T_a() +  Point::TriangleArea(v1b,v2b,c1->BC_b)*c1->T_b() + Point::TriangleArea(v1a,v2a,c2->BC_a)*c2->T_a() +  Point::TriangleArea(v1b,v2b,c2->BC_b)*c2->T_b();
//        
//        // FORCE FROM VOLUME PRESSURE
//        // save some values to avoid the repeated calculation
//        Point b_v1a = Point::dV_dp(b,v1a);
//        Point b_v1b = Point::dV_dp(b,v1b);
//        Point b_v2a = Point::dV_dp(b,v2a);
//        Point b_v2b = Point::dV_dp(b,v2b);
//        
//        if(c1)
//        {
//            Point ForceFromLateralBarycenter_Volume = dV_dBC*(-c1->Pressure);
//            
//            (c1->ForceFromApicalBarycenter_Volume + ForceFromLateralBarycenter_Volume).Print();
//          
//            
//            
//            //// FALSE:     ... Fa_v1 += c1->dAa_dBCa + ...   !?
//            
//            
//            Fa_v1 += c1->dAa_dBCa + c1->ForceFromApicalBarycenter_Volume + ForceFromLateralBarycenter_Volume
//            - Point::TriangleAreaChange(v1a,v2a,c1->BC_a)*c1->T_a() - (b_v2a-b_v1b+Point::dV_dp(v2a,c1->BC_a))*c1->Pressure;
//            
//            Fa_v2 += ForceFromLateralBarycenter_Volume - Point::TriangleAreaChange(v2a,v1a,c1->BC_a)*c1->T_a() - (b_v2b-b_v1a+Point::dV_dp(c1->BC_a,v1a))*c1->Pressure;
//            
//            Fb_v1 += c1->dAb_dBCb + c1->ForceFromBasalBarycenter_Volume + ForceFromLateralBarycenter_Volume
//            - Point::TriangleAreaChange(v1b,v2b,c1->BC_b)*c1->T_b() - (b_v1a-b_v2b+Point::dV_dp(c1->BC_b,v2b))*c1->Pressure;
//            
//            
//             //// FALSE:     ... Fa_v1 += c1->dAa_dBCa + ...   !?
//            
//            
//            Fb_v2 += ForceFromLateralBarycenter_Volume - Point::TriangleAreaChange(v2b,v1b,c1->BC_b)*c1->T_b() - (b_v1b-b_v2a+Point::dV_dp(v1b,c1->BC_b))*c1->Pressure;
//            
//            T->Energy_Surface +=  Point::TriangleArea(v1a,v2a,c1->BC_a)*c1->T_a() +  Point::TriangleArea(v1b,v2b,c1->BC_b)*c1->T_b();
//            
//        }
//        
//        if(c2)
//        {
//            Point ForceFromLateralBarycenter_Volume = dV_dBC*(c2->Pressure);
//            
//            Fa_v1 +=  ForceFromLateralBarycenter_Volume - Point::TriangleAreaChange(v1a,v2a,c1->BC_a)*c2->T_a() + c2->ForceFromApicalBarycenter_Volume +  (b_v2a-b_v1b+Point::dV_dp(v2a,c2->BC_a))*c2->Pressure;
//            
//            Fa_v2 +=  c2->dAa_dBCa + c2->ForceFromApicalBarycenter_Volume + ForceFromLateralBarycenter_Volume
//            - Point::TriangleAreaChange(v2a,v1a,c1->BC_a)*c2->T_a() + (b_v2b-b_v1a+Point::dV_dp(c2->BC_a,v1a))*c2->Pressure;
//            
//            Fb_v1 +=  ForceFromLateralBarycenter_Volume - Point::TriangleAreaChange(v1b,v2b,c2->BC_b)*c2->T_b() + (b_v1a-b_v2b+Point::dV_dp(c2->BC_b,v2b))*c2->Pressure;
//            
//            Fb_v2 += c2->dAb_dBCb + c2->ForceFromBasalBarycenter_Volume  + ForceFromLateralBarycenter_Volume
//            - Point::TriangleAreaChange(v2b,v1b,c2->BC_b)*c2->T_b() + (b_v1b-b_v2a+Point::dV_dp(v1b,c2->BC_b))*c2->Pressure;
//            
//            T->Energy_Surface += Point::TriangleArea(v1a,v2a,c2->BC_a)*c2->T_a() +  Point::TriangleArea(v1b,v2b,c2->BC_b)*c2->T_b();
//        }
//        
//        // ADD UP THE FORCE CONTRIBUTIONS IN THE 4 VERTEX POSITIONS THAT BELONG TO THE EDGE
//        v1->dW_dva += Fa_v1;
//        v2->dW_dva += Fa_v2;
//        v1->dW_dvb  += Fb_v1;
//        v2->dW_dvb  += Fb_v2;
//    }
//    
//    else { // periodic tissue
//        
//        // FORCE CONTRIBUTIONS FROM THE EDGE ON EVERY OF ITS VERTICES V1a,V1b,V2a,V2b
//        Point Fa_v1, Fb_v1, Fa_v2, Fb_v2, Fa_v1_b, Fb_v1_b, Fa_v2_b, Fb_v2_b;
//        
//        la = l_a();
//        lb = l_b();
//        
//        // positions
//        Point v1a = v1->coord;
//        Point v1b = v1->coord_b;
//        Point v2a = v2->coord;
//        Point v2b = v2->coord_b;
//        Point b = BC_l;
//        
//        //
//        Point dVc1_dV1a, dVc2_dV1a, dVc1_dV1b, dVc2_dV1b, dVc1_dV2a, dVc2_dV2a, dVc1_dV2b, dVc2_dV2b;
//        
//        // pure contributions from apical and basal barycenter, and from the lateral barycenter
//        int n1 = c1->ListEdges.size();
//        int n2 = c2->ListEdges.size();
//        
//        // NOTE: every vertex is once vertex one with respect to a cell
//        // therefore add the "force contribution" coming from the movement in the apical and basal barycenters of v1 wrt c1 and of v2 wrt to c2 (since v1 is first in direction of c1 and v2 to c2)
//        v1->dW_dva -= c1->dAa_dBCa*(c1->T_a() / c1->ListEdges.size());
//        v2->dW_dva -= c2->dAa_dBCa*(c2->T_a() / c2->ListEdges.size());
//        v1->dW_dvb -= c1->dAb_dBCb*(c1->T_b() / c1->ListEdges.size());
//        v2->dW_dvb -= c2->dAb_dBCb*(c2->T_b() / c2->ListEdges.size());
//        
//        dVc1_dV1a += c1->dV_dBCa/n1;
//        dVc1_dV1b += c1->dV_dBCb/n1;
//        dVc2_dV2a += c2->dV_dBCa/n2;
//        dVc2_dV2b += c2->dV_dBCb/n2;
//        
//        // systemSize
//        Quadrant SystemSize = T->SystemSize;
//        
//        bool crossBoundary_c1 = c1->CrossBoundary();
//        bool crossBoundary_c2 = c2->CrossBoundary();
//        
//        if(!CrossBoundary()) { // edge doesnt cross the boundary
//    
//            // APICAL AND BASAL LINE TENSION
//            // force on vertices
//            Fa_v1 += (v2a-v1a)*(Lambda_a()/la);
//            Fa_v2 += Fa_v1*(-1);
//
//            Fb_v1 += (v2b-v1b)*(Lambda_b()/lb);
//            Fb_v2 += Fb_v1*(-1);
//
//        
//            // force on the system size
//            T->dW_dL -= Point::lineStress(v1a,v2a)*Lambda_a() + Point::lineStress(v1b, v2b)*Lambda_b();
//            // energy contribution            
//            if(MIN_LENGTH<la) T->Energy_Line += Lambda_a()*la;
//            if(MIN_LENGTH<lb) T->Energy_Line += Lambda_b()*lb;
//            
//            
//            // LATEREAL SURFACE TENSION
//            // force on vertices
//            Fa_v1 -= (Point::TriangleAreaChange2(v1a,v2a,b,4)+Point::TriangleAreaChange2(v1a,v1b,b,4)+Point::TriangleAreaChange(b,v1b,v2b)/4+Point::TriangleAreaChange(b,v2a,v2b)/4)*T_l;
//            Fa_v2 -= (Point::TriangleAreaChange2(v2a,v1a,b,4)+Point::TriangleAreaChange2(v2a,v2b,b,4)+Point::TriangleAreaChange(b,v1a,v1b)/4+Point::TriangleAreaChange(b,v1b,v2b)/4)*T_l;
//            Fb_v1 -= (Point::TriangleAreaChange2(v1b,v1a,b,4)+Point::TriangleAreaChange2(v1b,v2b,b,4)+Point::TriangleAreaChange(b,v1a,v2a)/4+Point::TriangleAreaChange(b,v2a,v2b)/4)*T_l;
//            Fb_v2 -= (Point::TriangleAreaChange2(v2b,v2a,b,4)+Point::TriangleAreaChange2(v2b,v1b,b,4)+Point::TriangleAreaChange(b,v1a,v2a)/4+Point::TriangleAreaChange(b,v1a,v1b)/4)*T_l;
//            // force on the system size
//            T->dW_dL -= (Point::surfaceStress(v1a,v2a,b) + Point::surfaceStress(v2a,v2b,b ) + Point::surfaceStress(v2b,v1b,b ) + Point::surfaceStress(v1b,v1a,b ))* T_l;
//            // energy contribution
//            T->Energy_Surface += T_l*(Point::TriangleArea(v1a,v2a,b)+Point::TriangleArea(v2a,v2b,b)+Point::TriangleArea(v2b,v1b,b)+Point::TriangleArea(v1b,v1a,b));
//        
//        } else { // edge crosses boundary
//            // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE FIRST VERTEX OF THE EDGE
//            Point v1a_rel = v1a;
//            Point v1b_rel = shift( v1b, SystemSize*v1->q_b);
//            Point v2a_rel = shift( v2a, SystemSize*q_v2);
//            Point v2b_rel = shift( v2b, SystemSize*(q_v2+v2->q_b));
//            Point b_rel =   shift( b, SystemSize*q_BC_l);
//            
//            // APICAL AND BASAL LINE TENSION
//            // force on vertices
//            Fa_v1 += (v2a_rel-v1a_rel)*(Lambda_a()/la);
//            Fa_v2 += Fa_v1*(-1);
//            
//            Fb_v1 += (v2b_rel-v1b_rel)*(Lambda_b()/lb);
//            Fb_v2 += Fb_v1*(-1);
//
//            // force on the system size
//            T->dW_dL -= Point::lineStress(v1a_rel,v2a_rel)*Lambda_a() + Point::lineStress(v1b_rel, v2b_rel)*G_b;
//            // energy contribution
//            if(MIN_LENGTH<la) T->Energy_Line += Lambda_a()*la;
//            if(MIN_LENGTH<lb) T->Energy_Line += G_b*lb;
//        
//            // LATERAL SURFACE TENSION
//            // force on vertices
//            Fa_v1 -= (Point::TriangleAreaChange2(v1a_rel,v2a_rel,b_rel,4)+Point::TriangleAreaChange2(v1a_rel,v1b_rel,b_rel,4)+Point::TriangleAreaChange(b_rel,v1b_rel,v2b_rel)/4+Point::TriangleAreaChange(b_rel,v2a_rel,v2b_rel)/4)*T_l;
//            Fa_v2 -= (Point::TriangleAreaChange2(v2a_rel,v1a_rel,b_rel,4)+Point::TriangleAreaChange2(v2a_rel,v2b_rel,b_rel,4)+Point::TriangleAreaChange(b_rel,v1a_rel,v1b_rel)/4+Point::TriangleAreaChange(b_rel,v1b_rel,v2b_rel)/4)*T_l;
//            Fb_v1 -= (Point::TriangleAreaChange2(v1b_rel,v1a_rel,b_rel,4)+Point::TriangleAreaChange2(v1b_rel,v2b_rel,b_rel,4)+Point::TriangleAreaChange(b_rel,v1a_rel,v2a_rel)/4+Point::TriangleAreaChange(b_rel,v2a_rel,v2b_rel)/4)*T_l;
//            Fb_v2 -= (Point::TriangleAreaChange2(v2b_rel,v2a_rel,b_rel,4)+Point::TriangleAreaChange2(v2b_rel,v1b_rel,b_rel,4)+Point::TriangleAreaChange(b_rel,v1a_rel,v2a_rel)/4+Point::TriangleAreaChange(b_rel,v1a_rel,v1b_rel)/4)*T_l;
//            
//            // force on the system size
//            T->dW_dL -= (Point::surfaceStress(v1a_rel,v2a_rel,b_rel ) + Point::surfaceStress(v2a_rel,v2b_rel,b_rel ) + Point::surfaceStress(v2b_rel,v1b_rel,b_rel ) + Point::surfaceStress(v1b_rel,v1a_rel,b_rel ))*  T_l;
//            // energy contribution
//            T->Energy_Surface += T_l*(Point::TriangleArea(v1a_rel,v2a_rel,b_rel)+Point::TriangleArea(v2a_rel,v2b_rel,b_rel)+Point::TriangleArea(v2b_rel,v1b_rel,b_rel)+Point::TriangleArea(v1b_rel,v1a_rel,b_rel));
//        }
//    
//        if(!crossBoundary_c1)  {
//            // APICAL AND BASAL SURFACE TENSIONS
//            Fa_v1 -= Point::TriangleAreaChange3(v1a,v2a,c1->BC_a,n1) * c1->T_a();
//            Fb_v1 -= Point::TriangleAreaChange3(v1b,v2b,c1->BC_b,n1) * c1->T_b();
//            Fa_v2 -= Point::TriangleAreaChange3(v2a,v1a,c1->BC_a,n1) * c1->T_a();
//            Fb_v2 -= Point::TriangleAreaChange3(v2b,v1b,c1->BC_b,n1) * c1->T_b();
//            // force on system size
//            T->dW_dL -= (Point::surfaceStress(v1a, c1->BC_a, v2a)) * c1->T_a() + (Point::surfaceStress(v1b, c1->BC_b, v2b)) * c1->T_b();
//            // energy contribution
//            T->Energy_Surface +=  Point::TriangleArea(v1a,v2a,c1->BC_a)*c1->T_a() +  Point::TriangleArea(v1b,v2b,c1->BC_b)*c1->T_b();
//
//            // lateral contributions to volume changes
//            dVc1_dV1a += (Point::dV_dp(BC_l, v2->coord_b, v2->coord)+Point::dV_dp(BC_l, v1->coord_b, v2->coord_b))/4 + Point::dV_dp_b(v1->coord, v1->coord_b, BC_l, 4) - Point::dV_dp_b(v1->coord, v2->coord, BC_l,4);
//            dVc1_dV1b += (Point::dV_dp(BC_l, v2->coord, v1->coord)+Point::dV_dp(BC_l, v2->coord_b, v2->coord))/4 + Point::dV_dp_b(v1->coord_b, v2->coord_b, BC_l, 4) - Point::dV_dp_b(v1->coord_b, v1->coord, BC_l,4);
//            dVc1_dV2a += (Point::dV_dp(BC_l, v1->coord_b, v2->coord_b)+Point::dV_dp(BC_l, v1->coord, v1->coord_b))/4 + Point::dV_dp_b(v2->coord, v1->coord, BC_l, 4) - Point::dV_dp_b(v2->coord, v2->coord_b, BC_l,4);
//            dVc1_dV2b += (Point::dV_dp(BC_l, v1->coord, v1->coord_b)+Point::dV_dp(BC_l, v2->coord, v1->coord))/4 + Point::dV_dp_b(v2->coord_b, v2->coord, BC_l, 4) - Point::dV_dp_b(v2->coord_b, v1->coord_b, BC_l,4);
//            
//            // apical and basal contributions to volume changes
//            dVc1_dV1a += Point::dV_dp_b_1(v1->coord, v2->coord, c1->BC_a, n1);
//            dVc1_dV1b -= Point::dV_dp_b_1(v1->coord_b, v2->coord_b, c1->BC_b, n1);
//            dVc1_dV2a -= Point::dV_dp_b_1(v2->coord, v1->coord, c1->BC_a, n1);
//            dVc1_dV2b += Point::dV_dp_b_1(v2->coord_b, v1->coord_b, c1->BC_b, n1);
//            
//        } else {
//            
//            // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE APICAL BARYCENTER OF CELL 1
//            Point BC_l_rel = shift( b, SystemSize*( q_BC_l -  q_c1));
//            Point v1a_rel = shift( v1a, SystemSize* q_c1*(-1));
//            Point v1b_rel = shift( v1b, SystemSize*(v1->q_b -  q_c1));
//            Point v2a_rel = shift( v2a, SystemSize*( q_v2 -  q_c1));
//            Point v2b_rel = shift( v2b, SystemSize*( q_v2 + v2->q_b -  q_c1));
//            Point BC_b_rel = shift( c1->BC_b, SystemSize* c1->q_BC_b);
//            
//            // APICAL AND BASAL SURFACE TENSION
//            // force on vertices
//            Fa_v1 -= Point::TriangleAreaChange3(v1a_rel,v2a_rel,c1->BC_a,n1) * c1->T_a();
//            Fb_v1 -= Point::TriangleAreaChange3(v1b_rel,v2b_rel,BC_b_rel,n1) * c1->T_b();
//            Fa_v2 -= Point::TriangleAreaChange3(v2a_rel,v1a_rel,c1->BC_a,n1) * c1->T_a();
//            Fb_v2 -= Point::TriangleAreaChange3(v2b_rel,v1b_rel,BC_b_rel,n1) * c1->T_b();
//            // force on system size
//            T->dW_dL -= (Point::surfaceStress(v1a_rel, c1->BC_a, v2a_rel)) * c1->T_a() + (Point::surfaceStress(v1b_rel, BC_b_rel, v2b_rel)) * c1->T_b();
//            // energy contribution
//            T->Energy_Surface +=  Point::TriangleArea(v1a_rel,v2a_rel,c1->BC_a)*c1->T_a() +  Point::TriangleArea(v1b_rel,v2b_rel,BC_b_rel)*c1->T_b();
//
//            // lateral contributions to volume changes
//            dVc1_dV1a += (Point::dV_dp(BC_l_rel, v2b_rel, v2a_rel)+Point::dV_dp(BC_l_rel, v1b_rel, v2b_rel))/4 + Point::dV_dp_b(v1a_rel, v1b_rel, BC_l_rel, 4) - Point::dV_dp_b(v1a_rel, v2a_rel, BC_l_rel,4);
//            dVc1_dV1b += (Point::dV_dp(BC_l_rel, v2a_rel, v1a_rel)+Point::dV_dp(BC_l_rel, v2b_rel, v2a_rel))/4 + Point::dV_dp_b(v1b_rel, v2b_rel, BC_l_rel, 4) - Point::dV_dp_b(v1b_rel, v1a_rel, BC_l_rel,4);
//            dVc1_dV2a += (Point::dV_dp(BC_l_rel, v1b_rel, v2b_rel)+Point::dV_dp(BC_l_rel, v1a_rel, v1b_rel))/4 + Point::dV_dp_b(v2a_rel, v1a_rel, BC_l_rel, 4) - Point::dV_dp_b(v2a_rel, v2b_rel, BC_l_rel,4);
//            dVc1_dV2b += (Point::dV_dp(BC_l_rel, v1a_rel, v1b_rel)+Point::dV_dp(BC_l_rel, v2a_rel, v1a_rel))/4 + Point::dV_dp_b(v2b_rel, v2a_rel, BC_l_rel, 4) - Point::dV_dp_b(v2b_rel, v1b_rel, BC_l_rel,4);
//
//            // apical and basal contributions to volume changes
//            dVc1_dV1a += Point::dV_dp_b_1(v1a_rel, v2a_rel, c1->BC_a, n1);
//            dVc1_dV1b -= Point::dV_dp_b_1(v1b_rel, v2b_rel, BC_b_rel, n1);
//            dVc1_dV2a -= Point::dV_dp_b_1(v2a_rel, v1a_rel, c1->BC_a, n1);
//            dVc1_dV2b += Point::dV_dp_b_1(v2b_rel, v1b_rel, BC_b_rel, n1);
//
//        }
//        
//        if(!crossBoundary_c2) {
//            
//            // APICAL AND BASAL SURFACE TENSION
//            // force on vertices
//            Fa_v1 -= Point::TriangleAreaChange3(v1a,v2a,c2->BC_a,n2) * c2->T_a();
//            Fb_v1 -= Point::TriangleAreaChange3(v1b,v2b,c2->BC_b,n2) * c2->T_b();
//            Fa_v2 -= Point::TriangleAreaChange3(v2a,v1a,c2->BC_a,n2) * c2->T_a();
//            Fb_v2 -= Point::TriangleAreaChange3(v2b,v1b,c2->BC_b,n2) * c2->T_b();
//            // force on system size
//            T->dW_dL -= (Point::surfaceStress(v1a, c2->BC_a, v2a)) * c2->T_a() + (Point::surfaceStress(v1b, c2->BC_b, v2b))*c2->T_b();
//            // energy contribution
//            T->Energy_Surface +=  Point::TriangleArea(v1a,v2a,c2->BC_a)*c2->T_a() +  Point::TriangleArea(v1b,v2b,c2->BC_b)*c2->T_b();
//            
//            // lateral contributions to volume changes
//            dVc2_dV1a += ((Point::dV_dp(BC_l, v2->coord_b, v2->coord)+Point::dV_dp(BC_l, v1->coord_b, v2->coord_b))/4 + Point::dV_dp_b(v1->coord, v1->coord_b, BC_l, 4) - Point::dV_dp_b(v1->coord, v2->coord, BC_l,4))*(-1);
//            dVc2_dV1b += ((Point::dV_dp(BC_l, v2->coord, v1->coord)+Point::dV_dp(BC_l, v2->coord_b, v2->coord))/4 + Point::dV_dp_b(v1->coord_b, v2->coord_b, BC_l, 4) - Point::dV_dp_b(v1->coord_b, v1->coord, BC_l,4))*(-1);
//            dVc2_dV2a += ((Point::dV_dp(BC_l, v1->coord_b, v2->coord_b)+Point::dV_dp(BC_l, v1->coord, v1->coord_b))/4 + Point::dV_dp_b(v2->coord, v1->coord, BC_l, 4) - Point::dV_dp_b(v2->coord, v2->coord_b, BC_l,4))*(-1);
//            dVc2_dV2b += ((Point::dV_dp(BC_l, v1->coord, v1->coord_b)+Point::dV_dp(BC_l, v2->coord, v1->coord))/4 + Point::dV_dp_b(v2->coord_b, v2->coord, BC_l, 4) - Point::dV_dp_b(v2->coord_b, v1->coord_b, BC_l,4))*(-1);
//            
//            // apical and basal contributions to volume changes
//            dVc2_dV1a -= Point::dV_dp_b_1(v1->coord, v2->coord, c2->BC_a, n2);
//            dVc2_dV1b += Point::dV_dp_b_1(v1->coord_b, v2->coord_b, c2->BC_b, n2);
//            dVc2_dV2a += Point::dV_dp_b_1(v2->coord, v1->coord, c2->BC_a, n2);
//            dVc2_dV2b -= Point::dV_dp_b_1(v2->coord_b, v1->coord_b, c2->BC_b, n2);
//
//        } else {
//            
//            // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE APICAL BARYCENTER OF CELL 2
//            Point BC_l_rel = shift( b, SystemSize*( q_BC_l -  q_c2 -  q_v2));
//            Point v1a_rel = shift( v1a, SystemSize*(  q_c2 +  q_v2 ) *(-1));
//            Point v1b_rel = shift( v1b, SystemSize*(v1->q_b -  q_c2 -  q_v2));
//            Point v2a_rel = shift( v2a, SystemSize* q_c2*(-1));
//            Point v2b_rel = shift( v2b, SystemSize*(  q_c2 - v2->q_b)*(-1));
//            Point BC_b_rel = shift( c2->BC_b, SystemSize*c2->q_BC_b );
//            
//            // APICAL AND BASAL SURFACE TENSION
//            // force on vertices
//            Fa_v1 -= Point::TriangleAreaChange3(v1a_rel,v2a_rel,c2->BC_a,n2) * c2->T_a();
//            Fb_v1 -= Point::TriangleAreaChange3(v1b_rel,v2b_rel,BC_b_rel,n2) * c2->T_b();
//            Fa_v2 -= Point::TriangleAreaChange3(v2a_rel,v1a_rel,c2->BC_a,n2) * c2->T_a();
//            Fb_v2 -= Point::TriangleAreaChange3(v2b_rel,v1b_rel,BC_b_rel,n2) * c2->T_b();
//            // force on system size
//            T->dW_dL -= (Point::surfaceStress(v1a_rel, c2->BC_a, v2a_rel)) * c2->T_a() + (Point::surfaceStress(v1b_rel, BC_b_rel, v2b_rel)) *c2->T_b();
//            // energy contribution
//            T->Energy_Surface +=  Point::TriangleArea(v1a_rel,v2a_rel,c2->BC_a)*c2->T_a() +  Point::TriangleArea(v1b_rel,v2b_rel,BC_b_rel)*c2->T_b();
//
//            // lateral contributions to volume changes
//            dVc2_dV1a += ((Point::dV_dp(BC_l_rel, v2b_rel, v2a_rel)+Point::dV_dp(BC_l_rel, v1b_rel, v2b_rel))/4 + Point::dV_dp_b(v1a_rel, v1b_rel, BC_l_rel, 4) - Point::dV_dp_b(v1a_rel, v2a_rel, BC_l_rel,4))*(-1);
//            dVc2_dV1b += ((Point::dV_dp(BC_l_rel, v2a_rel, v1a_rel)+Point::dV_dp(BC_l_rel, v2b_rel, v2a_rel))/4 + Point::dV_dp_b(v1b_rel, v2b_rel, BC_l_rel, 4) - Point::dV_dp_b(v1b_rel, v1a_rel, BC_l_rel,4))*(-1);
//            dVc2_dV2a += ((Point::dV_dp(BC_l_rel, v1b_rel, v2b_rel)+Point::dV_dp(BC_l_rel, v1a_rel, v1b_rel))/4 + Point::dV_dp_b(v2a_rel, v1a_rel, BC_l_rel, 4) - Point::dV_dp_b(v2a_rel, v2b_rel, BC_l_rel,4))*(-1);
//            dVc2_dV2b += ((Point::dV_dp(BC_l_rel, v1a_rel, v1b_rel)+Point::dV_dp(BC_l_rel, v2a_rel, v1a_rel))/4 + Point::dV_dp_b(v2b_rel, v2a_rel, BC_l_rel, 4) - Point::dV_dp_b(v2b_rel, v1b_rel, BC_l_rel,4))*(-1);
//
//            // apical and basal contributions to volume changes
//            dVc2_dV1a -= Point::dV_dp_b_1(v1a_rel, v2a_rel, c2->BC_a, n2);
//            dVc2_dV1b += Point::dV_dp_b_1(v1b_rel, v2b_rel, BC_b_rel, n2);
//            dVc2_dV2a += Point::dV_dp_b_1(v2a_rel, v1a_rel, c2->BC_a, n2);
//            dVc2_dV2b -= Point::dV_dp_b_1(v2b_rel, v1b_rel, BC_b_rel, n2);
//        }
//        
//        // use the volume derivatives to calculate the forces acting on the vertices
//        Fa_v1 += dVc1_dV1a*c1->Pressure + dVc2_dV1a*c2->Pressure;
//        Fa_v2 += dVc1_dV2a*c1->Pressure + dVc2_dV2a*c2->Pressure;
//        Fb_v1 += dVc1_dV1b*c1->Pressure + dVc2_dV1b*c2->Pressure;
//        Fb_v2 += dVc1_dV2b*c1->Pressure + dVc2_dV2b*c2->Pressure;
//        
//        // ADD UP THE FORCE CONTRIBUTIONS IN THE 4 VERTEX POSITIONS THAT BELONG TO THE EDGE
//        v1->dW_dva += Fa_v1;
//        v2->dW_dva += Fa_v2;
//        v1->dW_dvb += Fb_v1;
//        v2->dW_dvb += Fb_v2;
//    }
//}


//void Edge::UpdateForcesAndEnergy()
//{
//    if(!T->isPeriodic) {
//        // SET THE FORCE CONTRIBUTIONS FROM THE EDGE ON EVERY OF ITS VERTICES
//        Point Fa_v1;
//        Point Fb_v1;
//        Point Fa_v2;
//        Point Fb_v2;
//        
//        // positions
//        Point v1a = v1->coord;
//        Point v1b = v1->coord_b;
//        Point v2a = v2->coord;
//        Point v2b = v2->coord_b;
//        Point b = BC_l;
//        
//        // APICAL AND BASAL LINE TENSION
//        // force on vertices
//        Fa_v1 = (v2a-v1a)*(Lambda_a()/la);
//        Fa_v2 = Fa_v1*(-1);
//
//        Fb_v1 = (v2b-v1b)*(Lambda_b()/lb);
//        Fb_v2 = Fb_v1*(-1);
//        
//        // energy contribution
//        T->Energy_Line+=Lambda_a()*la+Lambda_b()*lb;
//        
//        
//        // LATEREAL SURFACE TENSION
//        // force on vertices
//        Fa_v1 += ForceFromBarycenter_Area - (Point::TriangleAreaChange(v1a,v2a,b)+Point::TriangleAreaChange(v1a,v1b,b))*T_l;
//        Fa_v2 += ForceFromBarycenter_Area - (Point::TriangleAreaChange(v2a,v1a,b)+Point::TriangleAreaChange(v2a,v2b,b))*T_l;
//        Fb_v1 += ForceFromBarycenter_Area - (Point::TriangleAreaChange(v1b,v2b,b)+Point::TriangleAreaChange(v1b,v1a,b))*T_l;
//        Fb_v2 += ForceFromBarycenter_Area - (Point::TriangleAreaChange(v2b,v1b,b)+Point::TriangleAreaChange(v2b,v2a,b))*T_l;
//        
//        // energy contribution
//        T->Energy_Surface +=  Point::TriangleArea(v1a,v2a,c1->BC_a)*c1->T_a() +  Point::TriangleArea(v1b,v2b,c1->BC_b)*c1->T_b() + Point::TriangleArea(v1a,v2a,c2->BC_a)*c2->T_a() +  Point::TriangleArea(v1b,v2b,c2->BC_b)*c2->T_b();
//        
//        // FORCE FROM VOLUME PRESSURE
//        // save some values to avoid the repeated calculation
//        Point b_v1a = Point::dV_dp(b,v1a);
//        Point b_v1b = Point::dV_dp(b,v1b);
//        Point b_v2a = Point::dV_dp(b,v2a);
//        Point b_v2b = Point::dV_dp(b,v2b);
//        
//        if(c1)
//        {
//            Point ForceFromLateralBarycenter_Volume = dV_dBC*(-c1->Pressure);
//            
//            (c1->ForceFromApicalBarycenter_Volume + ForceFromLateralBarycenter_Volume).Print();
//            
//            Fa_v1 += c1->dAa_dBCa + c1->ForceFromApicalBarycenter_Volume + ForceFromLateralBarycenter_Volume
//            - Point::TriangleAreaChange(v1a,v2a,c1->BC_a)*c1->T_a() - (b_v2a-b_v1b+Point::dV_dp(v2a,c1->BC_a))*c1->Pressure;
//            
//            Fa_v2 += ForceFromLateralBarycenter_Volume - Point::TriangleAreaChange(v2a,v1a,c1->BC_a)*c1->T_a() - (b_v2b-b_v1a+Point::dV_dp(c1->BC_a,v1a))*c1->Pressure;
//            
//            Fb_v1 += c1->dAb_dBCb + c1->ForceFromBasalBarycenter_Volume + ForceFromLateralBarycenter_Volume
//            - Point::TriangleAreaChange(v1b,v2b,c1->BC_b)*c1->T_b() - (b_v1a-b_v2b+Point::dV_dp(c1->BC_b,v2b))*c1->Pressure;
//            
//            Fb_v2 += ForceFromLateralBarycenter_Volume - Point::TriangleAreaChange(v2b,v1b,c1->BC_b)*c1->T_b() - (b_v1b-b_v2a+Point::dV_dp(v1b,c1->BC_b))*c1->Pressure;
//            
//            T->Energy_Surface +=  Point::TriangleArea(v1a,v2a,c1->BC_a)*c1->T_a() +  Point::TriangleArea(v1b,v2b,c1->BC_b)*c1->T_b();
//            
//        }
//        
//        if(c2)
//        {
//            Point ForceFromLateralBarycenter_Volume = dV_dBC*(c2->Pressure);
//            
//            Fa_v1 +=  ForceFromLateralBarycenter_Volume - Point::TriangleAreaChange(v1a,v2a,c1->BC_a)*c2->T_a() + c2->ForceFromApicalBarycenter_Volume +  (b_v2a-b_v1b+Point::dV_dp(v2a,c2->BC_a))*c2->Pressure;
//            
//            Fa_v2 +=  c2->dAa_dBCa + c2->ForceFromApicalBarycenter_Volume + ForceFromLateralBarycenter_Volume
//            - Point::TriangleAreaChange(v2a,v1a,c1->BC_a)*c2->T_a() + (b_v2b-b_v1a+Point::dV_dp(c2->BC_a,v1a))*c2->Pressure;
//            
//            Fb_v1 +=  ForceFromLateralBarycenter_Volume - Point::TriangleAreaChange(v1b,v2b,c2->BC_b)*c2->T_b() + (b_v1a-b_v2b+Point::dV_dp(c2->BC_b,v2b))*c2->Pressure;
//            
//            Fb_v2 += c2->dAb_dBCb + c2->ForceFromBasalBarycenter_Volume  + ForceFromLateralBarycenter_Volume
//            - Point::TriangleAreaChange(v2b,v1b,c2->BC_b)*c2->T_b() + (b_v1b-b_v2a+Point::dV_dp(v1b,c2->BC_b))*c2->Pressure;
//            
//            T->Energy_Surface += Point::TriangleArea(v1a,v2a,c2->BC_a)*c2->T_a() +  Point::TriangleArea(v1b,v2b,c2->BC_b)*c2->T_b();
//        }
//        
//        // ADD UP THE FORCE CONTRIBUTIONS IN THE 4 VERTEX POSITIONS THAT BELONG TO THE EDGE
//        v1->dW_dva += Fa_v1;
//        v2->dW_dva += Fa_v2;
//        v1->dW_dvb  += Fb_v1;
//        v2->dW_dvb  += Fb_v2;
//    }
//    
//    else { // periodic tissue
//        
//        // FORCE CONTRIBUTIONS FROM THE EDGE ON EVERY OF ITS VERTICES V1a,V1b,V2a,V2b
//        Point Fa_v1, Fb_v1, Fa_v2, Fb_v2, Fa_v1_b, Fb_v1_b, Fa_v2_b, Fb_v2_b;
//        
//        // positions
//        Point v1a = v1->coord;
//        Point v1b = v1->coord_b;
//        Point v2a = v2->coord;
//        Point v2b = v2->coord_b;
//        Point b = BC_l;
//        
//        //
//        Point dVc1_dV1a, dVc2_dV1a, dVc1_dV1b, dVc2_dV1b, dVc1_dV2a, dVc2_dV2a, dVc1_dV2b, dVc2_dV2b;
//        
//        // pure contributions from apical and basal barycenter, and from the lateral barycenter
//        int n1 = c1->ListEdges.size();
//        int n2 = c2->ListEdges.size();
//       
//        
//        double c1_Ta = c1->T_a();
//        double c2_Ta = c2->T_a();
//        double c1_Tb = c1->T_b();
//        double c2_Tb = c2->T_b();
//        // NOTE: every vertex is once vertex one with respect to a cell
//        // therefore add the "force contribution" coming from the movement in the apical and basal barycenters of v1 wrt c1 and of v2 wrt to c2 (since v1 is first in direction of c1 and v2 to c2)
//        v1->dW_dva -= c1->dAa_dBCa*(c1_Ta / c1->ListEdges.size());
//        v2->dW_dva -= c2->dAa_dBCa*(c2_Ta / c2->ListEdges.size());
//        v1->dW_dvb -= c1->dAb_dBCb*(c1_Tb / c1->ListEdges.size());
//        v2->dW_dvb -= c2->dAb_dBCb*(c2_Tb / c2->ListEdges.size());
//        
//        dVc1_dV1a += c1->dV_dBCa/n1;
//        dVc1_dV1b += c1->dV_dBCb/n1;
//        dVc2_dV2a += c2->dV_dBCa/n2;
//        dVc2_dV2b += c2->dV_dBCb/n2;
//        
//        // systemSize
//        Quadrant SystemSize = T->SystemSize;
//        
//        bool crossBoundary_c1 = c1->CrossBoundary();
//        bool crossBoundary_c2 = c2->CrossBoundary();
//        
//        if(!CrossBoundary()) { // edge doesnt cross the boundary
//            // APICAL AND BASAL LINE TENSION
//            // force on vertices
//            Fa_v1 += (v2a-v1a)*(Lambda_a()/la);
//            Fa_v2 += Fa_v1*(-1);
//
//            Fb_v1 += (v2b-v1b)*(Lambda_b()/lb);
//            Fb_v2 += Fb_v1*(-1);
//            
//            // force on the system size
//            T->dW_dL -= Point::lineStress(v1a,v2a)*Lambda_a() + Point::lineStress(v1b, v2b)*Lambda_b();
//            // energy contribution
//            T->Energy_Line += Lambda_a()*la+Lambda_b()*lb;
//            
//            // LATEREAL SURFACE TENSION
//            // force on vertices
//            Fa_v1 -= (Point::TriangleAreaChange2(v1a,v2a,b,4)+Point::TriangleAreaChange2(v1a,v1b,b,4)+Point::TriangleAreaChange(b,v1b,v2b)/4+Point::TriangleAreaChange(b,v2a,v2b)/4)*T_l;
//            Fa_v2 -= (Point::TriangleAreaChange2(v2a,v1a,b,4)+Point::TriangleAreaChange2(v2a,v2b,b,4)+Point::TriangleAreaChange(b,v1a,v1b)/4+Point::TriangleAreaChange(b,v1b,v2b)/4)*T_l;
//            Fb_v1 -= (Point::TriangleAreaChange2(v1b,v1a,b,4)+Point::TriangleAreaChange2(v1b,v2b,b,4)+Point::TriangleAreaChange(b,v1a,v2a)/4+Point::TriangleAreaChange(b,v2a,v2b)/4)*T_l;
//            Fb_v2 -= (Point::TriangleAreaChange2(v2b,v2a,b,4)+Point::TriangleAreaChange2(v2b,v1b,b,4)+Point::TriangleAreaChange(b,v1a,v2a)/4+Point::TriangleAreaChange(b,v1a,v1b)/4)*T_l;
//            // force on the system size
//            T->dW_dL -= (Point::surfaceStress(v1a,v2a,b) + Point::surfaceStress(v2a,v2b,b ) + Point::surfaceStress(v2b,v1b,b ) + Point::surfaceStress(v1b,v1a,b ))* T_l;
//            // energy contribution
//            T->Energy_Surface += T_l*(Point::TriangleArea(v1a,v2a,b)+Point::TriangleArea(v2a,v2b,b)+Point::TriangleArea(v2b,v1b,b)+Point::TriangleArea(v1b,v1a,b));
//            
//        } else { // edge crosses boundary
//            // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE FIRST VERTEX OF THE EDGE
//            Point v1a_rel = v1a;
//            Point v1b_rel = shift( v1b, SystemSize*v1->q_b);
//            Point v2a_rel = shift( v2a, SystemSize*q_v2);
//            Point v2b_rel = shift( v2b, SystemSize*(q_v2+v2->q_b));
//            Point b_rel =   shift( b, SystemSize*q_BC_l);
//            
//            // APICAL AND BASAL LINE TENSION
//            // force on vertices
//            Fa_v1 += (v2a_rel-v1a_rel)*(Lambda_a()/la);
//            Fa_v2 += Fa_v1*(-1);
//            
//            Fb_v1 += (v2b_rel-v1b_rel)*(Lambda_b()/lb);
//            Fb_v2 += Fb_v1*(-1);
//                        
//            // force on the system size
//            T->dW_dL -= Point::lineStress(v1a_rel,v2a_rel)*Lambda_a() + Point::lineStress(v1b_rel, v2b_rel)*Lambda_b();
//            // energy contribution
//            T->Energy_Line+=Lambda_a()*la+Lambda_b()*lb;
//            
//            // LATERAL SURFACE TENSION
//            // force on vertices
//            Fa_v1 -= (Point::TriangleAreaChange2(v1a_rel,v2a_rel,b_rel,4)+Point::TriangleAreaChange2(v1a_rel,v1b_rel,b_rel,4)+Point::TriangleAreaChange(b_rel,v1b_rel,v2b_rel)/4+Point::TriangleAreaChange(b_rel,v2a_rel,v2b_rel)/4)*T_l;
//            Fa_v2 -= (Point::TriangleAreaChange2(v2a_rel,v1a_rel,b_rel,4)+Point::TriangleAreaChange2(v2a_rel,v2b_rel,b_rel,4)+Point::TriangleAreaChange(b_rel,v1a_rel,v1b_rel)/4+Point::TriangleAreaChange(b_rel,v1b_rel,v2b_rel)/4)*T_l;
//            Fb_v1 -= (Point::TriangleAreaChange2(v1b_rel,v1a_rel,b_rel,4)+Point::TriangleAreaChange2(v1b_rel,v2b_rel,b_rel,4)+Point::TriangleAreaChange(b_rel,v1a_rel,v2a_rel)/4+Point::TriangleAreaChange(b_rel,v2a_rel,v2b_rel)/4)*T_l;
//            Fb_v2 -= (Point::TriangleAreaChange2(v2b_rel,v2a_rel,b_rel,4)+Point::TriangleAreaChange2(v2b_rel,v1b_rel,b_rel,4)+Point::TriangleAreaChange(b_rel,v1a_rel,v2a_rel)/4+Point::TriangleAreaChange(b_rel,v1a_rel,v1b_rel)/4)*T_l;
//            
//            // force on the system size
//            T->dW_dL -= (Point::surfaceStress(v1a_rel,v2a_rel,b_rel ) + Point::surfaceStress(v2a_rel,v2b_rel,b_rel ) + Point::surfaceStress(v2b_rel,v1b_rel,b_rel ) + Point::surfaceStress(v1b_rel,v1a_rel,b_rel ))*  T_l;
//            // energy contribution
//            T->Energy_Surface += T_l*(Point::TriangleArea(v1a_rel,v2a_rel,b_rel)+Point::TriangleArea(v2a_rel,v2b_rel,b_rel)+Point::TriangleArea(v2b_rel,v1b_rel,b_rel)+Point::TriangleArea(v1b_rel,v1a_rel,b_rel));
//        }
//        
//        if(!crossBoundary_c1)  {
//            
//            //actually: double T_a = c1->T_a+c1->K_2D*(c1->CalculateApicalSurfaceArea()-T->cells_A0)
//            double T_a = c1_Ta;
//            double T_b = c1_Tb;
//            
//            // APICAL AND BASAL SURFACE TENSIONS
//            Fa_v1 -= Point::TriangleAreaChange3(v1a,v2a,c1->BC_a,n1) * T_a;
//            Fb_v1 -= Point::TriangleAreaChange3(v1b,v2b,c1->BC_b,n1) * T_b;
//            Fa_v2 -= Point::TriangleAreaChange3(v2a,v1a,c1->BC_a,n1) * T_a;
//            Fb_v2 -= Point::TriangleAreaChange3(v2b,v1b,c1->BC_b,n1) * T_b;
//            // force on system size
//            T->dW_dL -= (Point::surfaceStress(v1a, c1->BC_a, v2a)) * T_a + (Point::surfaceStress(v1b, c1->BC_b, v2b)) * T_b;
//            // energy contribution
//            T->Energy_Surface +=  Point::TriangleArea(v1a,v2a,c1->BC_a)*T_a +  Point::TriangleArea(v1b,v2b,c1->BC_b)*T_b;
//            
//            // lateral contributions to volume changes
//            dVc1_dV1a += (Point::dV_dp(BC_l, v2->coord_b, v2->coord)+Point::dV_dp(BC_l, v1->coord_b, v2->coord_b))/4 + Point::dV_dp_b(v1->coord, v1->coord_b, BC_l, 4) - Point::dV_dp_b(v1->coord, v2->coord, BC_l,4);
//            dVc1_dV1b += (Point::dV_dp(BC_l, v2->coord, v1->coord)+Point::dV_dp(BC_l, v2->coord_b, v2->coord))/4 + Point::dV_dp_b(v1->coord_b, v2->coord_b, BC_l, 4) - Point::dV_dp_b(v1->coord_b, v1->coord, BC_l,4);
//            dVc1_dV2a += (Point::dV_dp(BC_l, v1->coord_b, v2->coord_b)+Point::dV_dp(BC_l, v1->coord, v1->coord_b))/4 + Point::dV_dp_b(v2->coord, v1->coord, BC_l, 4) - Point::dV_dp_b(v2->coord, v2->coord_b, BC_l,4);
//            dVc1_dV2b += (Point::dV_dp(BC_l, v1->coord, v1->coord_b)+Point::dV_dp(BC_l, v2->coord, v1->coord))/4 + Point::dV_dp_b(v2->coord_b, v2->coord, BC_l, 4) - Point::dV_dp_b(v2->coord_b, v1->coord_b, BC_l,4);
//            
//            // apical and basal contributions to volume changes
//            dVc1_dV1a += Point::dV_dp_b_1(v1->coord, v2->coord, c1->BC_a, n1);
//            dVc1_dV1b -= Point::dV_dp_b_1(v1->coord_b, v2->coord_b, c1->BC_b, n1);
//            dVc1_dV2a -= Point::dV_dp_b_1(v2->coord, v1->coord, c1->BC_a, n1);
//            dVc1_dV2b += Point::dV_dp_b_1(v2->coord_b, v1->coord_b, c1->BC_b, n1);
//            
//        } else {
//           
//            //actually: double T_a = c1->T_a+c1->K_2D*(c1->CalculateApicalSurfaceArea()-T->cells_A0)
//            double T_a = c1_Ta;
//            double T_b = c1_Tb;
//            
//            // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE APICAL BARYCENTER OF CELL 1
//            Point BC_l_rel = shift( b, SystemSize*( q_BC_l -  q_c1));
//            Point v1a_rel = shift( v1a, SystemSize* q_c1*(-1));
//            Point v1b_rel = shift( v1b, SystemSize*(v1->q_b -  q_c1));
//            Point v2a_rel = shift( v2a, SystemSize*( q_v2 -  q_c1));
//            Point v2b_rel = shift( v2b, SystemSize*( q_v2 + v2->q_b -  q_c1));
//            Point BC_b_rel = shift( c1->BC_b, SystemSize* c1->q_BC_b);
//            
//            // APICAL AND BASAL SURFACE TENSION
//            // force on vertices
//            Fa_v1 -= Point::TriangleAreaChange3(v1a_rel,v2a_rel,c1->BC_a,n1) * T_a;
//            Fb_v1 -= Point::TriangleAreaChange3(v1b_rel,v2b_rel,BC_b_rel,n1) * T_b;
//            Fa_v2 -= Point::TriangleAreaChange3(v2a_rel,v1a_rel,c1->BC_a,n1) * T_a;
//            Fb_v2 -= Point::TriangleAreaChange3(v2b_rel,v1b_rel,BC_b_rel,n1) * T_b;
//            // force on system size
//            T->dW_dL -= (Point::surfaceStress(v1a_rel, c1->BC_a, v2a_rel)) * T_a + (Point::surfaceStress(v1b_rel, BC_b_rel, v2b_rel)) * T_b;
//            // energy contribution
//            T->Energy_Surface +=  Point::TriangleArea(v1a_rel,v2a_rel,c1->BC_a) * T_a +  Point::TriangleArea(v1b_rel,v2b_rel,BC_b_rel)*T_b;
//            
//            // lateral contributions to volume changes
//            dVc1_dV1a += (Point::dV_dp(BC_l_rel, v2b_rel, v2a_rel)+Point::dV_dp(BC_l_rel, v1b_rel, v2b_rel))/4 + Point::dV_dp_b(v1a_rel, v1b_rel, BC_l_rel, 4) - Point::dV_dp_b(v1a_rel, v2a_rel, BC_l_rel,4);
//            dVc1_dV1b += (Point::dV_dp(BC_l_rel, v2a_rel, v1a_rel)+Point::dV_dp(BC_l_rel, v2b_rel, v2a_rel))/4 + Point::dV_dp_b(v1b_rel, v2b_rel, BC_l_rel, 4) - Point::dV_dp_b(v1b_rel, v1a_rel, BC_l_rel,4);
//            dVc1_dV2a += (Point::dV_dp(BC_l_rel, v1b_rel, v2b_rel)+Point::dV_dp(BC_l_rel, v1a_rel, v1b_rel))/4 + Point::dV_dp_b(v2a_rel, v1a_rel, BC_l_rel, 4) - Point::dV_dp_b(v2a_rel, v2b_rel, BC_l_rel,4);
//            dVc1_dV2b += (Point::dV_dp(BC_l_rel, v1a_rel, v1b_rel)+Point::dV_dp(BC_l_rel, v2a_rel, v1a_rel))/4 + Point::dV_dp_b(v2b_rel, v2a_rel, BC_l_rel, 4) - Point::dV_dp_b(v2b_rel, v1b_rel, BC_l_rel,4);
//            
//            // apical and basal contributions to volume changes
//            dVc1_dV1a += Point::dV_dp_b_1(v1a_rel, v2a_rel, c1->BC_a, n1);
//            dVc1_dV1b -= Point::dV_dp_b_1(v1b_rel, v2b_rel, BC_b_rel, n1);
//            dVc1_dV2a -= Point::dV_dp_b_1(v2a_rel, v1a_rel, c1->BC_a, n1);
//            dVc1_dV2b += Point::dV_dp_b_1(v2b_rel, v1b_rel, BC_b_rel, n1);
//            
//        }
//        
//        //if(false)
//        if(!crossBoundary_c2) {
//            
//            //save time by computing the current surface tension only once
//            double T_a = c2_Ta;
//            double T_b = c2_Tb;
//
//            // APICAL AND BASAL SURFACE TENSION
//            // force on vertices
//            Fa_v1 -= Point::TriangleAreaChange3(v1a,v2a,c2->BC_a,n2) * T_a;
//            Fb_v1 -= Point::TriangleAreaChange3(v1b,v2b,c2->BC_b,n2) * T_b;
//            Fa_v2 -= Point::TriangleAreaChange3(v2a,v1a,c2->BC_a,n2) * T_a;
//            Fb_v2 -= Point::TriangleAreaChange3(v2b,v1b,c2->BC_b,n2) * T_b;
//            // force on system size
//            T->dW_dL -= (Point::surfaceStress(v1a, c2->BC_a, v2a)) * T_a + (Point::surfaceStress(v1b, c2->BC_b, v2b)) * T_b;
//            // energy contribution
//            T->Energy_Surface +=  Point::TriangleArea(v1a,v2a,c2->BC_a) * T_a +  Point::TriangleArea(v1b,v2b,c2->BC_b) * T_b;
//            
//            // lateral contributions to volume changes
//            dVc2_dV1a += ((Point::dV_dp(BC_l, v2->coord_b, v2->coord)+Point::dV_dp(BC_l, v1->coord_b, v2->coord_b))/4 + Point::dV_dp_b(v1->coord, v1->coord_b, BC_l, 4) - Point::dV_dp_b(v1->coord, v2->coord, BC_l,4))*(-1);
//            dVc2_dV1b += ((Point::dV_dp(BC_l, v2->coord, v1->coord)+Point::dV_dp(BC_l, v2->coord_b, v2->coord))/4 + Point::dV_dp_b(v1->coord_b, v2->coord_b, BC_l, 4) - Point::dV_dp_b(v1->coord_b, v1->coord, BC_l,4))*(-1);
//            dVc2_dV2a += ((Point::dV_dp(BC_l, v1->coord_b, v2->coord_b)+Point::dV_dp(BC_l, v1->coord, v1->coord_b))/4 + Point::dV_dp_b(v2->coord, v1->coord, BC_l, 4) - Point::dV_dp_b(v2->coord, v2->coord_b, BC_l,4))*(-1);
//            dVc2_dV2b += ((Point::dV_dp(BC_l, v1->coord, v1->coord_b)+Point::dV_dp(BC_l, v2->coord, v1->coord))/4 + Point::dV_dp_b(v2->coord_b, v2->coord, BC_l, 4) - Point::dV_dp_b(v2->coord_b, v1->coord_b, BC_l,4))*(-1);
//            
//            // apical and basal contributions to volume changes
//            dVc2_dV1a -= Point::dV_dp_b_1(v1->coord, v2->coord, c2->BC_a, n2);
//            dVc2_dV1b += Point::dV_dp_b_1(v1->coord_b, v2->coord_b, c2->BC_b, n2);
//            dVc2_dV2a += Point::dV_dp_b_1(v2->coord, v1->coord, c2->BC_a, n2);
//            dVc2_dV2b -= Point::dV_dp_b_1(v2->coord_b, v1->coord_b, c2->BC_b, n2);
//            
//        } else {
//            double T_a = c2_Ta;
//            double T_b = c2_Tb;
//
//            // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE APICAL BARYCENTER OF CELL 2
//            Point BC_l_rel = shift( b, SystemSize*( q_BC_l -  q_c2 -  q_v2));
//            Point v1a_rel = shift( v1a, SystemSize*(  q_c2 +  q_v2 ) *(-1));
//            Point v1b_rel = shift( v1b, SystemSize*(v1->q_b -  q_c2 -  q_v2));
//            Point v2a_rel = shift( v2a, SystemSize* q_c2*(-1));
//            Point v2b_rel = shift( v2b, SystemSize*(  q_c2 - v2->q_b)*(-1));
//            Point BC_b_rel = shift( c2->BC_b, SystemSize*c2->q_BC_b );
//            
//            // APICAL AND BASAL SURFACE TENSION
//            // force on vertices
//            Fa_v1 -= Point::TriangleAreaChange3(v1a_rel,v2a_rel,c2->BC_a,n2) * T_a;
//            Fb_v1 -= Point::TriangleAreaChange3(v1b_rel,v2b_rel,BC_b_rel,n2) * T_b;
//            Fa_v2 -= Point::TriangleAreaChange3(v2a_rel,v1a_rel,c2->BC_a,n2) * T_a;
//            Fb_v2 -= Point::TriangleAreaChange3(v2b_rel,v1b_rel,BC_b_rel,n2) * T_b;
//            // force on system size
//            T->dW_dL -= (Point::surfaceStress(v1a_rel, c2->BC_a, v2a_rel)) * T_a + (Point::surfaceStress(v1b_rel, BC_b_rel, v2b_rel)) * T_b;
//            // energy contribution
//            T->Energy_Surface +=  Point::TriangleArea(v1a_rel,v2a_rel,c2->BC_a)*T_a +  Point::TriangleArea(v1b_rel,v2b_rel,BC_b_rel)* T_b;
//            
//            // lateral contributions to volume changes
//            dVc2_dV1a += ((Point::dV_dp(BC_l_rel, v2b_rel, v2a_rel)+Point::dV_dp(BC_l_rel, v1b_rel, v2b_rel))/4 + Point::dV_dp_b(v1a_rel, v1b_rel, BC_l_rel, 4) - Point::dV_dp_b(v1a_rel, v2a_rel, BC_l_rel,4))*(-1);
//            dVc2_dV1b += ((Point::dV_dp(BC_l_rel, v2a_rel, v1a_rel)+Point::dV_dp(BC_l_rel, v2b_rel, v2a_rel))/4 + Point::dV_dp_b(v1b_rel, v2b_rel, BC_l_rel, 4) - Point::dV_dp_b(v1b_rel, v1a_rel, BC_l_rel,4))*(-1);
//            dVc2_dV2a += ((Point::dV_dp(BC_l_rel, v1b_rel, v2b_rel)+Point::dV_dp(BC_l_rel, v1a_rel, v1b_rel))/4 + Point::dV_dp_b(v2a_rel, v1a_rel, BC_l_rel, 4) - Point::dV_dp_b(v2a_rel, v2b_rel, BC_l_rel,4))*(-1);
//            dVc2_dV2b += ((Point::dV_dp(BC_l_rel, v1a_rel, v1b_rel)+Point::dV_dp(BC_l_rel, v2a_rel, v1a_rel))/4 + Point::dV_dp_b(v2b_rel, v2a_rel, BC_l_rel, 4) - Point::dV_dp_b(v2b_rel, v1b_rel, BC_l_rel,4))*(-1);
//            
//            // apical and basal contributions to volume changes
//            dVc2_dV1a -= Point::dV_dp_b_1(v1a_rel, v2a_rel, c2->BC_a, n2);
//            dVc2_dV1b += Point::dV_dp_b_1(v1b_rel, v2b_rel, BC_b_rel, n2);
//            dVc2_dV2a += Point::dV_dp_b_1(v2a_rel, v1a_rel, c2->BC_a, n2);
//            dVc2_dV2b -= Point::dV_dp_b_1(v2b_rel, v1b_rel, BC_b_rel, n2);
//        }
//        
//        // use the volume derivatives to calculate the forces acting on the vertices
//        Fa_v1 += dVc1_dV1a*c1->Pressure + dVc2_dV1a*c2->Pressure;
//        Fa_v2 += dVc1_dV2a*c1->Pressure + dVc2_dV2a*c2->Pressure;
//        Fb_v1 += dVc1_dV1b*c1->Pressure + dVc2_dV1b*c2->Pressure;
//        Fb_v2 += dVc1_dV2b*c1->Pressure + dVc2_dV2b*c2->Pressure;
//        
//        // ADD UP THE FORCE CONTRIBUTIONS IN THE 4 VERTEX POSITIONS THAT BELONG TO THE EDGE
//        v1->dW_dva += Fa_v1;
//        v2->dW_dva += Fa_v2;
//        v1->dW_dvb += Fb_v1;
//        v2->dW_dvb += Fb_v2;
//    }
//}


void Edge::update_forces()
{
    // important: if only apical contributions save computational time, since all other terms are zero anyways!
    if(T->onlyApicalContributions && !T->isPeriodic) {
        
        // positions
        Point v1a = v1->coord;
        Point v2a = v2->coord;

        // forces
        Point Fa_v1, Fa_v2;

        // APICAL LINE TENSION
        // force on vertices and system size from line tensions
        if(la>l0_a) {
            Fa_v1 += (v2a-v1a)*(G_a/la);
            Fa_v2 += Fa_v1*(-1);
        } else {
            Point d = v2a-v1a;
            Fa_v1 += d*(G_a/l0_a);
            Fa_v2 += Fa_v1*(-1);
        }
        
        // anisotropic apical line tension (along the z-axis)
        if(G_a_ai>0){
            if(v2a.z>v1a.z)
            {
                Fa_v1.z += G_a_ai;
                Fa_v2.z -= G_a_ai;
            } else {
                Fa_v1.z -= G_a_ai;
                Fa_v2.z += G_a_ai;
            }
        }

        // FORCES FROM PRESSURE AND APICAL AND BASAL SURFACE TENSIONS
        Point dVa_dV1a, dVa_dV2a;
        
        if(c1){
            int n1 = c1->ListEdges.size();
            double T_a = c1->T_a();
            
            // NOTE: every vertex is once vertex one with respect to a cell
            // therefore add the "force contribution" coming from the movement in the apical and basal barycenters of v1 wrt c1 and of v2 wrt to c2 (since v1 is first in direction of c1 and v2 to c2)
            v1->dW_dva -= c1->dAa_dBCa*(T_a / c1->ListEdges.size());
            // volume change from barycenter movement
            dVa_dV1a += c1->dV_dBCa/n1;

            Fa_v1 -= Point::TriangleAreaChange3(v1a,v2a,c1->BC_a,n1) * T_a;
            Fa_v2 -= Point::TriangleAreaChange3(v2a,v1a,c1->BC_a,n1) * T_a;
            
            dVa_dV1a += Point::dV_dp_b_1(v1->coord, v2->coord, c1->BC_a, n1);
            dVa_dV2a -= Point::dV_dp_b_1(v2->coord, v1->coord, c1->BC_a, n1);
            
        }
        
        if(c2)
        {
            // pure contributions from apical and basal barycenter, and from the lateral barycenter
            int n2 = c2->ListEdges.size();
            double T_a = c2->T_a();
            
            // NOTE: every vertex is once vertex one with respect to a cell
            // therefore add the "force contribution" coming from the movement in the apical and basal barycenters of v1 wrt c1 and of v2 wrt to c2 (since v1 is first in direction of c1 and v2 to c2)
            v2->dW_dva -= c2->dAa_dBCa*(T_a / c2->ListEdges.size());
            // volume change from barycenter movement
            dVa_dV2a += c2->dV_dBCa/n2;

            // APICAL AND BASAL SURFACE TENSION
            Fa_v1 -= Point::TriangleAreaChange3(v1a,v2a,c2->BC_a,n2) * T_a;
            Fa_v2 -= Point::TriangleAreaChange3(v2a,v1a,c2->BC_a,n2)* T_a;

            
            dVa_dV1a -= Point::dV_dp_b_1(v1->coord, v2->coord, c2->BC_a, n2);
            dVa_dV2a += Point::dV_dp_b_1(v2->coord, v1->coord, c2->BC_a, n2);
            
        }

        Fa_v1 += dVa_dV1a*T->P_Lumen_a;
        Fa_v2 += dVa_dV2a*T->P_Lumen_a;
        
        // ADD UP THE FORCE CONTRIBUTIONS
        v1->dW_dva += Fa_v1;
        v2->dW_dva += Fa_v2;
        
        return;
    }
    
    // FORCE CONTRIBUTIONS FROM THE EDGE ON EVERY OF ITS VERTICES V1a,V1b,V2a,V2b
    Point Fa_v1, Fb_v1, Fa_v2, Fb_v2, Fa_v1_b, Fb_v1_b, Fa_v2_b, Fb_v2_b;
    
    // positions
    Point v1a = v1->coord;
    Point v1b = v1->coord_b;
    Point v2a = v2->coord;
    Point v2b = v2->coord_b;
    Point b = BC_l;
    
    UpdateLengths();
    
    if(!T->isPeriodic) {
        
        // APICAL AND BASAL LINE TENSION
        // force on vertices and system size from line tensions
        if(la>l0_a) {
            Fa_v1 += (v2a-v1a)*(G_a/la);
            Fa_v2 += Fa_v1*(-1);
            T->dW_dL -= Point::lineStress(v1a,v2a)*Lambda_a();
        } else {
            Point d = v2a-v1a;
            Fa_v1 += d*(G_a/l0_a);
            Fa_v2 += Fa_v1*(-1);
            T->dW_dL -= Point(pow(d.x,2), pow(d.y,2),0)*(G_a/l0_a);
        }
        
        // anisotropic apical line tension (along the z-axis)
        if(G_a_ai>0){
            if(v2a.z>v1a.z)
            {
                Fa_v1.z += G_a_ai;
                Fa_v2.z -= G_a_ai;
            } else {
                Fa_v1.z -= G_a_ai;
                Fa_v2.z += G_a_ai;
            }
        }
        
        if(lb>l0_b) {
            Fb_v1 += (v2b-v1b)*(G_b/lb);
            Fb_v2 += Fb_v1*(-1);
            T->dW_dL -= Point::lineStress(v1b,v2b)*Lambda_b();
        } else {
            Point d = v2b-v1b;
            Fb_v1 += d*(G_b/l0_b);
            Fb_v2 += Fb_v1*(-1);
            T->dW_dL -= Point(pow(d.x,2), pow(d.y,2),0)*(G_b/l0_b);
        }


        
        double Tl = T_l;//lateralSurfaceTension();

        // LATERAL SURFACE TENSION
        // force on vertices
        Fa_v1 -= (Point::TriangleAreaChange2(v1a,v2a,b,4)+Point::TriangleAreaChange2(v1a,v1b,b,4)+Point::TriangleAreaChange(b,v1b,v2b)/4+Point::TriangleAreaChange(b,v2a,v2b)/4)*Tl;
        Fa_v2 -= (Point::TriangleAreaChange2(v2a,v1a,b,4)+Point::TriangleAreaChange2(v2a,v2b,b,4)+Point::TriangleAreaChange(b,v1a,v1b)/4+Point::TriangleAreaChange(b,v1b,v2b)/4)*Tl;
        Fb_v1 -= (Point::TriangleAreaChange2(v1b,v1a,b,4)+Point::TriangleAreaChange2(v1b,v2b,b,4)+Point::TriangleAreaChange(b,v1a,v2a)/4+Point::TriangleAreaChange(b,v2a,v2b)/4)*Tl;
        Fb_v2 -= (Point::TriangleAreaChange2(v2b,v2a,b,4)+Point::TriangleAreaChange2(v2b,v1b,b,4)+Point::TriangleAreaChange(b,v1a,v2a)/4+Point::TriangleAreaChange(b,v1a,v1b)/4)*Tl;

        if(c1) {
            
            //
            Point dVc1_dV1a, dVc2_dV1a, dVc1_dV1b, dVc2_dV1b, dVc1_dV2a, dVc2_dV2a, dVc1_dV2b, dVc2_dV2b;
            
            // pure contributions from apical and basal barycenter, and from the lateral barycenter
            int n1 = c1->ListEdges.size();
            
            double T_a = c1->T_a();
            double T_b = c1->T_b();
            
            // NOTE: every vertex is once vertex one with respect to a cell
            // therefore add the "force contribution" coming from the movement in the apical and basal barycenters of v1 wrt c1 and of v2 wrt to c2 (since v1 is first in direction of c1 and v2 to c2)
            v1->dW_dva -= c1->dAa_dBCa*(T_a / c1->ListEdges.size());
            v1->dW_dvb -= c1->dAb_dBCb*(T_b / c1->ListEdges.size());
            
            dVc1_dV1a += c1->dV_dBCa/n1;
            dVc1_dV1b += c1->dV_dBCb/n1;
            
            // APICAL AND BASAL SURFACE TENSIONS
            Fa_v1 -= Point::TriangleAreaChange3(v1a,v2a,c1->BC_a,n1) * T_a;
            Fb_v1 -= Point::TriangleAreaChange3(v1b,v2b,c1->BC_b,n1) * T_b;
            Fa_v2 -= Point::TriangleAreaChange3(v2a,v1a,c1->BC_a,n1) * T_a;
            Fb_v2 -= Point::TriangleAreaChange3(v2b,v1b,c1->BC_b,n1) * T_b;
            
            // lateral contributions to volume changes
            dVc1_dV1a += (Point::dV_dp(BC_l, v2->coord_b, v2->coord)+Point::dV_dp(BC_l, v1->coord_b, v2->coord_b))/4 + Point::dV_dp_b(v1->coord, v1->coord_b, BC_l, 4) - Point::dV_dp_b(v1->coord, v2->coord, BC_l,4);
            dVc1_dV1b += (Point::dV_dp(BC_l, v2->coord, v1->coord)+Point::dV_dp(BC_l, v2->coord_b, v2->coord))/4 + Point::dV_dp_b(v1->coord_b, v2->coord_b, BC_l, 4) - Point::dV_dp_b(v1->coord_b, v1->coord, BC_l,4);
            dVc1_dV2a += (Point::dV_dp(BC_l, v1->coord_b, v2->coord_b)+Point::dV_dp(BC_l, v1->coord, v1->coord_b))/4 + Point::dV_dp_b(v2->coord, v1->coord, BC_l, 4) - Point::dV_dp_b(v2->coord, v2->coord_b, BC_l,4);
            dVc1_dV2b += (Point::dV_dp(BC_l, v1->coord, v1->coord_b)+Point::dV_dp(BC_l, v2->coord, v1->coord))/4 + Point::dV_dp_b(v2->coord_b, v2->coord, BC_l, 4) - Point::dV_dp_b(v2->coord_b, v1->coord_b, BC_l,4);
            
            // apical and basal contributions to volume changes
            dVc1_dV1a += Point::dV_dp_b_1(v1->coord, v2->coord, c1->BC_a, n1);
            dVc1_dV1b -= Point::dV_dp_b_1(v1->coord_b, v2->coord_b, c1->BC_b, n1);
            dVc1_dV2a -= Point::dV_dp_b_1(v2->coord, v1->coord, c1->BC_a, n1);
            dVc1_dV2b += Point::dV_dp_b_1(v2->coord_b, v1->coord_b, c1->BC_b, n1);
            
            // use the volume derivatives to calculate the forces acting on the vertices
            if(T->K_Lumen_a==0 && T->P0_Lumen_a==0) {// only the pressure from the cells is taken into account
                Fa_v1 += dVc1_dV1a*c1->Pressure;
                Fa_v2 += dVc1_dV2a*c1->Pressure;
            } else {
                Fa_v1 += dVc1_dV1a*(c1->Pressure+T->P_Lumen_a);
                Fa_v2 += dVc1_dV2a*(c1->Pressure+T->P_Lumen_a);
            }
            
            if(T->K_Lumen_b==0 && T->P0_Lumen_b==0) {// only the pressure from the cells is taken into account
                Fb_v1 += dVc1_dV1b*c1->Pressure;
                Fb_v2 += dVc1_dV2b*c1->Pressure;
            } else {
                Fb_v1 += dVc1_dV1b*(c1->Pressure-T->P_Lumen_b);
                Fb_v2 += dVc1_dV2b*(c1->Pressure-T->P_Lumen_b);
            }
        }
        
        if(c2) {
            
            Point dVc2_dV1a, dVc2_dV1b, dVc2_dV2a, dVc2_dV2b;
            
            // pure contributions from apical and basal barycenter, and from the lateral barycenter
            int n2 = c2->ListEdges.size();
            
            double T_a = c2->T_a();
            double T_b = c2->T_b();
            
            // NOTE: every vertex is once vertex one with respect to a cell
            // therefore add the "force contribution" coming from the movement in the apical and basal barycenters of v1 wrt c1 and of v2 wrt to c2 (since v1 is first in direction of c1 and v2 to c2)
            v2->dW_dva -= c2->dAa_dBCa*(T_a / c2->ListEdges.size());
            v2->dW_dvb -= c2->dAb_dBCb*(T_b / c2->ListEdges.size());
            
            dVc2_dV2a += c2->dV_dBCa/n2;
            dVc2_dV2b += c2->dV_dBCb/n2;
            
            // APICAL AND BASAL SURFACE TENSION
            // force on vertices
            Fa_v1 -= Point::TriangleAreaChange3(v1a,v2a,c2->BC_a,n2) * T_a;
            Fb_v1 -= Point::TriangleAreaChange3(v1b,v2b,c2->BC_b,n2) * T_b;
            Fa_v2 -= Point::TriangleAreaChange3(v2a,v1a,c2->BC_a,n2) * T_a;
            Fb_v2 -= Point::TriangleAreaChange3(v2b,v1b,c2->BC_b,n2) * T_b;
            
            // lateral contributions to volume changes
            dVc2_dV1a += ((Point::dV_dp(BC_l, v2->coord_b, v2->coord)+Point::dV_dp(BC_l, v1->coord_b, v2->coord_b))/4 + Point::dV_dp_b(v1->coord, v1->coord_b, BC_l, 4) - Point::dV_dp_b(v1->coord, v2->coord, BC_l,4))*(-1);
            dVc2_dV1b += ((Point::dV_dp(BC_l, v2->coord, v1->coord)+Point::dV_dp(BC_l, v2->coord_b, v2->coord))/4 + Point::dV_dp_b(v1->coord_b, v2->coord_b, BC_l, 4) - Point::dV_dp_b(v1->coord_b, v1->coord, BC_l,4))*(-1);
            dVc2_dV2a += ((Point::dV_dp(BC_l, v1->coord_b, v2->coord_b)+Point::dV_dp(BC_l, v1->coord, v1->coord_b))/4 + Point::dV_dp_b(v2->coord, v1->coord, BC_l, 4) - Point::dV_dp_b(v2->coord, v2->coord_b, BC_l,4))*(-1);
            dVc2_dV2b += ((Point::dV_dp(BC_l, v1->coord, v1->coord_b)+Point::dV_dp(BC_l, v2->coord, v1->coord))/4 + Point::dV_dp_b(v2->coord_b, v2->coord, BC_l, 4) - Point::dV_dp_b(v2->coord_b, v1->coord_b, BC_l,4))*(-1);
            
            // apical and basal contributions to volume changes
            dVc2_dV1a -= Point::dV_dp_b_1(v1->coord, v2->coord, c2->BC_a, n2);
            dVc2_dV1b += Point::dV_dp_b_1(v1->coord_b, v2->coord_b, c2->BC_b, n2);
            dVc2_dV2a += Point::dV_dp_b_1(v2->coord, v1->coord, c2->BC_a, n2);
            dVc2_dV2b -= Point::dV_dp_b_1(v2->coord_b, v1->coord_b, c2->BC_b, n2);
            
            // use the volume derivatives to calculate the forces acting on the vertices
            if(T->K_Lumen_a==0 && T->P0_Lumen_a==0) {// only the pressure from the cells is taken into account
                Fa_v1 += dVc2_dV1a*c2->Pressure;
                Fa_v2 += dVc2_dV2a*c2->Pressure;
            } else {
                Fa_v1 += dVc2_dV1a*(c2->Pressure+T->P_Lumen_a);
                Fa_v2 += dVc2_dV2a*(c2->Pressure+T->P_Lumen_a);
            }
            
            if(T->K_Lumen_b==0 && T->P0_Lumen_b==0) {// only the pressure from the cells is taken into account
                Fb_v1 += dVc2_dV1b*c2->Pressure;
                Fb_v2 += dVc2_dV2b*c2->Pressure;
            } else {
                Fb_v1 += dVc2_dV1b*(c2->Pressure-T->P_Lumen_b);
                Fb_v2 += dVc2_dV2b*(c2->Pressure-T->P_Lumen_b);
            }

        }
        
        // ADD UP THE FORCE CONTRIBUTIONS IN THE 4 VERTEX POSITIONS THAT BELONG TO THE EDGE
        v1->dW_dva += Fa_v1;
        v2->dW_dva += Fa_v2;
        v1->dW_dvb += Fb_v1;
        v2->dW_dvb += Fb_v2;
    }
    
    else { // periodic tissue
        
        Point dVc1_dV1a, dVc2_dV1a, dVc1_dV1b, dVc2_dV1b, dVc1_dV2a, dVc2_dV2a, dVc1_dV2b, dVc2_dV2b;
        
        // pure contributions from apical and basal barycenter, and from the lateral barycenter
        int n1 = c1->ListEdges.size();
        int n2 = c2->ListEdges.size();
        
        
        double Tl = T_l;
        // calculate the effective current tensions dW/dAa and dW/dAb for both cells
        double c1_Ta = c1->T_a();
        double c2_Ta = c2->T_a();
        double c1_Tb = c1->T_b();
        double c2_Tb = c2->T_b();
    
        // NOTE: every vertex is once vertex one with respect to a cell
        // therefore add the "force contribution" coming from the movement in the apical and basal barycenters of v1 wrt c1 and of v2 wrt to c2 (since v1 is first in direction of c1 and v2 to c2)
        v1->dW_dva -= c1->dAa_dBCa*(c1_Ta / c1->ListEdges.size());
        v2->dW_dva -= c2->dAa_dBCa*(c2_Ta / c2->ListEdges.size());
        v1->dW_dvb -= c1->dAb_dBCb*(c1_Tb / c1->ListEdges.size());
        v2->dW_dvb -= c2->dAb_dBCb*(c2_Tb / c2->ListEdges.size());
        
        dVc1_dV1a += c1->dV_dBCa/n1;
        dVc1_dV1b += c1->dV_dBCb/n1;
        dVc2_dV2a += c2->dV_dBCa/n2;
        dVc2_dV2b += c2->dV_dBCb/n2;

        // systemSize
        Quadrant SystemSize = T->SystemSize;
        
        bool crossBoundary_c1 = c1->CrossBoundary();
        bool crossBoundary_c2 = c2->CrossBoundary();
        
        if(!CrossBoundary()) { // edge doesnt cross the boundary
            // APICAL AND BASAL LINE TENSION
            // force on vertices and system size
            if(la>l0_a) {
                Fa_v1 += (v2a-v1a)*(G_a/la);
                Fa_v2 += Fa_v1*(-1);
                T->dW_dL -= Point::lineStress(v1a,v2a)*G_a;
            } else {
                Point d = v2a-v1a;
                Fa_v1 += d*(G_a/l0_a);
                Fa_v2 += Fa_v1*(-1);
                T->dW_dL -= Point(pow(d.x,2), pow(d.y,2),0)*(G_a/l0_a);
            }
            
            if(lb>l0_b) {
                Fb_v1 += (v2b-v1b)*(G_b/lb);
                Fb_v2 += Fb_v1*(-1);
                T->dW_dL -= Point::lineStress(v1b,v2b)*G_a;
            } else {
                Point d = v2b-v1b;
                Fb_v1 += d*(G_b/l0_b);
                Fb_v2 += Fb_v1*(-1);
                T->dW_dL -= Point(pow(d.x,2), pow(d.y,2),0)*(G_b/l0_b);
            }
           
            // LATERAL SURFACE TENSION
            // force on vertices
            Fa_v1 -= (Point::TriangleAreaChange2(v1a,v2a,b,4)+Point::TriangleAreaChange2(v1a,v1b,b,4)+Point::TriangleAreaChange(b,v1b,v2b)/4+Point::TriangleAreaChange(b,v2a,v2b)/4)*Tl;
            Fa_v2 -= (Point::TriangleAreaChange2(v2a,v1a,b,4)+Point::TriangleAreaChange2(v2a,v2b,b,4)+Point::TriangleAreaChange(b,v1a,v1b)/4+Point::TriangleAreaChange(b,v1b,v2b)/4)*Tl;
            Fb_v1 -= (Point::TriangleAreaChange2(v1b,v1a,b,4)+Point::TriangleAreaChange2(v1b,v2b,b,4)+Point::TriangleAreaChange(b,v1a,v2a)/4+Point::TriangleAreaChange(b,v2a,v2b)/4)*Tl;
            Fb_v2 -= (Point::TriangleAreaChange2(v2b,v2a,b,4)+Point::TriangleAreaChange2(v2b,v1b,b,4)+Point::TriangleAreaChange(b,v1a,v2a)/4+Point::TriangleAreaChange(b,v1a,v1b)/4)*Tl;
            // force on the system size
            T->dW_dL -= (Point::surfaceStress(v1a,v2a,b) + Point::surfaceStress(v2a,v2b,b ) + Point::surfaceStress(v2b,v1b,b ) + Point::surfaceStress(v1b,v1a,b ))* T_l;
            
        } else { // edge crosses boundary
            // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE FIRST VERTEX OF THE EDGE
            Point v1a_rel = v1a;
            Point v2a_rel = shift( v2a, SystemSize*q_v2);

            Point v1b_rel = shift( v1b, SystemSize*v1->q_b);
            Point v2b_rel = shift( v2b, SystemSize*(q_v2+v2->q_b));
            Point b_rel =   shift( b, SystemSize*q_BC_l);
            
            // APICAL AND BASAL LINE TENSION
            // force on vertices
            if(la>l0_a) {
                Fa_v1 += (v2a_rel-v1a_rel)*(G_a/la);
                Fa_v2 += Fa_v1*(-1);
                T->dW_dL -= Point::lineStress(v1a_rel,v2a_rel)*Lambda_a();
            } else {
                Point d = v2a_rel-v1a_rel;
                Fa_v1 += d*(G_a/l0_a);
                Fa_v2 += Fa_v1*(-1);
                T->dW_dL -= Point(pow(d.x,2), pow(d.y,2),0)*(G_a/l0_a);
            }
            
            if(lb>l0_b) {
                Fb_v1 += (v2b_rel-v1b_rel)*(G_b/lb);
                Fb_v2 += Fb_v1*(-1);
                T->dW_dL -= Point::lineStress(v1b_rel,v2b_rel)*Lambda_b();
            } else {
                Point d = v2b_rel-v1b_rel;
                Fb_v1 += d*(G_b/l0_b);
                Fb_v2 += Fb_v1*(-1);
                T->dW_dL -= Point(pow(d.x,2), pow(d.y,2),0)*(G_b/l0_b);
            }

            
            // LATERAL SURFACE TENSION
            // force on vertices
            Fa_v1 -= (Point::TriangleAreaChange2(v1a_rel,v2a_rel,b_rel,4)+Point::TriangleAreaChange2(v1a_rel,v1b_rel,b_rel,4)+Point::TriangleAreaChange(b_rel,v1b_rel,v2b_rel)/4+Point::TriangleAreaChange(b_rel,v2a_rel,v2b_rel)/4)*Tl;
            Fa_v2 -= (Point::TriangleAreaChange2(v2a_rel,v1a_rel,b_rel,4)+Point::TriangleAreaChange2(v2a_rel,v2b_rel,b_rel,4)+Point::TriangleAreaChange(b_rel,v1a_rel,v1b_rel)/4+Point::TriangleAreaChange(b_rel,v1b_rel,v2b_rel)/4)*Tl;
            Fb_v1 -= (Point::TriangleAreaChange2(v1b_rel,v1a_rel,b_rel,4)+Point::TriangleAreaChange2(v1b_rel,v2b_rel,b_rel,4)+Point::TriangleAreaChange(b_rel,v1a_rel,v2a_rel)/4+Point::TriangleAreaChange(b_rel,v2a_rel,v2b_rel)/4)*Tl;
            Fb_v2 -= (Point::TriangleAreaChange2(v2b_rel,v2a_rel,b_rel,4)+Point::TriangleAreaChange2(v2b_rel,v1b_rel,b_rel,4)+Point::TriangleAreaChange(b_rel,v1a_rel,v2a_rel)/4+Point::TriangleAreaChange(b_rel,v1a_rel,v1b_rel)/4)*Tl;
            
            // force on the system size
            T->dW_dL -= (Point::surfaceStress(v1a_rel,v2a_rel,b_rel ) + Point::surfaceStress(v2a_rel,v2b_rel,b_rel ) + Point::surfaceStress(v2b_rel,v1b_rel,b_rel ) + Point::surfaceStress(v1b_rel,v1a_rel,b_rel ))*  T_l;
        }
        
        if(!crossBoundary_c1)  {
            double T_a = c1_Ta;
            double T_b = c1_Tb;
            
            // APICAL AND BASAL SURFACE TENSIONS
            Fa_v1 -= Point::TriangleAreaChange3(v1a,v2a,c1->BC_a,n1) * T_a;
            Fb_v1 -= Point::TriangleAreaChange3(v1b,v2b,c1->BC_b,n1) * T_b;
            Fa_v2 -= Point::TriangleAreaChange3(v2a,v1a,c1->BC_a,n1) * T_a;
            Fb_v2 -= Point::TriangleAreaChange3(v2b,v1b,c1->BC_b,n1) * T_b;
            // force on system size
            T->dW_dL -= (Point::surfaceStress(v1a, c1->BC_a, v2a)) * T_a + (Point::surfaceStress(v1b, c1->BC_b, v2b)) * T_b;
            
            // lateral contributions to volume changes
            dVc1_dV1a += (Point::dV_dp(BC_l, v2->coord_b, v2->coord)+Point::dV_dp(BC_l, v1->coord_b, v2->coord_b))/4 + Point::dV_dp_b(v1->coord, v1->coord_b, BC_l, 4) - Point::dV_dp_b(v1->coord, v2->coord, BC_l,4);
            dVc1_dV1b += (Point::dV_dp(BC_l, v2->coord, v1->coord)+Point::dV_dp(BC_l, v2->coord_b, v2->coord))/4 + Point::dV_dp_b(v1->coord_b, v2->coord_b, BC_l, 4) - Point::dV_dp_b(v1->coord_b, v1->coord, BC_l,4);
            dVc1_dV2a += (Point::dV_dp(BC_l, v1->coord_b, v2->coord_b)+Point::dV_dp(BC_l, v1->coord, v1->coord_b))/4 + Point::dV_dp_b(v2->coord, v1->coord, BC_l, 4) - Point::dV_dp_b(v2->coord, v2->coord_b, BC_l,4);
            dVc1_dV2b += (Point::dV_dp(BC_l, v1->coord, v1->coord_b)+Point::dV_dp(BC_l, v2->coord, v1->coord))/4 + Point::dV_dp_b(v2->coord_b, v2->coord, BC_l, 4) - Point::dV_dp_b(v2->coord_b, v1->coord_b, BC_l,4);
            
            // apical and basal contributions to volume changes
            dVc1_dV1a += Point::dV_dp_b_1(v1->coord, v2->coord, c1->BC_a, n1);
            dVc1_dV1b -= Point::dV_dp_b_1(v1->coord_b, v2->coord_b, c1->BC_b, n1);
            dVc1_dV2a -= Point::dV_dp_b_1(v2->coord, v1->coord, c1->BC_a, n1);
            dVc1_dV2b += Point::dV_dp_b_1(v2->coord_b, v1->coord_b, c1->BC_b, n1);
            
        } else {
            
            double T_a = c1_Ta;
            double T_b = c1_Tb;
            
            
            // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE APICAL BARYCENTER OF CELL 1
            Point BC_l_rel = shift( b, SystemSize*( q_BC_l -  q_c1));
            Point v1a_rel = shift( v1a, SystemSize* q_c1*(-1));
            Point v2a_rel = shift( v2a, SystemSize*( q_v2 -  q_c1));
            
            Point v1b_rel = shift( v1b, SystemSize*(v1->q_b -  q_c1));
            Point v2b_rel = shift( v2b, SystemSize*( q_v2 + v2->q_b -  q_c1));
            Point BC_b_rel = shift( c1->BC_b, SystemSize* c1->q_BC_b);
            
            // APICAL AND BASAL SURFACE TENSION
            // force on vertices
            Fa_v1 -= Point::TriangleAreaChange3(v1a_rel,v2a_rel,c1->BC_a,n1) * T_a;
            Fb_v1 -= Point::TriangleAreaChange3(v1b_rel,v2b_rel,BC_b_rel,n1) * T_b;
            Fa_v2 -= Point::TriangleAreaChange3(v2a_rel,v1a_rel,c1->BC_a,n1) * T_a;
            Fb_v2 -= Point::TriangleAreaChange3(v2b_rel,v1b_rel,BC_b_rel,n1) * T_b;
            // force on system size
            T->dW_dL -= (Point::surfaceStress(v1a_rel, c1->BC_a, v2a_rel)) * T_a + (Point::surfaceStress(v1b_rel, BC_b_rel, v2b_rel)) * T_b;
            
            // lateral contributions to volume changes
            dVc1_dV1a += (Point::dV_dp(BC_l_rel, v2b_rel, v2a_rel)+Point::dV_dp(BC_l_rel, v1b_rel, v2b_rel))/4 + Point::dV_dp_b(v1a_rel, v1b_rel, BC_l_rel, 4) - Point::dV_dp_b(v1a_rel, v2a_rel, BC_l_rel,4);
            dVc1_dV1b += (Point::dV_dp(BC_l_rel, v2a_rel, v1a_rel)+Point::dV_dp(BC_l_rel, v2b_rel, v2a_rel))/4 + Point::dV_dp_b(v1b_rel, v2b_rel, BC_l_rel, 4) - Point::dV_dp_b(v1b_rel, v1a_rel, BC_l_rel,4);
            dVc1_dV2a += (Point::dV_dp(BC_l_rel, v1b_rel, v2b_rel)+Point::dV_dp(BC_l_rel, v1a_rel, v1b_rel))/4 + Point::dV_dp_b(v2a_rel, v1a_rel, BC_l_rel, 4) - Point::dV_dp_b(v2a_rel, v2b_rel, BC_l_rel,4);
            dVc1_dV2b += (Point::dV_dp(BC_l_rel, v1a_rel, v1b_rel)+Point::dV_dp(BC_l_rel, v2a_rel, v1a_rel))/4 + Point::dV_dp_b(v2b_rel, v2a_rel, BC_l_rel, 4) - Point::dV_dp_b(v2b_rel, v1b_rel, BC_l_rel,4);
            
            // apical and basal contributions to volume changes
            dVc1_dV1a += Point::dV_dp_b_1(v1a_rel, v2a_rel, c1->BC_a, n1);
            dVc1_dV1b -= Point::dV_dp_b_1(v1b_rel, v2b_rel, BC_b_rel, n1);
            dVc1_dV2a -= Point::dV_dp_b_1(v2a_rel, v1a_rel, c1->BC_a, n1);
            dVc1_dV2b += Point::dV_dp_b_1(v2b_rel, v1b_rel, BC_b_rel, n1);
            
        }
        
        if(!crossBoundary_c2) {
            //save time by computing the current surface tension only once
            double T_a = c2_Ta;
            double T_b = c2_Tb;
            
            // APICAL AND BASAL SURFACE TENSION
            // force on vertices
            Fa_v1 -= Point::TriangleAreaChange3(v1a,v2a,c2->BC_a,n2) * T_a;
            Fb_v1 -= Point::TriangleAreaChange3(v1b,v2b,c2->BC_b,n2) * T_b;
            Fa_v2 -= Point::TriangleAreaChange3(v2a,v1a,c2->BC_a,n2) * T_a;
            Fb_v2 -= Point::TriangleAreaChange3(v2b,v1b,c2->BC_b,n2) * T_b;
            // force on system size
            T->dW_dL -= (Point::surfaceStress(v1a, c2->BC_a, v2a)) * T_a + (Point::surfaceStress(v1b, c2->BC_b, v2b)) * T_b;
            
            // lateral contributions to volume changes
            dVc2_dV1a += ((Point::dV_dp(BC_l, v2->coord_b, v2->coord)+Point::dV_dp(BC_l, v1->coord_b, v2->coord_b))/4 + Point::dV_dp_b(v1->coord, v1->coord_b, BC_l, 4) - Point::dV_dp_b(v1->coord, v2->coord, BC_l,4))*(-1);
            dVc2_dV1b += ((Point::dV_dp(BC_l, v2->coord, v1->coord)+Point::dV_dp(BC_l, v2->coord_b, v2->coord))/4 + Point::dV_dp_b(v1->coord_b, v2->coord_b, BC_l, 4) - Point::dV_dp_b(v1->coord_b, v1->coord, BC_l,4))*(-1);
            dVc2_dV2a += ((Point::dV_dp(BC_l, v1->coord_b, v2->coord_b)+Point::dV_dp(BC_l, v1->coord, v1->coord_b))/4 + Point::dV_dp_b(v2->coord, v1->coord, BC_l, 4) - Point::dV_dp_b(v2->coord, v2->coord_b, BC_l,4))*(-1);
            dVc2_dV2b += ((Point::dV_dp(BC_l, v1->coord, v1->coord_b)+Point::dV_dp(BC_l, v2->coord, v1->coord))/4 + Point::dV_dp_b(v2->coord_b, v2->coord, BC_l, 4) - Point::dV_dp_b(v2->coord_b, v1->coord_b, BC_l,4))*(-1);
            
            // apical and basal contributions to volume changes
            dVc2_dV1a -= Point::dV_dp_b_1(v1->coord, v2->coord, c2->BC_a, n2);
            dVc2_dV1b += Point::dV_dp_b_1(v1->coord_b, v2->coord_b, c2->BC_b, n2);
            dVc2_dV2a += Point::dV_dp_b_1(v2->coord, v1->coord, c2->BC_a, n2);
            dVc2_dV2b -= Point::dV_dp_b_1(v2->coord_b, v1->coord_b, c2->BC_b, n2);
            
        } else {
            double T_a = c2_Ta;
            double T_b = c2_Tb;            
            
            // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE APICAL BARYCENTER OF CELL 2
            Point BC_l_rel = shift( b, SystemSize*( q_BC_l -  q_c2 -  q_v2));
            
            Point v1a_rel = shift( v1a, SystemSize*(  q_c2 +  q_v2 ) *(-1));
            Point v2a_rel = shift( v2a, SystemSize* q_c2*(-1));

            Point v1b_rel = shift( v1b, SystemSize*(v1->q_b -  q_c2 -  q_v2));
            Point v2b_rel = shift( v2b, SystemSize*(  q_c2 - v2->q_b)*(-1));
            Point BC_b_rel = shift( c2->BC_b, SystemSize*c2->q_BC_b );
                
            // APICAL AND BASAL SURFACE TENSION
            // force on vertices
            Fa_v1 -= Point::TriangleAreaChange3(v1a_rel,v2a_rel,c2->BC_a,n2) * T_a;
            Fb_v1 -= Point::TriangleAreaChange3(v1b_rel,v2b_rel,BC_b_rel,n2) * T_b;
            Fa_v2 -= Point::TriangleAreaChange3(v2a_rel,v1a_rel,c2->BC_a,n2) * T_a;
            Fb_v2 -= Point::TriangleAreaChange3(v2b_rel,v1b_rel,BC_b_rel,n2) * T_b;
            // force on system size
            T->dW_dL -= (Point::surfaceStress(v1a_rel, c2->BC_a, v2a_rel)) * T_a + (Point::surfaceStress(v1b_rel, BC_b_rel, v2b_rel)) * T_b;
            
            // lateral contributions to volume changes
            dVc2_dV1a += ((Point::dV_dp(BC_l_rel, v2b_rel, v2a_rel)+Point::dV_dp(BC_l_rel, v1b_rel, v2b_rel))/4 + Point::dV_dp_b(v1a_rel, v1b_rel, BC_l_rel, 4) - Point::dV_dp_b(v1a_rel, v2a_rel, BC_l_rel,4))*(-1);
            dVc2_dV1b += ((Point::dV_dp(BC_l_rel, v2a_rel, v1a_rel)+Point::dV_dp(BC_l_rel, v2b_rel, v2a_rel))/4 + Point::dV_dp_b(v1b_rel, v2b_rel, BC_l_rel, 4) - Point::dV_dp_b(v1b_rel, v1a_rel, BC_l_rel,4))*(-1);
            dVc2_dV2a += ((Point::dV_dp(BC_l_rel, v1b_rel, v2b_rel)+Point::dV_dp(BC_l_rel, v1a_rel, v1b_rel))/4 + Point::dV_dp_b(v2a_rel, v1a_rel, BC_l_rel, 4) - Point::dV_dp_b(v2a_rel, v2b_rel, BC_l_rel,4))*(-1);
            dVc2_dV2b += ((Point::dV_dp(BC_l_rel, v1a_rel, v1b_rel)+Point::dV_dp(BC_l_rel, v2a_rel, v1a_rel))/4 + Point::dV_dp_b(v2b_rel, v2a_rel, BC_l_rel, 4) - Point::dV_dp_b(v2b_rel, v1b_rel, BC_l_rel,4))*(-1);
            
            // apical and basal contributions to volume changes
            dVc2_dV1a -= Point::dV_dp_b_1(v1a_rel, v2a_rel, c2->BC_a, n2);
            dVc2_dV1b += Point::dV_dp_b_1(v1b_rel, v2b_rel, BC_b_rel, n2);
            dVc2_dV2a += Point::dV_dp_b_1(v2a_rel, v1a_rel, c2->BC_a, n2);
            dVc2_dV2b -= Point::dV_dp_b_1(v2b_rel, v1b_rel, BC_b_rel, n2);
        }
        
        // use the volume derivatives to calculate the forces acting on the vertices
        Fa_v1 += dVc1_dV1a*c1->Pressure + dVc2_dV1a*c2->Pressure;
        Fa_v2 += dVc1_dV2a*c1->Pressure + dVc2_dV2a*c2->Pressure;
        Fb_v1 += dVc1_dV1b*c1->Pressure + dVc2_dV1b*c2->Pressure;
        Fb_v2 += dVc1_dV2b*c1->Pressure + dVc2_dV2b*c2->Pressure;
        
        // ADD UP THE FORCE CONTRIBUTIONS IN THE 4 VERTEX POSITIONS THAT BELONG TO THE EDGE
        v1->dW_dva += Fa_v1;
        v2->dW_dva += Fa_v2;
        v1->dW_dvb += Fb_v1;
        v2->dW_dvb += Fb_v2;
    }
}


//void Edge::update_energy()
//{
//    if(!T->isPeriodic)
//    {
//        // positions
//        Point v1a = v1->coord;
//        Point v1b = v1->coord_b;
//        Point v2a = v2->coord;
//        Point v2b = v2->coord_b;
//        Point b = BC_l;
//        
//        la = Point::Norm(v1a-v2a);
//        lb = Point::Norm(v1b-v2b);
//        
//        T->Energy_Surface += T_l*(Point::TriangleArea(v1a,v2a,b)+Point::TriangleArea(v2a,v2b,b)+Point::TriangleArea(v2b,v1b,b)+Point::TriangleArea(v1b,v1a,b));
//        T->Energy_Line += Lambda_a()*la + Lambda_b()*lb;
//        
//        if(c1) T->Energy_Surface += c1->T_a()*(Point::TriangleArea(v1a,v2a,c1->BC_a)) + c1->T_b()*(Point::TriangleArea(v1b,v2b,c1->BC_b));
//        if(c2) T->Energy_Surface += c2->T_a()*(Point::TriangleArea(v1a,v2a,c2->BC_a)) + c2->T_b()*(Point::TriangleArea(v1b,v2b,c2->BC_b));
//
//    }
//    else
//    {
//        Quadrant SystemSize = T->SystemSize;
//        
//        // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE APICAL BARYCENTER OF CELL 1
//        Point v1_rel = shift( v1->coord, SystemSize* q_c1*(-1));
//        Point v1b_rel = shift( v1->coord_b, SystemSize*(v1->q_b -  q_c1));
//        Point v2_rel = shift( v2->coord, SystemSize*( q_v2 -  q_c1));
//        Point v2b_rel = shift( v2->coord_b, SystemSize*( q_v2 + v2->q_b -  q_c1));
//        Point BC_b_rel = shift( c1->BC_b, SystemSize* c1->q_BC_b);
//        Point BC_l_rel = shift( BC_l, SystemSize*( q_BC_l -  q_c1));
//
//        la = Point::Norm(v1_rel-v2_rel);
//        lb = Point::Norm(v1b_rel-v2b_rel);
//        
//        // energy contribution from lateral surface and lines
//        T->Energy_Line += Lambda_a()*la + Lambda_b()*lb;
//        T->Energy_Surface += T_l*(Point::TriangleArea(v1_rel,v2_rel,BC_l_rel)+Point::TriangleArea(v2_rel,v2b_rel,BC_l_rel)+Point::TriangleArea(v2b_rel,v1b_rel,BC_l_rel)+Point::TriangleArea(v1b_rel,v1_rel,BC_l_rel));
//        
//        // energy contribution from the apical and basal surface of cell 1
//        T->Energy_Surface += c1->T_a()*(Point::TriangleArea(v1_rel,v2_rel,c1->BC_a)) + c1->T_b()*(Point::TriangleArea(v1b_rel,v2b_rel,BC_b_rel));
//        
//        // define all the points in the same relation to the center of cell number 2
//        // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE APICAL BARYCENTER OF CELL 2
//        v1_rel = shift( v1->coord, SystemSize*(  q_c2 +  q_v2 ) *(-1));
//        v1b_rel = shift( v1->coord_b, SystemSize*(v1->q_b -  q_c2 -  q_v2));
//        v2_rel = shift( v2->coord, SystemSize* q_c2*(-1));
//        v2b_rel = shift( v2->coord_b, SystemSize*(  q_c2 - v2->q_b)*(-1));
//        BC_b_rel = shift( c2->BC_b, SystemSize*c2->q_BC_b );
//
//        T->Energy_Surface += c2->T_a()*(Point::TriangleArea(v1_rel,v2_rel,c2->BC_a)) + c2->T_b()*(Point::TriangleArea(v1b_rel,v2b_rel,BC_b_rel));
//    }
//}

double Edge::W()
{
    double energy = 0;
    
    if(!T->isPeriodic) {
        // positions
        Point v1a = v1->coord;
        Point v1b = v1->coord_b;
        Point v2a = v2->coord;
        Point v2b = v2->coord_b;
        Point b = BC_l;
        
        la = Point::Norm(v1a-v2a);
        lb = Point::Norm(v1b-v2b);
        
        // energy contribution from apical and basal line tensions
        if(l0_a==0) {
            energy += la * G_a;
        }
        else if(la>l0_a) {
            energy += l0_a*G_a/2. + (la-l0_a)*G_a;
        } else {
            energy += la*la/2. * G_a/l0_a;
        };
        
        if(l0_b==0) {
            energy += lb * G_b;
        }
        else if(lb>l0_b) {
            energy += l0_b*G_b/2. + (lb-l0_b)*G_b;
        } else {
            energy += lb*lb/2. * G_b/l0_b;
        }
        
        if(G_a_ai>0) {
            double d = v1a.z-v2a.z;
            if(d>0) energy += G_a_ai*d;
            else energy -= G_a_ai*d;
        }
        
        energy += T_l*(Point::TriangleArea(v1a,v2a,b)+Point::TriangleArea(v2a,v2b,b)+Point::TriangleArea(v2b,v1b,b)+Point::TriangleArea(v1b,v1a,b));
        
//        if(anisotropicTension)
//        {
//            if(c1->Type == 2 || c2->Type == 2)
//            {
//                
//                
//                
//            }
//        }
    }
    else  {
        Quadrant SystemSize = T->SystemSize;
        
        // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE APICAL BARYCENTER OF CELL 1
        Point v1_rel = shift( v1->coord, SystemSize* q_c1*(-1));
        Point v1b_rel = shift( v1->coord_b, SystemSize*(v1->q_b -  q_c1));
        Point v2_rel = shift( v2->coord, SystemSize*( q_v2 -  q_c1));
        Point v2b_rel = shift( v2->coord_b, SystemSize*( q_v2 + v2->q_b -  q_c1));
        Point BC_l_rel = shift( BC_l, SystemSize*( q_BC_l -  q_c1));
        
        la = Point::Norm(v1_rel-v2_rel);
        lb = Point::Norm(v1b_rel-v2b_rel);
        
        // energy contribution from apical and basal line tensions
        if(la>l0_a) {
            energy += l0_a*G_a/2. + (la-l0_a)*G_a;
        } else {
            energy += la*la/2. * G_a/l0_a;
        };
        
        if(lb>l0_b) {
            energy += l0_b*G_b/2. + (lb-l0_b)*G_b;
        } else {
            energy += lb*lb/2. * G_b/l0_b;
        }

        // energy contribution from lateral surface tension
        double Al = Point::TriangleArea(v1_rel,v2_rel,BC_l_rel)+Point::TriangleArea(v2_rel,v2b_rel,BC_l_rel)+Point::TriangleArea(v2b_rel,v1b_rel,BC_l_rel)+Point::TriangleArea(v1b_rel,v1_rel,BC_l_rel); // lateral surface area
        energy += T_l*Al;
    }
    
    return energy;
}

double Edge::energy(Cell* c)
{
    UpdateLengths();
    
    double E = 0;
    
    if(!T->isPeriodic)
    {
        // positions
        Point v1a = v1->coord;
        Point v1b = v1->coord_b;
        Point v2a = v2->coord;
        Point v2b = v2->coord_b;
        Point b = BC_l;
        
        la = Point::Norm(v1a-v2a);
        lb = Point::Norm(v1b-v2b);
        
        E += T_l/2.*(Point::TriangleArea(v1a,v2a,b)+Point::TriangleArea(v2a,v2b,b)+Point::TriangleArea(v2b,v1b,b)+Point::TriangleArea(v1b,v1a,b));
        E += Lambda_a()*la + Lambda_b()*lb;
        
        if(c1==c) E += c1->T_a()*(Point::TriangleArea(v1a,v2a,c1->BC_a)) + c1->T_b()*(Point::TriangleArea(v1b,v2b,c1->BC_b));
        if(c2==c) E += c2->T_a()*(Point::TriangleArea(v1a,v2a,c2->BC_a)) + c2->T_b()*(Point::TriangleArea(v1b,v2b,c2->BC_b));
        
    }
    else
    {
        Quadrant SystemSize = T->SystemSize;
        
        // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE APICAL BARYCENTER OF CELL 1
        Point v1_rel = shift( v1->coord, SystemSize* q_c1*(-1));
        Point v1b_rel = shift( v1->coord_b, SystemSize*(v1->q_b -  q_c1));
        Point v2_rel = shift( v2->coord, SystemSize*( q_v2 -  q_c1));
        Point v2b_rel = shift( v2->coord_b, SystemSize*( q_v2 + v2->q_b -  q_c1));
        Point BC_b_rel = shift( c1->BC_b, SystemSize* c1->q_BC_b);
        Point BC_l_rel = shift( BC_l, SystemSize*( q_BC_l -  q_c1));
        
        la = Point::Norm(v1_rel-v2_rel);
        lb = Point::Norm(v1b_rel-v2b_rel);
        
        // energy contribution from lateral surface and lines
        E += Lambda_a()*la + Lambda_b()*lb;
        E += T_l/2.*(Point::TriangleArea(v1_rel,v2_rel,BC_l_rel)+Point::TriangleArea(v2_rel,v2b_rel,BC_l_rel)+Point::TriangleArea(v2b_rel,v1b_rel,BC_l_rel)+Point::TriangleArea(v1b_rel,v1_rel,BC_l_rel));
        
        
        // energy contribution from the apical and basal surface of cell 1
        if(c1==c) E += c1->T_a()*(Point::TriangleArea(v1_rel,v2_rel,c1->BC_a)) + c1->T_b()*(Point::TriangleArea(v1b_rel,v2b_rel,BC_b_rel));
        
        // define all the points in the same relation to the center of cell number 2
        // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE APICAL BARYCENTER OF CELL 2
        v1_rel = shift( v1->coord, SystemSize*(  q_c2 +  q_v2 ) *(-1));
        v1b_rel = shift( v1->coord_b, SystemSize*(v1->q_b -  q_c2 -  q_v2));
        v2_rel = shift( v2->coord, SystemSize* q_c2*(-1));
        v2b_rel = shift( v2->coord_b, SystemSize*(  q_c2 - v2->q_b)*(-1));
        BC_b_rel = shift( c2->BC_b, SystemSize*c2->q_BC_b );
        
        if(c2==c) E += c2->T_a()*(Point::TriangleArea(v1_rel,v2_rel,c2->BC_a)) + c2->T_b()*(Point::TriangleArea(v1b_rel,v2b_rel,BC_b_rel));
    }
    
    return E;
}

double Edge::SurroundingEnergy()
{
//    std::cout << "v1->SurroundingEnergy() = " << v1->SurroundingEnergy() << std::endl;
//    std::cout << "v2->SurroundingEnergy() = " << v2->SurroundingEnergy() << std::endl;
//    std::cout << "c1->energy() = " << c1->energy() << std::endl;
//    std::cout << "c2->energy() = " << c2->energy() << std::endl;
    
    return v1->SurroundingEnergy()+v2->SurroundingEnergy() - c1->energy() - c2->energy();
}

//Point Edge::ForcesFromLateralSurfaceTension_numerically()
//{
//    // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO THE FIRST VERTEX OF THE EDGE
//    Point v1_rel = v1->coord;
//    Point v1b_rel = shift( v1->coord_b, T->SystemSize*v1->q_b);
//    Point v2_rel = shift( v2->coord, T->SystemSize*q_v2);
//    Point v2b_rel = shift( v2->coord_b, T->SystemSize*(q_v2+v2->q_b));
//    Point BC_l_rel =   shift( BC_l, T->SystemSize*q_BC_l);
//    
//    double E0 = T_l*(Point::TriangleArea(v1_rel,v2_rel,BC_l_rel)+Point::TriangleArea(v2_rel,v2b_rel,BC_l_rel)+Point::TriangleArea(v2b_rel,v1b_rel,BC_l_rel)+Point::TriangleArea(v1b_rel,v1_rel,BC_l_rel));
//    
//    Point F1a;
//    
//    double d = 10e-5;
//    double d_2 = d/2;
//    
//    Point P_x_d_2(d_2,0,0),P_x_d(-d,0,0),P_y_d_2(0,d_2,0),P_y_d(0,-d,0),P_z_d_2(0,0,d_2),P_z_d(0,0,-d),P_null(0,0,0);
//    
//    v1_rel += P_x_d_2;
//    BC_l_rel += P_x_d_2/4;
//    double E_p1 = T_l*(Point::TriangleArea(v1_rel,v2_rel,BC_l_rel)+Point::TriangleArea(v2_rel,v2b_rel,BC_l_rel)+Point::TriangleArea(v2b_rel,v1b_rel,BC_l_rel)+Point::TriangleArea(v1b_rel,v1_rel,BC_l_rel));
//    v1_rel += P_x_d;
//    BC_l_rel += P_x_d/4;
//    double E_m1 = T_l*(Point::TriangleArea(v1_rel,v2_rel,BC_l_rel)+Point::TriangleArea(v2_rel,v2b_rel,BC_l_rel)+Point::TriangleArea(v2b_rel,v1b_rel,BC_l_rel)+Point::TriangleArea(v1b_rel,v1_rel,BC_l_rel));
//    v1_rel += P_x_d_2;
//    BC_l_rel += P_x_d_2/4;
//    F1a.x = (E_p1-E_m1)/d;
//    
//    v1_rel += P_y_d_2;
//    BC_l_rel += P_y_d_2/4;
//    E_p1 = T_l*(Point::TriangleArea(v1_rel,v2_rel,BC_l_rel)+Point::TriangleArea(v2_rel,v2b_rel,BC_l_rel)+Point::TriangleArea(v2b_rel,v1b_rel,BC_l_rel)+Point::TriangleArea(v1b_rel,v1_rel,BC_l_rel));
//    v1_rel += P_y_d;
//    BC_l_rel += P_y_d/4;
//    E_m1 = T_l*(Point::TriangleArea(v1_rel,v2_rel,BC_l_rel)+Point::TriangleArea(v2_rel,v2b_rel,BC_l_rel)+Point::TriangleArea(v2b_rel,v1b_rel,BC_l_rel)+Point::TriangleArea(v1b_rel,v1_rel,BC_l_rel));
//    v1_rel += P_y_d_2;
//    BC_l_rel += P_y_d_2/4;
//    F1a.y = (E_p1-E_m1)/d;
//
//    v1_rel += P_z_d_2;
//    BC_l_rel += P_z_d_2/4;
//    E_p1 = T_l*(Point::TriangleArea(v1_rel,v2_rel,BC_l_rel)+Point::TriangleArea(v2_rel,v2b_rel,BC_l_rel)+Point::TriangleArea(v2b_rel,v1b_rel,BC_l_rel)+Point::TriangleArea(v1b_rel,v1_rel,BC_l_rel));
//    v1_rel += P_z_d;
//    BC_l_rel += P_z_d/4;
//    E_m1 = T_l*(Point::TriangleArea(v1_rel,v2_rel,BC_l_rel)+Point::TriangleArea(v2_rel,v2b_rel,BC_l_rel)+Point::TriangleArea(v2b_rel,v1b_rel,BC_l_rel)+Point::TriangleArea(v1b_rel,v1_rel,BC_l_rel));
//    v1_rel += P_z_d_2;
//    BC_l_rel += P_z_d_2/4;
//    F1a.z = (E_p1-E_m1)/d;
//    
////    v1_rel += P_x_d_2;
////    BC_l_rel += P_x_d_2/4;
////    double E_p1 = T_l*Point::TriangleArea(v1_rel,v2_rel,BC_l_rel)+T_l*Point::TriangleArea(v1_rel,v1b_rel,BC_l_rel);
////    v1_rel += P_x_d;
////    BC_l_rel += P_x_d/4;
////    double E_m1 = T_l*Point::TriangleArea(v1_rel,v2_rel,BC_l_rel)+T_l*Point::TriangleArea(v1_rel,v1b_rel,BC_l_rel);;
////    v1_rel += P_x_d_2;
////    BC_l_rel += P_x_d_2/4;
////    F1a.x = (E_p1-E_m1)/d;
////    
////    v1_rel += P_y_d_2;
////    BC_l_rel += P_y_d_2/4;
////    E_p1 = T_l*Point::TriangleArea(v1_rel,v2_rel,BC_l_rel)+T_l*Point::TriangleArea(v1_rel,v1b_rel,BC_l_rel);;
////    v1_rel += P_y_d;
////    BC_l_rel += P_y_d/4;
////    E_m1 = T_l*Point::TriangleArea(v1_rel,v2_rel,BC_l_rel)+T_l*Point::TriangleArea(v1_rel,v1b_rel,BC_l_rel);;
////    v1_rel += P_y_d_2;
////    BC_l_rel += P_y_d_2/4;
////    F1a.y = (E_p1-E_m1)/d;
////    
////    v1_rel += P_z_d_2;
////    BC_l_rel += P_z_d_2/4;
////    E_p1 = T_l*Point::TriangleArea(v1_rel,v2_rel,BC_l_rel)+T_l*Point::TriangleArea(v1_rel,v1b_rel,BC_l_rel);;
////    v1_rel += P_z_d;
////    BC_l_rel += P_z_d/4;
////    E_m1 = T_l*Point::TriangleArea(v1_rel,v2_rel,BC_l_rel)+T_l*Point::TriangleArea(v1_rel,v1b_rel,BC_l_rel);;
////    v1_rel += P_z_d_2;
////    BC_l_rel += P_z_d_2/4;
////    F1a.z = (E_p1-E_m1)/d;
//
//    
//    return F1a;
//}
//    


double Edge::area()
{
    UpdateBarycenter();

    if(!T->isPeriodic || !CrossBoundary()) { // tissue is not periodic
        double area = 0;
        area += Point::TriangleArea(v1->coord,v2->coord,BC_l);
        area += Point::TriangleArea(v2->coord,v2->coord_b,BC_l);
        area += Point::TriangleArea(v2->coord_b,v1->coord_b,BC_l);
        area += Point::TriangleArea(v1->coord_b,v1->coord,BC_l);
        
        return area;
    }
    else { // edge crosses the boundary
        
        Quadrant SystemSize = T->SystemSize;
        
        // CALCULATE THE POSITIONS OF THE VERTICES WITH RESPECT TO VERTEX 1
        Point v1b_rel = shift( v1->coord_b, SystemSize*v1->q_b );
        Point v2_rel = shift( v2->coord, SystemSize*q_v2);
        Point v2b_rel = shift( v2->coord_b, SystemSize* (q_v2 + v2->q_b));
        Point BC_l_rel = shift( BC_l, SystemSize* q_BC_l);
        
        double area = 0;
        area += Point::TriangleArea(v1->coord,v2_rel,BC_l_rel);
        area += Point::TriangleArea(v2_rel,v2b_rel,BC_l_rel);
        area += Point::TriangleArea(v2b_rel,v1b_rel,BC_l_rel);
        area += Point::TriangleArea(v1b_rel,v1->coord,BC_l_rel);
        
        return area;
    }
}

void Edge::setMechanicalProperties()
{
    int type_b = std::max(c1->Type,c2->Type);
    int type_s = std::min(c1->Type,c2->Type);
    
    T_l = T->MechProp.Intercell_SurfaceTension[type_b][type_s];
    G_a = T->MechProp.Intercell_ApicalLineTension[type_b][type_s];
    G_b = T->MechProp.Intercell_BasalLineTension[type_b][type_s];
}