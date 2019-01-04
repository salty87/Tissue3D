/*
 *  vertex.cpp
 *  evolver1
 *  last modified: 07/03/2012 
 *
 */


#include "Vertex.h"

Vertex::Vertex(Tissue *T) : T(T)
{
    Number=T->maxVertexNumber()+1;

	coord=Point();
	coord_b=Point();

	// no fixation by default
	xafixed=false;
	yafixed=false;
	zafixed=false;
	
	xbfixed=false;
	ybfixed=false;
	zbfixed=false;
		  
    q_b = Quadrant();
}

// has been changed
Vertex::Vertex(Tissue *T, Point coord, Point coord_b) : coord(coord),coord_b(coord_b), T(T)
{
    Number=T->maxVertexNumber()+1;

	xafixed=false;
	yafixed=false;
	zafixed=false;
	
	xbfixed=false;
	ybfixed=false;
	zbfixed=false;
    
    q_b = Quadrant();
}

Vertex::Vertex(Tissue *T, Point coord, Point coord_b,std::list<Edge*> NeighbouringEdges_sorted) : coord(coord),coord_b(coord_b), T(T), NeighbouringEdges_sorted(NeighbouringEdges_sorted)
{
    Number=T->maxVertexNumber()+1;
	
	// no fixation by default
	xafixed=false;
	yafixed=false;
	zafixed=false;
	
	xbfixed=false;
	ybfixed=false;
	zbfixed=false;
	
    q_b = Quadrant();
}

// has been changed
void Vertex::Print()
{
	coord.Print();
	coord_b.Print();	
}

// has been adapted
bool Vertex::operator==(const Vertex& ve)
{
	// should involve some epsilon term, that stands for the accuracy!
	return((coord.x==ve.coord.x)&&(coord.y==ve.coord.y)&&(coord.z==ve.coord.z)&&(coord_b.x==ve.coord_b.x)&&(coord_b.y==ve.coord_b.y)&&(coord_b.z==ve.coord_b.z));
	
}

Point Vertex::lineForceOnApicalPoint()
{
	if(!CrossBoundary())
	{
		double dist = Point::Norm(coord_b-coord);
		if(dist==0)
		{
			return Point(0,0,0);
		}
		return(((coord_b-coord)/Point::Norm(coord_b-coord))*G_l);
	}
	else 
	{
		double dist = Point::Norm(coord_b-(shift(coord, T->SystemSize*q_b)));
		if(dist==0)
		{
			return Point(0,0,0);
		}
		return(((coord_b-coord)/Point::Norm(coord_b-coord))*G_l);
	}

}

// update the stress tensor induced by the line tension of the line inbetween the apical and the basal side
Point Vertex::UpdateLineStress()
{
	if(!CrossBoundary())
	{
		Stress=Point::lineStress(coord, coord_b)*G_l;
	}
	else 
	{
		Point vb = shift( coord_b, T->SystemSize*q_b);
		Stress=Point::lineStress(coord, vb)*G_l;
	}

	return Stress;
}


bool Vertex::CrossBoundary()
{
	return !(q_b.isZero());
}

/// iterates the vertex by step size mobility according to the forces acting on the points
double Vertex::movePoints(double mobility, double &Lx_new, double &Ly_new)
{
	if(!T->isPeriodic)
	{
		coord = coord + (dW_dva+ExternalForce_a)*mobility;
		coord_b = coord_b + (dW_dvb+ExternalForce_b)*mobility;
	}
	// if the tissue is periodic one has to check if the points crossed the boundary
	else 
	{
		// update the position inside the old system size
		coord = coord + (dW_dva+ExternalForce_a)*mobility;
		coord_b = coord_b + (dW_dvb+ExternalForce_b)*mobility;
        
		// check if the
		coord.x = coord.x/T->SystemSize.x;
		coord.y = coord.y/T->SystemSize.y;
		coord_b.x = coord_b.x/T->SystemSize.x;
		coord_b.y = coord_b.y/T->SystemSize.y;		
		
		double quadrantCount_a_x = floor(coord.x);
		double quadrantCount_b_x = floor(coord_b.x);
		double quadrantCount_a_y = floor(coord.y);
		double quadrantCount_b_y = floor(coord_b.y);
		
		coord.x = (coord.x - quadrantCount_a_x)*Lx_new;
		coord_b.x = (coord_b.x - quadrantCount_b_x)*Lx_new;
		coord.y = (coord.y - quadrantCount_a_y)*Ly_new;
		coord_b.y = (coord_b.y - quadrantCount_b_y)*Ly_new;
		
		
		
		// commit changes if the vertex actually has left one quadrant and entered another one
		if(quadrantCount_a_x!=0 || quadrantCount_a_y!=0 || quadrantCount_b_x!=0 ||quadrantCount_b_y!=0)
		{
			// update the basal quadrant of this vertex
			q_b.x += quadrantCount_b_x - quadrantCount_a_x;
			q_b.y += quadrantCount_b_y - quadrantCount_a_y;
			
			// iterate through all the edges and cells that are connected to this vertex and tell them that their quadrant has changed
			for(std::list<Edge*>::iterator it_e = NeighbouringEdges_sorted.begin();it_e!=NeighbouringEdges_sorted.end();it_e++)
			{
				// if this vertex is number one to the edge...
				if((*it_e)->v1 == this)
				{
					(*it_e)->q_v2.x-=quadrantCount_a_x;
					(*it_e)->q_v2.y-=quadrantCount_a_y;
				}
				else
				{
					(*it_e)->q_v2.x+=quadrantCount_a_x;
					(*it_e)->q_v2.y+=quadrantCount_a_y;
				}
				
				(*it_e)->crossBoundary = (*it_e)->q_v2.x!=0 || (*it_e)->q_v2.y!=0 || q_b.x!=0 || q_b.y!=0;
				
				// now, if the edge crosses the boundary, inform the two adjacent cells
				if((*it_e)->crossBoundary)
				{
					if((*it_e)->c1) (*it_e)->c1->crossBoundary = true;
					if((*it_e)->c2) (*it_e)->c2->crossBoundary = true;
				}
				// if it doesnt anymore check if the cells still cross the boundary
				else
				{
					if((*it_e)->c1) (*it_e)->c1->CrossBoundary();
					if((*it_e)->c2) (*it_e)->c2->CrossBoundary();

				}

					
				
			}


		}
		
	}

    return std::max(Point::Norm(dW_dva),Point::Norm(dW_dvb));
	
}

/// takes an edge, draws a plane with normal direction through the vertex, checks which edges lie on which side and then creates a new edge at (coord-direction/2, coord + direction/2)
int Vertex::UnpopEdge(Point direction, int* flag, Edge* e)
{
	// if there are less than 4 edges connected, the transition cannot take place
	if(NeighbouringEdges_sorted.size()<4)
	{
		*flag = 2;
		return -1;
	}
	
	// first: create the new vertex in the tissue
	Vertex vnew(T);T->ListVertex.push_back(vnew);
	Vertex* v = &(T->ListVertex.back());
	
	
	
	// second step: evaluate which edges of the tissue lie on which side of the plane, that goes through the apical point of the vertex and has the normal direction	
	std::list<Edge*> EdgesInDirection;
	for(std::list<Edge*>::iterator it_edge = NeighbouringEdges_sorted.begin();it_edge!=NeighbouringEdges_sorted.end();it_edge++)
	{
		int whichVertex = (*it_edge)->VertexIsPart(this);
		Point EdgeDirection;
		
		if(whichVertex == 1)
			EdgeDirection = (*it_edge)->getVertex2()->getCoord()-coord;
		else 
			EdgeDirection = (*it_edge)->getVertex1()->getCoord()-coord;

		
				
		if(EdgeDirection.dot(direction)<0) // if the angle between the two vectors is smaller than 0, the edge lies on the other side of the plane than the direction vector and hence it should point to the new vertex
		{
			//inDirection.push_front(&(*it_edge));
			if(whichVertex==1)
				(*it_edge)->setVertex1(v);
			else 
				(*it_edge)->setVertex2(v);
		}
		else 
		{
			EdgesInDirection.push_back((*it_edge));
		}

	}
	
	// third step: set the new coordinates of the two vertices of the edge
	v->setCoord(coord - direction*0.5);
	v->setBasalCoord(coord_b - direction*0.5);
	coord = coord + direction*0.5;
	coord_b = coord_b + direction*0.5;
	
	// fourth step: define a new edge
	Edge enew(T); T->ListEdge.push_back(enew);
	e = &(*(T->ListEdge.end()));
	e->setVertex1(this);
	e->setVertex2(v);
	
	// what are the cells e (ca, cb) is connected with? answer: the cells that the two edges towards direction do not have in common
	Cell *ca(NULL), *cb(NULL);
	if(EdgesInDirection.size()==2)
	{
		Edge* e1 = (*EdgesInDirection.begin());
		Edge* e2 = (*(++EdgesInDirection.begin()));

		if(e1->c1 == e2->c1)
		{
			ca = e1->c2; cb = e2->c2;
		}
		else if(e1->c1 == e2->c2)
		{
			ca = e1->c2; cb = e2->c1;
		}
		else if(e1->c2 == e2->c1)
		{
			ca = e1->c1; cb = e2->c2;
		}
		else if(e1->c2 == e2->c2)
		{
			ca = e1->c1; cb = e2->c1;
		}
		else
		{
			// just for testing!! :
			throw;
		}
	}
	
	

	// need to set the directionality of the edge with respect to the two cells	
	// now iterate through cell ca and find the existing edge that is connected to this vertex
	std::list<Edge*>::iterator it_edge = ca->ListEdges.begin();
	for(;it_edge!=ca->ListEdges.end();it_edge++)
	{
		int isPart = (*it_edge)->VertexIsPart(this);
		
		if(isPart==1)
		{
			if((*it_edge)->getCell1() == ca)
			{
				e->setCell1(cb); e->setCell2(ca);
			}
			else
			{
				e->setCell1(ca); e->setCell2(cb);
			}
			
			break;
			
		}
		else if(isPart==2)
		{
			if((*it_edge)->getCell1() == ca)
			{
				e->setCell1(ca); e->setCell2(cb);
			}
			else
			{
				e->setCell1(cb); e->setCell2(ca);
			}
			
			break;
		}
	}
	
	return 1;
}

// c1 and c2 are the cells that are going to be connected afterwards, popVector is the "force" that tries to separate the vertex and the two cell groups give the separation of the remaining cells in two groups 
//int Vertex::ForcesFromPop(Cell* c1, Cell* c2, Point* popVector, std::list<Cell*> *CellGroup1,std::list<Cell*> *CellGroup2, double SurfaceTensionNewEdge, double LineTensionNewEdge)
//{
//	CellGroup1->clear();
//	CellGroup2->clear();
//	
//	// vertex should at least be fourfold vertex
//	if(NeighbouringEdges_sorted.size()<=3)
//		return -1;
//	
//	// check if the two cells are actually adjacent to the vertex
//	//if(std::find(ListCells.begin(),ListCells.end(),c1) == ListCells.end() || std::find(ListCells.begin(),ListCells.end(),c2) == ListCells.end())
//	if(c1==c2)
//		return -2;
//	
//	// check if cells c1 and c2 are neighbours 
//	Edge* e = new Edge();
//	if(c1->hasCommonEdgeWith(c2, &e))
//		return -3;
//	
//	if(!T->CheckForPop) // the forces should have been updated (during the calculation of the current forces) before
//		return -5;
//	
//	// now collect the two groups of cells which will finally be on different sides of the edge
//	Point ForceFromCellGroup1=Point(), ForceFromCellGroup2=Point();
//	
//	// the edges that separate the two groups from cell 1 (GroupX_edge1) and from cell 2 (GroupX_edge2)
//	Edge *Group1_edge1, *Group1_edge2, *Group2_edge1, *Group2_edge2;
//	
//	// find neighbour to c1
//	Cell* lastFoundCell = c1;
////	std::list<Cell*>::iterator itc = ListCells.begin();
////	bool reachedC2 = false;
////	for(;itc!=ListCells.end();itc++)
////	{
////		if(*itc!=lastFoundCell && lastFoundCell->hasCommonEdgeWith(*itc,&e))
////		{
////			if(lastFoundCell!=c1)
////			{
////				int direction;
////				if(e->VertexIsPart(this)==1)
////					direction = 1;
////				else
////					direction = -1;
////				
////				// substract every inbetween edge once, since otherwise they are added twice later on
////				ForceFromCellGroup1 = ForceFromCellGroup1 - (e->ApicalLineForce()+e->BasalLineForce())*direction;
////			}
////			else
////			{
////				// safe the edge that belongs to cell 1
////				Group1_edge1 = e;
////			}
////				
////			if(*itc == c2)
////			{
////				reachedC2 = true;
////				Group1_edge2 = e;
////				break;
////			}
////			
////			CellGroup1->push_back((*itc));
////			c1 = lastFoundCell;
////			lastFoundCell = (*itc);
////			itc = ListCells.begin();
////
////		}
////	}
//	
//	
//	if(!reachedC2) // something went wrong
//		return -4;
//	
//	
//	// now that CellGroup1 has been created we can get CellGroup2 by taking the remaining cells
//	*CellGroup2 = ListCells;
//	CellGroup2->remove(c1);
//	CellGroup2->remove(c2);
//	for(itc=CellGroup1->begin();itc!=CellGroup1->end();itc++)
//	{
//		CellGroup2->remove(*itc);
//	}
//	
//	
//	// now also substract the edges inbetween group 2 from the forces once
//	for(itc=CellGroup2->begin();itc!=CellGroup2->end();itc++)
//	{
//		for(std::list<Cell*>::iterator itc2 = itc; itc2!=CellGroup2->end();itc2++)
//		{
//			(*itc)->hasCommonEdgeWith((*itc2),&e);
//		
//			int direction;
//			if(e->VertexIsPart(this)==1)
//				direction = 1;
//			else
//				direction = -1;
//		
//			// substract every inbetween edge once, since otherwise they are added twice later on
//			ForceFromCellGroup2 = ForceFromCellGroup2 - (e->ApicalLineForce()+e->BasalLineForce())*direction;
//		}
//	}
//	
//	// the edge of c1 that is connected with *this and that is not Group1_edge1, is Group2_edge1
//	for(std::list<Edge*>::iterator ite = c1->ListEdges.begin(); ite != c1->ListEdges.end();ite++)
//	{
//		if((*ite)->VertexIsPart(this)!=0 && (*ite) != Group1_edge1)
//		{
//			Group2_edge1 = (*ite);
//			break;
//		}
//	}
//
//	// the edge of c1 that is connected with *this and that is not Group1_edge1, is Group2_edge1
//	for(std::list<Edge*>::iterator ite = c2->ListEdges.begin(); ite != c2->ListEdges.end();ite++)
//	{
//		if((*ite)->VertexIsPart(this)!=0 && (*ite) != Group1_edge2)
//		{
//			Group2_edge2 = *ite;
//			break;
//		}
//	}
//	
//	
//	// now that the two groups have been separated we can calculate the forces/tensions exerted by the cells on the vertex
//	std::map<Cell*,Point>::iterator it_map;
//	for(itc = CellGroup1->begin();itc!=CellGroup1->end();itc++)
//	{
//		// find the cell in group 1 in the map that carries the information about all the forces
//		it_map = ForcesFromCells.find(*itc);
//		ForceFromCellGroup1 = ForceFromCellGroup1 + (*it_map).second;
//	}
//	
//	for(itc = CellGroup2->begin();itc!=CellGroup2->end();itc++)
//	{
//		// find the cell in group 2 in the map ForcesFromCells that carries the information about all the forces
//		it_map = ForcesFromCells.find(*itc);
//		ForceFromCellGroup2 = ForceFromCellGroup2 + (*it_map).second;
//	}
//	
//	// at that point, the forces exerted by the cells that are not c1 and c2 on the vertices are calculated, that is for the cells that do not change their topology
//	// now clearly the forces exerted by c1 and c2 have to be included as well:	
//	Point ForcesFromPressure = Group1_edge1->Force_PressureT2Transition(c1, this) + Group1_edge2->Force_PressureT2Transition(c2, this);
//	ForceFromCellGroup1 = ForceFromCellGroup1 + ForcesFromPressure;
//	
//	ForcesFromPressure = Group2_edge1->Force_PressureT2Transition(c1, this) + Group2_edge2->Force_PressureT2Transition(c2, this);
//	ForceFromCellGroup2 = ForceFromCellGroup2 + ForcesFromPressure;
//	
//	// now project the contributions coming from the two groups to the perpendicular plane
//	Point normal = (coord - coord_b)/Point::Norm(coord-coord_b);
//	ForceFromCellGroup1 = ForceFromCellGroup1 - normal*(ForceFromCellGroup1.dot(normal));
//	ForceFromCellGroup2 = ForceFromCellGroup2 - normal*(ForceFromCellGroup2.dot(normal));
//	
//	// finally subtract the forces coming from the surface increase of cell 1 and cell 2 
//	Point ForcesFromSurface = Group1_edge1->Force_SurfaceTensionT2Transition(c1, this, ForceFromCellGroup1) + Group1_edge2->Force_SurfaceTensionT2Transition(c2, this, ForceFromCellGroup1);
//	ForceFromCellGroup1 = ForceFromCellGroup1 + ForcesFromSurface;
//	ForcesFromSurface = Group2_edge1->Force_SurfaceTensionT2Transition(c1, this, ForceFromCellGroup2) + Group2_edge2->Force_SurfaceTensionT2Transition(c2, this, ForceFromCellGroup2);
//	ForceFromCellGroup2 = ForceFromCellGroup2 + ForcesFromSurface;	
//	
//	// the vector that characterizes the preferred shift in the direction of the group 1 (and opposite to the one of group 2) is now given by:
//	(*popVector) = (ForceFromCellGroup1-ForceFromCellGroup2)/2;
//	
//	// is this force bigger than the force coming from the surface and line tension of the new edge?
//	double test = Point::Norm(*popVector) - SurfaceTensionNewEdge * Point::Norm(coord - coord_b) - LineTensionNewEdge;
//	
//	if(test>0) 
//	{
//		(*popVector) = (*popVector) * (test/Point::Norm(*popVector));
//		return 1;
//	}
//	else
//	{
//		// pop cannot decrease the energy
//		return 0;
//	}
//	
//	
//	
//}

//int Vertex::T1_pop(double edgeLength, Edge* e)
//{
//
//    // ListEdge has always be maintained sorted!
//    
//    // no pop possible if the vertex has less than 4 neighbouring edges == number of cells!
//    if(ListEdge.size()<4)  return 0;
//
//        
//
//    
//}


bool Vertex::checkAndPop_energy(double edgeLength)
{
    updateNeighbouringEdges();
    std::list<Edge*> NeighbouringEdges_sorted2 = NeighbouringEdges_sorted;

    if(NeighbouringEdges_sorted.size()<4) // no pop possible if the vertex is less than 4 fold
        return 0;
    
    double minimalEnergy = T->W();
    unsigned minimalConfiguration = 0;
    unsigned counter = 1;
    
    //
    for(std::list<Edge*>::iterator it_e1 = NeighbouringEdges_sorted2.begin();it_e1 != NeighbouringEdges_sorted2.end(); it_e1++ )  { // if any of the neighbouring cells crosses the boundary do not pop
        if((*it_e1)->c1->CrossBoundary()) return 0;
    }
    
    std::cout << "Initial Energy before pop = " << minimalEnergy << std::endl;
    
    double testLength = 0.1;
    
    // iterate through the edges
    for(std::list<Edge*>::iterator it_e1 = NeighbouringEdges_sorted2.begin();; it_e1++ )
    {
        // go to the next element in the list and check if this is the last one
        std::list<Edge*>::iterator it_e2 = (++it_e1);
        --it_e1;
        if(it_e2==NeighbouringEdges_sorted2.end()) break;
        
        Cell* c1;
        if(this == (*it_e1)->v1)
        {
            c1 = (*it_e1)->c1;
        }
        else
        {
            c1 = (*it_e1)->c2;
        }
        
        
        // iterate with the 2nd iterator til the end
        ++it_e2;
        for(;it_e2!=NeighbouringEdges_sorted2.end();it_e2++)
        {
            
            // first and last edge together cannot pop since they are neighbours (only possible neighbour condition that has to be captured!)
            //if(it_e1== NeighbouringEdges_sorted2.begin() && (++it_e2)==NeighbouringEdges_sorted2.end()) {it_e2--;continue;}
            if((*it_e1)->c1 == (*it_e2)->c1 || (*it_e1)->c2 == (*it_e2)->c2 || (*it_e1)->c1 == (*it_e2)->c2 || (*it_e1)->c2 == (*it_e2)->c1) {continue;}
            
            Cell* c2;
            if(this == (*it_e2)->v1)
            {
                c2 = (*it_e2)->c1;
            }
            else
            {
                c2 = (*it_e2)->c2;
            }
                        
            // create new elements (v and e)
            T->ListVertex.push_back(*(new Vertex(T)));
            Vertex *v = &T->ListVertex.back();
            v->T = T;
            v->G_l = T->MechProp.Vertex_LineTension;
            Edge *e = new Edge(T);
            T->ListEdge.push_back(*e);
            e = &T->ListEdge.back();
            e->T = T;

            e->c1 = c1;
            e->c2 = c2;
            
            int type1, type2;
            if(c1->Type>=c2->Type) {type1=c1->Type;type2=c2->Type;}
            else {type1=c2->Type;type2=c1->Type;}
            e->T_l = T->MechProp.Intercell_SurfaceTension[type1][type2]; // mechanical properties of the edge are set by the types of the constituting cells
            e->G_a = T->MechProp.Intercell_ApicalLineTension[type1][type2];
            e->G_b = T->MechProp.Intercell_BasalLineTension[type1][type2];
            
            
            // connect the edges in (it_e1+1,...,it_e2) to the new vertex
            std::list<Edge*>::iterator it_e = it_e1;
            do
            {
                it_e++;
                if((*it_e)->v1==this) (*it_e)->v1=v;
                else (*it_e)->v2=v;
            }
            while(it_e!=it_e2);
            
            // assign the new edge to the cells c1 and c2 at the right position
            for(it_e = c1->ListEdges.begin(); it_e!=c1->ListEdges.end();it_e++)
            {
                if(*it_e==*it_e1)
                {
                    it_e++;
                    c1->ListEdges.insert(it_e,e);
                    break;
                }
            }

            // assign the new edge to the cells c1 and c2 at the right position
            for(it_e = c2->ListEdges.begin(); it_e != c2->ListEdges.end();it_e++)
            {
                if(*it_e==*it_e2)
                {
                    c2->ListEdges.insert(it_e,e);
                    break;
                }
            }
                        
            // calculate the apical/basal displacement: orthogonal to the plane through the basal/apical points and the apical/basal barycenters with length edgeLength
            Point displacement_a =  (c1->BC_a - coord_b)*(c2->BC_a - coord_b);  displacement_a = displacement_a*(testLength/(2*Point::Norm(displacement_a)));
            Point displacement_b =  (c1->BC_b - coord)*(c2->BC_b - coord);   displacement_b = displacement_b*(testLength/(2*Point::Norm(displacement_b)));
            if(T->STIFF_BASAL_MEMBRANE) displacement_b.z = 0;
            
            v->coord = coord-displacement_a; coord += displacement_a;
            v->coord_b = coord_b+displacement_b; coord_b -= displacement_b;
            
            T->W();
            
            if(T->Energy<minimalEnergy)
            {
                minimalEnergy=T->Energy;
                minimalConfiguration = counter;
            }
                
            std::cout << "New energy from counter " << counter << " equals: " << T->Energy << std::endl;
            
            // and reestablish the state before the pop to test for other pops
            coord = (coord+v->coord)/2;
            coord_b = (coord_b+v->coord_b)/2;
            
            // reconnect the edges in (it_e1+1,...,it_e2) to this
            it_e = it_e1;
            do
            {
                it_e++;
                if((*it_e)->v1==v) (*it_e)->v1=this;
                else (*it_e)->v2=this;
            }
            while(it_e!=it_e2);
            
            // remove the new edge from cells c1 and c2
            c1->ListEdges.remove(e);
            c2->ListEdges.remove(e);
            
            // remove the new cell and the new vertex from T
            T->ListEdge.pop_back();
            T->ListVertex.pop_back();
            
            counter++;
        }
    }
    
    unsigned counter2=1;
    // now go to the found minimal energy configuration and apply it
    if(minimalConfiguration!=0)
    {
        for(std::list<Edge*>::iterator it_e1 = NeighbouringEdges_sorted2.begin();; it_e1++ )
        {
            // go to the next element in the list and check if this is the last one
            std::list<Edge*>::iterator it_e2 = it_e1;
                       
            // iterate with the 2nd iterator til the end
            ++(++it_e2);
            for(;it_e2!=NeighbouringEdges_sorted2.end();it_e2++)
            {
                // first and last edge together cannot pop since they are neighbours (only possible neighbour condition that has to be captured!)
                if((*it_e1)->c1 == (*it_e2)->c1 || (*it_e1)->c2 == (*it_e2)->c2 || (*it_e1)->c1 == (*it_e2)->c2 || (*it_e1)->c2 == (*it_e2)->c1) {continue;}
                
                if(counter2==minimalConfiguration)
                {
                    Cell* c1;
                    if(this == (*it_e1)->v1)
                    {
                        c1 = (*it_e1)->c1;
                    }
                    else
                    {
                        c1 = (*it_e1)->c2;
                    }
                    
                    Cell* c2;
                    if(this == (*it_e2)->v1)
                    {
                        c2 = (*it_e2)->c1;
                    }
                    else
                    {
                        c2 = (*it_e2)->c2;
                    }
                    
                    // create new elements (v and e)
                    T->ListVertex.push_back(*(new Vertex(T)));
                    Vertex *v = &T->ListVertex.back();
                    v->T = T;
                    v->G_l = T->MechProp.Vertex_LineTension;

                    T->ListEdge.push_back(*(new Edge(T,v,this)));
                    Edge* e = &T->ListEdge.back();
                    e->T = T;
                    
                    e->c1 = c1;
                    e->c2 = c2;
                    
                    int type1, type2;
                    if(c1->Type<=c2->Type) {type1=c1->Type;type2=c2->Type;}
                    else {type1=c2->Type;type2=c1->Type;}
                    e->T_l = T->MechProp.Intercell_SurfaceTension[type1][type2]; // mechanical properties of the edge are set by the types of the constituting cells
                    e->G_a = T->MechProp.Intercell_ApicalLineTension[type1][type2];
                    e->G_b = T->MechProp.Intercell_BasalLineTension[type1][type2];
                    
                    
                    // connect the edges in (it_e1+1,...,it_e2) to the new vertex
                    std::list<Edge*>::iterator it_e = it_e1;
                    do
                    {
                        it_e++;
                        if((*it_e)->v1==this) (*it_e)->v1=v;
                        else (*it_e)->v2=v;
                    }
                    while(it_e!=it_e2);
                    
                    // assign the new edge to the cells c1 and c2 at the right position
                    for(it_e = c1->ListEdges.begin(); it_e!=c1->ListEdges.end();it_e++)
                    {
                        if(*it_e==*it_e1)
                        {
                            it_e++;
                            c1->ListEdges.insert(it_e,e);
                            break;
                        }
                    }
                    
                    // assign the new edge to the cells c1 and c2 at the right position
                    for(it_e = c2->ListEdges.begin(); it_e != c2->ListEdges.end();it_e++)
                    {
                        if(*it_e==*it_e2)
                        {
                            c2->ListEdges.insert(it_e,e);
                            break;
                        }
                    }
                    
                    // calculate the apical/basal displacement: orthogonal to the plane through the basal/apical points and the apical/basal barycenters with length edgeLength
                    Point displacement_a =  (c1->BC_a - coord_b)*(c2->BC_a - coord_b);  displacement_a = displacement_a*(edgeLength/(2*Point::Norm(displacement_a)));
                    Point displacement_b =  (c1->BC_b - coord)*(c2->BC_b - coord);   displacement_b = displacement_b*(edgeLength/(2*Point::Norm(displacement_b)));
                    if(T->STIFF_BASAL_MEMBRANE) displacement_b.z = 0;
                    v->coord = coord-displacement_a; coord += displacement_a;
                    v->coord_b = coord_b+displacement_b; coord_b -= displacement_b;

                    T->update();
                    
                    std::cout << "New energy = " << T->Energy << std::endl << "----------------------------\n";
                    
                    return 1;
                }
                
                counter2++;
            }
        }
        
    }
    else
    {
        return 0; // no pop has taken place
    }
}


// routine that checks for pop by calculating the forces that drive the vertices apart after the pop
//bool Vertex::checkAndPop_forces(double edgeLength)
//{
//    updateNeighbouringEdges();
//    std::list<Edge*> NeighbouringEdges_sorted2 = NeighbouringEdges_sorted;
//    
//    if(NeighbouringEdges_sorted.size()<4) // no pop possible if the vertex is less than 4 fold
//        return 0;
//    
//    double maximalSeparationForce = 0;
//    unsigned maximalSeparationForceConfiguration = 0;
//    unsigned counter = 1;
//    
//    for(std::list<Edge*>::iterator it_e1 = NeighbouringEdges_sorted2.begin();it_e1 != NeighbouringEdges_sorted2.end(); it_e1++ )  { // if any of the neighbouring cells crosses the boundary do not pop
//        if((*it_e1)->c1->CrossBoundary()) return 0;
//    }
//            
//    // iterate through the edges
//    for(std::list<Edge*>::iterator it_e1 = NeighbouringEdges_sorted2.begin();; it_e1++ )
//    {
//        // go to the next element in the list and check if this is the last one
//        std::list<Edge*>::iterator it_e2 = (++it_e1);
//        --it_e1;
//        if(it_e2==NeighbouringEdges_sorted2.end()) break;
//        
//        Cell* c1;
//        if(this == (*it_e1)->v1)
//        {
//            c1 = (*it_e1)->c1;
//        }
//        else
//        {
//            c1 = (*it_e1)->c2;
//        }
//        
//        
//        // iterate with the 2nd iterator til the end
//        ++it_e2;
//        for(;it_e2!=NeighbouringEdges_sorted2.end();it_e2++)
//        {
//            // first and last edge together cannot pop since they are neighbours (only possible neighbour condition that has to be captured!)
//            //if(it_e1== NeighbouringEdges_sorted2.begin() && (++it_e2)==NeighbouringEdges_sorted2.end()) {it_e2--;continue;}
//            if((*it_e1)->c1 == (*it_e2)->c1 || (*it_e1)->c2 == (*it_e2)->c2 || (*it_e1)->c1 == (*it_e2)->c2 || (*it_e1)->c2 == (*it_e2)->c1) {continue;}
//            
//            Cell* c2;
//            if(this == (*it_e2)->v1)
//            {
//                c2 = (*it_e2)->c1;
//            }
//            else
//            {
//                c2 = (*it_e2)->c2;
//            }
//            
//            // create new elements (v and e)
//            T->ListVertex.push_back(*(new Vertex()));
//            Vertex *v = &T->ListVertex.back();
//            v->T = T;
//            v->G_l = T->MechProp.Vertex_LineTension;
//            Edge *e = new Edge(v,this);
//            T->ListEdge.push_back(*e);
//            e = &T->ListEdge.back();
//            e->T = T;
//            
//            e->c1 = c1;
//            e->c2 = c2;
//            
//            int type1, type2;
//            if(c1->Type>=c2->Type) {type1=c1->Type;type2=c2->Type;}
//            else {type1=c2->Type;type2=c1->Type;}
//            e->T_l = T->MechProp.Intercell_SurfaceTension[type1][type2]; // mechanical properties of the edge are set by the types of the constituting cells
//            e->G_a = T->MechProp.Intercell_ApicalLineTension[type1][type2];
//            e->G_b = T->MechProp.Intercell_BasalLineTension[type1][type2];
//            
//            
//            // connect the edges in (it_e1+1,...,it_e2) to the new vertex
//            std::list<Edge*>::iterator it_e = it_e1;
//            do {
//                it_e++;
//                if((*it_e)->v1==this) (*it_e)->v1=v;
//                else (*it_e)->v2=v;
//            } while(it_e!=it_e2);
//            
//            for(it_e = c1->ListEdges.begin(); it_e!=c1->ListEdges.end();it_e++) {  // assign the new edge to the cells c1 and c2 at the right position
//
//                if(*it_e==*it_e1)
//                {
//                    it_e++;
//                    c1->ListEdges.insert(it_e,e);
//                    break;
//                }
//            }
//            
//
//            for(it_e = c2->ListEdges.begin(); it_e != c2->ListEdges.end();it_e++)  {  // assign the new edge to the cells c1 and c2 at the right position
//                if(*it_e==*it_e2)
//                {
//                    c2->ListEdges.insert(it_e,e);
//                    break;
//                }
//            }
//            
//            // both vertices are at the same position
//            v->coord = coord;
//            v->coord_b = coord_b;
//            T->update();
//            
//            // updating the tissues gives forces acting on the vertices after the pop. now clearly the forces from the edges area and the apical and basal lines is not well defined since there is no area or length of the edges. this problem is solved, by ignoring forces from edge with zero area
//            
//            // extract the apical and basal forces driving the vertices apart
//            Point ApicalForce = v->dW_dva - dW_dva;
//            Point BasalForce = v->dW_dvb - dW_dvb;
//            
//            // can those forces overcome the surface and line tensions
//            double norm_apicalForce = Point::Norm(ApicalForce);// - (e->G_a + Point::Norm(coord-coord_b)*e->T_l);
//            double norm_basalForce = Point::Norm(BasalForce);// - (e->G_b + Point::Norm(coord-coord_b)*e->T_l);
//            double separationForce = (norm_apicalForce+norm_basalForce)/2.;
//            
//            std::cout << "separationForce for counter " << counter << " is " << separationForce << std::endl;
//            
//            if(separationForce>maximalSeparationForce)  {
//                maximalSeparationForce = separationForce;
//                maximalSeparationForceConfiguration = counter;
//            }
//            
//            // now rebuild the tissue as it was before the pop:
//            
//            // reconnect the edges in (it_e1+1,...,it_e2) to this
//            it_e = it_e1;
//            do
//            {
//                it_e++;
//                if((*it_e)->v1==v) (*it_e)->v1=this;
//                else (*it_e)->v2=this;
//            }
//            while(it_e!=it_e2);
//            
//            // remove the new edge from cells c1 and c2
//            c1->ListEdges.remove(e);
//            c2->ListEdges.remove(e);
//            
//            // remove the new cell and the new vertex from T
//            T->ListEdge.pop_back();
//            T->ListVertex.pop_back();
//            
//            counter++;
//        }
//    }
//    
//    unsigned counter2=1;
//    // now go to the found minimal energy configuration and apply it
//    if(maximalSeparationForceConfiguration!=0)
//    {
//        std::list<Edge*> NeighbouringEdges_sorted2 = NeighbouringEdges_sorted;
//
//        for(std::list<Edge*>::iterator it_e1 = NeighbouringEdges_sorted2.begin();; it_e1++ )
//        {
//            // go to the next element in the list and check if this is the last one
//            std::list<Edge*>::iterator it_e2 = (++it_e1);
//            --it_e1;
//            
//            // iterate with the 2nd iterator til the end
//            for(;it_e2!=NeighbouringEdges_sorted2.end();it_e2++)
//            {
//                // first and last edge together cannot pop since they are neighbours (only possible neighbour condition that has to be captured!)
//                if((*it_e1)->c1 == (*it_e2)->c1 || (*it_e1)->c2 == (*it_e2)->c2 || (*it_e1)->c1 == (*it_e2)->c2 || (*it_e1)->c2 == (*it_e2)->c1) {continue;}
//                
//                if(counter2==maximalSeparationForceConfiguration)
//                {
//                    Cell* c1;
//                    if(this == (*it_e1)->v1)
//                    {
//                        c1 = (*it_e1)->c1;
//                    }
//                    else
//                    {
//                        c1 = (*it_e1)->c2;
//                    }
//                    
//                    Cell* c2;
//                    if(this == (*it_e2)->v1)
//                    {
//                        c2 = (*it_e2)->c1;
//                    }
//                    else
//                    {
//                        c2 = (*it_e2)->c2;
//                    }
//                    
//                    // create new elements (v and e)
//                    T->ListVertex.push_back(*(new Vertex()));
//                    Vertex *v = &T->ListVertex.back();
//                    v->T = T;
//                    v->G_l = T->MechProp.Vertex_LineTension;
//                    
//                    T->ListEdge.push_back(*(new Edge(v,this)));
//                    Edge* e = &T->ListEdge.back();
//                    e->T = T;
//                    
//                    e->c1 = c1;
//                    e->c2 = c2;
//                    
//                    int type1, type2;
//                    if(c1->Type<=c2->Type) {type1=c1->Type;type2=c2->Type;}
//                    else {type1=c2->Type;type2=c1->Type;}
//                    e->T_l = T->MechProp.Intercell_SurfaceTension[type1][type2]; // mechanical properties of the edge are set by the types of the constituting cells
//                    e->G_a = T->MechProp.Intercell_ApicalLineTension[type1][type2];
//                    e->G_b = T->MechProp.Intercell_BasalLineTension[type1][type2];
//                    
//                    
//                    // connect the edges in (it_e1+1,...,it_e2) to the new vertex
//                    std::list<Edge*>::iterator it_e = it_e1;
//                    do {
//                        it_e++;
//                        if((*it_e)->v1==this) (*it_e)->v1=v;
//                        else (*it_e)->v2=v;
//                    } while(it_e!=it_e2);
//                    
//                    // assign the new edge to the cells c1 and c2 at the right position
//                    for(it_e = c1->ListEdges.begin(); it_e!=c1->ListEdges.end();it_e++) {
//                        if(*it_e==*it_e1)
//                        {
//                            it_e++;
//                            c1->ListEdges.insert(it_e,e);
//                            break;
//                        }
//                    }
//                    
//                    // assign the new edge to the cells c1 and c2 at the right position
//                    for(it_e = c2->ListEdges.begin(); it_e != c2->ListEdges.end();it_e++) {
//                        if(*it_e==*it_e2)
//                        {
//                            c2->ListEdges.insert(it_e,e);
//                            break;
//                        }
//                    }
//                    
//                    // first update to calculate the forces ignoring the contributions from surface and line tension of the newly created edge
//                    T->update();
//
//                    Point displacement_a = dW_dva/(Point::Norm(dW_dva)/edgeLength);//(((Point::Norm(dW_dva) - (e->G_a + Point::Norm(coord-coord_b)*e->T_l))/Point::Norm(dW_dva)));
//                    Point displacement_b = dW_dvb/(Point::Norm(dW_dvb)/edgeLength);//(((Point::Norm(dW_dvb) - (e->G_b + Point::Norm(coord-coord_b)*e->T_l))/Point::Norm(dW_dvb)));
//
//                    v->coord = coord+displacement_a; coord -= displacement_a;
//                    v->coord_b = coord_b+displacement_b; coord_b -= displacement_b;
//
//                    //T->update();
//                    
//                    // apical and basal forces after pop
//                    dW_dva.Print();
//                    dW_dvb.Print();
//                    
//                    std::cout << "Separation force after Pop= " << Point::Norm(dW_dva-v->dW_dva) << std::endl;
//                    
//                    // calculate displacements
////
////                    std::cout << "Apical and Basal Displacement for Pop: ";
////                    displacement_a.Print();
////                    displacement_b.Print();
////                    
////                    if(T->STIFF_BASAL_MEMBRANE) displacement_b.z = 0;
////
////                    
////                    T->update();
////                    
////                    std::cout << "New energy = " << T->Energy << std::endl << "----------------------------\n";
//                    
//                    return 1;
//                }
//                
//                counter2++;
//            }
//        }
//        
//    }
//    else
//    {
//        return 0; // no pop has taken place
//    }
//}

bool Vertex::checkAndPop_forces_original(double edgeLength)
{
    updateNeighbouringEdges();
    std::list<Edge*> NeighbouringEdges_sorted2 = NeighbouringEdges_sorted;
    
    if(NeighbouringEdges_sorted.size()<4) // no pop possible if the vertex is less than 4 fold
        return 0;
    
    double maximalSeparationForce = 0;
    unsigned maximalSeparationForceConfiguration = 0;
    unsigned counter = 1;
    
    for(std::list<Edge*>::iterator it_e1 = NeighbouringEdges_sorted2.begin();it_e1 != NeighbouringEdges_sorted2.end(); it_e1++ )  { // if any of the neighbouring cells crosses the boundary do not pop
        if((*it_e1)->c1->CrossBoundary()) return 0;
    }
    
    // iterate through the edges
    for(std::list<Edge*>::iterator it_e1 = NeighbouringEdges_sorted2.begin();; it_e1++ )
    {
        // go to the next element in the list and check if this is the last one
        std::list<Edge*>::iterator it_e2 = (++it_e1);
        --it_e1;
        if(it_e2==NeighbouringEdges_sorted2.end()) break;
        
        Cell* c1;
        if(this == (*it_e1)->v1)
        {
            c1 = (*it_e1)->c1;
        }
        else
        {
            c1 = (*it_e1)->c2;
        }
        
        
        // iterate with the 2nd iterator til the end
        ++it_e2;
        for(;it_e2!=NeighbouringEdges_sorted2.end();it_e2++)
        {
            // first and last edge together cannot pop since they are neighbours (only possible neighbour condition that has to be captured!)
            //if(it_e1== NeighbouringEdges_sorted2.begin() && (++it_e2)==NeighbouringEdges_sorted2.end()) {it_e2--;continue;}
            if((*it_e1)->c1 == (*it_e2)->c1 || (*it_e1)->c2 == (*it_e2)->c2 || (*it_e1)->c1 == (*it_e2)->c2 || (*it_e1)->c2 == (*it_e2)->c1) {continue;}
            
            Cell* c2;
            
            if(this == (*it_e2)->v1)
            {
                c2 = (*it_e2)->c1;
            }
            else
            {
                c2 = (*it_e2)->c2;
            }
            
            // create new elements (v and e)
            T->ListVertex.push_back(*(new Vertex(T)));
            Vertex *v = &T->ListVertex.back();
            v->T = T;
            v->G_l = T->MechProp.Vertex_LineTension;
            Edge *e = new Edge(T,v,this);
            T->ListEdge.push_back(*e);
            e = &T->ListEdge.back();
            e->T = T;
            
            e->c1 = c1;
            e->c2 = c2;
            
            int type1, type2;
            if(c1->Type>=c2->Type) {type1=c1->Type;type2=c2->Type;}
            else {type1=c2->Type;type2=c1->Type;}
            e->T_l = T->MechProp.Intercell_SurfaceTension[type1][type2]; // mechanical properties of the edge are set by the types of the constituting cells
            e->G_a = T->MechProp.Intercell_ApicalLineTension[type1][type2];
            e->G_b = T->MechProp.Intercell_BasalLineTension[type1][type2];
            
            
            // connect the edges in (it_e1+1,...,it_e2) to the new vertex
            std::list<Edge*>::iterator it_e = it_e1;
            do {
                it_e++;
                if((*it_e)->v1==this) (*it_e)->v1=v;
                else (*it_e)->v2=v;
            } while(it_e!=it_e2);
            
            for(it_e = c1->ListEdges.begin(); it_e!=c1->ListEdges.end();it_e++) {  // assign the new edge to the cells c1 and c2 at the right position
                
                if(*it_e==*it_e1)
                {
                    it_e++;
                    c1->ListEdges.insert(it_e,e);
                    break;
                }
            }
            
            
            for(it_e = c2->ListEdges.begin(); it_e != c2->ListEdges.end();it_e++)  {  // assign the new edge to the cells c1 and c2 at the right position
                if(*it_e==*it_e2)
                {
                    c2->ListEdges.insert(it_e,e);
                    break;
                }
            }
            
            // both vertices are at the same position
            v->coord = coord;
            v->coord_b = coord_b;
            //T->update();
            //T->update_wo_pressure();
            // updating the tissues gives forces acting on the vertices after the pop. now clearly the forces from the edges area and the apical and basal lines is not well defined since there is no area or length of the edges. this problem is solved, by ignoring forces from edge with zero area
            
            // extract the apical and basal forces driving the vertices apart
            Point ApicalForce = v->dW_dva - dW_dva;
            Point BasalForce = v->dW_dvb - dW_dvb;
            
            // can those forces overcome the surface and line tensions
            double norm_apicalForce = Point::Norm(ApicalForce);// - (e->G_a + Point::Norm(coord-coord_b)*e->T_l);
            double norm_basalForce = Point::Norm(BasalForce);// - (e->G_b + Point::Norm(coord-coord_b)*e->T_l);
            //double norm_apicalForce = Point::Norm(v->dW_dva) + Point::Norm(dW_dva);
            //double norm_basalForce = Point::Norm(v->dW_dvb) + Point::Norm(dW_dvb);
            double separationForce = (norm_apicalForce+norm_basalForce)/2.;
            
            std::cout << "separationForce for counter " << counter << " is " << separationForce << std::endl;
            
            if(separationForce>maximalSeparationForce)  {
                maximalSeparationForce = separationForce;
                maximalSeparationForceConfiguration = counter;
            }
            
            // now rebuild the tissue as it was before the pop:
            
            // reconnect the edges in (it_e1+1,...,it_e2) to this
            it_e = it_e1;
            do
            {
                it_e++;
                if((*it_e)->v1==v) (*it_e)->v1=this;
                else (*it_e)->v2=this;
            }
            while(it_e!=it_e2);
            
            // remove the new edge from cells c1 and c2
            c1->ListEdges.remove(e);
            c2->ListEdges.remove(e);
            
            // remove the new cell and the new vertex from T
            T->ListEdge.pop_back();
            T->ListVertex.pop_back();
            
            counter++;
        }
    }
    
    unsigned counter2=1;
    
    // now go to the found minimal energy configuration and apply it
    if(maximalSeparationForceConfiguration!=0)
    {
        std::list<Edge*> NeighbouringEdges_sorted2 = NeighbouringEdges_sorted;
        
        // iterate through the edges
        for(std::list<Edge*>::iterator it_e1 = NeighbouringEdges_sorted2.begin();; it_e1++ )
        {
            // go to the next element in the list and check if this is the last one
            std::list<Edge*>::iterator it_e2 = (++it_e1);
            --it_e1;
            if(it_e2==NeighbouringEdges_sorted2.end()) break;
            
            Cell* c1;
            if(this == (*it_e1)->v1)
            {
                c1 = (*it_e1)->c1;
            }
            else
            {
                c1 = (*it_e1)->c2;
            }
            
            
            // iterate with the 2nd iterator til the end
            ++it_e2;
            for(;it_e2!=NeighbouringEdges_sorted2.end();it_e2++)
            {
                // first and last edge together cannot pop since they are neighbours (only possible neighbour condition that has to be captured!)
                //if(it_e1== NeighbouringEdges_sorted2.begin() && (++it_e2)==NeighbouringEdges_sorted2.end()) {it_e2--;continue;}
                if((*it_e1)->c1 == (*it_e2)->c1 || (*it_e1)->c2 == (*it_e2)->c2 || (*it_e1)->c1 == (*it_e2)->c2 || (*it_e1)->c2 == (*it_e2)->c1) {continue;}
                
                if(maximalSeparationForceConfiguration != counter2) { // not yet the right configuration
                    counter2++;
                } else { // right configuration
                    Cell* c2;
                    if(this == (*it_e2)->v1)
                    {
                        c2 = (*it_e2)->c1;
                    }
                    else
                    {
                        c2 = (*it_e2)->c2;
                    }
                    
                    // create new elements (v and e)
                    T->ListVertex.push_back(*(new Vertex(T)));
                    Vertex *v = &T->ListVertex.back();
                    v->T = T;
                    v->G_l = T->MechProp.Vertex_LineTension;
                    Edge *e = new Edge(T,v,this);
                    T->ListEdge.push_back(*e);
                    e = &T->ListEdge.back();
                    e->T = T;
                    
                    e->c1 = c1;
                    e->c2 = c2;
                    
                    int type1, type2;
                    if(c1->Type>=c2->Type) {type1=c1->Type;type2=c2->Type;}
                    else {type1=c2->Type;type2=c1->Type;}
                    e->T_l = T->MechProp.Intercell_SurfaceTension[type1][type2]; // mechanical properties of the edge are set by the types of the constituting cells
                    e->G_a = T->MechProp.Intercell_ApicalLineTension[type1][type2];
                    e->G_b = T->MechProp.Intercell_BasalLineTension[type1][type2];
                    
                    
                    // connect the edges in (it_e1+1,...,it_e2) to the new vertex
                    std::list<Edge*>::iterator it_e = it_e1;
                    do {
                        it_e++;
                        if((*it_e)->v1==this) (*it_e)->v1=v;
                        else (*it_e)->v2=v;
                    } while(it_e!=it_e2);
                    
                    for(it_e = c1->ListEdges.begin(); it_e!=c1->ListEdges.end();it_e++) {  // assign the new edge to the cells c1 and c2 at the right position
                        
                        if(*it_e==*it_e1)
                        {
                            it_e++;
                            c1->ListEdges.insert(it_e,e);
                            break;
                        }
                    }
                    
                    
                    for(it_e = c2->ListEdges.begin(); it_e != c2->ListEdges.end();it_e++)  {  // assign the new edge to the cells c1 and c2 at the right position
                        if(*it_e==*it_e2)
                        {
                            c2->ListEdges.insert(it_e,e);
                            break;
                        }
                    }
                    
                    // both vertices are at the same position
                    v->coord = coord;
                    v->coord_b = coord_b;
                    T->update();
                    
                    // updating the tissues gives forces acting on the vertices after the pop. now clearly the forces from the edges area and the apical and basal lines is not well defined since there is no area or length of the edges. this problem is solved, by ignoring forces from edge with zero area
                    
                    // extract the apical and basal forces driving the vertices apart
                    Point ApicalForce = v->dW_dva - dW_dva;
                    Point BasalForce = v->dW_dvb - dW_dvb;
                    
                    // can those forces overcome the surface and line tensions
                    double norm_apicalForce = Point::Norm(ApicalForce) - (e->Lambda_a() + Point::Norm(coord-coord_b)*e->T_l);
                    double norm_basalForce = Point::Norm(BasalForce) - (e->Lambda_b() + Point::Norm(coord-coord_b)*e->T_l);
                    double separationForce = (norm_apicalForce+norm_basalForce)/2.;
                    
                    std::cout << "separationForce after pop for counter " << maximalSeparationForceConfiguration << " is " << separationForce << std::endl;
                    
                    Point displacement_a = v->dW_dva/(Point::Norm(v->dW_dva)/edgeLength);
                    Point displacement_b = v->dW_dvb/(Point::Norm(v->dW_dvb)/edgeLength);
                    v->coord = coord+displacement_a;
                    v->coord_b = coord_b+displacement_b;
                    
                    displacement_a = dW_dva/(Point::Norm(dW_dva)/edgeLength);
                    displacement_b = dW_dvb/(Point::Norm(dW_dvb)/edgeLength);
                    
                    coord += displacement_a;
                    coord_b += displacement_b;
                    T->update();
                    
                    return 1;
                }
            }
        }
    }

    return 0; // no pop has taken place
}

bool Vertex::checkAndPop_forces(double edgeLength)
{
    updateNeighbouringEdges();
    std::list<Edge*> NeighbouringEdges_sorted2 = NeighbouringEdges_sorted;
    
    if(NeighbouringEdges_sorted.size()<4) // no pop possible if the vertex is less than 4 fold
        return 0;
    
    double maximalSeparationForce = 0;
    unsigned maximalSeparationForceConfiguration = 0;
    unsigned counter = 1;
    
    for(std::list<Edge*>::iterator it_e1 = NeighbouringEdges_sorted2.begin();it_e1 != NeighbouringEdges_sorted2.end(); it_e1++ )  { // if any of the neighbouring cells crosses the boundary do not pop
        if((*it_e1)->c1->CrossBoundary() || (*it_e1)->c2->CrossBoundary()) return 0;
    }
    
    // iterate over the edges
    for(std::list<Edge*>::iterator it_e1 = NeighbouringEdges_sorted2.begin();; it_e1++ )
    {
        // go to the next element in the list and check if this is the last one
        std::list<Edge*>::iterator it_e2 = (++it_e1);
        --it_e1;
        if(it_e2==NeighbouringEdges_sorted2.end()) break;
        
        Cell* c1;
        if(this == (*it_e1)->v1)
        {
            c1 = (*it_e1)->c1;
        }
        else
        {
            c1 = (*it_e1)->c2;
        }
        
        
        // iterate with the 2nd iterator til the end
        ++it_e2;
        for(;it_e2!=NeighbouringEdges_sorted2.end();it_e2++)
        {
            // first and last edge together cannot pop since they are neighbours (only possible neighbour condition that has to be captured!)
            //if(it_e1== NeighbouringEdges_sorted2.begin() && (++it_e2)==NeighbouringEdges_sorted2.end()) {it_e2--;continue;}
            if((*it_e1)->c1 == (*it_e2)->c1 || (*it_e1)->c2 == (*it_e2)->c2 || (*it_e1)->c1 == (*it_e2)->c2 || (*it_e1)->c2 == (*it_e2)->c1) {continue;}
            
            Cell* c2;
            
            if(this == (*it_e2)->v1)
            {
                c2 = (*it_e2)->c1;
            }
            else
            {
                c2 = (*it_e2)->c2;
            }
            
            // create new elements (v and e)
            T->ListVertex.push_back(*(new Vertex(T)));
            Vertex *v = &T->ListVertex.back();
            v->T = T;
            v->G_l = T->MechProp.Vertex_LineTension;
            Edge *e = new Edge(T,v,this);
            T->ListEdge.push_back(*e);
            e = &T->ListEdge.back();
            e->T = T;
            
            e->c1 = c1;
            e->c2 = c2;
            
            int type1, type2;
            if(c1->Type>=c2->Type) {type1=c1->Type;type2=c2->Type;}
            else {type1=c2->Type;type2=c1->Type;}
            e->T_l = T->MechProp.Intercell_SurfaceTension[type1][type2]; // mechanical properties of the edge are set by the types of the constituting cells
            e->G_a = T->MechProp.Intercell_ApicalLineTension[type1][type2];
            e->G_b = T->MechProp.Intercell_BasalLineTension[type1][type2];
            e->q_BC_l = Quadrant();
            
            if(T->isPeriodic)
            {
                e->q_v2 = Quadrant();
                e->q_c1 = Quadrant();
                e->q_c2 = Quadrant();
            }
            
            // connect the edges in (it_e1+1,...,it_e2) to the new vertex
            std::list<Edge*>::iterator it_e = it_e1;
            do {
                it_e++;
                if((*it_e)->v1==this) (*it_e)->v1=v;
                else (*it_e)->v2=v;
            } while(it_e!=it_e2);
            
            for(it_e = c1->ListEdges.begin(); it_e!=c1->ListEdges.end();it_e++) {  // assign the new edge to the cells c1 and c2 at the right position
                
                if(*it_e==*it_e1)
                {
                    it_e++;
                    c1->ListEdges.insert(it_e,e);
                    break;
                }
            }
            
            
            for(it_e = c2->ListEdges.begin(); it_e != c2->ListEdges.end();it_e++)  {  // assign the new edge to the cells c1 and c2 at the right position
                if(*it_e==*it_e2)
                {
                    c2->ListEdges.insert(it_e,e);
                    break;
                }
            }
            
            // both vertices are at the same position
            if(T->isPeriodic) v->q_b = q_b;
            v->coord = coord;
            v->coord_b = coord_b;
            
            // calculate the forces acting on the vertices
            T->update();
            //T->update_wo_pressure();
            // updating the tissues gives forces acting on the vertices after the pop. now clearly the forces from the edges area and the apical and basal lines is not well defined since there is no area or length of the edges. this problem is solved, by ignoring forces from edge with zero area
            
            // move the 2 vertices according to the current forces, and recalculate the forces projected back on their movement
//            Point displacement_a1 = v->dW_dva/(Point::Norm(v->dW_dva)/edgeLength_test);
//            Point displacement_b1 = v->dW_dvb/(Point::Norm(v->dW_dvb)/edgeLength_test);
//            v->moveVertex(displacement_a1, displacement_b1);
//            Point displacement_a2 = dW_dva/(Point::Norm(dW_dva)/edgeLength_test);
//            Point displacement_b2 = dW_dvb/(Point::Norm(dW_dvb)/edgeLength_test);
//            moveVertex(displacement_a2, displacement_b2);
//
//            Point displacementDirection_a = (displacement_a1-displacement_a2)/2;
//            Point displacementDirection_b = (displacement_b1-displacement_b2)/2;
//
//            T->update();
            
            // normalized vectors connecting centers of c1 and c2
            Point v_c1c2_a = c2->BC_a-c1->BC_a; v_c1c2_a = v_c1c2_a/Point::Norm(v_c1c2_a);
            Point v_c1c2_b = c2->BC_b-c1->BC_b; v_c1c2_b = v_c1c2_b/Point::Norm(v_c1c2_b);
            
            // extract the apical and basal forces driving the vertices apart
            Point force_a = v->dW_dva - dW_dva;
            Point force_b = v->dW_dvb - dW_dvb;
            // ... and project it on the plane orthogonal to the connections of the barycenters of the cells that are connected there after (v_c1c2_a and v_c1c2_b)
            // project the forces on the vector connecting the two points
            Point projectedForce_a = force_a;// - v_c1c2_a*force_a.dot(v_c1c2_a);
            Point projectedForce_b = force_b;// - v_c1c2_b*force_b.dot(v_c1c2_b);
            
            double e_T_l = T->MechProp.Intercell_SurfaceTension[type1][type2]; // mechanical properties of the edge are set by the types of the constituting cells
            double e_G_a = T->MechProp.Intercell_ApicalLineTension[type1][type2];
            double e_G_b = T->MechProp.Intercell_BasalLineTension[type1][type2];

            // can those forces overcome the surface and line tensions
            double norm_openingForce = Point::Norm(projectedForce_a)+Point::Norm(projectedForce_b) - e_G_a - 2*Point::Norm(coord-coord_b)*e_T_l- e_G_b;
            
            std::cout << "norm of separationForce for counter " << counter << " is " << norm_openingForce << std::endl;
            
            if(norm_openingForce>maximalSeparationForce)  {
                maximalSeparationForce = norm_openingForce;
                maximalSeparationForceConfiguration = counter;
            }
            
            // now rebuild the tissue as it was before the pop:
            
//            v->coord -= displacement_a1;
//            v->coord_b -= displacement_b1;
//            coord -= displacement_a2;
//            coord_b -= displacement_b2;
//
            // reconnect the edges in (it_e1+1,...,it_e2) to this
            it_e = it_e1;
            do
            {
                it_e++;
                if((*it_e)->v1==v) (*it_e)->v1=this;
                else (*it_e)->v2=this;
            }
            while(it_e!=it_e2);
            
            // remove the new edge from cells c1 and c2
            c1->ListEdges.remove(e);
            c2->ListEdges.remove(e);
            
            // remove the new cell and the new vertex from T
            T->ListEdge.pop_back();
            T->ListVertex.pop_back();
            
            std::list<Cell*> listCells = getNeighbouringCells();
            for(std::list<Cell*>::iterator itc = listCells.begin(); itc != listCells.end(); itc++) (*itc)->orderEdges();

            counter++;
        }
    }
    
    
    std::cout << "The maximalSeparationForceConfiguration is " << maximalSeparationForceConfiguration << std::endl;
    
    unsigned counter2=1;
    // now go to the found minimal energy configuration and apply it
    if(maximalSeparationForceConfiguration!=0)
    {
        std::list<Edge*> NeighbouringEdges_sorted2 = NeighbouringEdges_sorted;
        
        // iterate through the edges
        for(std::list<Edge*>::iterator it_e1 = NeighbouringEdges_sorted2.begin();; it_e1++ )
        {
            // go to the next element in the list and check if this is the last one
            std::list<Edge*>::iterator it_e2 = (++it_e1);
            --it_e1;
            if(it_e2==NeighbouringEdges_sorted2.end()) break;
            
            Cell* c1;
            if(this == (*it_e1)->v1)
            {
                c1 = (*it_e1)->c1;
            }
            else
            {
                c1 = (*it_e1)->c2;
            }
            
            
            // iterate with the 2nd iterator til the end
            ++it_e2;
            for(;it_e2!=NeighbouringEdges_sorted2.end();it_e2++)
            {
                // first and last edge together cannot pop since they are neighbours (only possible neighbour condition that has to be captured!)
                //if(it_e1== NeighbouringEdges_sorted2.begin() && (++it_e2)==NeighbouringEdges_sorted2.end()) {it_e2--;continue;}
                if((*it_e1)->c1 == (*it_e2)->c1 || (*it_e1)->c2 == (*it_e2)->c2 || (*it_e1)->c1 == (*it_e2)->c2 || (*it_e1)->c2 == (*it_e2)->c1) {continue;}
                
                if(maximalSeparationForceConfiguration != counter2) { // not yet the right configuration
                    counter2++;
                } else { // right configuration
                    Cell* c2;
                    if(this == (*it_e2)->v1)
                    {
                        c2 = (*it_e2)->c1;
                    }
                    else
                    {
                        c2 = (*it_e2)->c2;
                    }
                    
                    // create new elements (v and e)
                    T->ListVertex.push_back(*(new Vertex(T)));
                    Vertex *v = &T->ListVertex.back();
                    v->T = T;
                    v->G_l = T->MechProp.Vertex_LineTension;
                    Edge *e = new Edge(T,v,this);
                    T->ListEdge.push_back(*e);
                    e = &T->ListEdge.back();
                    e->T = T;
                    
                    e->c1 = c1;
                    e->c2 = c2;
                    
                    int type1, type2;
                    if(c1->Type>=c2->Type) {type1=c1->Type;type2=c2->Type;}
                    else {type1=c2->Type;type2=c1->Type;}
                    e->T_l = T->MechProp.Intercell_SurfaceTension[type1][type2]; // mechanical properties of the edge are set by the types of the constituting cells
                    e->G_a = T->MechProp.Intercell_ApicalLineTension[type1][type2];
                    e->G_b = T->MechProp.Intercell_BasalLineTension[type1][type2];
                    
                    if(T->isPeriodic)
                    {
                        e->q_BC_l = Quadrant();
                        e->q_v2 = Quadrant();
                        e->q_c1 = Quadrant();
                        e->q_c2 = Quadrant();
                    }
                    
                    // connect the edges in (it_e1+1,...,it_e2) to the new vertex
                    std::list<Edge*>::iterator it_e = it_e1;
                    do {
                        it_e++;
                        if((*it_e)->v1==this) (*it_e)->v1=v;
                        else (*it_e)->v2=v;
                    } while(it_e!=it_e2);
                    
                    for(it_e = c1->ListEdges.begin(); it_e!=c1->ListEdges.end();it_e++) {  // assign the new edge to the cells c1 and c2 at the right position
                        
                        if(*it_e==*it_e1)
                        {
                            it_e++;
                            c1->ListEdges.insert(it_e,e);
                            break;
                        }
                    }
                    
                    
                    for(it_e = c2->ListEdges.begin(); it_e != c2->ListEdges.end();it_e++)  {  // assign the new edge to the cells c1 and c2 at the right position
                        if(*it_e==*it_e2)
                        {
                            c2->ListEdges.insert(it_e,e);
                            break;
                        }
                    }
                    
                    // both vertices are at the same position
                    v->coord = coord;
                    v->coord_b = coord_b;
                    
                    if(T->isPeriodic)
                    {
                         v->q_b = Quadrant();
                    }
                    
                    std::list<Cell*> listCells = getNeighbouringCells();
                    for(std::list<Cell*>::iterator itc = listCells.begin(); itc != listCells.end(); itc++) (*itc)->orderEdges();

                    T->update();
                    
                    // updating the tissues gives forces acting on the vertices after the pop. now clearly the forces from the edges area and the apical and basal lines is not well defined since there is no area or length of the edges. this problem is solved, by ignoring forces from edge with zero area
                    
                    // extract the apical and basal forces driving the vertices apart
                    //Point ApicalForce = v->dW_dva - dW_dva;
                    //Point BasalForce = v->dW_dvb - dW_dvb;
                    
                    // can those forces overcome the surface and line tensions
                    //double norm_apicalForce = Point::Norm(ApicalForce) - (e->Lambda_a() + Point::Norm(coord-coord_b)*e->T_l);
                    //double norm_basalForce = Point::Norm(BasalForce) - (e->Lambda_b() + Point::Norm(coord-coord_b)*e->T_l);
                    //double separationForce = (norm_apicalForce+norm_basalForce)/2.;
                    std::cout << "T1 transition took place!!\n";
                    
                    Point displacement_a = v->dW_dva/(Point::Norm(v->dW_dva)/edgeLength);
                    Point displacement_b = v->dW_dvb/(Point::Norm(v->dW_dvb)/edgeLength);
                    v->coord = coord+displacement_a;
                    v->coord_b = coord_b+displacement_b;
                    
                    displacement_a = dW_dva/(Point::Norm(dW_dva)/edgeLength);
                    displacement_b = dW_dvb/(Point::Norm(dW_dvb)/edgeLength);
                    
                    coord += displacement_a;
                    coord_b += displacement_b;
                    T->update();

                    const char *fileName_char = "/Users/silvanus/MPI/EpithelialMechanics/SimulationResults/playground/T1Tests/saveAfterT1_noOptimization";
                    
                    T->writeFileToTxt(fileName_char);
                    
                    return 1; // pop took place
                    
                }
            }
        }
    }

    return 0; // no pop has taken place
}




//	Point bestPopVector, newPopVector;
//	std::list<Cell*> bestCellGroup1, bestCellGroup2, newCellGroup1, newCellGroup2;
//
//	//int flag;
//	//bool found = false;
//	Cell *c1, *c2; 
//	
//	
//	// take all pairs of cells that are neighbours to this vertex
//	for(;it_c1!=ListCells.end();it_c1++)
//	{
//		// only take neighbours so 
//		it_c2=it_c1;
//		it_c2++;
//
//		for(;it_c2!=ListCells.end();it_c2++)
//		{
//			
//			flag = ForcesFromPop(*it_c1,*it_c2,&newPopVector,&newCellGroup1,&newCellGroup2,LineTensionNewEdge,SurfaceTensionNewEdge);
//			
//			if(flag==1)
//			{
//				found=true;
//				// if the new tension vector is bigger than the old one
//				if(Point::Norm(newPopVector)>Point::Norm(bestPopVector))
//				{
//					bestPopVector = newPopVector;
//					
//					bestCellGroup1 = newCellGroup1;
//					bestCellGroup2 = newCellGroup2;
//					
//					c1 = *it_c1;
//					c2 = *it_c2;
//				}
//			}
//		}
//	}
//		
//	if(!found)
//		return 0;
//		
//	/***** initial test of the subroutine ***********/
//	
//	// for now asssume that all the poppable vertices are 4 fold and pop them somehow randomly
//	int flag;
//	bool found = false;
//	Cell *c1, *c2; 
//	
//	std::list<Cell*> cellList1, cellList2;
//	
//	for(;it_c1!=ListCells.end();it_c1++)
//	{
//		for(it_c2=ListCells.begin();it_c2!=ListCells.end();it_c2++)
//		{
//			c1 = (*it_c1); c2 = (*it_c2);
//			
//			if((*it_c1)!=(*it_c2) && !(*it_c1)->hasCommonEdgeWith(*it_c2,&e))
//			{
//				
//				Point popVector = (*it_c2)->BC_a - (*it_c1)->BC_a;
//				// calculate the new pop vectors
//				popVector = popVector/10;//*(edgeLength/Point::Norm(popVector));
//				
//				bool found1(false);
//				
//				for(std::list<Cell*>::iterator itc = ListCells.begin();itc!=ListCells.end();itc++)
//				{
//					if(!found1 && (*itc)!=c1 && (*itc)!=c2)
//					{
//						cellList1.push_back(*itc);
//						popVector = (*itc)->BC_a;
//						found1 = true;
//					}
//					else if(found1 && (*itc)!=c1 && (*itc)!=c2)
//					{
//						cellList2.push_back(*itc);
//						popVector = popVector - (*itc)->BC_a;
//						break;
//					}
//				}
//				
//				popVector = popVector / 100;
//				
//				//cellList2.push_back(*it_c1);
//				//cellList1.push_back(*it_c2);
//				
//				flag = pop(*it_c1,*it_c2,cellList1,cellList2,popVector,e);
//				
//				found = true;
//				break;
//
//			}
//			
//		}
//		
//		if(found) break;
//	}
//	
//	// now pop the edge inbetween the found cells c1 and c2
//	//flag = pop(c1,c2,bestCellGroup1,bestCellGroup2,bestPopVector*edgeLength/Point::Norm(bestPopVector),e);
//	
//	return flag;
//}



//int Vertex::pop(Cell* c1, Cell* c2, std::list<Cell*> CellGroup1, std::list<Cell*> CellGroup2, Point popVector, Edge* e )
//{
//	
//	if(c1==c2)
//	{
//		int i = 1+2;
//	}
//	
//	// first: create the new vertex in the tissue
//	Vertex vnew;T->ListVertex.push_back(vnew);
//	Vertex* v = &(T->ListVertex.back());
//	
//	// second: define a new edge
//	Edge enew; T->ListEdge.push_back(enew);
//	e = &(*(--T->ListEdge.end()));
//	e->setVertex1(this);
//	e->setVertex2(v);
//	e->setTissue(T);
//	e->setCell1(c1);
//	e->setCell2(c2);
//	e->setT_l(60);
//	e->setG_b(60);
//	e->setG_b(60);
//	
//	v->addCell(c1);
//	v->addCell(c2);
//	v->setG_l(G_l);
//	v->setTissue(T);
//
//    
//	// second step: set the new coordinates of the two vertices of the edge
//	v->setCoord(coord - popVector*0.5); // to the new vertex the second group is attached to
//	v->setBasalCoord(coord_b - popVector*0.5);
//	coord = coord + popVector*0.5; // here the first group is attached to
//	coord_b = coord_b + popVector*0.5;
//	
//	// attach the second group of cells to the new vertex
//	
//	// iterate through the group and replace this vertex everywhere by the new one
//	std::list<Cell*>::iterator it_cell = CellGroup2.begin();
//	
//	// take all the cells, remove the pointer to this vertex and replace it by the new vertex
//	for(;it_cell!=CellGroup2.end();it_cell++)
//	{
//		// add the cells to the new vertex
//		v->addCell(*it_cell);
//	}
//	
//	// take all the cells, find the edges that are connected to this vertex, remove it and push back the new one if it is not already included
//	for(it_cell = CellGroup2.begin();it_cell!=CellGroup2.end();it_cell++)
//	{
//		for(std::list<Edge*>::iterator it_edge = (*it_cell)->ListEdges.begin();it_edge != (*it_cell)->ListEdges.end();it_edge++)
//		{
//			if((*it_edge)->getVertex1() == this)
//			{
//				// set the first vertex of the edge
//				(*it_edge)->setVertex1(v);
//				
//				// now remove this edge from this as it is not connected anymore
//				ListEdges.remove(*it_edge);
//				
//				// add the edge to the new vertex
//				v->addEdge(*it_edge);
//			}
//			
//			else if ((*it_edge)->getVertex2() == this)
//			{
//				// set the second vertex of the edge
//				(*it_edge)->setVertex2(v);
//				
//				// now remove this edge from this as it is not connected anymore
//				ListEdges.remove(*it_edge);
//				
//				// add the edge to the new vertex
//				v->addEdge(*it_edge);
//			}
//				
//		}
//	}
//	
//	
//	// need to set the directionality of the edge with respect to the two cells	
//	// now iterate through the edges of cell c1 and find the existing edge that is connected to this vertex
//	std::list<Edge*>::iterator it_edge = c1->ListEdges.begin();
//
//	for(;it_edge!=c1->ListEdges.end();it_edge++)
//	{
//		int isPart = (*it_edge)->VertexIsPart(this);
//		
//		if(isPart==1)
//		{
//			if((*it_edge)->getCell1() == c1)
//			{
//				e->setCell1(c2); e->setCell2(c1);
//			}
//			else
//			{
//				e->setCell1(c1); e->setCell2(c2);
//			}
//			
//			break;
//			
//		}
//		else if(isPart==2)
//		{
//			if((*it_edge)->getCell1() == c1)
//			{
//				e->setCell1(c1); e->setCell2(c2);
//			}
//			else
//			{
//				e->setCell1(c2); e->setCell2(c1);
//			}
//			
//			break;
//		}
//	}
//	
//	v->addEdge(e);
//	addEdge(e);
//	
//	c1->ListEdges.push_back(e);
//	c2->ListEdges.push_back(e);		
//	
//	
//	ListCells.clear();
//	addCell(c1);addCell(c2);
//	for(it_cell = CellGroup1.begin();it_cell != CellGroup1.end();it_cell++)
//		addCell(*it_cell);
//	
//	
//	int size1 = ListEdges.size();
//	int size2 = ListCells.size();
//	
//	// remove the cells from cellGroup2 from this vertex
////	for(it_cell = CellGroup2.begin();it_cell != CellGroup2.end();it_cell++)
////	{
////		std::list<Cell*>::iterator found = std::find(ListCells.begin(),ListCells.end(),*it_cell);
////		
////		if(found!=ListCells.end())
////			ListCells.erase(found++);
////	}
//	
//	// connect the two 
//	
//	//		// now find out, which cell is c1 for the edge and which is c2...
//	//		ca->UpdateBarycenters();
//	//		cb->UpdateBarycenters();
//	//		
//	//		// edge vector
//	//		Point edgeVector = e->getVertex2()->getCoord() - e->getVertex1()->getCoord();
//	//		// vector to barycenters of c1 and c2
//	//		Point barycenterCaVector = ca->barycenter - e->getVertex1()->getCoord();
//	//		Point barycenterCbVector = cb->barycenter - e->getVertex1()->getCoord();
//	//		
//	//		double orientationCa = barycenterCaVector*edgeVector;
//	//		double orientationCb = barycenterCbVector*edgeVector;
//	//		
//	//		if (orientationCa>=0 && orientationCb<=0)
//	//		{
//	//			e->setCell1(ca);
//	//			e->setCell2(cb);
//	//		}
//	//		else if (orientationCa<=0 && orientationCb>=0)
//	//		{
//	//			e->setCell1(cb);
//	//			e->setCell2(ca);
//	//		}
//	//		else
//	//		{
//	//			// just for testing purposes to handle exceptions!!
//	//			throw;
//	//		}
//	
//	//T->update();
//	
//	return 1;
//	
//	
//	
//}

// algorithm that returns if the vertex is at the boundary
bool Vertex::get_isAtBoundary()
{
	// if any of the connected edges is at the boundary, i.e. has only one neighbouring cell, then the vertex is at the boundary
	for(std::list<Edge*>::iterator it_edge = NeighbouringEdges_sorted.begin();it_edge!=NeighbouringEdges_sorted.end();it_edge++)
	{
		if((*it_edge)->getCell1()==NULL || (*it_edge)->getCell2()==NULL)
			return true;
	}
	
	return false;
}

void Vertex::update_forces()
{
    // add the line force to the apical force and since forces are equal and opposite substract same term from the basal force
    Point p = lineForceOnApicalPoint();
    dW_dva += p;
    dW_dvb -= p;
    
    // contribution from elasticity of ellipsoidal basement membrane
    if(T->ellipsoid_ex>0 && T->ellipsoid_ez>0 && T->ellipsoid_K>0)
    {
        // if the apical point left the ellipsoid project it back on the surface of the ellipsoid
        if(std::pow(coord.x/T->ellipsoid_ex,2)+std::pow(coord.y/T->ellipsoid_ex,2)+std::pow(coord.z/T->ellipsoid_ez,2)>=1) {
            
            Point forceDirection = coord - coord.projectOnEllipse(T->ellipsoid_ex,T->ellipsoid_ez);
            dW_dva -= forceDirection*T->ellipsoid_K;
        }
        
        // if the basal point left the ellipsoid project it back on the surface of the ellipsoid
        if(!T->onlyApicalContributions && std::pow(coord_b.x/T->ellipsoid_ex,2)+std::pow(coord_b.y/T->ellipsoid_ex,2)+std::pow(coord_b.z/T->ellipsoid_ez,2)>=1) {
            Point forceDirection = coord_b - coord_b.projectOnEllipse(T->ellipsoid_ex,T->ellipsoid_ez);
            dW_dvb -= forceDirection*T->ellipsoid_K;
        }
    }
    
    // contribution from elasticity of basement membrane
    if(T->flatECM_stiffness>0)
    {
        // apical force on vertex
        if(coord.z > T->flatECM_apicalPosition)  {
            dW_dva += Point(0,0,T->flatECM_apicalPosition - coord.z)*T->flatECM_stiffness;
        }
        
        if(!T->onlyApicalContributions && (T->flatECM_basalPosition_bothways || coord_b.z < T->flatECM_basalPosition))
        {
            // the contribution depends on the kind of vertex - if not all connected connected vertices are of type 1, the stiffness is changed by prefactor flatECM_stiffness_alpha_nonCell1
            if(allNeighboursType1()) dW_dvb += Point(0,0,T->flatECM_basalPosition - coord_b.z)*T->flatECM_stiffness;
            else dW_dvb += Point(0,0,T->flatECM_basalPosition - coord_b.z)*(T->flatECM_stiffness*T->flatECM_stiffness_alpha_nonCell1);
        }

        // basal force on vertex
        if(!T->onlyApicalContributions && (T->flatECM_basalPosition_bothways || coord_b.z < T->flatECM_basalPosition)) {
            dW_dvb += Point(0,0,T->flatECM_basalPosition - coord_b.z)*T->flatECM_stiffness;
        }
    }

    // force on system size
    T->dW_dL -= UpdateLineStress();
}

// RETURNS THE POTENTIAL CONTRIBUTIONS FROM THE VERTEX
double Vertex::W()
{
    // lateral line tension
    double W = G_l*getLength();
    
    
    // contributions from confinement in an ellipsoid
    if(T->ellipsoid_ex>0 && T->ellipsoid_ez>0 && T->ellipsoid_K>0)
    {
        // if the apical point left the ellipsoid calculate the distance and add the energy
        if(std::pow(coord.x/T->ellipsoid_ex,2)+std::pow(coord.y/T->ellipsoid_ex,2)+std::pow(coord.z/T->ellipsoid_ez,2)>=1) {
            
            Point forceDirection = coord - coord.projectOnEllipse(T->ellipsoid_ex,T->ellipsoid_ez);

            W += T->ellipsoid_K/2 * (std::pow(forceDirection.x,2)+std::pow(forceDirection.y,2)+std::pow(forceDirection.z,2));
        }
        
        // if the basal point left the ellipsoid project it back on the surface of the ellipsoid
        if(!T->onlyApicalContributions && std::pow(coord_b.x/T->ellipsoid_ex,2)+std::pow(coord_b.y/T->ellipsoid_ex,2)+std::pow(coord_b.z/T->ellipsoid_ez,2)>=1) {
            
            Point forceDirection = coord_b - coord_b.projectOnEllipse(T->ellipsoid_ex,T->ellipsoid_ez);
            
            W += T->ellipsoid_K/2 * (std::pow(forceDirection.x,2)+std::pow(forceDirection.y,2)+std::pow(forceDirection.z,2));
        }
    }
    
    // if we consider the case of a constrained tissue between the ECM on the basal side and in fact the peripodial membrane on the apical side
    if(T->flatECM_stiffness>0)
    {
        // quadratic potential for apical side
        if(coord.z > T->flatECM_apicalPosition)  W += 0.5*T->flatECM_stiffness*pow(coord.z - T->flatECM_apicalPosition,2);
        // quadratic potential for basal side
        if(!T->onlyApicalContributions && (T->flatECM_basalPosition_bothways || coord_b.z < T->flatECM_basalPosition))
        {
            // the contribution depends on the kind of vertex - if not all connected connected vertices are of type 1, the stiffness is changed by prefactor flatECM_stiffness_alpha_nonCell1
            
            if(allNeighboursType1()) W += 0.5*T->flatECM_stiffness*pow(T->flatECM_basalPosition-coord_b.z,2);
            else W += T->flatECM_stiffness_alpha_nonCell1*0.5*T->flatECM_stiffness*pow(T->flatECM_basalPosition-coord_b.z,2);
        }
    }

    
    return W;
}

//void Vertex::updateNeighbouringEdges()
//{
//    
//    // create the list that will contain all the neighbouring edges
//    std::list<Edge*> NEdges = std::list<Edge*>();
//    
//    // iterator through all the edges of the tissue and add them to the list if they have one vertex in common
//    for(std::list<Edge>::iterator it_edge = T->ListEdge.begin();it_edge!=T->ListEdge.end();it_edge++)
//    {
//
//        Edge *e = &(*it_edge);
//        if(e->v1 == this || e->v2 == this) NEdges.push_back(e);
//    }
//    
//    // now that we have all edges sort them
//    NeighbouringEdges_sorted = std::list<Edge*>();
//    
//    // add the first edge to the sorted list and remove it from the other
//    NeighbouringEdges_sorted.push_back(NEdges.front());
//    NEdges.pop_front();
//    
//    
//    // iterate through the unsorted list and find the one, which has a same neighbour with the last one added and remove it from the old list
//    while(!NEdges.empty())
//    {
//        Edge* e = NeighbouringEdges_sorted.back();
//        
//        for(std::list<Edge*>::iterator it_e = NEdges.begin(); it_e!=NEdges.end();it_e++)
//        {
//            // if the edge has a neighbouring cell in common with the last one
//            if((*it_e)->c1 == e->c1 || (*it_e)->c1 == e->c2 || (*it_e)->c2 == e->c1 || (*it_e)->c2 == e->c2)
//            {
//                // add the edge to the sorted list
//                NeighbouringEdges_sorted.push_back(*it_e);
//                // remove it from the other one
//                NEdges.erase(it_e);
//                
//                break;
//            }
//        }
//    }
//    
////    for(std::list<Edge*>::iterator it_e = NeighbouringEdges_sorted.begin(); it_e != NeighbouringEdges_sorted.end(); ++it_e)
////    {
////        std::cout << (*it_e)->Number << " ";
////    }
//    
//    //std::cout << std::endl;
//    
//}

void Vertex::updateNeighbouringEdges()
{
    
    // create the list that will contain all the neighbouring edges
    std::list<Edge*> NEdges = std::list<Edge*>();
    
    // iterator through all the edges of the tissue and add them to the list if they have one vertex in common
    for(std::list<Edge>::iterator it_edge = T->ListEdge.begin();it_edge!=T->ListEdge.end();it_edge++)
    {
        
        Edge *e = &(*it_edge);
        if(e->v1 == this || e->v2 == this) NEdges.push_back(e);
    }
    
    // now that we have all edges sort them
    NeighbouringEdges_sorted = std::list<Edge*>();
    
    // add the first edge to the sorted list and remove it from the other
    NeighbouringEdges_sorted.push_back(NEdges.front());
    NEdges.pop_front();
    
    
    // iterate through the unsorted list and find the one, which has a same neighbour with the last one added and remove it from the old list
    while(!NEdges.empty())
    {
        Edge* e = NeighbouringEdges_sorted.back();
        bool found=false;
        for(std::list<Edge*>::iterator it_e = NEdges.begin(); it_e!=NEdges.end();it_e++)
        {
            // if the edge has the cell in ccw direction in common
            if(this == e->v1)
            {
                if((*it_e)->c1 == e->c1 || (*it_e)->c2 == e->c1)
                {
                    // add the edge to the sorted list
                    NeighbouringEdges_sorted.push_back(*it_e);
                    // remove it from the other one
                    NEdges.erase(it_e);
                    found=true;
                    break;
                }
            }
            else
            {
                if((*it_e)->c1 == e->c2 || (*it_e)->c2 == e->c2)
                {
                    // add the edge to the sorted list
                    NeighbouringEdges_sorted.push_back(*it_e);
                    // remove it from the other one
                    NEdges.erase(it_e);
                    found=true;                    
                    break;
                }
            }
            
            
        }
        
        if(!found)
        {
            break;
        }
    }
    
    //    for(std::list<Edge*>::iterator it_e = NeighbouringEdges_sorted.begin(); it_e != NeighbouringEdges_sorted.end(); ++it_e)
    //    {
    //        std::cout << (*it_e)->Number << " ";
    //    }
    
    //std::cout << std::endl;
    
}

// moves the vertex apically and basally under consideration of the possible periodicity of the tissue
Quadrant Vertex::moveVertex(Point apicalMovement, Point basalMovement)
{
    // dont allow the basal points to move in z direction - represents coupling to basal membrane
    bool allow_basal_z_change = true;
    if(!allow_basal_z_change) basalMovement.z = 0;
        
    if(!T->isPeriodic)
	{
        coord = coord + apicalMovement;
        // if there is a wall at (x,y,position_basal_wall) the vertices are not allowed to cross it
        if(T->apical_wall && coord.z>T->apical_wall_position) coord.z=T->apical_wall_position;

		if(T->onlyApicalContributions) return Quadrant(0,0);
        
        coord_b = coord_b + basalMovement;
        if(T->basal_wall && coord_b.z<T->basal_wall_position) coord_b.z=T->basal_wall_position;

        return Quadrant(0,0);
	}
	// if the tissue is periodic one has to check if the points crossed the boundary
	else
	{
		coord = coord + apicalMovement;
		coord_b = coord_b + basalMovement;

        // nothing needs to be done if the vertex does not leave the boundary
        //if(inBoundaries(T->SystemSize)) return Quadrant(0,0);
        
        // if there is a wall at (x,y,position_basal_wall) the vertices are not allowed to cross it
        //if(T->apical_wall && coord.z>T->apical_wall_position) coord.z=T->apical_wall_position;
        //if(T->basal_wall && coord_b.z<T->basal_wall_position) coord_b.z=T->basal_wall_position;

        // normalize the coordinates
        coord.x = coord.x/T->SystemSize.x;
        coord.y = coord.y/T->SystemSize.y;
        coord_b.x = coord_b.x/T->SystemSize.x;
        coord_b.y = coord_b.y/T->SystemSize.y;
        
        double quadrantCount_a_x = floor(coord.x);
        double quadrantCount_b_x = floor(coord_b.x);
        double quadrantCount_a_y = floor(coord.y);
        double quadrantCount_b_y = floor(coord_b.y);
        
        coord.x = (coord.x - quadrantCount_a_x)*T->SystemSize.x;
        coord_b.x = (coord_b.x - quadrantCount_b_x)*T->SystemSize.x;
        coord.y = (coord.y - quadrantCount_a_y)*T->SystemSize.y;
        coord_b.y = (coord_b.y - quadrantCount_b_y)*T->SystemSize.y;
        
        
        
        // commit changes if the vertex actually has left one quadrant and entered another one
        if(quadrantCount_a_x!=0 || quadrantCount_a_y!=0 || quadrantCount_b_x!=0 ||quadrantCount_b_y!=0)
        {
            // update the basal quadrant of this vertex
            q_b.x += quadrantCount_b_x - quadrantCount_a_x;
            q_b.y += quadrantCount_b_y - quadrantCount_a_y;
            
            // iterate through all the edges and cells that are connected to this vertex and tell them that their quadrant has changed
            for(std::list<Edge*>::iterator it_e = NeighbouringEdges_sorted.begin();it_e!=NeighbouringEdges_sorted.end();it_e++)
            {
                // if this vertex is number one to the edge...
                if((*it_e)->v1 == this)
                {
                    (*it_e)->q_v2.x-=quadrantCount_a_x;
                    (*it_e)->q_v2.y-=quadrantCount_a_y;
                    
                    //(*it_e)->q_c1.x-=quadrantCount_a_x;
                    //(*it_e)->q_c1.y-=quadrantCount_a_y;
                }
                else
                {
                    (*it_e)->q_v2.x+=quadrantCount_a_x;
                    (*it_e)->q_v2.y+=quadrantCount_a_y;
                    
                    //(*it_e)->q_c2.x+=quadrantCount_a_x;
                    //(*it_e)->q_c2.y+=quadrantCount_a_y;

                    //(*it_e)->q_c1.x-=quadrantCount_a_x;
                    //(*it_e)->q_c1.y-=quadrantCount_a_y;
                    //(*it_e)->q_c2.x-=quadrantCount_a_x;
                    //(*it_e)->q_c2.y-=quadrantCount_a_y;

                }
                
                (*it_e)->crossBoundary = (*it_e)->q_v2.x!=0 || (*it_e)->q_v2.y!=0 || q_b.x!=0 || q_b.y!=0;
                
                // now, if the edge crosses the boundary, inform the two adjacent cells
                if((*it_e)->crossBoundary)
                {
                    if((*it_e)->c1) (*it_e)->c1->crossBoundary = true;
                    if((*it_e)->c2) (*it_e)->c2->crossBoundary = true;
                }
                // if it doesnt anymore check if the cells still cross the boundary
                else
                {
                    if((*it_e)->c1) (*it_e)->c1->CrossBoundary();
                    if((*it_e)->c2) (*it_e)->c2->CrossBoundary();
                    
                }
            }
        }
        
        // update the center information of the neighbouring cells
        std::list<Cell*> cellList = getNeighbouringCells();
        for(std::list<Cell*>::iterator itc = cellList.begin(); itc!=cellList.end();itc++)
            (*itc)->UpdateCenter();
        
    }
    
    return q_b;
}

std::list<Cell*> Vertex::getNeighbouringCells()
{
    updateNeighbouringEdges();

    std::list<Cell*> neighbouringCells;
    
    for(std::list<Edge*>::iterator it_e = NeighbouringEdges_sorted.begin(); it_e != NeighbouringEdges_sorted.end(); it_e++)
    {
        if(this == (*it_e)->v1)  neighbouringCells.push_back((*it_e)->c1);
        else neighbouringCells.push_back((*it_e)->c2);
    }

    return neighbouringCells;
}

bool Vertex::allNeighboursType1()
{
    std::list<Cell*> neighbouringCells = getNeighbouringCells();
    
    bool allOfType1 = true;
    
    for(std::list<Cell*>::iterator it_c = neighbouringCells.begin(); it_c!=neighbouringCells.end(); it_c++)
    {
        if((*it_c)->Type!=1)
        {
            allOfType1 = false;
            break;
        }
    }
    
    return allOfType1;
}


// energy of the surrounding cells and surfaces
double Vertex::SurroundingEnergy()
{
    double energy = 0;
    
    // iterate through the edges
    for(std::list<Edge*>::iterator it_e = NeighbouringEdges_sorted.begin(); it_e != NeighbouringEdges_sorted.end(); it_e++ )
    {
        Edge* e = *it_e;
        // get the corresponding cell (in ccw direction from the edge)
        Cell* c;
        if(this == e->v1)
        {
            c = e->c1;
        }
        else
        {
            c = e->c2;
        }
        
        // add the elastic energy from the cell
        c->UpdateElasticEnergy();
        energy += c->energy();
    }
    
    return energy;
}

bool Vertex::neighbouringEdgeCrossesBoundary()
{
    updateNeighbouringEdges();

    for(std::list<Edge*>::iterator it_e = NeighbouringEdges_sorted.begin(); it_e != NeighbouringEdges_sorted.end(); it_e++) {
        if((*it_e)->CrossBoundary()) return true;
    }
    
    return false;
}

bool Vertex::inBoundaries(Quadrant Q)
{
    return coord.x>=0 && coord.x<Q.x && coord.y>=0 && coord.y<Q.y && coord_b.x>=0 && coord_b.x<Q.x && coord_b.y>=0 && coord_b.y<Q.y;
}