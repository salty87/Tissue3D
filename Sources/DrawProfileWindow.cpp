#ifdef _USE_QT_

#include <QPainter>
#include <QApplication>
#include <QPen>
#include <QKeySequence>
#include <QPrinter>
#include "DrawProfileWindow.h"
#include <string>
#include <list>
#include <fstream>
#include <sstream>
#include <math.h>
#include <time.h>
#include <iostream>
#include <QtGui>

typedef std::list<Point> CellIntersection;
typedef std::list<CellIntersection> TissueIntersection;

//#define PI 3.1415926535897932384

DrawProfileWindow::DrawProfileWindow(QWidget *parent, Integrator *integrator)
: QWidget(parent), integrator(integrator)
{
	
    setPalette(QPalette(Qt::black));
    setAutoFillBackground(true);
	setUpdatesEnabled(true);
	setGeometry(0,0,500,500);
    
	matrix.setMatrix(-1,0,0,1,0,rect().height());//matrice qui transforme les coordonnees pour les placer en coordonnees mathematiques, commencent en bas a gauche.
	matrix.scale(1,1);
	scaleFactor=1;
    drawBasalSide = true;
    
	setFocusPolicy(Qt::ClickFocus);
    
	QGraphicsScene scene;
	scene.addText("Tissue Window");
	
	
	QGraphicsView view(&scene);
	
	view.show();
	
	CurrentImage=QImage(500,500,QImage::Format_RGB32);
	CurrentImage.fill(QColor(245,208,123).rgb());
	
	propertiesOutputFile = "/Users/silvanus/MPI/EpithelialMechanics/SimulationResults/playground/test.txt";
}

//void DrawProfileWindow::paintEvent(QPaintEvent *event)
//{
//	QPainter painter(this);
//
//	painter.setMatrix(matrix, false);
//	double xCross = tissueWindow->CurrentPoint.x;
//	double yCross = tissueWindow->CurrentPoint.y;
//
//	tissueWindow->T->updateXprofile(xCross);
//	tissueWindow->T->updateYprofile(yCross);
//
//	if(tissueWindow->T->yProfileApical.size()==0 || tissueWindow->T->yProfileBasal.size()==0)
//	{
//		return;
//	}
//
//	// )T->updateXprofile(T->mainWindow->CurrentPoint.x);
//	//T->updateYprofile(T->mainWindow->CurrentPoint.y);
//
//	// define pen for the apical side
//	QPen pen_a(Qt::red, 1, Qt::SolidLine);
//	// define the pen for the basal side
//
//	tissueWindow->T->updateXprofile2(xCross);
//
//	painter.setPen(pen_a);
//
//	std::list<CellIntersection>::iterator it_isec=tissueWindow->T->xIntersection.begin();
//	std::list<Point>::iterator it_cell;
//
//	QPainterPath cellPaths;
//
//	for(;it_isec!=tissueWindow->T->xIntersection.end();it_isec++)
//	{
//		it_cell = (*it_isec).begin();
//
//		cellPaths.moveTo((*it_cell).y,(*it_cell).z);
//		it_cell++;
//
//		for(;it_cell!=(*it_isec).end();it_cell++)
//		{
//			cellPaths.lineTo((*it_cell).y,(*it_cell).z);
//		}
//
//		// close the cell:
//		cellPaths.lineTo((*it_isec).begin()->y,(*it_isec).begin()->z);
//	}
//
//	painter.drawPath(cellPaths);
//
//}


void DrawProfileWindow::paintEvent(QPaintEvent *event)
{
    Tissue *T = integrator->T;
    
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing,true);
    painter.setMatrix(matrix, false);
    
    Point positionVector_plane=currentPoint;
    
    // obtain the normal vector from the two angles alpha and beta
    Point normalVector_plane = Point(sin(sectionPlane_beta)*cos(sectionPlane_alpha),sin(sectionPlane_beta)*sin(sectionPlane_alpha),cos(sectionPlane_beta));
    
    Point e1, e2;
    
    if(normalVector_plane.x==0 && normalVector_plane.y==0)
    {
        e1 = Point(1,0,0);
        e2 = Point(0,1,0);
    }
    
    else
    {
        e1 = normalVector_plane * Point(0,0,1);
        e1 = e1 / Point::Norm(e1);
        
        e2 = normalVector_plane * e1;
    }
    
    QPainterPath LateralPaths10,LateralPaths11, LateralPaths20, LateralPaths21, LateralPaths22, ApicalPaths1, ApicalPaths2, ApicalPathsDef, BasalPaths1, BasalPaths2, BasalPathsDef;
    
    // local copy of the system size
    Quadrant SySi = T->SystemSize;
    
    // allocate local variables (vertex and barycenter positions)
    Point v1a, v1b, v2a, v2b, b_e, b_a, b_b;
    // allocate pointers for easier handling
    Cell *c;
    Edge *e;
    Vertex *v1, *v2;
    
    Point P1, P2, dummy;
    bool liesInPlane, crossesAB;
    
    double x_max, y_max, x_min, y_min;
    
    // iterate over all cells of the tissue
    for(std::list<Cell>::iterator it_cell = T->ListCell.begin();it_cell != T->ListCell.end();++it_cell)
    {
        // name cell c
        c = &(*it_cell);
        
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
            b_b = shift( c->BC_b , SySi*(c->q_BC_b));
            
            if(c == e->c1)
            {
                v1a = shift( v1->getCoord() ,	   SySi*(e->q_c1)*(-1) ) ;
                v1b = shift( v1->getBasalCoord() , SySi*(e->q_c1 - v1->q_b)*(-1)  ) ;
                v2a = shift( v2->getCoord() ,	   SySi*(e->q_c1 - e->q_v2)*(-1) );
                v2b = shift( v2->getBasalCoord() , SySi*(e->q_c1 - e->q_v2 - v2->q_b)*(-1) );
                b_e = shift( e->BC_l, SySi*(e->q_c1 - e->q_BC_l)*(-1) );
            }
            else
            {
                v2a = shift( v2->getCoord() ,	   SySi*(e->q_c2)*(-1) ) ;
                v2b = shift( v2->getBasalCoord() , SySi*(e->q_c2 - v2->q_b)*(-1)  ) ;
                v1a = shift( v1->getCoord() ,	   SySi*(e->q_c2 + e->q_v2)*(-1) );
                v1b = shift( v1->getBasalCoord() , SySi*(e->q_c2 + e->q_v2 - v1->q_b)*(-1) );
                b_e = shift( e->BC_l, SySi*(e->q_c2 + e->q_v2 - e->q_BC_l)*(-1) );
            }
             
        
        // draw the cross section with the lateral surface of the edge
        
        // first determine the kind of the lateral edge (inbetween (1 and 1) -> 11, (1 and 0) ->10, (2 and 0)->20, (2 and 1)->21, or (2 and 2)->22 ?
        
//        if((e->c1 == NULL && e->c2->Type==1) || (e->c2 == NULL && e->c1->Type==1) ) edgeKind = 10;
//        else if((e->c1 == NULL && e->c2->Type==2) || (e->c2 == NULL && e->c1->Type==2) ) edgeKind = 20;
//        else if((e->c1->Type == 1 && e->c2->Type==2) || (e->c2->Type == 1  && e->c1->Type==2) ) edgeKind = 21;
//        else if((e->c1->Type == 2 && e->c2->Type==2) || (e->c2->Type == 2  && e->c1->Type==2) ) edgeKind = 22;
//        else if((e->c1->Type == 1 && e->c2->Type==1) || (e->c2->Type == 1  && e->c1->Type==1) ) edgeKind = 11;
//        else if((e->c1->Type == 2 && e->c2->Type==3) || (e->c2->Type == 3  && e->c1->Type==2) ) edgeKind = 21;
//
        unsigned edgeKind;
        if(e->c1->Type == e->c2->Type) edgeKind = 11;
        else edgeKind = 21;
        
        
        // get the intersections with the triangles inside:
        // triangle towards v1
        if(Point::TrianglePlaneSection(v1a,v1b,b_e,normalVector_plane,positionVector_plane, P1, P2, liesInPlane, crossesAB))
        {
            double x1 = P1.dot(e1); double y1 = P1.dot(e2);
            double x2 = P2.dot(e1);double y2 = P2.dot(e2);
            
            if(edgeKind == 10) {LateralPaths10.moveTo(x1,y1);LateralPaths10.lineTo(x2,y2);}
            else if(edgeKind == 11) {LateralPaths11.moveTo(x1,y1);LateralPaths11.lineTo(x2,y2);}
            else if(edgeKind == 20) {LateralPaths20.moveTo(x1,y1);LateralPaths20.lineTo(x2,y2);}
            else if(edgeKind == 21) {LateralPaths21.moveTo(x1,y1);LateralPaths21.lineTo(x2,y2);}
            else if(edgeKind == 22) {LateralPaths22.moveTo(x1,y1);LateralPaths22.lineTo(x2,y2);}
            
            // update boundary limits
            x_max = std::max(x1, x_max);y_max = std::max(y_max, y1); x_min = std::min(x1, x_min); y_min = std::min(y1, y_min);
            x_max = std::max(x2, x_max);y_max = std::max(y_max, y2); x_min = std::min(x2, x_min); y_min = std::min(y2, y_min);
        }
        
        // triangle towards v2
        if(Point::TrianglePlaneSection(v2a,v2b,b_e,normalVector_plane,positionVector_plane, P1, P2, liesInPlane, crossesAB))
        {
            double x1 = P1.dot(e1); double y1 = P1.dot(e2);
            double x2 = P2.dot(e1);double y2 = P2.dot(e2);
            
            if(edgeKind == 10) {LateralPaths10.moveTo(x1,y1);LateralPaths10.lineTo(x2,y2);}
            else if(edgeKind == 11) {LateralPaths11.moveTo(x1,y1);LateralPaths11.lineTo(x2,y2);}
            else if(edgeKind == 20) {LateralPaths20.moveTo(x1,y1);LateralPaths20.lineTo(x2,y2);}
            else if(edgeKind == 21) {LateralPaths21.moveTo(x1,y1);LateralPaths21.lineTo(x2,y2);}
            else if(edgeKind == 22) {LateralPaths22.moveTo(x1,y1);LateralPaths22.lineTo(x2,y2);}
            
            // update boundary limits
            x_max = std::max(x1, x_max);y_max = std::max(y_max, y1); x_min = std::min(x1, x_min); y_min = std::min(y1, y_min);
            x_max = std::max(x2, x_max);y_max = std::max(y_max, y2); x_min = std::min(x2, x_min); y_min = std::min(y2, y_min);					}
        
        // triangle towards apical side
        if(Point::TrianglePlaneSection(v1a,v2a,b_e, normalVector_plane, positionVector_plane,  P1,  P2,  liesInPlane,  crossesAB))
        {
            double x1 = P1.dot(e1); double y1 = P1.dot(e2);
            double x2 = P2.dot(e1);double y2 = P2.dot(e2);
            
            if(edgeKind == 10) {LateralPaths10.moveTo(x1,y1);LateralPaths10.lineTo(x2,y2);}
            else if(edgeKind == 11) {LateralPaths11.moveTo(x1,y1);LateralPaths11.lineTo(x2,y2);}
            else if(edgeKind == 20) {LateralPaths20.moveTo(x1,y1);LateralPaths20.lineTo(x2,y2);}
            else if(edgeKind == 21) {LateralPaths21.moveTo(x1,y1);LateralPaths21.lineTo(x2,y2);}
            else if(edgeKind == 22) {LateralPaths22.moveTo(x1,y1);LateralPaths22.lineTo(x2,y2);}

            // update boundary limits
            x_max = std::max(x1, x_max);y_max = std::max(y_max, y1); x_min = std::min(x1, x_min); y_min = std::min(y1, y_min);
            x_max = std::max(x2, x_max);y_max = std::max(y_max, y2); x_min = std::min(x2, x_min); y_min = std::min(y2, y_min);
        }
        
        // triangle towards basal side
        if(Point::TrianglePlaneSection( v1b, v2b, b_e, normalVector_plane, positionVector_plane,  P1,  P2,  liesInPlane,  crossesAB))
        {
            double x1 = P1.dot(e1); double y1 = P1.dot(e2);
            double x2 = P2.dot(e1);double y2 = P2.dot(e2);
            
            // draw the edge depending on its kind
            if(edgeKind == 10) {LateralPaths10.moveTo(x1,y1);LateralPaths10.lineTo(x2,y2);}
            else if(edgeKind == 11) {LateralPaths11.moveTo(x1,y1);LateralPaths11.lineTo(x2,y2);}
            else if(edgeKind == 20) {LateralPaths20.moveTo(x1,y1);LateralPaths20.lineTo(x2,y2);}
            else if(edgeKind == 21) {LateralPaths21.moveTo(x1,y1);LateralPaths21.lineTo(x2,y2);}
            else if(edgeKind == 22) {LateralPaths22.moveTo(x1,y1);LateralPaths22.lineTo(x2,y2);}

            
            // update boundary limits
            x_max = std::max(x1, x_max);y_max = std::max(y_max, y1); x_min = std::min(x1, x_min); y_min = std::min(y1, y_min);
            x_max = std::max(x2, x_max);y_max = std::max(y_max, y2); x_min = std::min(x2, x_min); y_min = std::min(y2, y_min);
        }
        
        // now draw the possible intersections with the triangles of the apical and basal surface
        // apical triangle
        if(Point::TrianglePlaneSection( v1a, v2a,b_a, normalVector_plane, positionVector_plane,  P1,  P2,  liesInPlane,  crossesAB))
        {
            double x1 = P1.dot(e1);double y1 = P1.dot(e2);
            double x2 = P2.dot(e1);double y2 = P2.dot(e2);
            
            if(c->Type == 1) {ApicalPaths1.moveTo(x1,y1);	ApicalPaths1.lineTo(x2,y2);}
            else if(c->Type == 2) {ApicalPaths2.moveTo(x1,y1);	ApicalPaths2.lineTo(x2,y2);}
            else {ApicalPathsDef.moveTo(x1,y1);	ApicalPathsDef.lineTo(x2,y2);}
            // update the limits
            x_max = std::max(x1, x_max);	y_max = std::max(y_max, y1); x_min = std::min(x1, x_min); y_min = std::min(y1, y_min);
            x_max = std::max(x2, x_max);	y_max = std::max(y_max, y2); x_min = std::min(x2, x_min); y_min = std::min(y2, y_min);
            
        }
        
        // basal triangle
        if(Point::TrianglePlaneSection( v1b, v2b,b_b, normalVector_plane, positionVector_plane,  P1,  P2,  liesInPlane,  crossesAB))
        {
            double x1 = P1.dot(e1);double y1 = P1.dot(e2);
            double x2 = P2.dot(e1);double y2 = P2.dot(e2);
            
            if(c->Type == 1) {BasalPaths1.moveTo(x1,y1);	BasalPaths1.lineTo(x2,y2);}
            else if(c->Type == 2) {BasalPaths2.moveTo(x1,y1);	BasalPaths2.lineTo(x2,y2);}
            else {BasalPathsDef.moveTo(x1,y1);	BasalPathsDef.lineTo(x2,y2);}
            
            // update the limits
            x_max = std::max(x1, x_max);	y_max = std::max(y_max, y1); x_min = std::min(x1, x_min); y_min = std::min(y1, y_min);
            x_max = std::max(x2, x_max);	y_max = std::max(y_max, y2); x_min = std::min(x2, x_min); y_min = std::min(y2, y_min);
            
        }
        
    }
    }
    
    
    // draw rectangle around the tissue with 10% of the longer scale space in x- and y-direction of the tissue
    painter.setBrush(Qt::white);
    double dist = std::max(x_max-x_min,y_max-y_min);
    painter.drawRect(x_min - 0.1*dist,y_min - 0.1*dist,x_max-x_min + 0.2*dist,y_max - y_min + 0.2*dist);
    
    double thicknessFactor = 1;

    // pens for the apical sides
    QPen pen_a1(QColor(Qt::red), 2*thicknessFactor, Qt::SolidLine);
    QPen pen_a2(QColor(Qt::darkMagenta), 3*thicknessFactor, Qt::SolidLine);
    QPen pen_aDef(QColor(180,26,26), 3*thicknessFactor, Qt::SolidLine);
    

    // pens for the lateral intersections
    QPen pen_l10,pen_l11,pen_l20, pen_l21, pen_l22;
    // pens for the basal side
    QPen pen_b1,pen_b2,pen_bDef;
    
    // if basal side shall be drawn, also draw the lateral edges
    if(drawBasalSide) {
        pen_b1 = QPen(QColor(Qt::blue), 2*thicknessFactor, Qt::SolidLine);
        pen_b2 = QPen(QColor(Qt::darkBlue), 3*thicknessFactor, Qt::SolidLine);
        pen_bDef = QPen(QColor(16,150,139), 3*thicknessFactor, Qt::SolidLine);
        pen_l10 = QPen(Qt::darkGray, 2*thicknessFactor, Qt::SolidLine);
        pen_l11 = QPen(Qt::darkGray, 1*thicknessFactor, Qt::SolidLine);
        pen_l20 = QPen(Qt::darkGray, 2*thicknessFactor, Qt::SolidLine);
        pen_l21 = QPen(Qt::darkGray, 4*thicknessFactor, Qt::SolidLine);
        pen_l22 = QPen(Qt::darkGray, 1*thicknessFactor, Qt::SolidLine);

    } else {
        pen_b1 = Qt::NoPen;
        pen_b2 = Qt::NoPen;
        pen_bDef = Qt::NoPen;
        pen_l10 = Qt::NoPen;
        pen_l11 = Qt::NoPen;
        pen_l20 = Qt::NoPen;
        pen_l21 = Qt::NoPen;
        pen_l22 = Qt::NoPen;
    }
    
    // draw lateral sections
    painter.setPen(pen_l10); painter.drawPath(LateralPaths10);
    painter.setPen(pen_l11); painter.drawPath(LateralPaths11);
    painter.setPen(pen_l20); painter.drawPath(LateralPaths20);
    painter.setPen(pen_l21); painter.drawPath(LateralPaths21);
    painter.setPen(pen_l22); painter.drawPath(LateralPaths22);
    
    // draw apical and basal sections
    painter.setPen(pen_a1);	painter.drawPath(ApicalPaths1);
    painter.setPen(pen_a2);	painter.drawPath(ApicalPaths2);
    painter.setPen(pen_aDef); painter.drawPath(ApicalPathsDef);
    painter.setPen(pen_b1); painter.drawPath(BasalPaths1);
    painter.setPen(pen_b2);	painter.drawPath(BasalPaths2);
    painter.setPen(pen_bDef); painter.drawPath(BasalPathsDef);
    
    render();
}


std::list<Cell*> DrawProfileWindow::extractCellsInPlane()
{
    // list of all cells that are cut through by the current plane
    std::list<Cell*> cellsInPlane;

    Point positionVector_plane=currentPoint;
    
    // obtain the normal vector from the two angles alpha and beta
    Point normalVector_plane = Point(sin(sectionPlane_beta)*cos(sectionPlane_alpha),sin(sectionPlane_beta)*sin(sectionPlane_alpha),cos(sectionPlane_beta));

    // local copy of the system size
    Quadrant SySi = integrator->T->SystemSize;
    
    // allocate local variables (vertex and barycenter positions)
    Point v1a, v1b, v2a, v2b, b_e, b_a, b_b;
    // allocate pointers for easier handling
    Edge *e;
    Vertex *v1, *v2;
    
    Point P1, P2, dummy;
    bool liesInPlane, crossesAB;
    
    // iterate over all cells of the tissue
    for(std::list<Cell>::iterator it_cell = integrator->T->ListCell.begin();it_cell != integrator->T->ListCell.end();++it_cell)
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
            b_b = shift( c->BC_b , SySi*(c->q_BC_b));
            
            if(c == e->c1)
            {
                v1a = shift( v1->getCoord() ,	   SySi*(e->q_c1)*(-1) ) ;
                v1b = shift( v1->getBasalCoord() , SySi*(e->q_c1 - v1->q_b)*(-1)  ) ;
                v2a = shift( v2->getCoord() ,	   SySi*(e->q_c1 - e->q_v2)*(-1) );
                v2b = shift( v2->getBasalCoord() , SySi*(e->q_c1 - e->q_v2 - v2->q_b)*(-1) );
                b_e = shift( e->BC_l, SySi*(e->q_c1 - e->q_BC_l)*(-1) );
            }
            else
            {
                v2a = shift( v2->getCoord() ,	   SySi*(e->q_c2)*(-1) ) ;
                v2b = shift( v2->getBasalCoord() , SySi*(e->q_c2 - v2->q_b)*(-1)  ) ;
                v1a = shift( v1->getCoord() ,	   SySi*(e->q_c2 + e->q_v2)*(-1) );
                v1b = shift( v1->getBasalCoord() , SySi*(e->q_c2 + e->q_v2 - v1->q_b)*(-1) );
                b_e = shift( e->BC_l, SySi*(e->q_c2 + e->q_v2 - e->q_BC_l)*(-1) );
            }
            
            // check if the plane cuts through the cell at any triangle
            // triangle towards v1
            if(Point::TrianglePlaneSection(v1a,v2a,b_e, normalVector_plane, positionVector_plane,  P1,  P2,  liesInPlane,  crossesAB)) {cellsInPlane.push_back(c);break;}
            else if(Point::TrianglePlaneSection( v1a, v2a,b_a, normalVector_plane, positionVector_plane,  P1,  P2,  liesInPlane,  crossesAB)) {cellsInPlane.push_back(c);break;}
        }
        
    }

    return cellsInPlane;
}

// write properties from cells in Plane
void DrawProfileWindow::extractPropertiesFromCellsInPlane()
{
    // get a list of all cells that are cut by the plane
    std::list<Cell*> cellsInPlane = integrator->T->extractCellsInPlane(currentPoint, Point(sin(sectionPlane_beta)*cos(sectionPlane_alpha),sin(sectionPlane_beta)*sin(sectionPlane_alpha),cos(sectionPlane_beta)));
    
    // extract the maximal z-distance of cell of type 1 and type 2 from all cells that lie in the plane
    double max_y_1 = -999999999999;
    double min_y_2 = 999999999999;

    bool type1_found(false), type2_found(false);
    
    for(std::list<Cell*>::iterator it_c = cellsInPlane.begin(); it_c != cellsInPlane.end(); it_c++)
    {
        Cell* c = (*it_c);
        
        if(c->Type==1 && c->BC_a.y>max_y_1)
        {
            type1_found = true;
            max_y_1 = c->BC_a.y;
        } else if(c->Type==2 && c->BC_a.x<min_y_2) {
            type2_found = true;
            min_y_2 = c->BC_a.y;
        }
    }
    
    // if the height is not well defined set it to -1
    double depth;
    if(!type1_found || !type2_found) depth = -1;
    else depth = max_y_1 - min_y_2;
    
	std::ofstream output;
    
    if (!propertiesOutputFile.isEmpty())
    {
        // transform the QString fileName to a const char *
        QByteArray ba = propertiesOutputFile.toLocal8Bit();
        const char *fileName_char = ba.data();
    
        output.open(fileName_char, std::ofstream::out | std::ofstream::app);// initial position is set to be at the end of the file

        if(output.is_open()) {
            output << "depth = " << depth << std::endl;
            output.close();
        }
    }
}




//void DrawProfileWindow::updateFrame(double *xmin_, double *xmax_, double *ymin_, double *ymax_)
//{
//	// find the frame, that surrounds the whole geometry of the tissue
//
//	std::list<Vertex>::iterator it_v = tissueWindow->T->ListVertex.begin();
//
//
//	double xmin = (*it_v).planeSection.x;
//	double xmax = (*it_v).planeSection.x;
//
//	double ymin = (*it_v).planeSection.y;
//	double ymax = (*it_v).planeSection.y;
//
//	it_v++;
//
//	for(;it_v!=tissueWindow->T->ListVertex.end();it_v++)
//	{
//		if((*it_v).planeCrosses)
//		{
//
//			if((*it_v).planeSection.x < xmin) xmin = (*it_v).planeSection.x;
//			else if((*it_v).planeSection.x > xmax) xmax = (*it_v).planeSection.x;
//
//			if((*it_v).planeSection.y < ymin) ymin = (*it_v).planeSection.y;
//			else if((*it_v).planeSection.y > ymax) ymax = (*it_v).planeSection.y;
//
//		}
//	}
//
//	for(std::list<Edge>::iterator it_edge = tissueWindow->T->ListEdge.begin();it_edge!=tissueWindow->T->ListEdge.end();it_edge++)
//	{
//
//		if((*it_edge).planeCrosses_apical || (*it_edge).planeCrosses_basal)
//		{
//			if((*it_edge).planeSection_apical.x < xmin) xmin = (*it_edge).planeSection_apical.x;
//			else if((*it_edge).planeSection_apical.x > xmax) xmax = (*it_edge).planeSection_apical.x;
//
//			if((*it_edge).planeSection_apical.y < ymin) ymin = (*it_edge).planeSection_apical.y;
//			else if((*it_edge).planeSection_apical.y > ymax) ymax = (*it_edge).planeSection_apical.y;
//
//			if((*it_edge).planeSection_basal.x < xmin) xmin = (*it_edge).planeSection_basal.x;
//			else if((*it_edge).planeSection_basal.x > xmax) xmax = (*it_edge).planeSection_basal.x;
//
//			if((*it_edge).planeSection_basal.y < ymin) ymin = (*it_edge).planeSection_basal.y;
//			else if((*it_edge).planeSection_basal.y > ymax) ymax = (*it_edge).planeSection_basal.y;
//		}
//	}
//
//	*xmin_ = xmin;
//	*xmax_ = xmax;
//
//	*ymin_ = ymin;
//	*ymax_ = ymax;
//	// increase the frame by 10% in both directions
//	//xmin = xmin - (xmax - xmin)/10;
////	xmax = xmax + (xmax - xmin)/10;
////
////	ymin = ymin - (ymax - ymin)/10;
////	ymax = ymax + (ymax - ymin)/10;
////
//	//QPainter painter(this);
////
////	painter.setMatrix(matrix,false);
////
////	// draw surrounding box
////	painter.setBrush(Qt::white);
////	//painter.drawRect(xLeft-0.1*distance.x,yTop+0.1*distance.y,(int)distance.x,-(int)distance.y);
////	painter.drawRect(xmin,ymin,xmax,ymax);
////
//
//}
void DrawProfileWindow::keyPressEvent(QKeyEvent *event)
{
	if(event->key()==Qt::Key_Equal) {matrix.scale(1.2,1.2);scaleFactor*=1.2;}
	else if(event->key()==Qt::Key_Minus) {matrix.scale(0.8,0.8);scaleFactor*=0.8;}
	else if(event->key()==Qt::Key_Right) {matrix.translate(+height()/scaleFactor/10,0);}
	else if(event->key()==Qt::Key_Left) {matrix.translate(-height()/scaleFactor/10,0);}
	else if(event->key()==Qt::Key_Up) {matrix.translate(0,+height()/scaleFactor/10);}
	else if(event->key()==Qt::Key_Down) {matrix.translate(0,-height()/scaleFactor/10);}
	else if(event->key()==Qt::Key_A) {	matrix.rotate(5);}//matrice qui transforme les coordonnees pour les placer en coordonnees mathematiques, commencent en bas a gauche.
	else if(event->key()==Qt::Key_Z) {	matrix.rotate(-5);}//matrice qui transforme les coordonnees pour les placer en coordonnees mathematiques, commencent en bas a gauche.
	else if(event->key()==Qt::Key_R) {	matrix.setMatrix(-1,0,0,1,0,rect().height());matrix.scale(1,1);		scaleFactor=1;}
	else if(event->key()==Qt::Key_S) { // save the current plot as PDF
        QFileInfo info = QFileInfo(QFileDialog::getSaveFileName(this));
        
        // obtain the normal vector from the two angles alpha and beta
        Point normalVector_plane = Point(sin(sectionPlane_beta)*cos(sectionPlane_alpha),sin(sectionPlane_beta)*sin(sectionPlane_alpha),cos(sectionPlane_beta));

        writeCurrentImage(currentPoint,normalVector_plane,info.baseName(), info.path());  }
    else if(event->key()==Qt::Key_B) { drawBasalSide = !drawBasalSide;}
    else if(event->key()==Qt::Key_H) { extractPropertiesFromCellsInPlane();}
    update();
}



void DrawProfileWindow::writeCurrentImage(Point positionVector_plane, Point normalVector_plane, QString imageName, QString folderName)
{
    Tissue *T = integrator->T;
    
    bool drawPDF = true;
    bool drawJPG = false;
    
    Point e1, e2;
    
    QString fileName = folderName + "/" + imageName;
    
    std::cout << fileName.toStdString() << std::endl;
        
    if(normalVector_plane.x==0 && normalVector_plane.y==0)
    {
        e1 = Point(1,0,0);
        e2 = Point(0,1,0);
    }
    
    else
    {
        e1 = normalVector_plane * Point(0,0,1);
        e1 = e1 / Point::Norm(e1);
        
        e2 = normalVector_plane * e1;
    }
    
    QPainterPath LateralPaths10,LateralPaths11, LateralPaths20, LateralPaths21, LateralPaths22, ApicalPaths1, ApicalPaths2, ApicalPathsDef, BasalPaths1, BasalPaths2, BasalPathsDef;
        
    // local copy of the system size
    Quadrant SySi = T->SystemSize;
    
    // allocate local variables (vertex and barycenter positions)
    Point v1a, v1b, v2a, v2b, b_e, b_a, b_b;
    // allocate pointers for easier handling
    Cell *c;
    Edge *e;
    Vertex *v1, *v2;
    
    Point P1, P2, dummy;
    bool liesInPlane, crossesAB;
    
    double x_max(-1000000), y_max(-1000000), x_min(1000000), y_min(1000000);//, x, y;
    
    // iterate over all cells of the tissue
    for(std::list<Cell>::iterator it_cell = T->ListCell.begin();it_cell != T->ListCell.end();++it_cell)
    {
        // name cell c
        c = &(*it_cell);
        
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
            b_b = shift( c->BC_b , SySi*(c->q_BC_b));
            
            if(c == e->c1)
            {
                v1a = shift( v1->getCoord() ,	   SySi*(e->q_c1)*(-1) ) ;
                v1b = shift( v1->getBasalCoord() , SySi*(e->q_c1 - v1->q_b)*(-1)  ) ;
                v2a = shift( v2->getCoord() ,	   SySi*(e->q_c1 - e->q_v2)*(-1) );
                v2b = shift( v2->getBasalCoord() , SySi*(e->q_c1 - e->q_v2 - v2->q_b)*(-1) );
                b_e = shift( e->BC_l, SySi*(e->q_c1 - e->q_BC_l)*(-1) );
            }
            else
            {
                v2a = shift( v2->getCoord() ,	   SySi*(e->q_c2)*(-1) ) ;
                v2b = shift( v2->getBasalCoord() , SySi*(e->q_c2 - v2->q_b)*(-1)  ) ;
                v1a = shift( v1->getCoord() ,	   SySi*(e->q_c2 + e->q_v2)*(-1) );
                v1b = shift( v1->getBasalCoord() , SySi*(e->q_c2 + e->q_v2 - v1->q_b)*(-1) );
                b_e = shift( e->BC_l, SySi*(e->q_c2 + e->q_v2 - e->q_BC_l)*(-1) );
            }
            
            // first draw the possible intersections with the triangles of the apical and basal surface
            // apical triangle
            if(Point::TrianglePlaneSection( v1a, v2a,b_a, normalVector_plane, positionVector_plane,  P1,  P2,  liesInPlane,  crossesAB))
            {
                double x1 = P1.dot(e1);double y1 = P1.dot(e2);
                double x2 = P2.dot(e1);double y2 = P2.dot(e2);
                
                // draw intersections with colour depending on the type of the cell
                //            switch(c->Type)
                //            {
                //                case 1: ApicalPaths1.moveTo(x1,y1);	ApicalPaths1.lineTo(x2,y2);
                //                case 2: ApicalPaths2.moveTo(x1,y1);	ApicalPaths2.lineTo(x2,y2);
                //                default: ApicalPathsDef.moveTo(x1,y1); ApicalPathsDef.lineTo(x2,y2);
                //            }
                
                if(c->Type == 1) {ApicalPaths1.moveTo(x1,y1);	ApicalPaths1.lineTo(x2,y2);}
                else if(c->Type == 2) {ApicalPaths2.moveTo(x1,y1);	ApicalPaths2.lineTo(x2,y2);}
                else {ApicalPathsDef.moveTo(x1,y1);	ApicalPathsDef.lineTo(x2,y2);}
                // update the limits
                x_max = std::max(x1, x_max);	y_max = std::max(y_max, y1); x_min = std::min(x1, x_min); y_min = std::min(y1, y_min);
                x_max = std::max(x2, x_max);	y_max = std::max(y_max, y2); x_min = std::min(x2, x_min); y_min = std::min(y2, y_min);
                
            }

            // if wanted draw the rest
            if(drawBasalSide)
            {
                // draw the cross section with the lateral surface of the edge
                
                // first determine the kind of the lateral edge (inbetween (1 and 1) -> 11, (1 and 0) ->10, (2 and 0)->20, (2 and 1)->21, or (2 and 2)->22 ?
                unsigned edgeKind;
                if((e->c1 == NULL && e->c2->Type==1) || (e->c2 == NULL && e->c1->Type==1) ) edgeKind = 10;
                else if((e->c1 == NULL && e->c2->Type==2) || (e->c2 == NULL && e->c1->Type==2) ) edgeKind = 20;
                else if((e->c1->Type == 1 && e->c2->Type==2) || (e->c2->Type == 1  && e->c1->Type==2) ) edgeKind = 21;
                else if((e->c1->Type == 2 && e->c2->Type==2) || (e->c2->Type == 2  && e->c1->Type==2) ) edgeKind = 22;
                else if((e->c1->Type == 1 && e->c2->Type==1) || (e->c2->Type == 1  && e->c1->Type==1) ) edgeKind = 11;
                
                
                // get the intersections with the triangles inside:
                // triangle towards v1
                if(Point::TrianglePlaneSection(v1a,v1b,b_e,normalVector_plane,positionVector_plane, P1, P2, liesInPlane, crossesAB))
                {
                    double x1 = P1.dot(e1); double y1 = P1.dot(e2);
                    double x2 = P2.dot(e1);double y2 = P2.dot(e2);
                    
                    if(edgeKind == 10) {LateralPaths10.moveTo(x1,y1);LateralPaths10.lineTo(x2,y2);}
                    else if(edgeKind == 11) {LateralPaths11.moveTo(x1,y1);LateralPaths11.lineTo(x2,y2);}
                    else if(edgeKind == 20) {LateralPaths20.moveTo(x1,y1);LateralPaths20.lineTo(x2,y2);}
                    else if(edgeKind == 21) {LateralPaths21.moveTo(x1,y1);LateralPaths21.lineTo(x2,y2);}
                    else if(edgeKind == 22) {LateralPaths22.moveTo(x1,y1);LateralPaths22.lineTo(x2,y2);}
                    
                    // update boundary limits
                    x_max = std::max(x1, x_max);y_max = std::max(y_max, y1); x_min = std::min(x1, x_min); y_min = std::min(y1, y_min);
                    x_max = std::max(x2, x_max);y_max = std::max(y_max, y2); x_min = std::min(x2, x_min); y_min = std::min(y2, y_min);
                }
                
                // triangle towards v2
                if(Point::TrianglePlaneSection(v2a,v2b,b_e,normalVector_plane,positionVector_plane, P1, P2, liesInPlane, crossesAB))
                {
                    double x1 = P1.dot(e1); double y1 = P1.dot(e2);
                    double x2 = P2.dot(e1);double y2 = P2.dot(e2);
                    
                    if(edgeKind == 10) {LateralPaths10.moveTo(x1,y1);LateralPaths10.lineTo(x2,y2);}
                    else if(edgeKind == 11) {LateralPaths11.moveTo(x1,y1);LateralPaths11.lineTo(x2,y2);}
                    else if(edgeKind == 20) {LateralPaths20.moveTo(x1,y1);LateralPaths20.lineTo(x2,y2);}
                    else if(edgeKind == 21) {LateralPaths21.moveTo(x1,y1);LateralPaths21.lineTo(x2,y2);}
                    else if(edgeKind == 22) {LateralPaths22.moveTo(x1,y1);LateralPaths22.lineTo(x2,y2);}
                    
                    // update boundary limits
                    x_max = std::max(x1, x_max);y_max = std::max(y_max, y1); x_min = std::min(x1, x_min); y_min = std::min(y1, y_min);
                    x_max = std::max(x2, x_max);y_max = std::max(y_max, y2); x_min = std::min(x2, x_min); y_min = std::min(y2, y_min);					}
                
                // triangle towards apical side
                if(Point::TrianglePlaneSection(v1a,v2a,b_e, normalVector_plane, positionVector_plane,  P1,  P2,  liesInPlane,  crossesAB))
                {
                    double x1 = P1.dot(e1); double y1 = P1.dot(e2);
                    double x2 = P2.dot(e1);double y2 = P2.dot(e2);
                    
                    if(edgeKind == 10) {LateralPaths10.moveTo(x1,y1);LateralPaths10.lineTo(x2,y2);}
                    else if(edgeKind == 11) {LateralPaths11.moveTo(x1,y1);LateralPaths11.lineTo(x2,y2);}
                    else if(edgeKind == 20) {LateralPaths20.moveTo(x1,y1);LateralPaths20.lineTo(x2,y2);}
                    else if(edgeKind == 21) {LateralPaths21.moveTo(x1,y1);LateralPaths21.lineTo(x2,y2);}
                    else if(edgeKind == 22) {LateralPaths22.moveTo(x1,y1);LateralPaths22.lineTo(x2,y2);}
                    
                    // update boundary limits
                    x_max = std::max(x1, x_max);y_max = std::max(y_max, y1); x_min = std::min(x1, x_min); y_min = std::min(y1, y_min);
                    x_max = std::max(x2, x_max);y_max = std::max(y_max, y2); x_min = std::min(x2, x_min); y_min = std::min(y2, y_min);
                }
                
                // triangle towards basal side
                if(Point::TrianglePlaneSection( v1b, v2b, b_e, normalVector_plane, positionVector_plane,  P1,  P2,  liesInPlane,  crossesAB))
                {
                    double x1 = P1.dot(e1); double y1 = P1.dot(e2);
                    double x2 = P2.dot(e1);double y2 = P2.dot(e2);
                    
                    // draw the edge depending on its kind
                    if(edgeKind == 10) {LateralPaths10.moveTo(x1,y1);LateralPaths10.lineTo(x2,y2);}
                    else if(edgeKind == 11) {LateralPaths11.moveTo(x1,y1);LateralPaths11.lineTo(x2,y2);}
                    else if(edgeKind == 20) {LateralPaths20.moveTo(x1,y1);LateralPaths20.lineTo(x2,y2);}
                    else if(edgeKind == 21) {LateralPaths21.moveTo(x1,y1);LateralPaths21.lineTo(x2,y2);}
                    else if(edgeKind == 22) {LateralPaths22.moveTo(x1,y1);LateralPaths22.lineTo(x2,y2);}
                    
                    
                    // update boundary limits
                    x_max = std::max(x1, x_max);y_max = std::max(y_max, y1); x_min = std::min(x1, x_min); y_min = std::min(y1, y_min);
                    x_max = std::max(x2, x_max);y_max = std::max(y_max, y2); x_min = std::min(x2, x_min); y_min = std::min(y2, y_min);
                }
            
                
                // basal triangle
                if(Point::TrianglePlaneSection( v1b, v2b,b_b, normalVector_plane, positionVector_plane,  P1,  P2,  liesInPlane,  crossesAB))
                {
                    double x1 = P1.dot(e1);double y1 = P1.dot(e2);
                    double x2 = P2.dot(e1);double y2 = P2.dot(e2);
                    
                    if(c->Type == 1) {BasalPaths1.moveTo(x1,y1);	BasalPaths1.lineTo(x2,y2);}
                    else if(c->Type == 2) {BasalPaths2.moveTo(x1,y1);	BasalPaths2.lineTo(x2,y2);}
                    else {BasalPathsDef.moveTo(x1,y1);	BasalPathsDef.lineTo(x2,y2);}
                    
                    // update the limits
                    x_max = std::max(x1, x_max);	y_max = std::max(y_max, y1); x_min = std::min(x1, x_min); y_min = std::min(y1, y_min);
                    x_max = std::max(x2, x_max);	y_max = std::max(y_max, y2); x_min = std::min(x2, x_min); y_min = std::min(y2, y_min);
                    
                }
            }
        }
        
    }
        
    double smootheningFactor = 1;
    
    CurrentImage=QImage(1.2*(x_max-x_min)*smootheningFactor,1.2*(y_max-y_min)*smootheningFactor,QImage::Format_RGB32);
    CurrentImage.fill(Qt::white);

    QPrinter CurrentPDFImage(QPrinter::HighResolution);
	CurrentPDFImage.setOutputFormat(QPrinter::PdfFormat);
	CurrentPDFImage.setOutputFileName(folderName+"/"+imageName+".pdf");
	CurrentPDFImage.setFullPage(true);//impose d'utiliser toute la page
	//CurrentPDFImage.setPageSize(QPrinter::Custom);//taille del a page qui peut etre imprimee
	int pdfWidth=1.2*(x_max-x_min)*smootheningFactor;
	int pdfHeight=1.2*(y_max-y_min)*smootheningFactor;
    
    //int pdfWidth=1;
    //int pdfHeight=1;
    CurrentPDFImage.setPaperSize(QSizeF(pdfWidth,pdfHeight),QPrinter::Millimeter);//taille du papier
	CurrentPDFImage.setPageMargins(0,0,0,0,QPrinter::Millimeter);
    
    QPainter painter(&CurrentImage); // define painter on CurrentImage
    QPainter painter_pdf(&CurrentPDFImage);
    
    
    painter.setRenderHint(QPainter::Antialiasing,true);
    QTransform TPiXEL = QTransform::fromScale(painter_pdf.device()->physicalDpiX() / 25.400, painter_pdf.device()->physicalDpiY() / 25.400);
    painter_pdf.setWorldTransform(TPiXEL, false);
	
    // draw rectangle around the tissue with 10% of the longer scale space in x- and y-direction of the tissue
    painter.setBrush(Qt::white);
        
    
    QPointF offset(smootheningFactor*(-x_min+0.1*(x_max-x_min)),smootheningFactor*(-y_min+0.1*(y_max-y_min)));
    //QPointF offset(-x_min,-y_min);
    //LateralPaths10.scale(smootheningFactor,smootheningFactor);
    
//    QPainterPathStroker stroker;
//    stroker.setWidth( 2 * smootheningFactor );
//    const QPainterPath stroked = stroker.createStroke( ApicalPaths1 );
//    ApicalPaths1 = stroked.united( ApicalPaths1 );
//    
    LateralPaths10.translate(offset);
    LateralPaths11.translate(offset);
    LateralPaths20.translate(offset);
    LateralPaths21.translate(offset);
    LateralPaths22.translate(offset);
    ApicalPaths1.translate(offset);
    ApicalPaths2.translate(offset);
    ApicalPathsDef.translate(offset);
    BasalPaths1.translate(offset);
    BasalPaths2.translate(offset);
    BasalPathsDef.translate(offset);
    
    double thicknessFactor = 1;
    // pens for the lateral intersections
    QPen pen_l10(Qt::darkGray, 2*thicknessFactor, Qt::SolidLine);
    QPen pen_l11(Qt::darkGray, 1*thicknessFactor, Qt::SolidLine);
    QPen pen_l20(Qt::darkGray, 2*thicknessFactor, Qt::SolidLine);
    QPen pen_l21(Qt::darkGray, 4*thicknessFactor, Qt::SolidLine);
    QPen pen_l22(Qt::darkGray, 1*thicknessFactor, Qt::SolidLine);
    // pens for the apical sides
    QPen pen_a1(QColor(Qt::red), 2*thicknessFactor, Qt::SolidLine);
    QPen pen_a2(QColor(Qt::darkMagenta), 3*thicknessFactor, Qt::SolidLine);
    QPen pen_aDef(QColor(180,26,26), 3*thicknessFactor, Qt::SolidLine);
    
    // pens for the basal side
    QPen pen_b1(QColor(Qt::blue), 2*thicknessFactor, Qt::SolidLine);
    QPen pen_b2(QColor(Qt::darkBlue), 3*thicknessFactor, Qt::SolidLine);
    QPen pen_bDef(QColor(16,150,139), 3*thicknessFactor, Qt::SolidLine);
    
    
    // draw lateral sections
    painter.setPen(pen_l10); painter.drawPath(LateralPaths10);
    painter.setPen(pen_l11); painter.drawPath(LateralPaths11);
    painter.setPen(pen_l20); painter.drawPath(LateralPaths20);
    painter.setPen(pen_l21); painter.drawPath(LateralPaths21);
    painter.setPen(pen_l22); painter.drawPath(LateralPaths22);
    
    // draw apical and basal sections
    painter.setPen(pen_a1);	painter.drawPath(ApicalPaths1);
    painter.setPen(pen_a2);	painter.drawPath(ApicalPaths2);
    painter.setPen(pen_aDef); painter.drawPath(ApicalPathsDef);
    painter.setPen(pen_b1); painter.drawPath(BasalPaths1);
    painter.setPen(pen_b2);	painter.drawPath(BasalPaths2);
    painter.setPen(pen_bDef); painter.drawPath(BasalPathsDef);

    // draw lateral sections
    painter_pdf.setPen(pen_l10); painter_pdf.drawPath(LateralPaths10);
    painter_pdf.setPen(pen_l11); painter_pdf.drawPath(LateralPaths11);
    painter_pdf.setPen(pen_l20); painter_pdf.drawPath(LateralPaths20);
    painter_pdf.setPen(pen_l21); painter_pdf.drawPath(LateralPaths21);
    painter_pdf.setPen(pen_l22); painter_pdf.drawPath(LateralPaths22);
    
    // draw apical and basal sections
    painter_pdf.setPen(pen_a1);	painter_pdf.drawPath(ApicalPaths1);
    painter_pdf.setPen(pen_a2);	painter_pdf.drawPath(ApicalPaths2);
    painter_pdf.setPen(pen_aDef); painter_pdf.drawPath(ApicalPathsDef);
    painter_pdf.setPen(pen_b1); painter_pdf.drawPath(BasalPaths1);
    painter_pdf.setPen(pen_b2);	painter_pdf.drawPath(BasalPaths2);
    painter_pdf.setPen(pen_bDef); painter_pdf.drawPath(BasalPathsDef);

    
    //QPainter painter2 = painter;
    
    
    //painter.drawImage( QRect(0,0,500,500), CurrentImage,QRect(x_min - 0.1*dist,y_min - 0.1*dist,x_max-x_min + 0.2*dist,y_max - y_min + 0.2*dist));
    
    
    //QString path = QDir::currentPath() + "/Tissue1.png";
    //std::cout << (folderName+"/"+imageName+".jpg").toStdString() << std::endl;
    
    if(drawJPG) CurrentImage.save(folderName+"/"+imageName+".jpg","JPG",-1);
}

void DrawProfileWindow::mutateCellsInHalf(){
    
    Point positionVector_plane=currentPoint;
    
    // obtain the normal vector from the two angles alpha and beta
    Point normalVector_plane = Point(sin(sectionPlane_beta)*cos(sectionPlane_alpha),sin(sectionPlane_beta)*sin(sectionPlane_alpha),cos(sectionPlane_beta));

    for(std::list<Cell>::iterator it_cell = integrator->T->ListCell.begin();it_cell != integrator->T->ListCell.end();++it_cell)
    {
        Cell* c = &(*it_cell);
        if((c->BC_a-positionVector_plane).dot(normalVector_plane)>0)
        {
            c->Type = 2;
        } else {
            c->Type = 1;
        }
    }
    
    integrator->T->setMechanicalProperties();
    
    integrator->T->setBasalTensionGradient(3,5);
    
}

void DrawProfileWindow::render(void){
//    QImage image(this->size(),QImage::Format_RGB444);
//    QPainter p(&image);
//    image.fill(255);
//    this->render(&p);
//    
//    QString dest = "/Programming/Tissue/" + QString::number(1).rightJustified(4,'0')+".jpg";
//    image.save(dest);
}


//void FOO::paintEvent(QPaintEvent *event){
//    QDeclarativeView::paintEvent(event);
//    render();
//}

#endif // _USE_QT_