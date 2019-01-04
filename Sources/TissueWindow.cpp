/****************************************************************************
**
** Copyright (C) 2008 Nokia Corporation and/or its subsidiary(-ies).
** Contact: Qt Software Information (qt-info@nokia.com)
**
** This file is part of the example classes of the Qt Toolkit.
**s
** Commercial Usage
** Licensees holding valid Qt Commercial licenses may use this file in
** accordance with the Qt Commercial License Agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and Nokia.
**
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License versions 2.0 or 3.0 as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file.  Please review the following information
** to ensure GNU General Public Licensing requirements will be met:
** http://www.fsf.org/licensing/licenses/info/GPLv2.html and
** http://www.gnu.org/copyleft/gpl.html.  In addition, as a special
** exception, Nokia gives you certain additional rights. These rights
** are described in the Nokia Qt GPL Exception version 1.3, included in
** the file GPL_EXCEPTION.txt in this package.
**
** Qt for Windows(R) Licensees
** As a special exception, Nokia, as the sole copyright holder for Qt
** Designer, grants users of the Qt/Eclipse Integration plug-in the
** right for the Qt/Eclipse Integration to link to functionality
** provided by Qt Designer and its related libraries.
**
** If you are unsure which license is appropriate for your use, please
** contact the sales department at qt-sales@nokia.com.
**
****************************************************************************/


#ifdef _USE_QT_

#include <QPainter>
#include <QPen>
#include <QKeySequence>
#include <QPrinter>
#include <QImage>
#include "TissueWindow.h"
#include <string>
#include <sstream>
#include <math.h>
#include <time.h>
#include <iostream>
#include <QtGui>

//#define PI 3.1415926535897932384

TissueWindow::TissueWindow(QWidget *parent, Integrator *integrator)
    : QWidget(parent), integrator(integrator)
{

    // INITIALIZE OPTIONS
    QColor lightGreen = Qt::white;
    lightGreen.setAlpha(60);
	setPalette(QPalette(Qt::white));
	setAutoFillBackground(true);
	setUpdatesEnabled(true);	
	setFocusPolicy(Qt::ClickFocus);
	setGeometry(1000,50,200,600);
	matrix.setMatrix(1,0,0,-1,0,rect().height());
	matrix.scale(1,1);
	scaleFactor=1;
	CurrentImage=QImage(1000,600,QImage::Format_RGB32);
	
	
	drawBarycenters=false;
	drawApicalSide=true;
	drawBasalSide=true;
	drawRectangleForPeriodicTissue=false;
    ProfileIsActivated = false;
    
    // define the pens for the apical side
	pen_a10 = QPen(Qt::red, 1, Qt::SolidLine);
	pen_a20 = QPen(Qt::red, 1, Qt::SolidLine);
	pen_a11 = QPen(Qt::gray, 1, Qt::SolidLine);
	pen_a22 = QPen(Qt::gray, 1, Qt::SolidLine);
	pen_a12 = QPen(Qt::red, 3, Qt::SolidLine);
    
	// define the pens for the basal side
	pen_b10 = QPen(Qt::blue, 3, Qt::SolidLine);
	pen_b20 = QPen(Qt::blue, 3, Qt::SolidLine);
	pen_b11 = QPen(Qt::blue, 1, Qt::SolidLine);
	pen_b22 = QPen(Qt::blue, 1, Qt::SolidLine);
	pen_b12 = QPen(Qt::darkGreen, 3, Qt::SolidLine);
	
	// and the pen for the lateral edges
	pen_l = QPen(Qt::white, 1, Qt::SolidLine);
}
	
	
void TissueWindow::highlightEdge(QPainter& painter, QPainter& painter_b, Edge* e)
{
	// set colours
	painter.setBrush(Qt::NoBrush);
	painter_b.setBrush(Qt::NoBrush);
	
	// set pens
	QPen pen(Qt::blue, 2, Qt::SolidLine);
	QPen pen_b(Qt::red, 2, Qt::SolidLine);

	painter_b.setPen(pen_b);
	painter.setPen(pen);
	
	QPainterPath CellPath;
	QPainterPath CellPath_b;

	// go to starting point
	CellPath.moveTo(e->getVertex1()->getCoord().x,e->getVertex1()->getCoord().y);
	CellPath.lineTo(e->getVertex2()->getCoord().x,e->getVertex2()->getCoord().y);

	// draw
	CellPath_b.moveTo(e->getVertex1()->getBasalCoord().x,e->getVertex1()->getBasalCoord().y);
	CellPath_b.lineTo(e->getVertex2()->getBasalCoord().x,e->getVertex2()->getBasalCoord().y);
	
	painter.drawPath(CellPath);
	painter_b.drawPath(CellPath_b);

}

void TissueWindow::highlightVertex(QPainter& painter, QPainter& painter_b, Vertex* v)
{
	painter.setBrush(Qt::NoBrush);
	QPen pen(Qt::blue, 1, Qt::SolidLine);
	painter.setPen(pen);
	painter.drawEllipse(QPointF(v->getCoord().x,v->getCoord().y),3,3);
	
	painter_b.setBrush(Qt::NoBrush);
	QPen pen_b(Qt::red, 1, Qt::SolidLine);
	painter_b.setPen(pen_b);
	painter_b.drawEllipse(QPointF(v->getBasalCoord().x,v->getBasalCoord().y),3,3);
	
}

void TissueWindow::paintEvent(QPaintEvent * /* event */)
{ 
    // if there is no tissue, return right away

	QPainter painter(this);
	
	painter.setWorldMatrix(matrix,false);
	
	// set coordinate system according to the position of the tissue
	integrator->T->UpdateExtPositions();
	Point bottomLeft=integrator->T->extPosition_bottomLeft;
   	Point topRight=integrator->T->extPosition_topRight;
	distance=(topRight-bottomLeft);

	bottomLeft = bottomLeft - distance*0.1;
	topRight = topRight + distance*0.1;
	
	// draw surrounding box
	painter.setBrush(Qt::white);
	//painter.drawRect(xLeft-0.1*distance.x,yTop+0.1*distance.y,(int)distance.x,-(int)distance.y);
	//if(!integrator->T->boundary_isSet) painter.drawRect(ExtPositions.BottomLeft.x,ExtPositions.BottomLeft.y,distance.x*1.2,distance.y*1.2);

	
	painter.setRenderHint(QPainter::Antialiasing, true);


	
	//painter.setViewport(xLeft,yTop,(int)distance.x,-(int)distance.y);
	//painter.setWindow(xLeft,yTop,(int)distance.x,-(int)distance.y);
	//int side = qMin(width(), height());
	//int x = (width() - side / 2);
	//int y = (height() - side / 2);
	
	//painter.setViewport(x, y, side, side);
	//painter.setWindow(x,y,side,side);
	//painter.translate(distance.x*0.1, -distance.y*0.1);
	if(!integrator->T->isPeriodic)	{
        
		std::list<Cell>::iterator it;
		
		// if there exists a circular boundary around the tissue, then draw it 
		if(integrator->T->boundary_isSet)
		{
			QPen pen_boundary(Qt::black, 3, Qt::SolidLine);
			painter.setPen(pen_boundary);
			painter.setBrush(Qt::white);

			QPoint center = QPoint(0,0);
			int radius = int(std::max((int)integrator->T->boundary_apicalRadius*1.1,(int)integrator->T->boundary_basalRadius*1.1));
			painter.drawEllipse( center, radius, radius);
		}	
		
		drawTissue(painter);

    } else { // periodic tissue...
	
		if(drawRectangleForPeriodicTissue) {
            // first draw the boundary of the tissue
            QPainterPath BoundaryPath;
            
            if(integrator->T->SystemSize.y>0)
            {
                
                BoundaryPath.moveTo(0,0);
                
                BoundaryPath.lineTo(integrator->T->SystemSize.x,0);
                BoundaryPath.lineTo(integrator->T->SystemSize.x,integrator->T->SystemSize.y);
                BoundaryPath.lineTo(0,integrator->T->SystemSize.y);
                BoundaryPath.lineTo(0,0);
            }
            else 
            {
                BoundaryPath.moveTo(0,-10000);
                BoundaryPath.lineTo(0,10000);
                
                BoundaryPath.moveTo(integrator->T->SystemSize.x,-10000);
                BoundaryPath.lineTo(0,10000);
            }
            
            
            QPen pen_boundary(Qt::black, 3, Qt::SolidLine);
            
            painter.setPen(pen_boundary);
            painter.drawPath(BoundaryPath);
		}
				
		drawTissue(painter);
		
		
	}
	
	
	QImage img(this->size(),QImage::Format_RGB32);
	QPainter painter2(&img);
	this->render(&painter2);
	img.save("/file1.jpg");

    // std::cout << "finished Painting TissueWindow" << std::endl;

	//painter.setWindow(QRect(-50, -50, 100, 100));

	
}

void TissueWindow::drawTissue(QPainter& painter)
{
    if(!integrator->T->isPeriodic) { // draw not periodic tissue
        
        std::list<Edge>::iterator it_edge;
        
        QPainterPath CellPath;
        QPainterPath CellPath_b;
        QPainterPath CellPath_i;
        QPainterPath CellPath_BC_a;
        bool fillCells_basal(false);
        bool fillCells_apical(true);
        bool fillLateralSides(true);
        
        
        // at first: draw the basal side of the tissue
        for(it_edge = integrator->T->ListEdge.begin();it_edge!=integrator->T->ListEdge.end();it_edge++)
        {
            Edge* e = &(*it_edge);
            
            // remove everything from CellPath_b
            CellPath_b = QPainterPath();
            
            CellPath_b.moveTo(e->v1->coord_b.x ,e->v1->coord_b.y);
            CellPath_b.lineTo(e->v2->coord_b.x ,e->v2->coord_b.y);
            

            fillCells_basal=false;
            // draw the apical line of the edge
            if(drawBasalSide)
            {
                
                if(fillCells_basal)
                {
                    //QPen noPen(Qt::NoPen);
                    //painter.setPen(noPen);
                    
                    QPen noBarycenterConnections(QColor(24,116,205),2);
                    painter.setPen(noBarycenterConnections);
                    
                    if(e->getCell1()!=NULL)
                    {
                        if(e->getCell1()->Type==1)	
                        {
                            noBarycenterConnections.setColor(QColor(24,116,205));
                            painter.setPen(noBarycenterConnections);
                            painter.setBrush(QColor(24,116,205));
                        }
                        else 
                        {
                            noBarycenterConnections.setColor(QColor(16,78,139));
                            painter.setPen(noBarycenterConnections);
                            painter.setBrush(QColor(16,78,139));
                        }
                        
                        QPolygon polygon;
                        polygon << QPoint(e->getVertex1()->getBasalCoord().x, e->getVertex1()->getBasalCoord().y) << QPoint(e->getVertex2()->getBasalCoord().x, e->getVertex2()->getBasalCoord().y) << QPoint(e->getCell1()->BC_b.x,e->getCell1()->BC_b.y);
                        //painter.drawPolygon(polygon);
                    }
                    
                    if(e->getCell2()!=NULL)
                    {
                        if(e->getCell2()->Type==1)	
                        {
                            noBarycenterConnections.setColor(QColor(24,116,205));
                            painter.setPen(noBarycenterConnections);
                            painter.setBrush(QColor(24,116,205));
                        }
                        else if(e->getCell2()->Type==2)
                        {
                            noBarycenterConnections.setColor(QColor(16,78,139));
                            painter.setPen(noBarycenterConnections);
                            painter.setBrush(QColor(16,78,139));
                        }
                        else
                        {
                            noBarycenterConnections.setColor(QColor(30,4,150));
                            painter.setPen(noBarycenterConnections);
                            painter.setBrush(QColor(12,178,39));
                        }
                        QPolygon polygon;
                        polygon << QPoint(e->getVertex1()->getBasalCoord().x, e->getVertex1()->getBasalCoord().y) << QPoint(e->getVertex2()->getBasalCoord().x, e->getVertex2()->getBasalCoord().y) << QPoint(e->getCell2()->BC_b.x,e->getCell2()->BC_b.y);
                        //painter.drawPolygon(polygon);
                    }
                    
                }
                
                // if the cell is at the boundary
                if(e->getCell1()==NULL)
                {
                    if(e->getCell2()->Type==1)	painter.setPen(pen_b10);
                    else painter.setPen(pen_b20);
                }
                else if(e->getCell2()==NULL)
                {
                    if(e->getCell1()->Type==1)	painter.setPen(pen_b10);
                    else painter.setPen(pen_b20);
                }
                // if the edge is in the middle, check which two cell type it connects
                else
                {
                    if((e->getCell1()->Type==1 && e->getCell2()->Type!=1) || (e->getCell2()->Type==1 && e->getCell1()->Type!=1) )	painter.setPen(pen_b12);
                    else if(e->getCell1()->Type==1 || e->getCell2()->Type==1) painter.setPen(pen_b11);
                    else painter.setPen(pen_b11);
                }
                
                
                painter.drawPath(CellPath_b);
                
            }
            
        }
        
        // now draw the lateral parts of the tissue
        if(drawApicalSide&&drawBasalSide)
        {
            
            for(it_edge = integrator->T->ListEdge.begin();it_edge!=integrator->T->ListEdge.end();it_edge++)
            {
                // clear path
                CellPath_i = QPainterPath();

                Edge* e = &(*it_edge);

                if(fillLateralSides && (e->getCell1()==NULL || e->getCell2()==NULL))
                {
                    painter.setPen(pen_l);
                    
                    if((e->getCell1()) && e->getCell1()->Type==1)	painter.setBrush(Qt::lightGray);
                    else painter.setBrush(QColor(216,191,216));
                    QPolygon polygon;
                    polygon << QPoint(e->getVertex1()->getCoord().x, e->getVertex1()->getCoord().y) << QPoint(e->getVertex2()->getCoord().x, e->getVertex2()->getCoord().y) << QPoint(e->getVertex2()->getBasalCoord().x, e->getVertex2()->getBasalCoord().y) << QPoint(e->getVertex1()->getBasalCoord().x, e->getVertex1()->getBasalCoord().y);
                    painter.drawPolygon(polygon);
                    
                }
            }
        }
        
        
        for(it_edge = integrator->T->ListEdge.begin();it_edge!=integrator->T->ListEdge.end();it_edge++)
        {
            // empty the cell path
            CellPath = QPainterPath();
            
            Edge* e = &(*it_edge);
            
            CellPath.moveTo(e->getVertex1()->getCoord().x ,e->getVertex1()->getCoord().y);
            CellPath.lineTo(e->getVertex2()->getCoord().x ,e->getVertex2()->getCoord().y);
        
            // draw the apical line of the edge	
            if(drawApicalSide)
            {
                
                if(fillCells_apical)
                {
                    QPen noPen(Qt::NoPen);
                    painter.setPen(noPen);
                    
                    if(e->getCell1()!=NULL)
                    {
                        if(e->getCell1()->Type==1)	painter.setBrush(QColor(205,38,38));
                        else painter.setBrush(QColor(139,26,26));
                        //painter.setOpacity(0.7);
                        QPolygon polygon;
                        //if(fillCells_basal) painter.setOpacity(50);
                        //painter.drawPolygon(polygon);
                        //painter.setOpacity(100);
                    }
                    
                    if(e->getCell2()!=NULL)
                    {
                        if(e->getCell2()->Type==1)	painter.setBrush(QColor(205,38,38));
                        else painter.setBrush(QColor(139,26,26));
                        QPolygon polygon;
                        polygon << QPoint(e->getVertex1()->getCoord().x, e->getVertex1()->getCoord().y) << QPoint(e->getVertex2()->getCoord().x, e->getVertex2()->getCoord().y) << QPoint(e->getCell2()->BC_a.x,e->getCell2()->BC_a.y);
                        //if(fillCells_basal) painter.setOpacity(50);
                        //painter.drawPolygon(polygon);
                        //painter.setOpacity(100);			
                    }
                    
                }
                
                // if the cell is at the boundary
                if(e->getCell1()==NULL)
                {
                    if(e->getCell2()->Type==1)	painter.setPen(pen_a10);
                    else painter.setPen(pen_a20);
                }
                else if(e->getCell2()==NULL)
                {
                    if(e->getCell1()->Type==1)	painter.setPen(pen_a10);
                    else painter.setPen(pen_a20);
                }
                // if the edge is in the middle, check which two cell type it connects
                else
                {
                    if(e->getCell1()->Type==1 && e->getCell2()->Type==1)	painter.setPen(pen_a11);
                    else if(e->getCell1()->Type==2 && e->getCell2()->Type==2) painter.setPen(pen_a22);
                    else painter.setPen(pen_a12);
                }
                
                painter.drawPath(CellPath);
                
            }
        }
        
        
        
    } else { // draw the tissue in the periodic case
        
        drawBasalSide = false;
        
        Quadrant SystemSize = integrator->T->SystemSize;
        
        std::list<Cell>::iterator it_cell = integrator->T->ListCell.begin();
        
        
        QPen redPen_thick(Qt::red);
        redPen_thick.setWidth(3.0);
        QPen redPen_thin(Qt::red);
        redPen_thin.setWidth(0.4);
        
        
        QPen whitePen_thick(Qt::white);
        whitePen_thick.setWidth(1);
        
        QPen whitePen_thin(Qt::white);
        
        QPen bluePen_thick(Qt::blue);
        bluePen_thick.setWidth(3.0);
        
        QPen magentaPen_thick(Qt::darkMagenta);
        magentaPen_thick.setWidth(3.0);

        QPen blackPen(Qt::black);
        blackPen.setWidth(0.6);
        
        
        
        QPen greenPen_thick(Qt::darkGreen);
        greenPen_thick.setWidth(2);
        QPen darkgreenPen(Qt::darkGreen);
        darkgreenPen.setWidth(1);
        
        //QPen bluePen_thin(Qt::blue);
        QPen bluePen_thin(Qt::NoPen);
        bluePen_thin.setWidthF(0.3);
        
        QPen noPen(Qt::NoPen);
        
        //QColor QColor_basal_type1(205,38,38);
        QColor QColor_basal_type1(Qt::darkGray);
        QColor QColor_basal_type2(Qt::darkGreen);
        QColor QColor_basal_type3(90,160,120);
        QColor QColor_basal_default(235,138,8);
        
        // iterate through all the cells twice: first draw the basal side, then draw the apical lines
        if(drawBasalSide) {
            
            for(it_cell = integrator->T->ListCell.begin();it_cell!=integrator->T->ListCell.end();it_cell++) {
                
                Cell* c = &(*it_cell);
                QPainterPath BasalLines, BasalLines12;
                
                // set the fill colour of the basal sides, depending on the kind of the cells
                if(c->Type==1) painter.setBrush(Qt::darkGray);
                else if(c->Type==2) painter.setBrush(Qt::darkGreen);
                else if(c->Type==3) painter.setBrush(QColor_basal_type3);
                else  painter.setBrush(QColor_basal_default);
                
                // the following way of drawing allows for the inclusion of the lines to the BC_a of the cells
                
                if(!c->crossBoundary) { // if the cell doesn't cross the boundary
                    // dont draw the basal of the triangles
                    painter.setPen(noPen);
                    
                    for(std::list<Edge*>::iterator it_edge = c->ListEdges.begin();it_edge != c->ListEdges.end();it_edge++) {
                        Edge* e = *it_edge;
                        
                        if(e->c1->Type == e->c2->Type) { // draw normal edges
                            BasalLines.moveTo(e->getVertex1()->getBasalCoord().x, e->getVertex1()->getBasalCoord().y);
                            BasalLines.lineTo(e->getVertex2()->getBasalCoord().x, e->getVertex2()->getBasalCoord().y);
                        } else { // draw edges between cells of types 1 and 2
                            BasalLines12.moveTo(e->getVertex1()->getBasalCoord().x, e->getVertex1()->getBasalCoord().y);
                            BasalLines12.lineTo(e->getVertex2()->getBasalCoord().x, e->getVertex2()->getBasalCoord().y);
                        }
                        
                        
                        QPolygon polygon;
                        polygon << QPoint(e->getVertex1()->getBasalCoord().x, e->getVertex1()->getBasalCoord().y) << QPoint(e->getVertex2()->getBasalCoord().x, e->getVertex2()->getBasalCoord().y) << QPoint(c->BC_b.x,c->BC_b.y);
                        painter.drawPolygon(polygon);
                    }
                    
                    painter.setPen(blackPen);
                    painter.drawPath(BasalLines);
                    painter.setPen(magentaPen_thick);
                    painter.drawPath(BasalLines12);
                } else { // if the cell crosses the boundary draw the phantom cell around the point where the actual BC_a of the cell is
                    
                    // dont draw the basal of the triangles
                    painter.setPen(noPen);
                    //painter.setPen(redPen_thick);
                    
                    //painter.setBrush(Qt::NoBrush);
                    
                    for(std::list<Edge*>::iterator it_edge = c->ListEdges.begin();it_edge != c->ListEdges.end();it_edge++)
                    {
                        Edge* e = *it_edge;
                        
                        if(e->c1 == c) {

                            // position the vertices such that they are all connected to the BC_a of the cell
                            Point BC_b_rel = shift( c->BC_b, SystemSize*c->q_BC_b );
                            Point v1b_rel = shift( e->v1->coord_b, SystemSize*(e->q_c1-e->v1->q_b)*(-1)   );
                            Point v2b_rel = shift( e->v2->coord_b, SystemSize*(e->q_c1-e->q_v2-e->v2->q_b)*(-1) );
                            
                            // draw filled polygons on the basal side
                            QPolygon polygon;
                            polygon << QPoint(v1b_rel.x, v1b_rel.y) << QPoint(v2b_rel.x, v2b_rel.y) << QPoint(BC_b_rel.x,BC_b_rel.y);
                            painter.drawPolygon(polygon);
                            
                            // draw thicker lines, at the connections of the vertices
                            
                            if(e->c1->Type == e->c2->Type) { // draw normal edges
                                BasalLines.moveTo(v1b_rel.x, v1b_rel.y);
                                BasalLines.lineTo(v2b_rel.x, v2b_rel.y);
                            } else { // draw edges between cells of types 1 and 2
                                BasalLines12.moveTo(v1b_rel.x, v1b_rel.y);
                                BasalLines12.lineTo(v2b_rel.x, v2b_rel.y);
                            }
                            
                        } else {
                            
                            Point BC_b_rel = shift( c->BC_b, SystemSize*c->q_BC_b );
                            Point v1b_rel = shift( e->v1->coord_b, SystemSize*(e->q_c2 + e->q_v2 - e->v1->q_b)*(-1)   );
                            Point v2b_rel = shift( e->v2->coord_b, SystemSize*(e->q_c2-e->v2->q_b)*(-1) );
                            
                            // draw filled polygons on the basal side
                            QPolygon polygon;
                            polygon << QPoint(v1b_rel.x, v1b_rel.y) << QPoint(v2b_rel.x, v2b_rel.y) << QPoint(BC_b_rel.x,BC_b_rel.y);
                            painter.drawPolygon(polygon);
                            
                            // draw thicker lines, at the connections of the vertices
                            if(e->c1->Type == e->c2->Type) { // draw normal edges
                                BasalLines.moveTo(v1b_rel.x, v1b_rel.y);
                                BasalLines.lineTo(v2b_rel.x, v2b_rel.y);
                            } else { // draw edges between cells of types 1 and 2
                                BasalLines12.moveTo(v1b_rel.x, v1b_rel.y);
                                BasalLines12.lineTo(v2b_rel.x, v2b_rel.y);
                            }

                            
                        }
                        
                        
                    }
                    
                    painter.setPen(blackPen);
                    painter.drawPath(BasalLines);
                    painter.setPen(magentaPen_thick);
                    painter.drawPath(BasalLines12);

                    
                }
                
            }
        }
        
        // iterate through all the cells twice: first draw the basal side, then draw the apical lines
        if(drawApicalSide)
        {
            // dont fill out the apical surfaces, such that we can still see the basal side
            //painter.setBrush(QColor(190,30,30));
            painter.setBrush(Qt::NoBrush);
            
            QPen apicalPen11 = redPen_thin;
            QPen apicalPen22 = QPen(Qt::darkMagenta);
            QPen apicalPen12 = QPen(Qt::magenta);;
            
            for(it_cell = integrator->T->ListCell.begin();it_cell!=integrator->T->ListCell.end();it_cell++)
            {
                Cell* c = &(*it_cell);
                QPainterPath ApicalLines11, ApicalLines22, ApicalLines12;
                
                // set the fill colour of the basal sides, depending on the kind of the cells
                QColor lightGreen = Qt::green;
                lightGreen = lightGreen.light();
                
                QColor lightRed = Qt::red;
                lightRed = lightRed.light();
                
                if(c->Type==1) painter.setBrush(Qt::lightGray);
                else if(c->Type==2) painter.setBrush(Qt::lightGray);
                //else if(c->Type==2) painter.setBrush(lightGreen);
                else painter.setBrush(lightRed);
                
                
                // the following way of drawing allows for the inclusion of the lines to the BC_a of the cells
                // if the cell doesnt cross the boundary draw it simply once
                if(!c->crossBoundary)
                {
                    painter.setPen(bluePen_thin);
                    
                    for(std::list<Edge*>::iterator it_edge = c->ListEdges.begin();it_edge != c->ListEdges.end();it_edge++)
                    {		
                        Edge* e = *it_edge;
                        if(e->c1->Type == 1 && e->c2->Type==1) { // draw 11 - edges
                            ApicalLines11.moveTo(e->getVertex1()->getCoord().x, e->getVertex1()->getCoord().y);
                            ApicalLines11.lineTo(e->getVertex2()->getCoord().x, e->getVertex2()->getCoord().y);
                        } else if(e->c1->Type == 2 && e->c2->Type==2) { // draw 12-edges
                                ApicalLines22.moveTo(e->getVertex1()->getCoord().x, e->getVertex1()->getCoord().y);
                                ApicalLines22.lineTo(e->getVertex2()->getCoord().x, e->getVertex2()->getCoord().y);
                        } else { // draw edges between cells of types 1 and 2
                            ApicalLines12.moveTo(e->getVertex1()->getCoord().x, e->getVertex1()->getCoord().y);
                            ApicalLines12.lineTo(e->getVertex2()->getCoord().x, e->getVertex2()->getCoord().y);
                        }
                        
                        // could draw filled polygons on the apical side as well
                        QPolygon polygon;
                        polygon << QPoint(e->getVertex1()->getCoord().x, e->getVertex1()->getCoord().y) << QPoint(e->getVertex2()->getCoord().x, e->getVertex2()->getCoord().y) << QPoint(c->BC_a.x,c->BC_a.y);
                                                
                        painter.drawPolygon(polygon);
                    }
                    
                    painter.setPen(pen_a11); painter.drawPath(ApicalLines11);
                    painter.setPen(pen_a12); painter.drawPath(ApicalLines12);
                    painter.setPen(pen_a22); painter.drawPath(ApicalLines22);
                    
                }
                
                // if the cell crosses the boundary draw the phantom cell around the point where the actual BC_a of the cell is
                else if(c->crossBoundary) {
                    //painter.setPen(bluePen_thin);
                    painter.setPen(bluePen_thin);
                    
                    for(std::list<Edge*>::iterator it_edge = c->ListEdges.begin();it_edge != c->ListEdges.end();it_edge++)
                    {						
                        Edge* e = *it_edge;
                        
                        if(e->c1 == c)
                        {
                            // position the vertices such that they are all connected to the BC_a of the cell
                            Point v1a_rel = shift( e->v1->coord, SystemSize * e->q_c1 * (-1));
                            Point v2a_rel = shift( e->v2->coord, SystemSize * (e->q_v2 - e->q_c1));
                            
                            if(e->c1->Type == 1 && e->c2->Type==1) { // draw normal edges
                                ApicalLines11.moveTo(v1a_rel.x, v1a_rel.y);
                                ApicalLines11.lineTo(v2a_rel.x, v2a_rel.y);
                            } else if(e->c1->Type == 2 && e->c2->Type==2) { // draw 12-edges
                                ApicalLines22.moveTo(v1a_rel.x, v1a_rel.y);
                                ApicalLines22.lineTo(v2a_rel.x, v2a_rel.y);
                            } else { // draw edges between cells of types 1 and 2
                                ApicalLines12.moveTo(v1a_rel.x, v1a_rel.y);
                                ApicalLines12.lineTo(v2a_rel.x, v2a_rel.y);
                            }
                            
                            // could draw filled polygons on the apical side as well
                            QPolygon polygon;
                            polygon << QPoint(v1a_rel.x, v1a_rel.y) << QPoint(v2a_rel.x, v2a_rel.y) << QPoint(c->BC_a.x,c->BC_a.y);
                            painter.drawPolygon(polygon);
                        }
                        
                        else
                        {
                            Point v1a_rel = shift( e->v1->coord, SystemSize * (e->q_c2 + e->q_v2) * (-1));
                            Point v2a_rel = shift( e->v2->coord, SystemSize * e->q_c2 * (-1));
                            
                            if(e->c1->Type == 1 && e->c2->Type==1) { // draw normal edges
                                ApicalLines11.moveTo(v1a_rel.x, v1a_rel.y);
                                ApicalLines11.lineTo(v2a_rel.x, v2a_rel.y);
                            } else if(e->c1->Type == 2 && e->c2->Type==2) { // draw 12-edges
                                ApicalLines22.moveTo(v1a_rel.x, v1a_rel.y);
                                ApicalLines22.lineTo(v2a_rel.x, v2a_rel.y);
                            } else { // draw edges between cells of types 1 and 2
                                ApicalLines12.moveTo(v1a_rel.x, v1a_rel.y);
                                ApicalLines12.lineTo(v2a_rel.x, v2a_rel.y);
                            }
                            
                            // could draw filled polygons on the apical side as well
                            QPolygon polygon;
                            polygon << QPoint(v1a_rel.x, v1a_rel.y) << QPoint(v2a_rel.x, v2a_rel.y) << QPoint(c->BC_a.x,c->BC_a.y);
                            painter.drawPolygon(polygon);
                            
                        }
                        
                    }
                    
                    painter.setPen(pen_a11); painter.drawPath(ApicalLines11);
                    painter.setPen(pen_a12); painter.drawPath(ApicalLines12);
                    painter.setPen(pen_a22); painter.drawPath(ApicalLines22);
                }
            }
        }
        
        // draw the box around the tissue
        QPainterPath boxPath;
        boxPath.moveTo(0,0);
        boxPath.lineTo(SystemSize.x,0);
        boxPath.lineTo(SystemSize.x,SystemSize.y);
        boxPath.lineTo(0,SystemSize.y);
        boxPath.lineTo(0,0);
        
        //painter.setPen(darkgreenPen_thick); painter.drawPath(boxPath);
        
    }
}

		
		
		
	//std::list<Edge>::iterator it_edge;
//	
//	QPainterPath CellPath;
//		
//	for(it_edge = integrator->T->ListEdge.begin();it_edge!=integrator->T->ListEdge.end();it_edge++)
//	{
//		Edge* e = &(*it_edge);
//		
//		// if the edge does not cross the boundary draw it normally
//		if(!e->crossBoundary)
//		{
//			double v1a_x = e->getVertex1()->getCoord().x;
//			double v1a_y = e->getVertex1()->getCoord().y;
//			double v2a_x = e->getVertex2()->getCoord().x;
//			double v2a_y = e->getVertex2()->getCoord().y;
//			
//			double v1b_x = e->getVertex1()->getBasalCoord().x;
//			double v1b_y = e->getVertex1()->getBasalCoord().y;
//			double v2b_x = e->getVertex2()->getBasalCoord().x;
//			double v2b_y = e->getVertex2()->getBasalCoord().y;
//			
//			// first draw the apical edges:
//			painter.setPen(pen_a11); //set the pen
//			CellPath = QPainterPath(); //remove all the paths from CellPath
//			
//			CellPath.moveTo(v1a_x, v1a_y);
//			CellPath.lineTo(v2a_x, v2a_y);
//			painter.drawPath(CellPath);
//			
//			// then draw the basal edges:
//			painter.setPen(pen_b11); //set the pen
//			CellPath = QPainterPath(); //remove all the paths from CellPath
//			
//			CellPath.moveTo(v1b_x, v1b_y);
//			CellPath.lineTo(v2b_x, v2b_y);
//			painter.drawPath(CellPath);
//			
//			
//		}
//		// if it crosses the boundary, then draw it twice	 
//		else
//		{	
//			double v1a_1_x = e->getVertex1()->getCoord().x;
//			double v1a_1_y = e->getVertex1()->getCoord().y;
//			double v2a_1_x = e->getVertex2()->getCoord().x + e->q_v2.x * integrator->T->SystemSize.x ;
//			double v2a_1_y = e->getVertex2()->getCoord().y + e->q_v2.y * integrator->T->SystemSize.y ;		
//			
//			double v1b_1_x = e->getVertex1()->getCoord().x + e->v1->q_b.x * integrator->T->SystemSize.x;
//			double v1b_1_y = e->getVertex1()->getCoord().y + e->v1->q_b.y * integrator->T->SystemSize.y;
//			double v2b_1_x = e->getVertex2()->getCoord().x + (e->v2->q_b.x + e->q_v2.x) * integrator->T->SystemSize.x;
//			double v2b_1_y = e->getVertex2()->getCoord().y + (e->v2->q_b.y + e->q_v2.y) * integrator->T->SystemSize.y;;		
//		
//			double v2a_2_x = e->getVertex2()->getCoord().x;
//			double v2a_2_y = e->getVertex2()->getCoord().y;
//			double v1a_2_x = e->getVertex1()->getCoord().x - e->q_v2.x * integrator->T->SystemSize.x ;
//			double v1a_2_y = e->getVertex1()->getCoord().y - e->q_v2.y * integrator->T->SystemSize.y ;		
//			
//			double v2b_2_x = e->getVertex2()->getCoord().x + e->v2->q_b.x * integrator->T->SystemSize.x;
//			double v2b_2_y = e->getVertex2()->getCoord().y + e->v2->q_b.y * integrator->T->SystemSize.y;
//			double v1b_2_x = e->getVertex1()->getCoord().x + (e->v1->q_b.x - e->q_v2.x) * integrator->T->SystemSize.x;
//			double v1b_2_y = e->getVertex1()->getCoord().y + (e->v1->q_b.y - e->q_v2.y) * integrator->T->SystemSize.y;;		
//			
//			
//			// now draw the apical lines:
//			
//			painter.setPen(pen_a11); //set the pen
//			CellPath = QPainterPath(); //remove all the paths from CellPath
//			
//			CellPath.moveTo(v1a_1_x, v1a_1_y);
//			CellPath.lineTo(v2a_1_x, v2a_1_y);
//			CellPath.moveTo(v1a_2_x, v1a_2_y);
//			CellPath.lineTo(v2a_2_x, v2a_2_y);
//			painter.drawPath(CellPath);
//
//			// now draw the apical lines:
//			painter.setPen(pen_b11); //set the pen
//			CellPath = QPainterPath(); //remove all the paths from CellPath
//					
//			CellPath.moveTo(v1b_1_x, v1b_1_y);
//			CellPath.lineTo(v2b_1_x, v2b_1_y);
//			CellPath.moveTo(v1b_2_x, v1b_2_y);
//			CellPath.lineTo(v2b_2_x, v2b_2_y);
//			painter.drawPath(CellPath);
//			
//		
//		}
//	}


void TissueWindow::mousePressEvent(QMouseEvent *event)
{
	// if no tissue has be defined dont do anything
	if(integrator->T->ListCell.empty())	return;
	
	bool invertible;
	QMatrix inverse=matrix.inverted(&invertible);
	int tx;
	int ty;
	inverse.map(event->x(),event->y(),&tx,&ty);
	
	// signal_currentPointChanged(tx, ty);
	
	
    // if left mouse button has been pressed
	if (event->button() == Qt::LeftButton)
    {
        if(!event->modifiers().testFlag(Qt::ShiftModifier)) // if shift also has been pressed
        {
            signal_currentPoint(tx,ty);
        }
        else
        {
            signal_addCurrentPoint(tx,ty);
        }
        
    }
    
	
}


void TissueWindow::keyPressEvent(QKeyEvent *event)
{
	if(event->key()==Qt::Key_Equal) {matrix.scale(1.25,1.25);scaleFactor*=1.25;}
	else if(event->key()==Qt::Key_Minus) {matrix.scale(0.8,0.8);scaleFactor*=0.8;}
	else if(event->key()==Qt::Key_Right) {matrix.translate(-height()/scaleFactor/10,0);}
	else if(event->key()==Qt::Key_Left) {matrix.translate(height()/scaleFactor/10,0);}
	else if(event->key()==Qt::Key_Up) {matrix.translate(0,-height()/scaleFactor/10);}
	else if(event->key()==Qt::Key_Down) {matrix.translate(0,height()/scaleFactor/10);}
	else if(event->key()==Qt::Key_Q) {	matrix.rotate(5);}
	else if(event->key()==Qt::Key_Z) {	matrix.rotate(-5);}
	else if(event->key()==Qt::Key_R) {	matrix.setMatrix(1,0,0,-1,0,rect().height());matrix.scale(1,1);scaleFactor=1;}
    else if(event->key()==Qt::Key_C) {	signal_printCellData();}
    else if(event->key()==Qt::Key_D) {  signal_addToCurrentPoint(0.02*distance.x,0);  }
    else if(event->key()==Qt::Key_A) {  signal_addToCurrentPoint(-0.02*distance.x,0); }
    else if(event->key()==Qt::Key_W) {  signal_addToCurrentPoint(0,0.02*distance.y);  }
    else if(event->key()==Qt::Key_S) {  signal_addToCurrentPoint(0,-0.02*distance.y); }
    else if(event->key()==Qt::Key_B) {  drawBasalSide = !drawBasalSide;}
	
    update();
}


void TissueWindow::slot_drawingModeChanged(int drawBarycenters_, int drawBasalSide_, int drawApicalSide_)
{
	drawBarycenters=drawBarycenters_;
	drawBasalSide = drawBasalSide_;
    drawApicalSide = drawApicalSide_;
    update();
	repaint();
	QApplication::processEvents();
	
}

//void TissueWindow::RecordImageFile()
//{
//	//QLOG_INFO()<<"saving image in folder"<<CurrentRecordPath<<"\n";
//	std::string s;
//	std::stringstream out;
//	//int numero=(tissueWindow->nbIterationsDone)/(tissueWindow->frameRate);
//	int number=nbImagesRecorded;
//	
//	//make sure to have zeros in the beginning to simplify the ordering
//	if(number<10) out<<0<<0<<number;
//	else
//    {
//		if(number<100) out<<0<<number;	else out<<number;
//	}
//	
//	s = out.str();
//    QString fileName = CurrentRecordPath+CurrentRecordPath.fromStdString(s)+".png";
//	CurrentImage.fill(QColor(245,208,123).rgb());//remet l'image a 0
//	writeCurrentImage();// save the current image
//    CurrentImage.save(fileName,"PNG",-1);
//}
//
//void TissueWindow::writeCurrentImage()
//{
//        
//	//double recordfactor=4.0;
//	double recordfactor=2.0;
//	if(T->xPeriodic&&T->yPeriodic)
//	{
//        CurrentImage=QImage(recordfactor*((int)T->SystemSize.x+1),recordfactor*((int)T->SystemSize.y+1),QImage::Format_RGB32);
//        CurrentImage.fill(QColor(245,208,123).rgb());//remet l'image a 0
//	}
//	QPainter painter(&CurrentImage);
//    
//	QPen pen(Qt::black, 2, Qt::SolidLine);
//    painter.setPen(pen);
//    painter.setBrush(Qt::white);
//	painter.setRenderHint(QPainter::Antialiasing,true);//proprietes du pinceau
//	if(T->xPeriodic&&!T->yPeriodic) painter.setClipRect(QRectF(0-1,-10000,T->Lx,20000), Qt::ReplaceClip);
//	if(T->yPeriodic&&!T->xPeriodic) painter.setClipRect(QRectF(-10000,0-1,20000,T->Ly), Qt::ReplaceClip);
//	if(T->xPeriodic&&T->yPeriodic)
//	{
//		QMatrix matrix2;
//		matrix2.setMatrix(recordfactor,0,0,-recordfactor,0,recordfactor*(T->Ly));//matrice qui transforme les coordonnees pour les placer en coordonnees mathematiques, commencent en bas a gauche.
//		//matrix2.scale(1,1);
//		//matrix2.scale(2,2);
//		painter.setWorldMatrix(matrix2,false);
//		//painter.setClipRect(QRectF(-1,-1,T->Lx+1,T->Ly+1),Qt::ReplaceClip);
//		painter.setClipRect(QRectF(-1,-1,recordfactor*(T->Lx+1),recordfactor*(T->Ly+1)),Qt::ReplaceClip);
//	}
//    
//	//painter.setClipRect(QRectF(0-1,0-1,10,10), Qt::ReplaceClip);
//    
//	std::list<Edge>::iterator it;//it itere sur la liste de pointeurs listant les edges d'un vertex
//	std::list<Cell>::iterator itCell;
//	for(itCell=T->ListCell.begin();itCell!=T->ListCell.end();itCell++)
//	{
//		drawCell(painter,&(*itCell));
//		painter.setPen(pen);
//        
//	}
//	//pour avoir une image 2 fois plus grande
//	//QImage ImageTemp=(CurrentImage->scaled(QSize(2*((int)T->Lx+2),2*((int)T->Ly+2)),Qt::KeepAspectRatio,Qt::SmoothTransformation));
//	//CurrentImage=&ImageTemp;
//}


#endif // _USE_QT_
