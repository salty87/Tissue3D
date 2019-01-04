#ifdef _USE_QT_


/****************************************************************************
 **
 ** Copyright (C) 2012 Digia Plc and/or its subsidiary(-ies).
 ** Contact: http://www.qt-project.org/legal
 **
 ** This file is part of the examples of the Qt Toolkit.
 **
 ** $QT_BEGIN_LICENSE:BSD$
 ** You may use this file under the terms of the BSD license as follows:
 **
 ** "Redistribution and use in source and binary forms, with or without
 ** modification, are permitted provided that the following conditions are
 ** met:
 **   * Redistributions of source code must retain the above copyright
 **     notice, this list of conditions and the following disclaimer.
 **   * Redistributions in binary form must reproduce the above copyright
 **     notice, this list of conditions and the following disclaimer in
 **     the documentation and/or other materials provided with the
 **     distribution.
 **   * Neither the name of Digia Plc and its Subsidiary(-ies) nor the names
 **     of its contributors may be used to endorse or promote products derived
 **     from this software without specific prior written permission.
 **
 **
 ** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 ** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 ** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 ** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 ** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 ** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 ** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 ** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 ** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 ** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 ** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
 **
 ** $QT_END_LICENSE$
 **
 ****************************************************************************/

#include <QGLWidget>
#include <QMatrix4x4>
#include <QVector3D>

#include <iostream>
#include <qmath.h>
#include "Tissue.h"

#include "GLTissueImage.h"

struct Geometry
{
    QVector<unsigned int> faces;
    QVector<unsigned int> lines;
    QVector<QVector3D> vertices;
    QVector<QVector3D> normals;
    void appendVertex(const QVector3D &a);
    void finalize();
    void preloadArrays() const;
};

// first step of drawing routine -> initialize the GLarrays
void Geometry::preloadArrays() const
{
    glVertexPointer(3, GL_FLOAT, 0, vertices.constData());
    glNormalPointer(GL_FLOAT, 0, normals.constData());
}

void Geometry::appendVertex(const QVector3D &a)
{
    int v = vertices.count();
    vertices.append(a);
    faces.append(v);
}

// either an apical, basal or lateral surface of a cell
class SurfacePatch
{
public:
    SurfacePatch(Geometry *);
    void draw() const;
    void addTri(const Point &pa, const Point &pb, const Point &c);
    /// add several triangles
    void addTri(const Point &pa, const Point &pb, const Point &pc, const int &xr, const int &yr, const double &Lx, const double &Ly);
    
    unsigned int start;
    unsigned int count;
    unsigned int initv;
    
    GLfloat faceColor[4]; //vector of the face colours (later set both ambient and diffuse)
    Geometry *geom;
};

static inline void qSetColor(float colorVec[], QColor c, float alpha)
{
    colorVec[0] = c.redF();
    colorVec[1] = c.greenF();
    colorVec[2] = c.blueF();
    colorVec[3] = alpha;
}

//static inline void qSetColor(float colorVec[], QColor c)
//{
//    colorVec[0] = c.redF();
//    colorVec[1] = c.greenF();
//    colorVec[2] = c.blueF();
//    colorVec[3] = c.alphaF();
//}


SurfacePatch::SurfacePatch(Geometry *g)
: start(g->faces.count()), count(0), initv(g->vertices.count()), geom(g)
{
    qSetColor(faceColor, QColor(Qt::white),0.0);
}

//! [2]
void SurfacePatch::draw() const
{
    //glDisable(GL_CULL_FACE);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, faceColor);
    //glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, faceColor);
    
    // enable for transparency
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, faceColor);
        
    const unsigned int *indices = geom->faces.constData();
    glDrawElements(GL_TRIANGLES, count, GL_UNSIGNED_INT, indices + start);
    
}
//! [2]

// add triangle to patch
void SurfacePatch::addTri(const Point &pa, const Point &pb, const Point &pc)
{
    QVector3D a((float)pa.x,(float)pa.y,(float)pa.z);
    QVector3D b((float)pb.x,(float)pb.y,(float)pb.z);
    QVector3D c((float)pc.x,(float)pc.y,(float)pc.z);
    
    geom->appendVertex(a);
    geom->appendVertex(b);
    geom->appendVertex(c);
    
    count += 3;
}

// add triangles into different boxes
void SurfacePatch::addTri(const Point &pa, const Point &pb, const Point &pc, const int &xr, const int &yr, const double &Lx, const double &Ly)
{
    for(int i=0;i<xr; i++) {
        for(int j=0;j<yr; j++) {
    
            QVector3D a((float)pa.x+i*Lx,(float)pa.y+j*Ly,(float)pa.z);
            QVector3D b((float)pb.x+i*Lx,(float)pb.y+j*Ly,(float)pb.z);
            QVector3D c((float)pc.x+i*Lx,(float)pc.y+j*Ly,(float)pc.z);
            
            geom->appendVertex(a);
            geom->appendVertex(b);
            geom->appendVertex(c);
            
            count += 3;
        }
    }
}

class LinePatch
{
public: 
    LinePatch(Geometry *);
    void draw() const;
    void addLine(const Point &pa, const Point &pb);
    void addLine(const Point &pa, const Point &pb, const int &xr, const int &yr, const double &Lx, const double &Ly);
    
    unsigned int start;
    unsigned int count;
    unsigned int initv;
    
    GLfloat faceColor[4]; //vector of the colours of the line
    Geometry *geom;
};

LinePatch::LinePatch(Geometry *g)
: start(g->faces.count()), count(0), initv(g->vertices.count()), geom(g)
{
    qSetColor(faceColor, QColor(Qt::red),0.0);
}

//! [2]
void LinePatch::draw() const
{
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, faceColor);
    //glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, faceColor);
    
    const unsigned int *indices = geom->faces.constData();
    glDrawElements(GL_LINES, count, GL_UNSIGNED_INT, indices + start);
    
}
//! [2]

// add line to patch
void LinePatch::addLine(const Point &pa, const Point &pb)
{
    QVector3D a((float)pa.x,(float)pa.y,(float)pa.z);
    QVector3D b((float)pb.x,(float)pb.y,(float)pb.z);
    
    geom->appendVertex(a);
    geom->appendVertex(b);
    
    count += 2;
}

// add line to patch
void LinePatch::addLine(const Point &pa, const Point &pb, const int &xr, const int &yr, const double &Lx, const double &Ly)
{
    
    for(int i=0;i<xr; i++) {
        for(int j=0;j<yr; j++) {
            
            QVector3D a((float)pa.x+i*Lx,(float)pa.y+j*Ly,(float)pa.z);
            QVector3D b((float)pb.x+i*Lx,(float)pb.y+j*Ly,(float)pb.z);
            
            geom->appendVertex(a);
            geom->appendVertex(b);
            
            count += 2;
        }
    }

}

GLTissueImage::GLTissueImage(QWidget *parent, Integrator *integrator)
: QObject(parent), integrator(integrator), tissueChanged(true)
, geom(new Geometry())
{
    x_reps = 2;
    y_reps = 2;
}

GLTissueImage::~GLTissueImage()
{
    qDeleteAll(AParts);
    qDeleteAll(AParts_c);
    qDeleteAll(BParts);
    qDeleteAll(LParts);
    qDeleteAll(ALineParts);
    qDeleteAll(BLineParts);
    qDeleteAll(LLineParts);
    qDeleteAll(LLineParts_cyst);
    delete geom;
}

//void GLTissueImage::loadParts()
//{
//    
//    QColor dummy;
//    //QColor apicalPatchColors [] = { QColor(Qt::gray), QColor(Qt::magenta), QColor(Qt::red), QColor(Qt::darkRed), QColor(Qt::darkYellow) };
//    
//    QColor apicalType1; apicalType1.setNamedColor("black");
//    QColor apicalType2(Qt::magenta);
//    QColor basalType1;  basalType1.setNamedColor("black");
//    QColor basalType2;  basalType2.setNamedColor("black");
//    
//    QColor apicalPatchColors [] = { QColor(Qt::gray), QColor(Qt::red), QColor(Qt::green), QColor(Qt::cyan), QColor(Qt::darkCyan) };//{ QColor(Qt::gray), QColor(Qt::red), QColor(Qt::darkGreen), QColor(Qt::blue), QColor(Qt::darkGray) };
//    double apicalPatchOpacities [] = { 0.0, 0.8, 0.8, 0.3, 0.5 };//{ 0.0, 0.1, 0.5, 0.3, 0.5 };
//    
//    QColor lightGreen = Qt::green;
//    lightGreen = lightGreen.light();
//    
//    QColor basalPatchColors  [] = { QColor(Qt::gray), QColor(Qt::lightGray), QColor(lightGreen), QColor(Qt::cyan), QColor(Qt::darkCyan) };
//    double basalPatchOpacities [] = { 0.0, 0.05, 0.5, 0.3, 0.5 };
//
//    QColor lateralPatchColor(Qt::darkGray);
//    QColor lateralCystColor(Qt::blue);
//    
//    
//    QColor apicalEdgeColor_sameCellType(Qt::black);
//    QColor apicalEdgeColor_differentCellType(Qt::magenta);
//    QColor basalEdgeColor_sameCellType(Qt::black);
//    QColor basalEdgeColor_differentCellType(Qt::magenta);
//    
//    int numberOfCellTypes = 4;
//    
//    // clear the list of parts
//    removeParts();
//    Tissue* T = integrator->T;
//    
//    // update the systems outermost points
//    T->UpdateExtPositions();
//    
//    // get the biggest size in the system
//    double Lz = T->extPosition_topRight.z-T->extPosition_bottomLeft.z;
//    double L = std::max(T->SystemSize.x, std::max(T->SystemSize.y,Lz));;
//
//    // draw the box around the tissue
//    bool drawBox = false;
//    
//    if(T->isPeriodic && drawBox)
//    {
//        
//        // system size and the box size, that gives the size of the box which is drawn around the tissue at relative distance exceed
//        double exceed = 0.1;
//        Point boxSize = Point(T->SystemSize.x,T->SystemSize.y,Lz*(1+2*exceed))/L;
//        Point box_origin = Point(0,0,T->extPosition_bottomLeft.z/L-Lz*(exceed/L));
//        
//        // draw the box outlines
//        LinePatch *box_lines = new LinePatch(geom);
//        
//        box_lines->addLine(box_origin,box_origin+Point(boxSize.x,0,0));
//        box_lines->addLine(box_origin,box_origin+Point(0,boxSize.y,0));
//        box_lines->addLine(box_origin,box_origin+Point(0,0,boxSize.z));
//        
//        box_lines->addLine(box_origin+Point(boxSize.x,0,0),box_origin+Point(boxSize.x,boxSize.y,0));
//        box_lines->addLine(box_origin+Point(boxSize.x,0,0),box_origin+Point(boxSize.x,0,boxSize.z));
//        
//        box_lines->addLine(box_origin+Point(0,boxSize.y,0),box_origin+Point(boxSize.x,boxSize.y,0));
//        box_lines->addLine(box_origin+Point(0,boxSize.y,0),box_origin+Point(0,boxSize.y,boxSize.z));
//        
//        box_lines->addLine(box_origin+Point(0,0,boxSize.z),box_origin+Point(0,boxSize.y,boxSize.z));
//        box_lines->addLine(box_origin+Point(0,0,boxSize.z),box_origin+Point(boxSize.x,0,boxSize.z));
//        
//        box_lines->addLine(box_origin+Point(boxSize.x,boxSize.y,0),box_origin+Point(boxSize.x,boxSize.y,boxSize.z));
//        box_lines->addLine(box_origin+Point(boxSize.x,0,boxSize.z),box_origin+Point(boxSize.x,boxSize.y,boxSize.z));
//        box_lines->addLine(box_origin+Point(0,boxSize.y,boxSize.z),box_origin+Point(boxSize.x,boxSize.y,boxSize.z));
//        
//        qSetColor(box_lines->faceColor, QColor(Qt::black),0.5);
//        BoxLineParts<<box_lines;
//        
//    }
//    
//    for(std::list<Cell>::iterator it_c=T->ListCell.begin(); it_c!=T->ListCell.end(); ++it_c)
//    {
//        if(!(*it_c).CrossBoundary())
//        {
//            
//            // normalize the variables
//            Cell* c = &(*it_c);
//            Point BC_a = c->BC_a;
//            Point BC_a_norm = BC_a/L;
//            
//            Point BC_b_norm = c->BC_b/L;
//            
//            // define surface patches for apical and basal surfaces
//            SurfacePatch *APatch = new SurfacePatch(geom);
//            
//            // add the apical side
//            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
//            {
//                Edge *e = (*it_e);
//                
//                //std::cout << "(" << BC_a.x/T->SystemSize.x << "," << BC_a.y/T->SystemSize.y << "," << BC_a.z << ")" << std::endl;
//                
//                // calculate the normalized coordinates
//                Point v1_a_norm = e->v1->coord/L;
//                Point v2_a_norm = e->v2->coord/L;
//                
//                // add the triangles (normal is defined by the direction towards the cell)
////                if(e->c2==c)
////                {
////                    APatch->addTri(v1_a_norm,BC_a_norm,v2_a_norm);
////                }
////                else
////                {
////                    APatch->addTri(BC_a_norm,v1_a_norm,v2_a_norm);
////                }
////                
//                APatch->addTri(v1_a_norm,BC_a_norm,v2_a_norm);
//                APatch->addTri(BC_a_norm,v1_a_norm,v2_a_norm);
//
//                
//            }
//            
//            
//            qSetColor(APatch->faceColor, apicalPatchColors[c->Type],apicalPatchOpacities[c->Type]);
//            AParts << APatch;
//            
//            
//            // add lines between same cell types
//            LinePatch *ALinePatch = new LinePatch(geom);
//            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
//            {
//               if((*it_e)->c1->Type == (*it_e)->c2->Type) ALinePatch->addLine((*it_e)->v1->coord/L,(*it_e)->v2->coord/L);
//            }
//            qSetColor(ALinePatch->faceColor, apicalEdgeColor_sameCellType,1.0);
//            ALineParts << ALinePatch;
//            
//            // add lines between different cell types
//            LinePatch *ALinePatch_d = new LinePatch(geom);
//            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
//            {
//                if((*it_e)->c1->Type != (*it_e)->c2->Type) ALinePatch_d->addLine((*it_e)->v1->coord/L,(*it_e)->v2->coord/L);
//            }
//            qSetColor(ALinePatch_d->faceColor, apicalEdgeColor_differentCellType,1.0);
//            ALineParts << ALinePatch_d;
//            
//            // add basal surface patches
//            SurfacePatch *BPatch = new SurfacePatch(geom);
//            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
//            {
//                Edge *e = (*it_e);
//                
//                // calculate the normalized coordinates
//                Point v1_b_norm = e->v1->coord_b/L;
//                Point v2_b_norm = e->v2->coord_b/L;
//
//                BPatch->addTri(BC_b_norm,v1_b_norm,v2_b_norm);
//                BPatch->addTri(v1_b_norm,BC_b_norm,v2_b_norm);
//            }
//            qSetColor(BPatch->faceColor, basalPatchColors[c->Type],basalPatchOpacities[c->Type]);
//            BParts << BPatch;
//            
//            // add basal lines
//            LinePatch *BLinePatch = new LinePatch(geom);
//            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
//            {
//                BLinePatch->addLine((*it_e)->v1->coord_b/L,(*it_e)->v2->coord_b/L);
//            }
//            qSetColor(ALinePatch->faceColor, basalEdgeColor_sameCellType, 0.0);
//            BLineParts << BLinePatch;
//        }
//        else { // if the cell crosses the boundary
//            
//            Cell* c = &(*it_c);
//            
//            // define surface patches for apical and basal surfaces
//            SurfacePatch *APatch = new SurfacePatch(geom);
//            
//            // add the apical side
//            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
//            {
//                Edge *e = (*it_e);
//                
//                // add the triangles (normal is defined by the direction towards the cell)
//                if(e->c1==c)
//                {
//                    // calculate the normalized coordinates
//                    Point v1_a_norm = shift(e->v1->coord,T->SystemSize*e->q_c1*(-1))/L;
//                    Point v2_a_norm = shift(e->v2->coord,T->SystemSize*(e->q_v2-e->q_c1))/L;
//                    Point BC_a_norm = c->BC_a/L;
//                    APatch->addTri(v1_a_norm,BC_a_norm,v2_a_norm);
//                    APatch->addTri(v1_a_norm,v2_a_norm,BC_a_norm);
//                }
//                else
//                {
//                    Point v1_a_norm = shift(e->v1->coord,T->SystemSize*(e->q_c2+e->q_v2)*(-1))/L;
//                    Point v2_a_norm = shift(e->v2->coord,T->SystemSize*(e->q_c2)*(-1))/L;
//                    Point BC_a_norm = c->BC_a/L;
//                    APatch->addTri(v1_a_norm,BC_a_norm,v2_a_norm);
//                    APatch->addTri(v1_a_norm,v2_a_norm,BC_a_norm);
//                }
//                
//            }
//            
//            qSetColor(APatch->faceColor, apicalPatchColors[c->Type],apicalPatchOpacities[c->Type]);
//
//            AParts << APatch;
//            
//            LinePatch *ALinePatch = new LinePatch(geom);
//            
//            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
//            {
//                Edge *e = (*it_e);
//                
//                if(e->c1==c)
//                {
//                    // calculate the normalized coordinates
//                    Point v1_a_norm = shift(e->v1->coord,T->SystemSize*e->q_c1*(-1))/L;
//                    Point v2_a_norm = shift(e->v2->coord,T->SystemSize*(e->q_v2-e->q_c1))/L;
//                    ALinePatch->addLine(v1_a_norm,v2_a_norm);
//                }
//                else
//                {
//                    Point v1_a_norm = shift(e->v1->coord,T->SystemSize*(e->q_c2+e->q_v2)*(-1))/L;
//                    Point v2_a_norm = shift(e->v2->coord,T->SystemSize*(e->q_c2)*(-1))/L;
//                    ALinePatch->addLine(v1_a_norm,v2_a_norm);
//                }
//            }
//            
//            qSetColor(ALinePatch->faceColor, apicalEdgeColor_sameCellType, 0.0);
//            ALineParts << ALinePatch;
//            
//            
//            SurfacePatch *BPatch = new SurfacePatch(geom);
//            // add the basal side
//            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
//            {
//                Edge *e = (*it_e);
//                
//                if(e->c1==c)
//                {
//                    // calculate the normalized coordinates
//                    Point v1_b_norm = shift(e->v1->coord_b,T->SystemSize*(e->v1->q_b-e->q_c1))/L;
//                    Point v2_b_norm = shift(e->v2->coord_b,T->SystemSize*(e->v2->q_b+e->q_v2-e->q_c1))/L;
//                    Point BC_b_norm = shift(c->BC_b,T->SystemSize*c->q_BC_b)/L;
//                    BPatch->addTri(v1_b_norm,BC_b_norm,v2_b_norm);
//                    BPatch->addTri(v1_b_norm,v2_b_norm,BC_b_norm);
//                }
//                else
//                {
//                    Point v1_b_norm = shift(e->v1->coord_b,T->SystemSize*(e->v1->q_b-e->q_c2-e->q_v2))/L;
//                    Point v2_b_norm = shift(e->v2->coord_b,T->SystemSize*(e->v2->q_b-e->q_c2))/L;
//                    Point BC_b_norm = shift(c->BC_b,T->SystemSize*c->q_BC_b)/L;
//                    BPatch->addTri(v1_b_norm,BC_b_norm,v2_b_norm);
//                    BPatch->addTri(v1_b_norm,v2_b_norm,BC_b_norm);
//                }
//                
//            }
//            
//            
//            qSetColor(BPatch->faceColor, basalPatchColors[c->Type],basalPatchOpacities[c->Type]);
//            BParts << BPatch;
//            
//            LinePatch *BLinePatch = new LinePatch(geom);
//            
//            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
//            {
//                Edge *e = (*it_e);
//                
//                if(e->c1==c)
//                {
//                    // calculate the normalized coordinates
//                    Point v1_b_norm = shift(e->v1->coord_b,T->SystemSize*(e->v1->q_b-e->q_c1))/L;
//                    Point v2_b_norm = shift(e->v2->coord_b,T->SystemSize*(e->q_v2+e->v2->q_b-e->q_c1))/L;
//                    BLinePatch->addLine(v1_b_norm,v2_b_norm);
//                }
//                else
//                {
//                    Point v1_b_norm = shift(e->v1->coord_b,T->SystemSize*(e->v1->q_b-e->q_c2-e->q_v2))/L;
//                    Point v2_b_norm = shift(e->v2->coord_b,T->SystemSize*(e->v2->q_b-e->q_c2))/L;
//                    BLinePatch->addLine(v1_b_norm,v2_b_norm);
//                }
//            }
//            
//            qSetColor(BLinePatch->faceColor, apicalEdgeColor_sameCellType,0.0);
//            BLineParts << BLinePatch;
//            
//        }
//        
//    }
//    
//    // iterate through the edges and draw all, that are at the boundary i.e. can be seen OR if they are at mut-wt interface
//    for(std::list<Edge>::iterator it_e = T->ListEdge.begin(); it_e!=T->ListEdge.end(); ++it_e)
//    {
//        
//        if(it_e->isAtBoundary() || it_e->c1->Type != it_e->c2->Type)
//        {
//            if(!T->isPeriodic)
//            {
//                // normalize the coordinates
//                Point v1_a_norm = it_e->v1->coord/L;
//                Point v2_a_norm = it_e->v2->coord/L;
//                Point v1_b_norm = it_e->v1->coord_b/L;
//                Point v2_b_norm = it_e->v2->coord_b/L;
//                Point b_norm = it_e->BC_l/L;
//                
//                SurfacePatch *LPatch = new SurfacePatch(geom);
//                
//                // add the four triangles to the patch
//                LPatch->addTri(v1_a_norm,v2_a_norm,b_norm);
//                LPatch->addTri(v2_a_norm,v2_b_norm,b_norm);
//                LPatch->addTri(v2_b_norm,v1_b_norm,b_norm);
//                LPatch->addTri(v1_b_norm,v1_a_norm,b_norm);
//                
//                qSetColor(LPatch->faceColor, QColor(Qt::black),0.0);
//                LParts << LPatch;
//                
//                // draw the vertex connection
//                LinePatch *LLinePatch = new LinePatch(geom);
//                LLinePatch->addLine(v1_a_norm,v1_b_norm);
//                LLinePatch->addLine(v2_a_norm,v2_b_norm);
//                
//                qSetColor(LLinePatch->faceColor, QColor(Qt::black),0.0);
//                LLineParts << LLinePatch;
//                
//            }
//            else {
//                // draw the edge on the side of cell 1
//                // normalize the coordinates
//                Point v1_a_norm = shift(it_e->v1->coord, T->SystemSize*it_e->q_c1*(-1))/L;
//                Point v2_a_norm = shift(it_e->v2->coord, T->SystemSize*(it_e->q_v2-it_e->q_c1))/L;
//                Point v1_b_norm = shift(it_e->v1->coord_b, T->SystemSize*(it_e->v1->q_b - it_e->q_c1))/L;
//                Point v2_b_norm = shift(it_e->v2->coord_b, T->SystemSize*(it_e->q_v2 + it_e->v2->q_b - it_e->q_c1))/L;
//                Point b_norm = shift(it_e->BC_l, T->SystemSize*(it_e->q_BC_l - it_e->q_c1))/L;
//                
//                // add the four triangles to the patch
//                SurfacePatch *LPatch_c1 = new SurfacePatch(geom);
//                LPatch_c1->addTri(v1_a_norm,v2_a_norm,b_norm);
//                LPatch_c1->addTri(v1_a_norm,b_norm,v2_a_norm);
//                LPatch_c1->addTri(v2_a_norm,v2_b_norm,b_norm);
//                LPatch_c1->addTri(v2_a_norm,b_norm,v2_b_norm);
//                LPatch_c1->addTri(v2_b_norm,v1_b_norm,b_norm);
//                LPatch_c1->addTri(v2_b_norm,b_norm,v1_b_norm);
//                LPatch_c1->addTri(v1_b_norm,v1_a_norm,b_norm);
//                LPatch_c1->addTri(v1_b_norm,b_norm,v1_a_norm);
//                
//                
//                if(it_e->isAtBoundary()) qSetColor(LPatch_c1->faceColor, lateralPatchColor, 0.5);
//                else  qSetColor(LPatch_c1->faceColor, lateralCystColor, 0.5);
//                
//                LParts << LPatch_c1;
//                
//                // draw the vertex connection
//                LinePatch *LLinePatch_c1 = new LinePatch(geom);
//                LLinePatch_c1->addLine(v1_a_norm,v1_b_norm);
//                LLinePatch_c1->addLine(v2_a_norm,v2_b_norm);
//                
//                qSetColor(LLinePatch_c1->faceColor, QColor(Qt::black), 0.0);
//                LLineParts << LLinePatch_c1;
//                
//                // draw the edge on the side of cell 2
//                // normalize the coordinates
//                v2_a_norm = shift(it_e->v2->coord, T->SystemSize*it_e->q_c2*(-1))/L;
//                v1_a_norm = shift(it_e->v1->coord, T->SystemSize*(it_e->q_v2+it_e->q_c2)*(-1))/L;
//                v2_b_norm = shift(it_e->v2->coord_b, T->SystemSize*(it_e->v2->q_b - it_e->q_c2))/L;
//                v1_b_norm = shift(it_e->v1->coord_b, T->SystemSize*(it_e->v1->q_b - it_e->q_v2 - it_e->q_c2))/L;
//                b_norm = shift(it_e->BC_l, T->SystemSize*(it_e->q_BC_l - it_e->q_c2 - it_e->q_v2))/L;
//                
//                // add the four triangles to the patch
//                SurfacePatch *LPatch_c2 = new SurfacePatch(geom);
//                LPatch_c2->addTri(v1_a_norm,v2_a_norm,b_norm);
//                LPatch_c2->addTri(v2_a_norm,v2_b_norm,b_norm);
//                LPatch_c2->addTri(v2_b_norm,v1_b_norm,b_norm);
//                LPatch_c2->addTri(v1_b_norm,v1_a_norm,b_norm);
//                
//                LPatch_c2->addTri(v1_a_norm,b_norm,v2_a_norm);
//                LPatch_c2->addTri(v2_a_norm,b_norm,v2_b_norm);
//                LPatch_c2->addTri(v2_b_norm,b_norm,v1_b_norm);
//                LPatch_c2->addTri(v1_b_norm,b_norm,v1_a_norm);
//
//                
//                if(it_e->isAtBoundary()) qSetColor(LPatch_c2->faceColor, lateralPatchColor, 0.5);
//                else  qSetColor(LPatch_c2->faceColor, lateralCystColor, 0.5);
////                    
////                    qSetColor(LPatch_c2->faceColor, lateralPatchColor, 0.0);
//                LParts << LPatch_c2;
//                
//                // draw the vertex connection
//                LinePatch *LLinePatch_c2 = new LinePatch(geom);
//                LLinePatch_c2->addLine(v1_a_norm,v1_b_norm);
//                LLinePatch_c2->addLine(v2_a_norm,v2_b_norm);
//                
//                qSetColor(LLinePatch_c2->faceColor, QColor(Qt::black), 0.0);
//                LLineParts << LLinePatch_c2;
//            }
//        }
//        
//        
//    }
//}

void GLTissueImage::loadPartsInRadius(double radius)
{
    // colors & opacities of the apical patches for each cell type
    QColor apicalPatchColors [4] =  {};
    double apicalPatchOpacities [4] = {};
    QColor apicalPatchColor_fixed = QColor(Qt::black); // if cells are fixed they get a different color

    // colors of the basal patches for each cell type
    QColor basalPatchColors [4] =  {};//{ QColor(Qt::gray), QColor(Qt::white), QColor(Qt::green), QColor(Qt::blue), QColor(Qt::darkGray) };
    double basalPatchOpacities [4] = {};//{ 0.0, 0.5, 0.7, 0.3, 0.5 };

    // colors of apical & basal lines
    QColor apicalLineColors [4][4] =  {};
    QColor basalLineColors [4][4] =  {};
    
    // color of lateral interfaces
    QColor lateralPatchColors [4][4] = {};
    double lateralPatchOpacities [4][4] = {};

    // set the colors depending on the simulation type
    if(integrator->T->onlyApicalContributions) { // gastrulation
        apicalPatchColors[1] = QColor(Qt::red);
        apicalPatchColors[2] = QColor(Qt::darkGreen);
        apicalPatchColors[3] = QColor(Qt::darkGreen);
        
        apicalPatchOpacities[1] = 0.4;
        apicalPatchOpacities[2] = 0.7;
        apicalPatchOpacities[3] = 0.3;

        apicalLineColors[1][1] = QColor(Qt::black);
        apicalLineColors[2][2] = QColor(Qt::black);
        apicalLineColors[2][1] = QColor(Qt::magenta);
        apicalLineColors[1][2] = QColor(Qt::magenta);
        
    } else { // cyst formation
        
        // apical side
        apicalPatchColors[1] = QColor(Qt::white);
        apicalPatchColors[2] = QColor(Qt::green);
        apicalPatchColors[3] = QColor(Qt::red);
        apicalPatchOpacities[1] = 0.5;
        apicalPatchOpacities[2] = 0.7;
        apicalPatchOpacities[3] = 0.3;
        apicalLineColors[1][1] = QColor(Qt::red);
        apicalLineColors[2][2] = QColor(Qt::red);
        apicalLineColors[2][1] = QColor(Qt::magenta);
        apicalLineColors[1][2] = QColor(Qt::magenta);
        
        // basal side
        basalPatchColors[1] = QColor(Qt::lightGray);
        basalPatchColors[2] = QColor(Qt::green); basalPatchColors[2] = basalPatchColors[2].light();
        basalPatchColors[3] = QColor(Qt::cyan);
        basalPatchOpacities[1] = 1.0;
        basalPatchOpacities[2] = 1.0;
        basalPatchOpacities[3] = 0.3;
        basalLineColors[1][1] = QColor(Qt::black);
        basalLineColors[2][2] = QColor(Qt::black);
        basalLineColors[2][1] = QColor(Qt::magenta);

        // lateral side
        lateralPatchColors[1][1] = QColor(Qt::darkGray);
        lateralPatchColors[2][2] = QColor(Qt::darkGray);
        lateralPatchColors[1][2] = QColor(Qt::blue);
        lateralPatchColors[2][1] = QColor(Qt::blue);
        lateralPatchOpacities[1][1] = 0.4;
        lateralPatchOpacities[2][2] = 0.4;
        lateralPatchOpacities[1][2] = 0.4;
        lateralPatchOpacities[2][1] = 0.4;

    }
   
    
    // clear the list of parts
    removeParts();
    Tissue* T = integrator->T;
    
    // update the systems outermost points
    T->UpdateExtPositions();
    
    // get the biggest size in the system
    double Lz = T->extPosition_topRight.z-T->extPosition_bottomLeft.z;
    double L = std::max(T->SystemSize.x, std::max(T->SystemSize.y,Lz));;
    
    // draw the box around the tissue
    bool drawBox = true;
    
    Point midPoint = Point(T->SystemSize.x/2,T->SystemSize.y/2,0);
    
    
    
    if(T->isPeriodic && drawBox)
    {
        // system size and the box size, that gives the size of the box which is drawn around the tissue at relative distance exceed
        double exceed = 0.1;
        Point boxSize = Point(T->SystemSize.x,T->SystemSize.y,Lz*(1+2*exceed))/L;
        Point box_origin = Point(0,0,T->extPosition_bottomLeft.z/L-Lz*(exceed/L));
        
        // draw the box outlines
        LinePatch *box_lines = new LinePatch(geom);
        
        box_lines->addLine(box_origin,box_origin+Point(boxSize.x,0,0));
        box_lines->addLine(box_origin,box_origin+Point(0,boxSize.y,0));
        box_lines->addLine(box_origin,box_origin+Point(0,0,boxSize.z));
        
        box_lines->addLine(box_origin+Point(boxSize.x,0,0),box_origin+Point(boxSize.x,boxSize.y,0));
        box_lines->addLine(box_origin+Point(boxSize.x,0,0),box_origin+Point(boxSize.x,0,boxSize.z));
        
        box_lines->addLine(box_origin+Point(0,boxSize.y,0),box_origin+Point(boxSize.x,boxSize.y,0));
        box_lines->addLine(box_origin+Point(0,boxSize.y,0),box_origin+Point(0,boxSize.y,boxSize.z));
        
        box_lines->addLine(box_origin+Point(0,0,boxSize.z),box_origin+Point(0,boxSize.y,boxSize.z));
        box_lines->addLine(box_origin+Point(0,0,boxSize.z),box_origin+Point(boxSize.x,0,boxSize.z));
        
        box_lines->addLine(box_origin+Point(boxSize.x,boxSize.y,0),box_origin+Point(boxSize.x,boxSize.y,boxSize.z));
        box_lines->addLine(box_origin+Point(boxSize.x,0,boxSize.z),box_origin+Point(boxSize.x,boxSize.y,boxSize.z));
        box_lines->addLine(box_origin+Point(0,boxSize.y,boxSize.z),box_origin+Point(boxSize.x,boxSize.y,boxSize.z));
        
        qSetColor(box_lines->faceColor, QColor(Qt::black),0.5);
        BoxLineParts<<box_lines;
        
        // draw the box surfaces
        SurfacePatch *BoxPatch = new SurfacePatch(geom); // add surfaces of the box with weak opacity
        QColor boxColor = QColor(Qt::black);
        double boxOpacity = 0.0;
        
        // define basic elements of box
        Point P = box_origin;
        Point vx = Point(boxSize.x,0,0);
        Point vy = Point(0,boxSize.y,0);
        Point vz = Point(0,0,boxSize.z);
        
        // front patch
        BoxPatch->addTri(P, P+vx, P+vx+vz);
        BoxPatch->addTri(P, P+vz, P+vx+vz);

        // back patch
        BoxPatch->addTri(P+vy, P+vx+vy, P+vx+vy+vz);
        BoxPatch->addTri(P+vy, P+vy+vz, P+vx+vy+vz);

        // left patch
        BoxPatch->addTri(P, P+vy, P+vy+vz);
        BoxPatch->addTri(P, P+vz, P+vy+vz);

        // right patch
        BoxPatch->addTri(P+vx, P+vx+vy, P+vx+vy+vz);
        BoxPatch->addTri(P+vx, P+vx+vz, P+vx+vy+vz);
        
        qSetColor(BoxPatch->faceColor, boxColor,boxOpacity);
        //BoxParts<<BoxPatch;

        
    }
    
    for(std::list<Cell>::iterator it_c=T->ListCell.begin(); it_c!=T->ListCell.end(); ++it_c)
    {
        // if the cell's midpoint is not close enough to the midpoint skip it (to show round pieces of tissue)
        double dist = sqrt(pow(it_c->BC_a.x-midPoint.x,2) + pow(it_c->BC_a.y-midPoint.y,2));
        if(dist>radius) continue;
        
        if(!(*it_c).CrossBoundary())
        {
            
            // normalize the variables
            Cell* c = &(*it_c);
            Point BC_a = c->BC_a;
            Point BC_a_norm = BC_a/L;
            
            Point BC_b_norm = c->BC_b/L;
            
            // define surface patches for apical and basal surfaces
            SurfacePatch *APatch = new SurfacePatch(geom);
            
            // add the apical side
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                Edge *e = (*it_e);
                
                // calculate the normalized coordinates
                Point v1_a_norm = e->v1->coord/L;
                Point v2_a_norm = e->v2->coord/L;
                
                APatch->addTri(v1_a_norm,BC_a_norm,v2_a_norm);
                APatch->addTri(BC_a_norm,v1_a_norm,v2_a_norm);
            }
            
            qSetColor(APatch->faceColor, apicalPatchColors[c->Type],apicalPatchOpacities[c->Type]);
            
            //if(!((*it_c).positionFixed))  qSetColor(APatch->faceColor, apicalPatchColors[c->Type],apicalPatchOpacities[c->Type]);
            //else qSetColor(APatch->faceColor, apicalPatchColor_fixed,apicalPatchOpacities[c->Type]);
            
            if(c->Type==1) AParts << APatch;
            else AParts_c << APatch;
            
            
            // add lines between same cell types
            LinePatch *ALinePatch = new LinePatch(geom);
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                if((*it_e)->c1->Type == (*it_e)->c2->Type) ALinePatch->addLine((*it_e)->v1->coord/L,(*it_e)->v2->coord/L);
            }
            qSetColor(ALinePatch->faceColor, apicalLineColors[1][1], 1.0);
            ALineParts << ALinePatch;
            
            // add lines between different cell types
            LinePatch *ALinePatch_d = new LinePatch(geom);
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                if((*it_e)->c1->Type != (*it_e)->c2->Type) ALinePatch_d->addLine((*it_e)->v1->coord/L,(*it_e)->v2->coord/L);
            }
            
            qSetColor(ALinePatch_d->faceColor, apicalLineColors[2][1], 1.0);
            ALineParts << ALinePatch_d;
            
            // add basal surface patches
            SurfacePatch *BPatch = new SurfacePatch(geom);
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                Edge *e = (*it_e);
                
                // calculate the normalized coordinates
                Point v1_b_norm = e->v1->coord_b/L;
                Point v2_b_norm = e->v2->coord_b/L;
                
                BPatch->addTri(BC_b_norm,v1_b_norm,v2_b_norm);
                BPatch->addTri(v1_b_norm,BC_b_norm,v2_b_norm);
            }
            qSetColor(BPatch->faceColor, basalPatchColors[c->Type],basalPatchOpacities[c->Type]);
            BParts << BPatch;
            
            // add basal lines
            LinePatch *BLinePatch = new LinePatch(geom);
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                BLinePatch->addLine((*it_e)->v1->coord_b/L,(*it_e)->v2->coord_b/L);
            }
            qSetColor(BLinePatch->faceColor, basalLineColors[1][1], 0.0);
            BLineParts << BLinePatch;
        }
        else { // if the cell crosses the boundary
            
            Cell* c = &(*it_c);
            
            // define surface patches for apical and basal surfaces
            SurfacePatch *APatch = new SurfacePatch(geom);
            
            // add the apical side
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                Edge *e = (*it_e);
                
                // add the triangles (normal is defined by the direction towards the cell)
                if(e->c1==c)
                {
                    // calculate the normalized coordinates
                    Point v1_a_norm = shift(e->v1->coord,T->SystemSize*e->q_c1*(-1))/L;
                    Point v2_a_norm = shift(e->v2->coord,T->SystemSize*(e->q_v2-e->q_c1))/L;
                    Point BC_a_norm = c->BC_a/L;
                    APatch->addTri(v1_a_norm,BC_a_norm,v2_a_norm);
                    APatch->addTri(v1_a_norm,v2_a_norm,BC_a_norm);
                }
                else
                {
                    Point v1_a_norm = shift(e->v1->coord,T->SystemSize*(e->q_c2+e->q_v2)*(-1))/L;
                    Point v2_a_norm = shift(e->v2->coord,T->SystemSize*(e->q_c2)*(-1))/L;
                    Point BC_a_norm = c->BC_a/L;
                    APatch->addTri(v1_a_norm,BC_a_norm,v2_a_norm);
                    APatch->addTri(v1_a_norm,v2_a_norm,BC_a_norm);
                }
                
            }
            
            qSetColor(APatch->faceColor, apicalPatchColors[c->Type],apicalPatchOpacities[c->Type]);
            
            AParts << APatch;
            
            LinePatch *ALinePatch = new LinePatch(geom);
            
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                Edge *e = (*it_e);
                
                if(e->c1==c)
                {
                    // calculate the normalized coordinates
                    Point v1_a_norm = shift(e->v1->coord,T->SystemSize*e->q_c1*(-1))/L;
                    Point v2_a_norm = shift(e->v2->coord,T->SystemSize*(e->q_v2-e->q_c1))/L;
                    ALinePatch->addLine(v1_a_norm,v2_a_norm);
                }
                else
                {
                    Point v1_a_norm = shift(e->v1->coord,T->SystemSize*(e->q_c2+e->q_v2)*(-1))/L;
                    Point v2_a_norm = shift(e->v2->coord,T->SystemSize*(e->q_c2)*(-1))/L;
                    ALinePatch->addLine(v1_a_norm,v2_a_norm);
                }
            }
            
            qSetColor(ALinePatch->faceColor, apicalLineColors[c->Type][c->Type], 0.0);
            ALineParts << ALinePatch;
            
            
            SurfacePatch *BPatch = new SurfacePatch(geom);
            // add the basal side
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                Edge *e = (*it_e);
                
                if(e->c1==c)
                {
                    // calculate the normalized coordinates
                    Point v1_b_norm = shift(e->v1->coord_b,T->SystemSize*(e->v1->q_b-e->q_c1))/L;
                    Point v2_b_norm = shift(e->v2->coord_b,T->SystemSize*(e->v2->q_b+e->q_v2-e->q_c1))/L;
                    Point BC_b_norm = shift(c->BC_b,T->SystemSize*c->q_BC_b)/L;
                    BPatch->addTri(v1_b_norm,BC_b_norm,v2_b_norm);
                    BPatch->addTri(v1_b_norm,v2_b_norm,BC_b_norm);
                }
                else
                {
                    Point v1_b_norm = shift(e->v1->coord_b,T->SystemSize*(e->v1->q_b-e->q_c2-e->q_v2))/L;
                    Point v2_b_norm = shift(e->v2->coord_b,T->SystemSize*(e->v2->q_b-e->q_c2))/L;
                    Point BC_b_norm = shift(c->BC_b,T->SystemSize*c->q_BC_b)/L;
                    BPatch->addTri(v1_b_norm,BC_b_norm,v2_b_norm);
                    BPatch->addTri(v1_b_norm,v2_b_norm,BC_b_norm);
                }
                
            }
            
            
            qSetColor(BPatch->faceColor, basalPatchColors[c->Type],basalPatchOpacities[c->Type]);
            BParts << BPatch;
            
            LinePatch *BLinePatch = new LinePatch(geom);
            
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                Edge *e = (*it_e);
                
                if(e->c1==c)
                {
                    // calculate the normalized coordinates
                    Point v1_b_norm = shift(e->v1->coord_b,T->SystemSize*(e->v1->q_b-e->q_c1))/L;
                    Point v2_b_norm = shift(e->v2->coord_b,T->SystemSize*(e->q_v2+e->v2->q_b-e->q_c1))/L;
                    BLinePatch->addLine(v1_b_norm,v2_b_norm);
                }
                else
                {
                    Point v1_b_norm = shift(e->v1->coord_b,T->SystemSize*(e->v1->q_b-e->q_c2-e->q_v2))/L;
                    Point v2_b_norm = shift(e->v2->coord_b,T->SystemSize*(e->v2->q_b-e->q_c2))/L;
                    BLinePatch->addLine(v1_b_norm,v2_b_norm);
                }
            }
            
            qSetColor(BLinePatch->faceColor, basalLineColors[c->Type][c->Type], 0.0);
            BLineParts << BLinePatch;
            
        }
        
    }
    
    
    // iterate through the edges deepest inside which are supposed to be visible i.e. the ones at mut-wt interface
    for(std::list<Edge>::iterator it_e = T->ListEdge.begin(); it_e!=T->ListEdge.end(); ++it_e)
    {
        if(it_e->c1->Type==it_e->c2->Type) continue;
        
        if(!T->isPeriodic)
        {
            // normalize the coordinates
            Point v1_a_norm = it_e->v1->coord/L;
            Point v2_a_norm = it_e->v2->coord/L;
            Point v1_b_norm = it_e->v1->coord_b/L;
            Point v2_b_norm = it_e->v2->coord_b/L;
            Point b_norm = it_e->BC_l/L;
            
            SurfacePatch *LPatch = new SurfacePatch(geom);
            
            // add the four triangles to the patch
            LPatch->addTri(v1_a_norm,v2_a_norm,b_norm);
            LPatch->addTri(v2_a_norm,v2_b_norm,b_norm);
            LPatch->addTri(v2_b_norm,v1_b_norm,b_norm);
            LPatch->addTri(v1_b_norm,v1_a_norm,b_norm);
            
            qSetColor(LPatch->faceColor, lateralPatchColors[it_e->c1->Type][it_e->c2->Type],lateralPatchOpacities[it_e->c1->Type][it_e->c2->Type]);
            LParts << LPatch;
            
            // draw the vertex connection
            LinePatch *LLinePatch = new LinePatch(geom);
            LLinePatch->addLine(v1_a_norm,v1_b_norm);
            LLinePatch->addLine(v2_a_norm,v2_b_norm);
            
            qSetColor(LLinePatch->faceColor, QColor(Qt::black),0.0);
            LLineParts_cyst << LLinePatch;
            
        }
        else {
            // draw the edge on the side of cell 1
            // normalize the coordinates
            Point v1_a_norm = shift(it_e->v1->coord, T->SystemSize*it_e->q_c1*(-1))/L;
            Point v2_a_norm = shift(it_e->v2->coord, T->SystemSize*(it_e->q_v2-it_e->q_c1))/L;
            Point v1_b_norm = shift(it_e->v1->coord_b, T->SystemSize*(it_e->v1->q_b - it_e->q_c1))/L;
            Point v2_b_norm = shift(it_e->v2->coord_b, T->SystemSize*(it_e->q_v2 + it_e->v2->q_b - it_e->q_c1))/L;
            Point b_norm = shift(it_e->BC_l, T->SystemSize*(it_e->q_BC_l - it_e->q_c1))/L;
            
            // add the four triangles to the patch
            SurfacePatch *LPatch_c1 = new SurfacePatch(geom);
            LPatch_c1->addTri(v1_a_norm,v2_a_norm,b_norm);
            LPatch_c1->addTri(v1_a_norm,b_norm,v2_a_norm);
            LPatch_c1->addTri(v2_a_norm,v2_b_norm,b_norm);
            LPatch_c1->addTri(v2_a_norm,b_norm,v2_b_norm);
            LPatch_c1->addTri(v2_b_norm,v1_b_norm,b_norm);
            LPatch_c1->addTri(v2_b_norm,b_norm,v1_b_norm);
            LPatch_c1->addTri(v1_b_norm,v1_a_norm,b_norm);
            LPatch_c1->addTri(v1_b_norm,b_norm,v1_a_norm);
            
            qSetColor(LPatch_c1->faceColor, lateralPatchColors[it_e->c1->Type][it_e->c2->Type],lateralPatchOpacities[it_e->c1->Type][it_e->c2->Type]);
            LParts << LPatch_c1;

            // draw the vertex connection
            LinePatch *LLinePatch_c1 = new LinePatch(geom);
            LLinePatch_c1->addLine(v1_a_norm,v1_b_norm);
            LLinePatch_c1->addLine(v2_a_norm,v2_b_norm);
            
            qSetColor(LLinePatch_c1->faceColor, QColor(Qt::black), 1);
            LLineParts_cyst << LLinePatch_c1;
            
            // draw the edge on the side of cell 2
            // normalize the coordinates
            v2_a_norm = shift(it_e->v2->coord, T->SystemSize*it_e->q_c2*(-1))/L;
            v1_a_norm = shift(it_e->v1->coord, T->SystemSize*(it_e->q_v2+it_e->q_c2)*(-1))/L;
            v2_b_norm = shift(it_e->v2->coord_b, T->SystemSize*(it_e->v2->q_b - it_e->q_c2))/L;
            v1_b_norm = shift(it_e->v1->coord_b, T->SystemSize*(it_e->v1->q_b - it_e->q_v2 - it_e->q_c2))/L;
            b_norm = shift(it_e->BC_l, T->SystemSize*(it_e->q_BC_l - it_e->q_c2 - it_e->q_v2))/L;
            
            // add the four triangles to the patch
            SurfacePatch *LPatch_c2 = new SurfacePatch(geom);
            LPatch_c2->addTri(v1_a_norm,v2_a_norm,b_norm);
            LPatch_c2->addTri(v2_a_norm,v2_b_norm,b_norm);
            LPatch_c2->addTri(v2_b_norm,v1_b_norm,b_norm);
            LPatch_c2->addTri(v1_b_norm,v1_a_norm,b_norm);
            
            LPatch_c2->addTri(v1_a_norm,b_norm,v2_a_norm);
            LPatch_c2->addTri(v2_a_norm,b_norm,v2_b_norm);
            LPatch_c2->addTri(v2_b_norm,b_norm,v1_b_norm);
            LPatch_c2->addTri(v1_b_norm,b_norm,v1_a_norm);
            
            
            qSetColor(LPatch_c2->faceColor, lateralPatchColors[it_e->c1->Type][it_e->c2->Type],lateralPatchOpacities[it_e->c1->Type][it_e->c2->Type]);
            LParts << LPatch_c2;
            
            // draw the vertex connection
            LinePatch *LLinePatch_c2 = new LinePatch(geom);
            LLinePatch_c2->addLine(v1_a_norm,v1_b_norm);
            LLinePatch_c2->addLine(v2_a_norm,v2_b_norm);
            
            qSetColor(LLinePatch_c2->faceColor, QColor(Qt::black), 1);
            LLineParts_cyst << LLinePatch_c2;
        }
    }
    
    // draw opaque edges at boundary
    for(std::list<Edge>::iterator it_e = T->ListEdge.begin(); it_e!=T->ListEdge.end(); ++it_e)
    {
        // check which of the neighbours are in the radius
        bool c1InRadius = sqrt(pow(it_e->c1->BC_a.x-midPoint.x,2) + pow(it_e->c1->BC_a.y-midPoint.y,2))<radius;
        bool c2InRadius = sqrt(pow(it_e->c2->BC_a.x-midPoint.x,2) + pow(it_e->c2->BC_a.y-midPoint.y,2))<radius;
        
        // draw edge only if it is at interface between visible and invisible cells or at the boundary of the tissue
        if(c1InRadius==c2InRadius && !it_e->isAtBoundary()) continue;
        
        if(!T->isPeriodic)
        {
            // normalize the coordinates
            Point v1_a_norm = it_e->v1->coord/L;
            Point v2_a_norm = it_e->v2->coord/L;
            Point v1_b_norm = it_e->v1->coord_b/L;
            Point v2_b_norm = it_e->v2->coord_b/L;
            Point b_norm = it_e->BC_l/L;
            
            SurfacePatch *LPatch = new SurfacePatch(geom);
            
            // add the four triangles to the patch
            LPatch->addTri(v1_a_norm,v2_a_norm,b_norm);
            LPatch->addTri(v2_a_norm,v2_b_norm,b_norm);
            LPatch->addTri(v2_b_norm,v1_b_norm,b_norm);
            LPatch->addTri(v1_b_norm,v1_a_norm,b_norm);
            
            qSetColor(LPatch->faceColor, QColor(Qt::black),0.0);
            LParts << LPatch;
            
            // draw the vertex connection
            LinePatch *LLinePatch = new LinePatch(geom);
            LLinePatch->addLine(v1_a_norm,v1_b_norm);
            LLinePatch->addLine(v2_a_norm,v2_b_norm);
            
            qSetColor(LLinePatch->faceColor, QColor(Qt::black),0.0);
            LLineParts << LLinePatch;
            
        }
        else {
            // draw the edge on the side of cell 1
            // normalize the coordinates
            Point v1_a_norm = shift(it_e->v1->coord, T->SystemSize*it_e->q_c1*(-1))/L;
            Point v2_a_norm = shift(it_e->v2->coord, T->SystemSize*(it_e->q_v2-it_e->q_c1))/L;
            Point v1_b_norm = shift(it_e->v1->coord_b, T->SystemSize*(it_e->v1->q_b - it_e->q_c1))/L;
            Point v2_b_norm = shift(it_e->v2->coord_b, T->SystemSize*(it_e->q_v2 + it_e->v2->q_b - it_e->q_c1))/L;
            Point b_norm = shift(it_e->BC_l, T->SystemSize*(it_e->q_BC_l - it_e->q_c1))/L;
            
            // add the four triangles to the patch
            SurfacePatch *LPatch_c1 = new SurfacePatch(geom);
            LPatch_c1->addTri(v1_a_norm,v2_a_norm,b_norm);
            LPatch_c1->addTri(v1_a_norm,b_norm,v2_a_norm);
            LPatch_c1->addTri(v2_a_norm,v2_b_norm,b_norm);
            LPatch_c1->addTri(v2_a_norm,b_norm,v2_b_norm);
            LPatch_c1->addTri(v2_b_norm,v1_b_norm,b_norm);
            LPatch_c1->addTri(v2_b_norm,b_norm,v1_b_norm);
            LPatch_c1->addTri(v1_b_norm,v1_a_norm,b_norm);
            LPatch_c1->addTri(v1_b_norm,b_norm,v1_a_norm);
            
            
            qSetColor(LPatch_c1->faceColor, lateralPatchColors[it_e->c1->Type][it_e->c1->Type], lateralPatchOpacities[it_e->c1->Type][it_e->c1->Type]);
            LParts << LPatch_c1;
            
            // draw the vertex connection
            LinePatch *LLinePatch_c1 = new LinePatch(geom);
            LLinePatch_c1->addLine(v1_a_norm,v1_b_norm);
            LLinePatch_c1->addLine(v2_a_norm,v2_b_norm);
            
            qSetColor(LLinePatch_c1->faceColor, QColor(Qt::black), 0.0);
            LLineParts << LLinePatch_c1;
            
            // draw the edge on the side of cell 2
            // normalize the coordinates
            v2_a_norm = shift(it_e->v2->coord, T->SystemSize*it_e->q_c2*(-1))/L;
            v1_a_norm = shift(it_e->v1->coord, T->SystemSize*(it_e->q_v2+it_e->q_c2)*(-1))/L;
            v2_b_norm = shift(it_e->v2->coord_b, T->SystemSize*(it_e->v2->q_b - it_e->q_c2))/L;
            v1_b_norm = shift(it_e->v1->coord_b, T->SystemSize*(it_e->v1->q_b - it_e->q_v2 - it_e->q_c2))/L;
            b_norm = shift(it_e->BC_l, T->SystemSize*(it_e->q_BC_l - it_e->q_c2 - it_e->q_v2))/L;
            
            // add the four triangles to the patch
            SurfacePatch *LPatch_c2 = new SurfacePatch(geom);
            LPatch_c2->addTri(v1_a_norm,v2_a_norm,b_norm);
            LPatch_c2->addTri(v2_a_norm,v2_b_norm,b_norm);
            LPatch_c2->addTri(v2_b_norm,v1_b_norm,b_norm);
            LPatch_c2->addTri(v1_b_norm,v1_a_norm,b_norm);
            
            LPatch_c2->addTri(v1_a_norm,b_norm,v2_a_norm);
            LPatch_c2->addTri(v2_a_norm,b_norm,v2_b_norm);
            LPatch_c2->addTri(v2_b_norm,b_norm,v1_b_norm);
            LPatch_c2->addTri(v1_b_norm,b_norm,v1_a_norm);
            
            qSetColor(LPatch_c2->faceColor, lateralPatchColors[it_e->c1->Type][it_e->c1->Type], lateralPatchOpacities[it_e->c1->Type][it_e->c1->Type]);
            LParts << LPatch_c2;
            
            // draw the vertex connection
            LinePatch *LLinePatch_c2 = new LinePatch(geom);
            LLinePatch_c2->addLine(v1_a_norm,v1_b_norm);
            LLinePatch_c2->addLine(v2_a_norm,v2_b_norm);
            
            qSetColor(LLinePatch_c2->faceColor, QColor(Qt::black), 0.0);
            LLineParts << LLinePatch_c2;
        }
    }
    
    
    
}


void GLTissueImage::removeParts()
{
    
    // delete all the patches
    AParts.clear();
    AParts_c.clear();
    BParts.clear();
    LParts.clear();
    ALineParts.clear();
    BLineParts.clear();
    LLineParts.clear();
    LLineParts_cyst.clear();

    BoxLineParts.clear();
    
    geom->vertices.clear();
    geom->faces.clear();
    geom->lines.clear();
    geom->normals.clear();
}

void GLTissueImage::loadPartsInRadius(double radius, double xr, double yr)
{
    // colors & opacities of the apical patches for each cell type
    QColor apicalPatchColors [4] =  {};
    double apicalPatchOpacities [4] = {};
    QColor apicalPatchColor_fixed = QColor(Qt::black); // if cells are fixed they get a different color
    
    // colors of the basal patches for each cell type
    QColor basalPatchColors [4] =  {};//{ QColor(Qt::gray), QColor(Qt::white), QColor(Qt::green), QColor(Qt::blue), QColor(Qt::darkGray) };
    double basalPatchOpacities [4] = {};//{ 0.0, 0.5, 0.7, 0.3, 0.5 };
    
    // colors of apical & basal lines
    QColor apicalLineColors [4][4] =  {};
    QColor basalLineColors [4][4] =  {};
    
    // color of lateral interfaces
    QColor lateralPatchColors [4][4] = {};
    double lateralPatchOpacities [4][4] = {};
    
    // set the colors depending on the simulation type
    if(integrator->T->onlyApicalContributions) { // gastrulation
        apicalPatchColors[1] = QColor(Qt::red);
        apicalPatchColors[2] = QColor(Qt::darkGreen);
        apicalPatchColors[3] = QColor(Qt::darkGreen);
        
        apicalPatchOpacities[1] = 0.4;
        apicalPatchOpacities[2] = 0.7;
        apicalPatchOpacities[3] = 0.3;
        
        apicalLineColors[1][1] = QColor(Qt::black);
        apicalLineColors[2][2] = QColor(Qt::black);
        apicalLineColors[2][1] = QColor(Qt::magenta);
        apicalLineColors[1][2] = QColor(Qt::magenta);
        
    } else { // cyst formation
        
        // apical side
        apicalPatchColors[1] = QColor(Qt::white);
        apicalPatchColors[2] = QColor(Qt::green);
        apicalPatchColors[3] = QColor(Qt::red);
        apicalPatchOpacities[1] = 0.5;
        apicalPatchOpacities[2] = 0.7;
        apicalPatchOpacities[3] = 0.3;
        apicalLineColors[1][1] = QColor(Qt::red);
        apicalLineColors[2][2] = QColor(Qt::red);
        apicalLineColors[2][1] = QColor(Qt::magenta);
        apicalLineColors[1][2] = QColor(Qt::magenta);
        
        // basal side
        basalPatchColors[1] = QColor(Qt::lightGray);
        basalPatchColors[2] = QColor(Qt::green); basalPatchColors[2] = basalPatchColors[2].light();
        basalPatchColors[3] = QColor(Qt::cyan);
        basalPatchOpacities[1] = 1.0;
        basalPatchOpacities[2] = 1.0;
        basalPatchOpacities[3] = 0.3;
        basalLineColors[1][1] = QColor(Qt::black);
        basalLineColors[2][2] = QColor(Qt::black);
        basalLineColors[2][1] = QColor(Qt::magenta);
        
        // lateral side
        lateralPatchColors[1][1] = QColor(Qt::darkGray);
        lateralPatchColors[2][2] = QColor(Qt::darkGray);
        lateralPatchColors[1][2] = QColor(Qt::blue);
        lateralPatchColors[2][1] = QColor(Qt::blue);
        lateralPatchOpacities[1][1] = 0.4;
        lateralPatchOpacities[2][2] = 0.4;
        lateralPatchOpacities[1][2] = 0.4;
        lateralPatchOpacities[2][1] = 0.4;
        
    }
    
    
    // clear the list of parts
    removeParts();
    Tissue* T = integrator->T;
    
    // update the systems outermost points
    T->UpdateExtPositions();
    
    // get the biggest size in the system
    double Lz = T->extPosition_topRight.z-T->extPosition_bottomLeft.z;
    double L = std::max(T->SystemSize.x, std::max(T->SystemSize.y,Lz));;
    
    // draw the box around the tissue
    bool drawBox = false;
    
    Point midPoint = Point(T->SystemSize.x/2,T->SystemSize.y/2,0);
    
    if(T->isPeriodic && drawBox)
    {
        // system size and the box size, that gives the size of the box which is drawn around the tissue at relative distance exceed
        double exceed = 0.0;
        Point boxSize = Point(T->SystemSize.x,T->SystemSize.y,Lz*(1+2*exceed))/L;
        Point box_origin = Point(0,0,T->extPosition_bottomLeft.z/L-Lz*(exceed/L));
        
        // draw the box outlines
        LinePatch *box_lines = new LinePatch(geom);
        
        box_lines->addLine(box_origin,box_origin+Point(boxSize.x,0,0), xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
        box_lines->addLine(box_origin,box_origin+Point(0,boxSize.y,0), xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
        box_lines->addLine(box_origin,box_origin+Point(0,0,boxSize.z), xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
        
        box_lines->addLine(box_origin+Point(boxSize.x,0,0),box_origin+Point(boxSize.x,boxSize.y,0), xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
        box_lines->addLine(box_origin+Point(boxSize.x,0,0),box_origin+Point(boxSize.x,0,boxSize.z), xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
        
        box_lines->addLine(box_origin+Point(0,boxSize.y,0),box_origin+Point(boxSize.x,boxSize.y,0), xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
        box_lines->addLine(box_origin+Point(0,boxSize.y,0),box_origin+Point(0,boxSize.y,boxSize.z), xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
        
        box_lines->addLine(box_origin+Point(0,0,boxSize.z),box_origin+Point(0,boxSize.y,boxSize.z), xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
        box_lines->addLine(box_origin+Point(0,0,boxSize.z),box_origin+Point(boxSize.x,0,boxSize.z), xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
        
        box_lines->addLine(box_origin+Point(boxSize.x,boxSize.y,0),box_origin+Point(boxSize.x,boxSize.y,boxSize.z), xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
        box_lines->addLine(box_origin+Point(boxSize.x,0,boxSize.z),box_origin+Point(boxSize.x,boxSize.y,boxSize.z), xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
        box_lines->addLine(box_origin+Point(0,boxSize.y,boxSize.z),box_origin+Point(boxSize.x,boxSize.y,boxSize.z), xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
        
        qSetColor(box_lines->faceColor, QColor(Qt::black),0.5);
        BoxLineParts<<box_lines;
        
        // draw the box surfaces
        SurfacePatch *BoxPatch = new SurfacePatch(geom); // add surfaces of the box with weak opacity
        QColor boxColor = QColor(Qt::black);
        double boxOpacity = 0.0;
        
        // define basic elements of box
        Point P = box_origin;
        Point vx = Point(boxSize.x,0,0);
        Point vy = Point(0,boxSize.y,0);
        Point vz = Point(0,0,boxSize.z);
        
        // front patch
        BoxPatch->addTri(P, P+vx, P+vx+vz);
        BoxPatch->addTri(P, P+vz, P+vx+vz);
        
        // back patch
        BoxPatch->addTri(P+vy, P+vx+vy, P+vx+vy+vz);
        BoxPatch->addTri(P+vy, P+vy+vz, P+vx+vy+vz);
        
        // left patch
        BoxPatch->addTri(P, P+vy, P+vy+vz);
        BoxPatch->addTri(P, P+vz, P+vy+vz);
        
        // right patch
        BoxPatch->addTri(P+vx, P+vx+vy, P+vx+vy+vz);
        BoxPatch->addTri(P+vx, P+vx+vz, P+vx+vy+vz);
        
        qSetColor(BoxPatch->faceColor, boxColor,boxOpacity);
        //BoxParts<<BoxPatch;
        
        
    }
    
    for(std::list<Cell>::iterator it_c=T->ListCell.begin(); it_c!=T->ListCell.end(); ++it_c)
    {
        // if the cell's midpoint is not close enough to the midpoint skip it (to show round pieces of tissue)
        double dist = sqrt(pow(it_c->BC_a.x-midPoint.x,2) + pow(it_c->BC_a.y-midPoint.y,2));
        if(dist>radius) continue;
        
        if(!(*it_c).CrossBoundary())
        {
            
            // normalize the variables
            Cell* c = &(*it_c);
            Point BC_a = c->BC_a;
            Point BC_a_norm = BC_a/L;
            
            Point BC_b_norm = c->BC_b/L;
            
            // define surface patches for apical and basal surfaces
            SurfacePatch *APatch = new SurfacePatch(geom);
            
            // add the apical side
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                Edge *e = (*it_e);
                
                // calculate the normalized coordinates
                Point v1_a_norm = e->v1->coord/L;
                Point v2_a_norm = e->v2->coord/L;
                
                APatch->addTri(v1_a_norm,BC_a_norm,v2_a_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
                APatch->addTri(BC_a_norm,v1_a_norm,v2_a_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            }
            
            qSetColor(APatch->faceColor, apicalPatchColors[c->Type],apicalPatchOpacities[c->Type]);
            
            //if(!((*it_c).positionFixed))  qSetColor(APatch->faceColor, apicalPatchColors[c->Type],apicalPatchOpacities[c->Type]);
            //else qSetColor(APatch->faceColor, apicalPatchColor_fixed,apicalPatchOpacities[c->Type]);
            
            if(c->Type==1) AParts << APatch;
            else AParts_c << APatch;
            
            
            // add lines between same cell types
            LinePatch *ALinePatch = new LinePatch(geom);
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                if((*it_e)->c1->Type == (*it_e)->c2->Type) ALinePatch->addLine((*it_e)->v1->coord/L,(*it_e)->v2->coord/L, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            }
            qSetColor(ALinePatch->faceColor, apicalLineColors[1][1], 1.0);
            ALineParts << ALinePatch;
            
            // add lines between different cell types
            LinePatch *ALinePatch_d = new LinePatch(geom);
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                if((*it_e)->c1->Type != (*it_e)->c2->Type) ALinePatch_d->addLine((*it_e)->v1->coord/L,(*it_e)->v2->coord/L, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            }
            
            qSetColor(ALinePatch_d->faceColor, apicalLineColors[2][1], 1.0);
            ALineParts << ALinePatch_d;
            
            // add basal surface patches
            SurfacePatch *BPatch = new SurfacePatch(geom);
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                Edge *e = (*it_e);
                
                // calculate the normalized coordinates
                Point v1_b_norm = e->v1->coord_b/L;
                Point v2_b_norm = e->v2->coord_b/L;
                
                BPatch->addTri(BC_b_norm,v1_b_norm,v2_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
                BPatch->addTri(v1_b_norm,BC_b_norm,v2_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            }
            qSetColor(BPatch->faceColor, basalPatchColors[c->Type],basalPatchOpacities[c->Type]);
            BParts << BPatch;
            
            // add basal lines
            LinePatch *BLinePatch = new LinePatch(geom);
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                BLinePatch->addLine((*it_e)->v1->coord_b/L,(*it_e)->v2->coord_b/L, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            }
            qSetColor(BLinePatch->faceColor, basalLineColors[1][1], 0.0);
            BLineParts << BLinePatch;
        }
        else { // if the cell crosses the boundary
            
            Cell* c = &(*it_c);
            
            // define surface patches for apical and basal surfaces
            SurfacePatch *APatch = new SurfacePatch(geom);
            
            // add the apical side
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                Edge *e = (*it_e);
                
                // add the triangles (normal is defined by the direction towards the cell)
                if(e->c1==c)
                {
                    // calculate the normalized coordinates
                    Point v1_a_norm = shift(e->v1->coord,T->SystemSize*e->q_c1*(-1))/L;
                    Point v2_a_norm = shift(e->v2->coord,T->SystemSize*(e->q_v2-e->q_c1))/L;
                    Point BC_a_norm = c->BC_a/L;
                    APatch->addTri(v1_a_norm,BC_a_norm,v2_a_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
                    APatch->addTri(v1_a_norm,v2_a_norm,BC_a_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
                }
                else
                {
                    Point v1_a_norm = shift(e->v1->coord,T->SystemSize*(e->q_c2+e->q_v2)*(-1))/L;
                    Point v2_a_norm = shift(e->v2->coord,T->SystemSize*(e->q_c2)*(-1))/L;
                    Point BC_a_norm = c->BC_a/L;
                    APatch->addTri(v1_a_norm,BC_a_norm,v2_a_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
                    APatch->addTri(v1_a_norm,v2_a_norm,BC_a_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
                }
                
            }
            
            qSetColor(APatch->faceColor, apicalPatchColors[c->Type],apicalPatchOpacities[c->Type]);
            
            AParts << APatch;
            
            LinePatch *ALinePatch = new LinePatch(geom);
            
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                Edge *e = (*it_e);
                
                if(e->c1==c)
                {
                    // calculate the normalized coordinates
                    Point v1_a_norm = shift(e->v1->coord,T->SystemSize*e->q_c1*(-1))/L;
                    Point v2_a_norm = shift(e->v2->coord,T->SystemSize*(e->q_v2-e->q_c1))/L;
                    ALinePatch->addLine(v1_a_norm,v2_a_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
                }
                else
                {
                    Point v1_a_norm = shift(e->v1->coord,T->SystemSize*(e->q_c2+e->q_v2)*(-1))/L;
                    Point v2_a_norm = shift(e->v2->coord,T->SystemSize*(e->q_c2)*(-1))/L;
                    ALinePatch->addLine(v1_a_norm,v2_a_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
                }
            }
            
            qSetColor(ALinePatch->faceColor, apicalLineColors[c->Type][c->Type], 0.0);
            ALineParts << ALinePatch;
            
            
            SurfacePatch *BPatch = new SurfacePatch(geom);
            // add the basal side
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                Edge *e = (*it_e);
                
                if(e->c1==c)
                {
                    // calculate the normalized coordinates
                    Point v1_b_norm = shift(e->v1->coord_b,T->SystemSize*(e->v1->q_b-e->q_c1))/L;
                    Point v2_b_norm = shift(e->v2->coord_b,T->SystemSize*(e->v2->q_b+e->q_v2-e->q_c1))/L;
                    Point BC_b_norm = shift(c->BC_b,T->SystemSize*c->q_BC_b)/L;
                    BPatch->addTri(v1_b_norm,BC_b_norm,v2_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
                    BPatch->addTri(v1_b_norm,v2_b_norm,BC_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
                }
                else
                {
                    Point v1_b_norm = shift(e->v1->coord_b,T->SystemSize*(e->v1->q_b-e->q_c2-e->q_v2))/L;
                    Point v2_b_norm = shift(e->v2->coord_b,T->SystemSize*(e->v2->q_b-e->q_c2))/L;
                    Point BC_b_norm = shift(c->BC_b,T->SystemSize*c->q_BC_b)/L;
                    BPatch->addTri(v1_b_norm,BC_b_norm,v2_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
                    BPatch->addTri(v1_b_norm,v2_b_norm,BC_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
                }
                
            }
            
            
            qSetColor(BPatch->faceColor, basalPatchColors[c->Type],basalPatchOpacities[c->Type]);
            BParts << BPatch;
            
            LinePatch *BLinePatch = new LinePatch(geom);
            
            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
            {
                Edge *e = (*it_e);
                
                if(e->c1==c)
                {
                    // calculate the normalized coordinates
                    Point v1_b_norm = shift(e->v1->coord_b,T->SystemSize*(e->v1->q_b-e->q_c1))/L;
                    Point v2_b_norm = shift(e->v2->coord_b,T->SystemSize*(e->q_v2+e->v2->q_b-e->q_c1))/L;
                    BLinePatch->addLine(v1_b_norm,v2_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
                }
                else
                {
                    Point v1_b_norm = shift(e->v1->coord_b,T->SystemSize*(e->v1->q_b-e->q_c2-e->q_v2))/L;
                    Point v2_b_norm = shift(e->v2->coord_b,T->SystemSize*(e->v2->q_b-e->q_c2))/L;
                    BLinePatch->addLine(v1_b_norm,v2_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
                }
            }
            
            qSetColor(BLinePatch->faceColor, basalLineColors[c->Type][c->Type], 0.0);
            BLineParts << BLinePatch;
            
        }
        
    }
    
    
    // iterate through the edges deepest inside which are supposed to be visible i.e. the ones at mut-wt interface
    for(std::list<Edge>::iterator it_e = T->ListEdge.begin(); it_e!=T->ListEdge.end(); ++it_e)
    {
        if(it_e->c1->Type==it_e->c2->Type) continue;
        
        if(!T->isPeriodic)
        {
            // normalize the coordinates
            Point v1_a_norm = it_e->v1->coord/L;
            Point v2_a_norm = it_e->v2->coord/L;
            Point v1_b_norm = it_e->v1->coord_b/L;
            Point v2_b_norm = it_e->v2->coord_b/L;
            Point b_norm = it_e->BC_l/L;
            
            SurfacePatch *LPatch = new SurfacePatch(geom);
            
            // add the four triangles to the patch
            LPatch->addTri(v1_a_norm,v2_a_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch->addTri(v2_a_norm,v2_b_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch->addTri(v2_b_norm,v1_b_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch->addTri(v1_b_norm,v1_a_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            
            qSetColor(LPatch->faceColor, lateralPatchColors[it_e->c1->Type][it_e->c2->Type],lateralPatchOpacities[it_e->c1->Type][it_e->c2->Type]);
            LParts << LPatch;
            
            // draw the vertex connection
            LinePatch *LLinePatch = new LinePatch(geom);
            LLinePatch->addLine(v1_a_norm,v1_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LLinePatch->addLine(v2_a_norm,v2_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            
            qSetColor(LLinePatch->faceColor, QColor(Qt::black),0.0);
            LLineParts_cyst << LLinePatch;
            
        }
        else {
            // draw the edge on the side of cell 1
            // normalize the coordinates
            Point v1_a_norm = shift(it_e->v1->coord, T->SystemSize*it_e->q_c1*(-1))/L;
            Point v2_a_norm = shift(it_e->v2->coord, T->SystemSize*(it_e->q_v2-it_e->q_c1))/L;
            Point v1_b_norm = shift(it_e->v1->coord_b, T->SystemSize*(it_e->v1->q_b - it_e->q_c1))/L;
            Point v2_b_norm = shift(it_e->v2->coord_b, T->SystemSize*(it_e->q_v2 + it_e->v2->q_b - it_e->q_c1))/L;
            Point b_norm = shift(it_e->BC_l, T->SystemSize*(it_e->q_BC_l - it_e->q_c1))/L;
            
            // add the four triangles to the patch
            SurfacePatch *LPatch_c1 = new SurfacePatch(geom);
            LPatch_c1->addTri(v1_a_norm,v2_a_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c1->addTri(v1_a_norm,b_norm,v2_a_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c1->addTri(v2_a_norm,v2_b_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c1->addTri(v2_a_norm,b_norm,v2_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c1->addTri(v2_b_norm,v1_b_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c1->addTri(v2_b_norm,b_norm,v1_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c1->addTri(v1_b_norm,v1_a_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c1->addTri(v1_b_norm,b_norm,v1_a_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            
            qSetColor(LPatch_c1->faceColor, lateralPatchColors[it_e->c1->Type][it_e->c2->Type],lateralPatchOpacities[it_e->c1->Type][it_e->c2->Type]);
            LParts << LPatch_c1;
            
            // draw the vertex connection
            LinePatch *LLinePatch_c1 = new LinePatch(geom);
            LLinePatch_c1->addLine(v1_a_norm,v1_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LLinePatch_c1->addLine(v2_a_norm,v2_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            
            qSetColor(LLinePatch_c1->faceColor, QColor(Qt::black), 1);
            LLineParts_cyst << LLinePatch_c1;
            
            // draw the edge on the side of cell 2
            // normalize the coordinates
            v2_a_norm = shift(it_e->v2->coord, T->SystemSize*it_e->q_c2*(-1))/L;
            v1_a_norm = shift(it_e->v1->coord, T->SystemSize*(it_e->q_v2+it_e->q_c2)*(-1))/L;
            v2_b_norm = shift(it_e->v2->coord_b, T->SystemSize*(it_e->v2->q_b - it_e->q_c2))/L;
            v1_b_norm = shift(it_e->v1->coord_b, T->SystemSize*(it_e->v1->q_b - it_e->q_v2 - it_e->q_c2))/L;
            b_norm = shift(it_e->BC_l, T->SystemSize*(it_e->q_BC_l - it_e->q_c2 - it_e->q_v2))/L;
            
            // add the four triangles to the patch
            SurfacePatch *LPatch_c2 = new SurfacePatch(geom);
            LPatch_c2->addTri(v1_a_norm,v2_a_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c2->addTri(v2_a_norm,v2_b_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c2->addTri(v2_b_norm,v1_b_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c2->addTri(v1_b_norm,v1_a_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            
            LPatch_c2->addTri(v1_a_norm,b_norm,v2_a_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c2->addTri(v2_a_norm,b_norm,v2_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c2->addTri(v2_b_norm,b_norm,v1_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c2->addTri(v1_b_norm,b_norm,v1_a_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            
            
            qSetColor(LPatch_c2->faceColor, lateralPatchColors[it_e->c1->Type][it_e->c2->Type],lateralPatchOpacities[it_e->c1->Type][it_e->c2->Type]);
            LParts << LPatch_c2;
            
            // draw the vertex connection
            LinePatch *LLinePatch_c2 = new LinePatch(geom);
            LLinePatch_c2->addLine(v1_a_norm,v1_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LLinePatch_c2->addLine(v2_a_norm,v2_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            
            qSetColor(LLinePatch_c2->faceColor, QColor(Qt::black), 1);
            LLineParts_cyst << LLinePatch_c2;
        }
    }
    
    // draw opaque edges at boundary
    for(std::list<Edge>::iterator it_e = T->ListEdge.begin(); it_e!=T->ListEdge.end(); ++it_e)
    {
        // check which of the neighbours are in the radius
        bool c1InRadius = sqrt(pow(it_e->c1->BC_a.x-midPoint.x,2) + pow(it_e->c1->BC_a.y-midPoint.y,2))<radius;
        bool c2InRadius = sqrt(pow(it_e->c2->BC_a.x-midPoint.x,2) + pow(it_e->c2->BC_a.y-midPoint.y,2))<radius;
        
        // draw edge only if it is at interface between visible and invisible cells or at the boundary of the tissue
        if(c1InRadius==c2InRadius && !it_e->isAtBoundary()) continue;
        
        if(!T->isPeriodic)
        {
            // normalize the coordinates
            Point v1_a_norm = it_e->v1->coord/L;
            Point v2_a_norm = it_e->v2->coord/L;
            Point v1_b_norm = it_e->v1->coord_b/L;
            Point v2_b_norm = it_e->v2->coord_b/L;
            Point b_norm = it_e->BC_l/L;
            
            SurfacePatch *LPatch = new SurfacePatch(geom);
            
            // add the four triangles to the patch
            LPatch->addTri(v1_a_norm,v2_a_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch->addTri(v2_a_norm,v2_b_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch->addTri(v2_b_norm,v1_b_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch->addTri(v1_b_norm,v1_a_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            
            qSetColor(LPatch->faceColor, QColor(Qt::black),0.0);
            LParts << LPatch;
            
            // draw the vertex connection
            LinePatch *LLinePatch = new LinePatch(geom);
            LLinePatch->addLine(v1_a_norm,v1_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LLinePatch->addLine(v2_a_norm,v2_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            
            qSetColor(LLinePatch->faceColor, QColor(Qt::black),0.0);
            LLineParts << LLinePatch;
            
        }
        else {
            // draw the edge on the side of cell 1
            // normalize the coordinates
            Point v1_a_norm = shift(it_e->v1->coord, T->SystemSize*it_e->q_c1*(-1))/L;
            Point v2_a_norm = shift(it_e->v2->coord, T->SystemSize*(it_e->q_v2-it_e->q_c1))/L;
            Point v1_b_norm = shift(it_e->v1->coord_b, T->SystemSize*(it_e->v1->q_b - it_e->q_c1))/L;
            Point v2_b_norm = shift(it_e->v2->coord_b, T->SystemSize*(it_e->q_v2 + it_e->v2->q_b - it_e->q_c1))/L;
            Point b_norm = shift(it_e->BC_l, T->SystemSize*(it_e->q_BC_l - it_e->q_c1))/L;
            
            // add the four triangles to the patch
            SurfacePatch *LPatch_c1 = new SurfacePatch(geom);
            LPatch_c1->addTri(v1_a_norm,v2_a_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c1->addTri(v1_a_norm,b_norm,v2_a_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c1->addTri(v2_a_norm,v2_b_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c1->addTri(v2_a_norm,b_norm,v2_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c1->addTri(v2_b_norm,v1_b_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c1->addTri(v2_b_norm,b_norm,v1_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c1->addTri(v1_b_norm,v1_a_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c1->addTri(v1_b_norm,b_norm,v1_a_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            
            
            qSetColor(LPatch_c1->faceColor, lateralPatchColors[it_e->c1->Type][it_e->c1->Type], lateralPatchOpacities[it_e->c1->Type][it_e->c1->Type]);
            LParts << LPatch_c1;
            
            // draw the vertex connection
            LinePatch *LLinePatch_c1 = new LinePatch(geom);
            LLinePatch_c1->addLine(v1_a_norm,v1_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LLinePatch_c1->addLine(v2_a_norm,v2_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            
            qSetColor(LLinePatch_c1->faceColor, QColor(Qt::black), 0.0);
            LLineParts << LLinePatch_c1;
            
            // draw the edge on the side of cell 2
            // normalize the coordinates
            v2_a_norm = shift(it_e->v2->coord, T->SystemSize*it_e->q_c2*(-1))/L;
            v1_a_norm = shift(it_e->v1->coord, T->SystemSize*(it_e->q_v2+it_e->q_c2)*(-1))/L;
            v2_b_norm = shift(it_e->v2->coord_b, T->SystemSize*(it_e->v2->q_b - it_e->q_c2))/L;
            v1_b_norm = shift(it_e->v1->coord_b, T->SystemSize*(it_e->v1->q_b - it_e->q_v2 - it_e->q_c2))/L;
            b_norm = shift(it_e->BC_l, T->SystemSize*(it_e->q_BC_l - it_e->q_c2 - it_e->q_v2))/L;
            
            // add the four triangles to the patch
            SurfacePatch *LPatch_c2 = new SurfacePatch(geom);
            LPatch_c2->addTri(v1_a_norm,v2_a_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c2->addTri(v2_a_norm,v2_b_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c2->addTri(v2_b_norm,v1_b_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c2->addTri(v1_b_norm,v1_a_norm,b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            
            LPatch_c2->addTri(v1_a_norm,b_norm,v2_a_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c2->addTri(v2_a_norm,b_norm,v2_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c2->addTri(v2_b_norm,b_norm,v1_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LPatch_c2->addTri(v1_b_norm,b_norm,v1_a_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            
            qSetColor(LPatch_c2->faceColor, lateralPatchColors[it_e->c1->Type][it_e->c1->Type], lateralPatchOpacities[it_e->c1->Type][it_e->c1->Type]);
            LParts << LPatch_c2;
            
            // draw the vertex connection
            LinePatch *LLinePatch_c2 = new LinePatch(geom);
            LLinePatch_c2->addLine(v1_a_norm,v1_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            LLinePatch_c2->addLine(v2_a_norm,v2_b_norm, xr, yr, T->SystemSize.x/L, T->SystemSize.y/L);
            
            qSetColor(LLinePatch_c2->faceColor, QColor(Qt::black), 0.0);
            LLineParts << LLinePatch_c2;
        }
    }
    
    
    
}



void GLTissueImage::draw()
{
    
    int x_r = 1;
    int y_r = 1;
        
    geom->preloadArrays();
    
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    
    if(tissueChanged) {
        if(!integrator->T->onlyApicalContributions) loadPartsInRadius(10000000,x_r,y_r);//loadPartsInRadius(integrator->T->SystemSize.x/3.0);
        else loadPartsInRadius(10000000);//loadParts();
        tissueChanged=false;
    }


    // apical lines (edges)
    for (int i = 0; i < ALineParts.count(); ++i)
        ALineParts[i]->draw();

    
    // basal lines (edges)
    for (int i = 0; i < BLineParts.count(); ++i)
        BLineParts[i]->draw();
    
    if(!integrator->T->onlyApicalContributions){
    // lateral lines of vertices visible at the boundary of the tissue
    for (int i = 0; i < LLineParts_cyst.count(); ++i)
        LLineParts_cyst[i]->draw();
       
    // lateral surfaces of cells at the boundary
    for (int i = 0; i < LParts.count(); ++i)
        LParts[i]->draw();
     
    
    // basal surface
    for (int i = 0; i < BParts.count(); ++i)
        BParts[i]->draw();
    }

    // apical surface - mutant cells
    for (int i = 0; i < AParts_c.count(); ++i)
        AParts_c[i]->draw();
    
    // wild type cells
    for (int i = 0; i < AParts.count(); ++i)
        AParts[i]->draw();
    
    for (int i = 0; i < BoxLineParts.count(); ++i)
        BoxLineParts[i]->draw();

    for (int i = 0; i < BoxParts.count(); ++i)
        BoxParts[i]->draw();

    if(!integrator->T->onlyApicalContributions){
    // lateral lines of vertices visible at the boundary of the tissue
    for (int i = 0; i < LLineParts.count(); ++i)
        LLineParts[i]->draw();
    }
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);    
}


//! [1]
#endif // _USE_QT_