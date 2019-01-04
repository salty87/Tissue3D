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
    QVector<GLushort> faces;
    QVector<GLushort> lines;
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
    
    GLushort start;
    GLushort count;
    GLushort initv;
    
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

static inline void qSetColor(float colorVec[], QColor c)
{
    colorVec[0] = c.redF();
    colorVec[1] = c.greenF();
    colorVec[2] = c.blueF();
    colorVec[3] = c.alphaF();
}


SurfacePatch::SurfacePatch(Geometry *g)
: start(g->faces.count()), count(0), initv(g->vertices.count()), geom(g)
{
    qSetColor(faceColor, QColor(Qt::white),0.0);
}

//! [2]
void SurfacePatch::draw() const
{
    glDisable(GL_CULL_FACE);
    //glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,1.);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, faceColor);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, faceColor);
    
    //glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, Qt::white);
    const GLushort *indices = geom->faces.constData();
    glDrawElements(GL_TRIANGLES, count, GL_UNSIGNED_SHORT, indices + start);
    
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

class LinePatch
{
public: 
    LinePatch(Geometry *);
    void draw() const;
    void addLine(const Point &pa, const Point &pb);
    
    GLushort start;
    GLushort count;
    GLushort initv;
    
    GLfloat faceColor[4]; //vector of the colours of the line
    Geometry *geom;
};

LinePatch::LinePatch(Geometry *g)
: start(g->faces.count()), count(0), initv(g->vertices.count()), geom(g)
{
    qSetColor(faceColor, QColor(Qt::white),0.0);
}

//! [2]
void LinePatch::draw() const
{
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, faceColor);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, faceColor);
    
    const GLushort *indices = geom->faces.constData();
    glDrawElements(GL_LINES, count, GL_UNSIGNED_SHORT, indices + start);
    
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

GLTissueImage::GLTissueImage(QObject *parent, Integrator *integrator)
: QObject(parent), integrator(integrator), tissueChanged(true)
, geom(new Geometry())
{
}

GLTissueImage::~GLTissueImage()
{
    qDeleteAll(AParts);
    qDeleteAll(BParts);
    qDeleteAll(LParts);
    qDeleteAll(ALineParts);
    qDeleteAll(BLineParts);
    qDeleteAll(LLineParts);
    delete geom;
}

void GLTissueImage::loadParts()
{
    
    QColor dummy;
    //QColor apicalPatchColors [] = { QColor(Qt::gray), QColor(Qt::magenta), QColor(Qt::red), QColor(Qt::darkRed), QColor(Qt::darkYellow) };
    
    QColor apicalType1; apicalType1.setNamedColor("red");
    QColor apicalType2; apicalType2.setNamedColor("blue");
    QColor basalType1;  basalType1.setNamedColor("blue");
    QColor basalType2;  basalType2.setNamedColor("purple");
    
    QColor apicalPatchColors [] = { QColor(Qt::gray), apicalType1, apicalType2, QColor(Qt::blue), QColor(Qt::darkGray), };
    QColor basalPatchColors  [] = { QColor(Qt::gray), basalType1 , basalType2, QColor(Qt::cyan), QColor(Qt::darkCyan) };
    QColor lateralPatchColor(Qt::gray);
    
    QColor apicalEdgeColor_sameCellType(Qt::black);
    QColor apicalEdgeColor_differentCellType(Qt::green);
    QColor basalEdgeColor_sameCellType(Qt::black);
    QColor basalEdgeColor_differentCellType(Qt::green);
    
    int numberOfCellTypes = 4;
    
    if(1){
    //if(integrator->T->isPeriodic) {
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
            
        }
        
        for(std::list<Cell>::iterator it_c=T->ListCell.begin(); it_c!=T->ListCell.end(); ++it_c)
        {
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
                    
                    //std::cout << "(" << BC_a.x/T->SystemSize.x << "," << BC_a.y/T->SystemSize.y << "," << BC_a.z << ")" << std::endl;
                    
                    // calculate the normalized coordinates
                    Point v1_a_norm = e->v1->coord/L;
                    Point v2_a_norm = e->v2->coord/L;
                    
                    // add the triangles (normal is defined by the direction towards the cell)
                    if(e->c1==c)
                    {
                        APatch->addTri(v1_a_norm,BC_a_norm,v2_a_norm);
                    }
                    else
                    {
                        APatch->addTri(BC_a_norm,v1_a_norm,v2_a_norm);
                    }
                    
                }
                
                
                if(c->Type<=numberOfCellTypes && c->Type > 0) qSetColor(APatch->faceColor, apicalPatchColors[c->Type],0.5);
                else                         qSetColor(APatch->faceColor, apicalPatchColors[0],0.5);
                
                AParts << APatch;
                
                LinePatch *ALinePatch = new LinePatch(geom);
                
                for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
                {
                    ALinePatch->addLine((*it_e)->v1->coord/L,(*it_e)->v2->coord/L);
                }
                
                qSetColor(ALinePatch->faceColor, apicalEdgeColor_sameCellType,0.0);
                ALineParts << ALinePatch;
                
                
                SurfacePatch *BPatch = new SurfacePatch(geom);
                // add the basal side
                for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
                {
                    Edge *e = (*it_e);
                    
                    // calculate the normalized coordinates
                    Point v1_b_norm = e->v1->coord_b/L;
                    Point v2_b_norm = e->v2->coord_b/L;
                    
                    // add the triangles (normal is defined by the direction towards the cell)
                    if(e->c1==c)
                    {
                        BPatch->addTri(BC_b_norm,v1_b_norm,v2_b_norm);
                    }
                    else
                    {
                        BPatch->addTri(v1_b_norm,BC_b_norm,v2_b_norm);
                    }
                    
                }
                
                if(c->Type<=numberOfCellTypes && c->Type > 0)    qSetColor(BPatch->faceColor, basalPatchColors[c->Type],0.5);
                else                                            qSetColor(BPatch->faceColor, basalPatchColors[0],0.5);
                
                BParts << BPatch;
                
                LinePatch *BLinePatch = new LinePatch(geom);
                
                for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
                {
                    BLinePatch->addLine((*it_e)->v1->coord_b/L,(*it_e)->v2->coord_b/L);
                }
                
                qSetColor(ALinePatch->faceColor, basalEdgeColor_sameCellType);
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
                    }
                    else
                    {
                        Point v1_a_norm = shift(e->v1->coord,T->SystemSize*(e->q_c2+e->q_v2)*(-1))/L;
                        Point v2_a_norm = shift(e->v2->coord,T->SystemSize*(e->q_c2)*(-1))/L;
                        Point BC_a_norm = c->BC_a/L;
                        APatch->addTri(v1_a_norm,BC_a_norm,v2_a_norm);
                    }
                    
                }
                
                if(c->Type<=numberOfCellTypes && c->Type > 0) qSetColor(APatch->faceColor, apicalPatchColors[c->Type],0.5);
                else                         qSetColor(APatch->faceColor, apicalPatchColors[0],0.5);

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
                
                qSetColor(ALinePatch->faceColor, apicalEdgeColor_sameCellType, 0.0);
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
                    }
                    else
                    {
                        Point v1_b_norm = shift(e->v1->coord_b,T->SystemSize*(e->v1->q_b-e->q_c2-e->q_v2))/L;
                        Point v2_b_norm = shift(e->v2->coord_b,T->SystemSize*(e->v2->q_b-e->q_c2))/L;
                        Point BC_b_norm = shift(c->BC_b,T->SystemSize*c->q_BC_b)/L;
                        BPatch->addTri(v1_b_norm,BC_b_norm,v2_b_norm);
                    }
                    
                }
                
                if(c->Type<=numberOfCellTypes && c->Type > 0) qSetColor(BPatch->faceColor, basalPatchColors[c->Type],0.0);
                else                                         qSetColor(BPatch->faceColor, basalPatchColors[0],0.0);
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
                
                qSetColor(BLinePatch->faceColor, apicalEdgeColor_sameCellType,0.0);
                BLineParts << BLinePatch;
                
            }
            
        }
        
        // iterate through the edges and draw all, that are at the boundary i.e. can be seen
        for(std::list<Edge>::iterator it_e = T->ListEdge.begin(); it_e!=T->ListEdge.end(); ++it_e)
        {
            
            if(it_e->isAtBoundary())
            {
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
                    
                    qSetColor(LPatch->faceColor, QColor(Qt::white),0.0);
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
                    LPatch_c1->addTri(v2_a_norm,v2_b_norm,b_norm);
                    LPatch_c1->addTri(v2_b_norm,v1_b_norm,b_norm);
                    LPatch_c1->addTri(v1_b_norm,v1_a_norm,b_norm);
                    
                    qSetColor(LPatch_c1->faceColor, lateralPatchColor, 0.3);
                    LParts << LPatch_c1;
                    
                    // draw the vertex connection
                    LinePatch *LLinePatch_c1 = new LinePatch(geom);
                    LLinePatch_c1->addLine(v1_a_norm,v1_b_norm);
                    LLinePatch_c1->addLine(v2_a_norm,v2_b_norm);
                    
                    qSetColor(LLinePatch_c1->faceColor, QColor(Qt::black), 0.3);
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
                    
                    qSetColor(LPatch_c2->faceColor, lateralPatchColor, 0.0);
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
    } else { // periodic tissue
        // clear the list of parts
        removeParts();
        Tissue* T = integrator->T;
        
        // update the systems outermost points
        T->UpdateExtPositions();
        
        double Lx = T->extPosition_topRight.x-T->extPosition_bottomLeft.x;
        double Ly = T->extPosition_topRight.y-T->extPosition_bottomLeft.y;
        double Lz = T->extPosition_topRight.z-T->extPosition_bottomLeft.z;
        double L = std::max(Lx, std::max(Ly,Lz)); // biggest size in the system
        
        for(std::list<Edge>::iterator it_e=T->ListEdge.begin(); it_e!=T->ListEdge.end(); ++it_e) { // iterate through all edges and draw them one after another
            
            Edge *e = &(*it_e);
            // normalized coordinates
            Point v1a = e->v1->coord/L;
            Point v2a = e->v2->coord/L;
            Point v1b = e->v1->coord_b/L;
            Point v2b = e->v2->coord_b/L;
            Point BC_l = e->BC_l/L;
            
            // define the patch that contains the lateral side
            SurfacePatch *LPatch = new SurfacePatch(geom); // patch for the lateral surfaces
            
            // add the four triangles to the patch
            LPatch->addTri(v1a,v2a,BC_l);
            LPatch->addTri(v2a,v2b,BC_l);
            LPatch->addTri(v2b,v1b,BC_l);
            LPatch->addTri(v1b,v1a,BC_l);
            
            // set color of lateral interface
            if(e->c1 && e->c2) { // if the cell has two neighbours
                if(e->c1->Type == e->c2->Type) { // both cell types are equal
                    qSetColor(LPatch->faceColor,lateralPatchColor,0.3);
                } else {
                    qSetColor(LPatch->faceColor, lateralPatchColor,0.3);
                }
            } else { // edge is a free boundary
                qSetColor(LPatch->faceColor, QColor(Qt::lightGray),0.3);
            }
            
            LParts << LPatch;
            LinePatch *LLinePatch = new LinePatch(geom);
            LLinePatch->addLine(v1a,v1b);
            LLinePatch->addLine(v2a,v2b);

            qSetColor(LLinePatch->faceColor, QColor(Qt::black),0.0);
            LLineParts << LLinePatch;

            // now draw apical and basal contributions of cell 1 if it exists
            if(e->c1) { // check if c1 exists
                // define the surface patches
                SurfacePatch *APatch = new SurfacePatch(geom);
                SurfacePatch *BPatch = new SurfacePatch(geom);

                Point BC_a = e->c1->BC_a/L;
                Point BC_b = e->c1->BC_b/L;
                
                // add the triangles (normal is defined by the direction towards the cell)
                APatch->addTri(v1a,BC_a,v2a);
                BPatch->addTri(BC_b,v1b,v2b);
                // set color of the apical side
                if(e->c1->Type==1) {
                    qSetColor(APatch->faceColor, QColor(Qt::red));
                    qSetColor(BPatch->faceColor, QColor(Qt::blue));
                } else { //
                    qSetColor(APatch->faceColor, QColor(Qt::darkRed));
                    qSetColor(BPatch->faceColor, QColor(Qt::darkBlue));
                }
                
                //AParts << APatch;
                //BParts << BPatch;
               
                // draw the lines
                LinePatch *ALinePatch = new LinePatch(geom);                
                ALinePatch->addLine(v1a,v2a);
                qSetColor(ALinePatch->faceColor, QColor(Qt::red));
                ALineParts << ALinePatch;
                
                LinePatch *BLinePatch = new LinePatch(geom);
                BLinePatch->addLine(v1b,v2b);
                qSetColor(BLinePatch->faceColor, QColor(Qt::blue));
                BLineParts << BLinePatch;

                ALinePatch = new LinePatch(geom);
                ALinePatch->addLine(v1a,BC_a);
                ALinePatch->addLine(v2a,BC_a);
                qSetColor(ALinePatch->faceColor, QColor(Qt::darkGray));
                
                BLinePatch = new LinePatch(geom);
                BLinePatch->addLine(v1b,BC_b);
                BLinePatch->addLine(v2b,BC_b);
                qSetColor(BLinePatch->faceColor, QColor(Qt::darkGray));
            }
            // ... and apical and basal contributions of cell 2 if it exists
            if(e->c2) { // check if c1 exists
                // define the surface patches
                SurfacePatch *APatch = new SurfacePatch(geom);
                SurfacePatch *BPatch = new SurfacePatch(geom);
                
                Point BC_a = e->c2->BC_a/L;
                Point BC_b = e->c2->BC_b/L;
                
                // add the triangles (normal is defined by the direction towards the cell)
                APatch->addTri(v1a,v2a,BC_a);
                BPatch->addTri(v1b,BC_b,v2b);
                
                // set color of the apical side
                if(e->c2->Type==1) {
                    qSetColor(APatch->faceColor, QColor(Qt::red));
                    qSetColor(BPatch->faceColor, QColor(Qt::blue));
                } else { //
                    qSetColor(APatch->faceColor, QColor(Qt::darkRed));
                    qSetColor(BPatch->faceColor, QColor(Qt::darkBlue));
                }
                
                //AParts << APatch;
                //BParts << BPatch;
                
                // draw the lines
                LinePatch *ALinePatch = new LinePatch(geom);                                
                ALinePatch->addLine(v1a,v2a);
                qSetColor(ALinePatch->faceColor, QColor(Qt::red));
                ALineParts << ALinePatch;

                LinePatch *BLinePatch = new LinePatch(geom);
                BLinePatch->addLine(v1b,v2b);
                qSetColor(BLinePatch->faceColor, QColor(Qt::blue));
                BLineParts << BLinePatch;
                
                ALinePatch = new LinePatch(geom);
                ALinePatch->addLine(v1a,BC_a);
                ALinePatch->addLine(v2a,BC_a);
                qSetColor(ALinePatch->faceColor, QColor(Qt::darkGray));
                //ALineParts << ALinePatch;
                
                BLinePatch = new LinePatch(geom);
                BLinePatch->addLine(v1b,BC_b);
                BLinePatch->addLine(v2b,BC_b);
                qSetColor(BLinePatch->faceColor, QColor(Qt::darkGray));
                //BLineParts << BLinePatch;
            }
        }
        
    }
}

//            
//        for(std::list<Cell>::iterator it_c=T->ListCell.begin(); it_c!=T->ListCell.end(); ++it_c)
//        {
//            // normalize the variables
//            Cell* c = &(*it_c);
//            Point BC_a = c->BC_a;
//            Point BC_a_norm = BC_a/L;
//            
//            Point BC_b = c->BC_b;
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
//                if(e->c1==c)
//                {
//                    APatch->addTri(v1_a_norm,BC_a_norm,v2_a_norm);
//                }
//                else
//                {
//                    APatch->addTri(BC_a_norm,v1_a_norm,v2_a_norm);
//                }
//                
//            }
//            
//            
//            if(c->Type==1)
//            {
//                qSetColor(APatch->faceColor, QColor(Qt::red));
//            }
//            else
//            {
//                qSetColor(APatch->faceColor, QColor(Qt::red).dark());
//            }
//            
//            //AParts << APatch;
//            
//            LinePatch *ALinePatch = new LinePatch(geom);
//            //LinePatch *ABCLinePatch = new LinePatch(geom);
//            //LinePatch *BBCLinePatch = new LinePatch(geom);
//            
//            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
//            {
//                ALinePatch->addLine((*it_e)->v1->coord/L,(*it_e)->v2->coord/L);
//                ALinePatch->addLine((*it_e)->v1->coord/L,c->BC_a/L);
//                ALinePatch->addLine((*it_e)->v2->coord/L,c->BC_a/L);
//                //BBCLinePatch->addLine((*it_e)->v1->coord_b/L,c->BC_b/L);
//                //BBCLinePatch->addLine((*it_e)->v2->coord_b/L,c->BC_b/L);
//            }
//            
//            qSetColor(ALinePatch->faceColor, QColor(Qt::red));
//            //qSetColor(ABCLinePatch->faceColor, QColor(Qt::black));
//            //qSetColor(BBCLinePatch->faceColor, QColor(Qt::black));
//            ALineParts << ALinePatch;
//            //ALineParts << ABCLinePatch;
//            
//            SurfacePatch *BPatch = new SurfacePatch(geom);
//            // add the basal side
//            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
//            {
//                Edge *e = (*it_e);
//                
//                // calculate the normalized coordinates
//                Point v1_b_norm = e->v1->coord_b/L;
//                Point v2_b_norm = e->v2->coord_b/L;
//                
//                // add the triangles (normal is defined by the direction towards the cell)
//                if(e->c1==c)
//                {
//                    BPatch->addTri(BC_b_norm,v1_b_norm,v2_b_norm);
//                }
//                else
//                {
//                    BPatch->addTri(v1_b_norm,BC_b_norm,v2_b_norm);
//                }
//                
//            }
//            
//            if(c->Type==1)
//            {
//                qSetColor(BPatch->faceColor, QColor(Qt::blue));
//            }
//            else
//            {
//                qSetColor(BPatch->faceColor, QColor(Qt::blue));
//            }
//            
//            //BParts << BPatch;
//            
//            LinePatch *BLinePatch = new LinePatch(geom);
//            
//            for(std::list<Edge*>::iterator it_e = c->ListEdges.begin(); it_e!=c->ListEdges.end();++it_e)
//            {
//                BLinePatch->addLine((*it_e)->v1->coord_b/L,(*it_e)->v2->coord_b/L);
//            }
//            
//            qSetColor(BLinePatch->faceColor, QColor(Qt::blue));
//            BLineParts << BLinePatch;
//        }
//    
//        
//        // iterate through the edges and draw all, that are at the boundary i.e. can be seen
//        for(std::list<Edge>::iterator it_e = T->ListEdge.begin(); it_e!=T->ListEdge.end(); ++it_e)
//        {
//            
//            Point v1_a_norm, v2_a_norm, v1_b_norm, v2_b_norm, b_norm;
//        
//            // normalize the coordinates
//            v1_a_norm = it_e->v1->coord/L;
//            v2_a_norm = it_e->v2->coord/L;
//            v1_b_norm = it_e->v1->coord_b/L;
//            v2_b_norm = it_e->v2->coord_b/L;
//            b_norm = it_e->BC_l/L;
//                
//            SurfacePatch *LPatch = new SurfacePatch(geom);
//            
//            // add the four triangles to the patch
//            LPatch->addTri(v1_a_norm,v2_a_norm,b_norm);
//            LPatch->addTri(v2_a_norm,v2_b_norm,b_norm);
//            LPatch->addTri(v2_b_norm,v1_b_norm,b_norm);
//            LPatch->addTri(v1_b_norm,v1_a_norm,b_norm);
//            
//            qSetColor(LPatch->faceColor, QColor(Qt::green));
//            LParts << LPatch;
//            
//            // draw the vertex connection
//            LinePatch *LLinePatch = new LinePatch(geom);
//            LLinePatch->addLine(v1_a_norm,v1_b_norm);
//            LLinePatch->addLine(v2_a_norm,v2_b_norm);
//            
//            qSetColor(LLinePatch->faceColor, QColor(Qt::black));
//            LLineParts << LLinePatch;
//            
//            // draw the edge on the side of cell 2
//            // normalize the coordinates
//            v2_a_norm = shift(it_e->v2->coord, T->SystemSize*it_e->q_c2*(-1))/L;
//            v1_a_norm = shift(it_e->v1->coord, T->SystemSize*(it_e->q_v2+it_e->q_c2)*(-1))/L;
//            v2_b_norm = shift(it_e->v2->coord_b, T->SystemSize*(it_e->v2->q_b - it_e->q_c2))/L;
//            v1_b_norm = shift(it_e->v1->coord_b, T->SystemSize*(it_e->v1->q_b - it_e->q_v2 - it_e->q_c2))/L;
//            b_norm = shift(it_e->BC_l, T->SystemSize*(it_e->q_BC_l - it_e->q_c2 - it_e->q_v2))/L;
//            
//            // add the four triangles to the patch
//            SurfacePatch *LPatch_c2 = new SurfacePatch(geom);
//            LPatch_c2->addTri(v1_a_norm,v2_a_norm,b_norm);
//            LPatch_c2->addTri(v2_a_norm,v2_b_norm,b_norm);
//            LPatch_c2->addTri(v2_b_norm,v1_b_norm,b_norm);
//            LPatch_c2->addTri(v1_b_norm,v1_a_norm,b_norm);
//            
//            qSetColor(LPatch_c2->faceColor, QColor(Qt::green).dark());
//            LParts << LPatch_c2;
//            
//            // draw the vertex connection
//            LinePatch *LLinePatch_c2 = new LinePatch(geom);
//            LLinePatch_c2->addLine(v1_a_norm,v1_b_norm);
//            LLinePatch_c2->addLine(v2_a_norm,v2_b_norm);
//            
//            qSetColor(LLinePatch_c2->faceColor, QColor(Qt::black));
//            LLineParts << LLinePatch_c2;
//        }
//    
//    }
//}

void GLTissueImage::removeParts()
{
    
    // delete all the patches
    AParts.clear();
    BParts.clear();
    LParts.clear();
    ALineParts.clear();
    BLineParts.clear();
    LLineParts.clear();
    
    BoxLineParts.clear();
    
    geom->vertices.clear();
    geom->faces.clear();
    geom->lines.clear();
    geom->normals.clear();
}


void GLTissueImage::draw()
{
        
    geom->preloadArrays();
    
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    
    if(tissueChanged) {loadParts(); tissueChanged=false;}
    
    glEnable(GL_BLEND);
    glBlendFunc(GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
    
    // apical surface
    for (int i = 0; i < AParts.count(); ++i)
        AParts[i]->draw();

    // basal surface
    for (int i = 0; i < BParts.count(); ++i)
        BParts[i]->draw();
    
    // lateral surfaces of cells at the boundary
    for (int i = 0; i < LParts.count(); ++i)
        LParts[i]->draw();
    
    // apical lines (edges)
    for (int i = 0; i < ALineParts.count(); ++i)
        ALineParts[i]->draw();
    
    // basal lines (edges)
    for (int i = 0; i < BLineParts.count(); ++i)
        BLineParts[i]->draw();
    
    // lateral lines of vertices visible at the boundary of the tissue
    for (int i = 0; i < LLineParts.count(); ++i)
        LLineParts[i]->draw();
    
    for (int i = 0; i < BoxLineParts.count(); ++i)
        BoxLineParts[i]->draw();
    
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    
    
}


//! [1]
#endif // _USE_QT_