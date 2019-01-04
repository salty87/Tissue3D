/****************************************************************************
**
** Copyright (C) 2008 Nokia Corporation and/or its subsidiary(-ies).
** Contact: Qt Software Information (qt-info@nokia.com)
**
** This file is part of the example classes of the Qt Toolkit.
**
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

#ifndef TISSUEWINDOW_H
#define TISSUEWINDOW_H

#include <QWidget>
#include <QObject>
#include <QGraphicsView>
#include <QString>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QImage>
#include <QColormap>
#include <QList>
#include <list>
#include "Tissue.h"
#include "Edge.h"
#include "Cell.h"


// struct that contains the maximum stretch in x,y,z direction of the tissue
struct ExtPositionsTissue 
{
	Point BottomLeft;
	
	Point TopRight;
};

class  MainWindow;
class TissueWindow : public QWidget
{
    Q_OBJECT
    
public:
    /*  CONSTRUCTOR */
    TissueWindow(QWidget* parent, Integrator *integrator);
    
protected:
    /* PAINT EVENT */
    void paintEvent(QPaintEvent *event);
	
private:
    /* FUNCTIONS TO DRAW THE TISSUE */
    void drawTissue(QPainter& painter);
    
    /* DRAWING OPTIONS */
    QPen pen_a11, pen_a12, pen_a10, pen_a20, pen_a22;
	QPen pen_b11, pen_b12, pen_b10, pen_b20, pen_b22;
	QPen pen_l, pen_barycenter;
	bool drawBarycenters;
	bool drawApicalSide;
	bool drawBasalSide;
    bool drawRectangleForPeriodicTissue;
	const static int NbTypes=5; // number of cell types
	std::string ApicalColors[NbTypes]; // colours for the different cell types
    
    /* VIEW OPTIONS */
	QTransform* tMatrix;
	QMatrix matrix;
	double scaleFactor;
	
    /* INTEGRATOR INFORMATION */
	Integrator *integrator;
    
    /* HIGHLIGHTING OF CHOSEN CELLS/VERTICES/EDGES */
    void highlightCell(QPainter& painter,QPainter& painter_b, Cell* c);
	void highlightEdge(QPainter& painter,QPainter& painter_b, Edge* e);
	void highlightVertex(QPainter& painter,QPainter& painter_b, Vertex* v);
    
    /* SAVE CURRENT IMAGE */
    QImage CurrentImage;
    void writeCurrentImage();
	void writeCurrentPDFImage(const QString &fileName);//dessine le tissu dans CurrentPDFImage
    void recordImageFile();
    int nbImagesRecorded;
    /* OTHER VARIABLES */
	bool ProfileIsActivated; // if the cross section is shown, then the line should also be dispayed
	Point ClickedPoint; // point that has currently been clicked at in the tissue window
	int mainWindow_tab;
    Point distance;
    
    
public slots:
	void mousePressEvent(QMouseEvent *event);
	void keyPressEvent(QKeyEvent *event);
    void slot_drawingModeChanged(int, int, int);
    void slot_update(){update(); /*repaint();*/ QApplication::processEvents();}
    void slot_tabChanged(int newTab) {mainWindow_tab=newTab;}
signals:
	void signal_currentPointChanged(double,double);
    void signal_addCurrentPoint(double,double);
    void signal_addToCurrentPoint(double,double);
    void signal_currentPoint(double,double);
    void signal_printCellData();
};
#include "MainWindow.h"
#endif

#endif // _USE_QT_