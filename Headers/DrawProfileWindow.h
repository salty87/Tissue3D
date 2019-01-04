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

#ifndef DrawProfileWindow_H
#define DrawProfileWindow_H


#include <list>
#include <QApplication>
#include <QWidget>
#include <QObject>
#include <QGraphicsView>
#include <QString>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QImage>
#include <QColormap>
#include <QList>
#include "Tissue.h"
#include "Integrator.h"
#include "Edge.h"
#include "Cell.h"

#define PI 3.141592653589793238462643383279

class  ProfileWindow;


class DrawProfileWindow : public QWidget
{
    Q_OBJECT
	
protected:
    void paintEvent(QPaintEvent *event);
	QMatrix matrix;
	
public:
    DrawProfileWindow(QWidget* parent, Integrator *integrator);
    QString propertiesOutputFile;
	void updateFrame(double *xmin_, double *xmax_, double *ymin_, double *ymax_);
    void extractPropertiesFromCellsInPlane();
    std::list<Cell*> extractCellsInPlane();

private:
	void keyPressEvent(QKeyEvent *event);
	double scaleFactor;
    double sectionPlane_alpha;
    double sectionPlane_beta;
    Point currentPoint;
	bool drawBasalSide;
	
    Integrator *integrator;
    
    void total_update(){update(); /*repaint()*/; QApplication::processEvents();}

    QImage CurrentImage;
    
public slots:
    //writeCurrentImage(Point positionVector_plane, double alpha, double beta, QString imageName, QString folderName)
    void writeCurrentImage(Point positionVector_plane, Point normalVector_plane, QString imageName, QString folderName = QDir::currentPath());
	void slot_sectionNormal_alpha(double alpha) {sectionPlane_alpha=alpha*PI/180.;total_update();};
	void slot_sectionNormal_beta(double beta) {sectionPlane_beta=beta*PI/180.;total_update();};
    void slot_sectionPosition_z(double z) {currentPoint.z = z;total_update();}
    void slot_currentPoint(double x, double y) {currentPoint.x = x; currentPoint.y=y;total_update();}
    void slot_addToCurrentPoint(double x, double y) {currentPoint.x += x; currentPoint.y += y; total_update();}
    void render();
    /// mutate all cells that are on one half of the cross section plane
    void mutateCellsInHalf();
    
};

#include "MainWindow.h"
#endif

#endif //_USE_QT_