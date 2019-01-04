#ifdef _USE_QT_

#ifndef ProfileWindow_H
#define ProfileWindow_H

#include <QApplication>
#include <QWidget>
#include <QDialog>
#include <QSpinBox>
#include "Tissue.h"
#include "TissueWindow.h"
#include "Edge.h"
#include "DrawProfileWindow.h"
#include <QGraphicsView>

//#define PI 3.14159265

QT_BEGIN_NAMESPACE
class QGroupBox;
class QLabel;
class QLineEdit;
class MainWindow;
QT_END_NAMESPACE

class DrawProfileWindow;

class ProfileWindow : public QWidget
{
    Q_OBJECT
public:
    ProfileWindow(Integrator *integrator);
    DrawProfileWindow* drawProfileWindow;
    
private:
    Integrator *integrator;
    QPushButton *QP_mutate_cells_in_half;
    void keyPressEvent(QKeyEvent *event);
    double scaleFactor;
    
    bool showProfile;
    Point SectionPlane_PositionVector;
    Point SectionPlane_NormalVector;
    
    double sectionPlane_alpha;
    double sectionPlane_beta;
	double sectionPosition_z;
    
public slots:
    void slot_sectionNormal_alpha(double a) {sectionPlane_alpha=a;};
    void slot_sectionNormal_beta(double b) {sectionPlane_beta=b;};
    void slot_sectionPosition_z(double z) {sectionPosition_z=z;};
    void slot_set_propertiesOutputFile(QString fileName) {drawProfileWindow->propertiesOutputFile = fileName;std::cout << drawProfileWindow->propertiesOutputFile.toStdString() << std::endl;}
    //    void slot_sectionPlane_PositionVector(Point pos) {SectionPlane_PositionVector = pos;};
    void slot_showProfile();
    void slot_updateProfile();
    void slot_update(){update(); /*repaint()*/; QApplication::processEvents();}
};
#include "MainWindow.h"
#endif

#endif // _USE_QT_