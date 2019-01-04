#ifdef _USE_OPENGL_
#ifdef _USE_QT_

#ifndef Tissue3DWindow_H
#define Tissue3DWindow_H

#include <QApplication>
#include <QWidget>
#include <QDialog>
#include <QSpinBox>
#include <QGraphicsView>
#include <QtOpenGL>

#include "Tissue.h"
#include "Edge.h"
#include "GLWidget.h"

QT_BEGIN_NAMESPACE
class QGroupBox;
class QLabel;
class QLineEdit;
class MainWindow;
QT_END_NAMESPACE

//class DrawTissue3DWindow;

class Tissue3DWindow : public QWidget
{
    Q_OBJECT
public:
    Tissue3DWindow(Integrator *integrator);
    GLWidget* glWidget;
    
public slots:
    void slot_show3DTissue();
    void slot_update(){update(); QApplication::processEvents();}
protected:
    void keyPressEvent(QKeyEvent *event);
    
private:
    QSlider *createSlider();

    Integrator* integrator;
 
    QSlider *xSlider;
    QSlider *ySlider;
    QSlider *zSlider;
    
};
#endif

#endif // _USE_OPENGL_
#endif //_USE_QT_