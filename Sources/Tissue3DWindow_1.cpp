#ifdef _USE_QT_
#ifdef _USE_OPENGL_

#include <QtGui>
#include <QtOpenGL>
#include <iostream>
#include "Tissue.h"
#include "Tissue3DWindow.h"

using namespace std;

Tissue3DWindow::Tissue3DWindow(MainWindow *mainWindow, Integrator *integrator)
: QWidget()
{
	setGeometry(40, 40, 500, 500);
	
    // give it a name
    setWindowTitle("3D view");
    
    // create the openGL widget
    glWidget = new GLWidget(this, integrator);

    // create and initialize the sliders
    xSlider = createSlider();
    ySlider = createSlider();
    zSlider = createSlider();
    
    // connect sliders
    connect(xSlider, SIGNAL(valueChanged(int)), glWidget, SLOT(setXRotation(int)));
    connect(glWidget, SIGNAL(xRotationChanged(int)), xSlider, SLOT(setValue(int)));
    connect(ySlider, SIGNAL(valueChanged(int)), glWidget, SLOT(setYRotation(int)));
    connect(glWidget, SIGNAL(yRotationChanged(int)), ySlider, SLOT(setValue(int)));
    connect(zSlider, SIGNAL(valueChanged(int)), glWidget, SLOT(setZRotation(int)));
    connect(glWidget, SIGNAL(zRotationChanged(int)), zSlider, SLOT(setValue(int)));
    
    // set layout
    QHBoxLayout *mainLayout = new QHBoxLayout;
    mainLayout->addWidget(glWidget);
    mainLayout->addWidget(xSlider);
    mainLayout->addWidget(ySlider);
    mainLayout->addWidget(zSlider);
    setLayout(mainLayout);
    
    xSlider->setValue(15 * 16);
    ySlider->setValue(345 * 16);
    zSlider->setValue(0 * 16);
    
    // show();
    // raise();
}

QSlider *Tissue3DWindow::createSlider()
{
    QSlider *slider = new QSlider(Qt::Vertical);
    slider->setRange(0, 360 * 16);
    slider->setSingleStep(16);
    slider->setPageStep(15 * 16);
    slider->setTickInterval(15 * 16);
    slider->setTickPosition(QSlider::TicksRight);
    return slider;
}

void Tissue3DWindow::keyPressEvent(QKeyEvent *e)
{
    if (e->key() == Qt::Key_Escape)
        close();
    else
        QWidget::keyPressEvent(e);
}

void Tissue3DWindow::slot_show3DTissue()
{    
    if(!isVisible()) { std::cout << "Pressed" << std::endl; show(); raise();}
    else { std::cout << "Pressed" << std::endl; hide();}
}

#endif // _USE_OPENGL_
#endif // _USE_QT_