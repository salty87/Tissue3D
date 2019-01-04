#ifdef _USE_QT_

#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>
#include <QApplication>
#include "Tissue.h"
// #include "qtlogo.h"
#include "GLTissueImage.h"



class GLWidget : public QGLWidget
{
    Q_OBJECT
    
public:
    GLWidget(QWidget* parent, Integrator *integrator);
    ~GLWidget();

    
    QSize minimumSizeHint() const;
    QSize sizeHint() const;
    
    GLTissueImage *GLImage;
    
public slots:
    void setXRotation(int angle);
    void setYRotation(int angle);
    void setZRotation(int angle);
    void slot_TissueChanged();
    
signals:
    void xRotationChanged(int angle);
    void yRotationChanged(int angle);
    void zRotationChanged(int angle);
    
    
protected:
    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void keyPressEvent(QKeyEvent *event);
    void draw();
    
private:
    bool tissueChanged;
    Integrator* integrator;
    
    int xRot;
    int yRot;
    int zRot;
    double scaleFactor;
    QPoint lastPos;
    QColor qtGreen;
    QColor qtPurple;
};

#endif // GLWIDGET_H
#endif // _USE_QT_