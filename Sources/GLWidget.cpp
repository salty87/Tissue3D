#ifdef _USE_QT_


#include <QResizeEvent>
#include "GLWidget.h"
//#include "qtlogo.h"

using namespace std;

GLWidget::GLWidget(QWidget *parent, Integrator *integrator)
: QGLWidget(QGLFormat(QGL::SampleBuffers), parent), integrator(integrator)
{
    // set the image and the initial angles to zero
    GLImage = 0;
    xRot = 0;
    yRot = 0;
    zRot = 0;
    
    tissueChanged = true;
    
    // define back ground colours
    qtGreen = QColor::fromCmykF(0.40, 0.0, 1.0, 0.0);
    qtPurple = QColor::fromCmykF(0.39, 0.39, 0.0, 0.0);
    
    setFocusPolicy(Qt::StrongFocus);
}

GLWidget::~GLWidget()
{
}

QSize GLWidget::minimumSizeHint() const
{
    return QSize(50, 50);
}

QSize GLWidget::sizeHint() const
{
    return QSize(800, 800);
}

static void qNormalizeAngle(int &angle)
{
    while (angle < 0)
        angle += 360 * 16;
    while (angle > 360 * 16)
        angle -= 360 * 16;
}

void GLWidget::setXRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != xRot) {
        xRot = angle;
        emit xRotationChanged(angle);
        updateGL();
    }
}

void GLWidget::setYRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != yRot) {
        yRot = angle;
        emit yRotationChanged(angle);
        updateGL();
    }
}

void GLWidget::setZRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != zRot) {
        zRot = angle;
        emit zRotationChanged(angle);
        updateGL();
    }
}

void GLWidget::initializeGL()
{
    // draw the background
    if(!integrator->T->onlyApicalContributions){
        qglClearColor(Qt::white);
    }  else {
        qglClearColor(Qt::black);
    }
    
    qglClearColor(Qt::white);
    
    GLImage = new GLTissueImage(this, integrator);

    scaleFactor = 1;
    
//    static float ambient[] =  {1.0, 1.0, 1.0, 1.0};
//    static float diffuse[] =  {1.0, 1.0, 1.0, 1.0};//{0.5, 1.0, 1.0, 1.0};
//    static float lmodel_ambient[] =  {1.0, 1.0, 1.0, 1.0};

    static float position0[] =  {0, 0, 0.5, 1};
    static float position1[] =  {0, 0, 0.5, 1};
    static float direction0[] =  {0, 0, -1};
    static float direction1[] =  {0, 0, 1};

    
    static float position2[] =  {0, 1000, 0, 0};
    static float position3[] =  {0, -1000, 0, 0};
    static float direction2[] =  {0, -1, 0};
    static float direction3[] =  {0, 1, 0};

    
    // use two light sources
    glLightfv(GL_LIGHT0, GL_POSITION, position0);
    glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, direction0);
    glLightfv(GL_LIGHT1, GL_POSITION, position1);
    glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, direction1);
    glLightfv(GL_LIGHT2, GL_POSITION, position2);
    glLightfv(GL_LIGHT2, GL_SPOT_DIRECTION, direction2);
    glLightfv(GL_LIGHT3, GL_POSITION, position3);
    glLightfv(GL_LIGHT3, GL_SPOT_DIRECTION, direction3);
    
    // create the lighting conditions
    glEnable(GL_LIGHTING);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LIGHT1);
    //glEnable(GL_LIGHT2);
    //glEnable(GL_LIGHT3);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
}

void GLWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0.0, 0.0, -10.0);
    glRotatef(xRot / 16.0, 1.0, 0.0, 0.0);
    glRotatef(yRot / 16.0, 0.0, 1.0, 0.0);
    glRotatef(zRot / 16.0, 0.0, 0.0, 1.0);
    
    int x_r = 1;
    int y_r = 10;

    if(tissueChanged) {
        
        if(!integrator->T->onlyApicalContributions) GLImage->loadPartsInRadius(100000,x_r,y_r);//integrator->T->SystemSize.x/3.0);
        else GLImage->loadPartsInRadius(100000);//GLImage->loadParts();
        
        tissueChanged=false;}
    
    GLImage->draw();
    
    QImage  img = grabFrameBuffer(false);
    img.save("/Users/silvanus/MPI/EpithelialMechanics/Tissue3D/imageSave.jpg",0,100);

}

void GLWidget::resizeGL(int width, int height)
{
    //glViewport((width - side) / 2, (height - side) / 2, side, side);
    glViewport(0, 0, (GLsizei) width, (GLsizei) height);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
#ifdef QT_OPENGL_ES_1
    glOrthof(-0.5, +0.5, -0.5, +0.5, 4.0, 15.0);
#else
    glOrtho(-0.5, +0.5, -0.5, +0.5, 4.0, 15.0);
#endif
    glMatrixMode(GL_MODELVIEW);
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    lastPos = event->pos();
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    // press and release changes the angle of view
    int dx = event->x() - lastPos.x();
    int dy = event->y() - lastPos.y();
    
    if (event->buttons() & Qt::LeftButton) {
        setXRotation(xRot + 8 * dy);
        setYRotation(yRot + 8 * dx);
    } else if (event->buttons() & Qt::RightButton) {
        setXRotation(xRot + 8 * dy);
        setZRotation(zRot + 8 * dx);
    }
    lastPos = event->pos();
}


void GLWidget::slot_TissueChanged()
{
    //std::cout << "GLWidget::slot_TissueChanged()" << std::endl;
    tissueChanged = true;
    if(isVisible())
    {
        paintGL();
        update();
        QApplication::processEvents();
     
        
        QImage  img = grabFrameBuffer(true);
        img.save("/Users/silvanus/MPI/EpithelialMechanics/Tissue3D/imageSave.png");        
        
        // save the current image
//        QImage img(sizeHint(),QImage::Format_RGB32);
//        QPainter painter(&img);
//        render(&painter);
//        img.save("/Users/silvanus/MPI/EpithelialMechanics/Tissue3D/imageSave.png");
        
    }
}

void GLWidget::keyPressEvent(QKeyEvent *event)
{
	//if(event->key()==Qt::Key_Equal) {matrix.scale(1.2,1.2);scaleFactor*=1.2;}
	if(event->key()==Qt::Key_Minus) {glMatrixMode(GL_PROJECTION);glScalef(0.8,0.8,0.8);scaleFactor*=0.8;glMatrixMode(GL_MODELVIEW);}
	else if(event->key()==Qt::Key_Equal) {glMatrixMode(GL_PROJECTION);glScalef(1.25,1.25,1.25);scaleFactor*=1.25;glMatrixMode(GL_MODELVIEW);}
	else if(event->key()==Qt::Key_Down) {glMatrixMode(GL_PROJECTION);glTranslatef(0,scaleFactor/10,0);glMatrixMode(GL_MODELVIEW);}
	else if(event->key()==Qt::Key_Up) {glMatrixMode(GL_PROJECTION);glTranslatef(0,-scaleFactor/10,0);glMatrixMode(GL_MODELVIEW);}
	else if(event->key()==Qt::Key_Right) {glMatrixMode(GL_PROJECTION);glTranslatef(-scaleFactor/10,0,0);glMatrixMode(GL_MODELVIEW);}
	else if(event->key()==Qt::Key_Left) {glMatrixMode(GL_PROJECTION);glTranslatef(scaleFactor/10,0,0);glMatrixMode(GL_MODELVIEW);}
    
	update();
}

#endif // _USE_QT_