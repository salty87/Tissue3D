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
    //qglClearColor(qtGreen.light());
    qglClearColor(Qt::black);
    GLImage = new GLTissueImage(this, integrator);

    //logo->setColor(qtGreen.dark());
    
    scaleFactor = 1;
    
    static float ambient[] =
    {0.1, 0.1, 0.1, 1.0};
    static float diffuse[] =
    {0.5, 1.0, 1.0, 1.0};
    static float position[] =
    {90.0, 90.0, 150.0, 0.0};
    
    static float front_mat_shininess[] =
    {60.0};
    static float front_mat_specular[] =
    {0.2, 0.2, 0.2, 1.0};
    static float front_mat_diffuse[] =
    {0.5, 0.5, 0.28, 1.0};
    static float back_mat_shininess[] =
    {60.0};
    static float back_mat_specular[] =
    {0.5, 0.5, 0.2, 1.0};
    static float back_mat_diffuse[] =
    {1.0, 0.9, 0.2, 1.0};
    
    static float lmodel_ambient[] =
    {1.0, 1.0, 1.0, 1.0};
    
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, position);
    
    glMaterialfv(GL_FRONT, GL_SHININESS, front_mat_shininess);
    glMaterialfv(GL_FRONT, GL_SPECULAR, front_mat_specular);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, front_mat_diffuse);
    glMaterialfv(GL_BACK, GL_SHININESS, back_mat_shininess);
    glMaterialfv(GL_BACK, GL_SPECULAR, back_mat_specular);
    glMaterialfv(GL_BACK, GL_DIFFUSE, back_mat_diffuse);
    
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
    //glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE, 1.0);
    
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glShadeModel(GL_SMOOTH);
}

void GLWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0.0, 0.0, -10.0);
    glRotatef(xRot / 16.0, 1.0, 0.0, 0.0);
    glRotatef(yRot / 16.0, 0.0, 1.0, 0.0);
    glRotatef(zRot / 16.0, 0.0, 0.0, 1.0);
    
    if(tissueChanged) {GLImage->loadParts(); tissueChanged=false;}
    
    //qglClearColor(Qt::darkGray);
    GLImage->draw();
    
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
        //update();
        QApplication::processEvents();
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
	//else if(event->key()==Qt::Key_Left) {matrix.translate(height()/scaleFactor/10,0);}
//	else if(event->key()==Qt::Key_Up) {matrix.translate(0,-height()/scaleFactor/10);}
//	else if(event->key()==Qt::Key_Down) {matrix.translate(0,height()/scaleFactor/10);}
//	else if(event->key()==Qt::Key_Q) {	matrix.rotate(5);}//matrice qui transforme les coordonnees pour les placer en coordonnees mathematiques, commencent en bas a gauche.
//	else if(event->key()==Qt::Key_Z) {	matrix.rotate(-5);}//matrice qui transforme les coordonnees pour les placer en coordonnees mathematiques, commencent en bas a gauche.
//	else if(event->key()==Qt::Key_R) {	matrix.setMatrix(1,0,0,-1,0,rect().height());matrix.scale(1,1);scaleFactor=1;}
    //	else if(event->key()==Qt::Key_D) {CurrentPoint.x = CurrentPoint.x + 0.02*distance.x;}
    //	else if(event->key()==Qt::Key_A) {CurrentPoint.x = CurrentPoint.x - 0.02*distance.x;}
    //	else if(event->key()==Qt::Key_W) {CurrentPoint.y = CurrentPoint.y + 0.02*distance.y;}
    //	else if(event->key()==Qt::Key_S) {CurrentPoint.y = CurrentPoint.y - 0.02*distance.y;}
    
	update();
}

#endif // _USE_QT_