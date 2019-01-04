#ifdef _USE_QT_


#include <QtGui>
#include "ProfileWindow.h"


ProfileWindow::ProfileWindow(Integrator *integrator)
: QWidget()
{
	
	// define Geometry
	setGeometry(40, 40, 500, 500);
    
	// the tissue buttons are a new widget of this tissue window
	QWidget *sectionProperties=new QWidget(this);

	QGridLayout *layout = new QGridLayout(sectionProperties);
	
	//setCentralWidget(sectionProperties);
	
	QLabel *labelAlpha= new QLabel("Angle in x-y plane");
	labelAlpha->setFont(QFont("Times", 12));
	
	QDoubleSpinBox *Normal_alpha = new QDoubleSpinBox(sectionProperties);
	Normal_alpha->setValue(0);
	Normal_alpha->setMaximum(10000);
	Normal_alpha->setMinimum(-10000);
	Normal_alpha->setFont(QFont("Times", 12));
	Normal_alpha->setDecimals(2);
	Normal_alpha->setSingleStep(1);
	
	QLabel *labelBeta= new QLabel("Angle in x-z plane");
	labelBeta->setFont(QFont("Times", 12));
	
	QDoubleSpinBox *Normal_beta = new QDoubleSpinBox(sectionProperties);
	Normal_beta->setValue(0);
	Normal_beta->setMaximum(10000);
	Normal_beta->setMinimum(-10000);
	Normal_beta->setFont(QFont("Times", 12));
	Normal_beta->setDecimals(2);
	Normal_beta->setSingleStep(1);
	
	QLabel *labelz= new QLabel("z-coordinate of starting point (0,0,z)");
	labelz->setFont(QFont("Times", 12));
	
	QDoubleSpinBox *Position_z = new QDoubleSpinBox(sectionProperties);
	Position_z->setValue(0);
	Position_z->setMaximum(10000);
	Position_z->setMinimum(-10000);
	Position_z->setFont(QFont("Times", 12));
	Position_z->setDecimals(2);
	Position_z->setSingleStep(1);
	
    QP_mutate_cells_in_half = new QPushButton("Mutate Half");
	QP_mutate_cells_in_half->setFont(QFont("Times", 12, QFont::Bold));
    
	//layout->setColumnStretch(1, 1);
    //layout->setColumnStretch(2, 1);
	layout->addWidget(labelAlpha,1,0,1,1);
	layout->addWidget(Normal_alpha,2,0,1,1);
	layout->addWidget(labelBeta,3,0,1,1);
	layout->addWidget(Normal_beta,4,0,1,1);
	//layout->addWidget(Normal_z,3,0,1,1);
	layout->addWidget(labelz,5,0,1,1);
	layout->addWidget(Position_z,6,0,1,1);
    layout->addWidget(QP_mutate_cells_in_half,7,0,1,1);

	setWindowTitle("Tissue Cross Section");

    drawProfileWindow = new DrawProfileWindow(this, integrator);

    layout->addWidget(drawProfileWindow,0,3,12,6);
    
	setLayout(layout);
    
	connect(Normal_alpha, SIGNAL(valueChanged(double)), drawProfileWindow, SLOT(slot_sectionNormal_alpha(double)));
	connect(Normal_beta, SIGNAL(valueChanged(double)), drawProfileWindow, SLOT(slot_sectionNormal_beta(double)));
	connect(Position_z, SIGNAL(valueChanged(double)), drawProfileWindow, SLOT(slot_sectionPosition_z(double)));
    connect(QP_mutate_cells_in_half, SIGNAL(pressed()),drawProfileWindow,SLOT(mutateCellsInHalf()));
    
    QApplication::processEvents();

	//show();
	//raise();
	
}

void ProfileWindow::slot_showProfile()
{

    if(isVisible())
	{
		hide();
	}
	else
	{
        //drawProfileWindow->paintEvent();
		raise();
		show();
	}
}


void ProfileWindow::slot_updateProfile()
{
	
	update();
    
	//show();
	//raise();
}

void ProfileWindow::keyPressEvent(QKeyEvent *event)
{
    if(event->key()==Qt::Key_H) { drawProfileWindow->extractPropertiesFromCellsInPlane();}

}

//void ProfileWindow::updateProfile()
//{
//	drawProfileWindow->updateProfile();
//}

//
//void ProfileWindow::setShowYprofile()
//{
//	showYprofile = !showYprofile;
//	raise();
//}

#endif // _USE_QT_