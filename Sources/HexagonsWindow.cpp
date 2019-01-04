#ifdef _USE_QT_


#include "HexagonsWindow.h"

HexagonsWindow::HexagonsWindow(MainWindow * mainWindow)
: QWidget()
{
    // set initial values
	height=100;
    length=10;
    width=10;
    noCols=20;
    noRows=20;

    // define Geometry
	setGeometry(40, 40, 500, 500);
    setWindowTitle("Create hexagonal packing");

    // set the layout
	QC_isPeriodic = new QCheckBox("Periodic Tissue");
    QC_isPeriodic->setChecked(true);
    
    QLabel *QL_noRows = new QLabel("Number of rows");
    QS_noRows = new SpecialQSpinBox();
    QS_noRows->setValue(noRows);
    QS_noRows->setMaximum(100);
    QS_noRows->setMinimum(1);
    
    QLabel *QL_noCols = new QLabel("Number of columns");
    QS_noCols = new SpecialQSpinBox();
    QS_noCols->setValue(noCols);
    QS_noCols->setMaximum(100);
    QS_noCols->setMinimum(1);

    QLabel *QL_radiusInCells = new QLabel("Number of cells in radius");
    QS_radiusInCells = new SpecialQSpinBox();
    QS_radiusInCells->setValue(10);
    QS_radiusInCells->setMaximum(100);
    QS_radiusInCells->setMinimum(1);
    
    QLabel *QL_height = new QLabel("cell height");
    QS_height = new SpecialQDoubleSpinBox();
    QS_height->setValue(height);
    QS_height->setMinimum(0);
    QS_height->setMaximum(5000);
    
    QLabel *QL_width = new QLabel("cell width");
    QS_width = new SpecialQDoubleSpinBox();
    QS_width->setMinimum(0);
    QS_width->setMaximum(5000);
    QS_width->setValue(width);
    
    QLabel *QL_length = new QLabel("cell length");
    QS_length = new SpecialQDoubleSpinBox();
    QS_length->setMinimum(0);
    QS_length->setMaximum(5000);
    QS_length->setValue(length);
    
    QP_generate= new QPushButton("Generate Tissue");
    QP_create_round_tissue= new QPushButton("Create round tissue");
    
    QGridLayout *layout = new QGridLayout;
	
    layout->addWidget(QC_isPeriodic, 0, 1);
    layout->addWidget(QL_noRows,     1, 1);
    layout->addWidget(QS_noRows,     1, 2);
    layout->addWidget(QL_noCols,     2, 1);
    layout->addWidget(QS_noCols,     2, 2);
    
    layout->addWidget(QL_radiusInCells,     3, 1);
    layout->addWidget(QS_radiusInCells,     3, 2);
    
    layout->addWidget(QL_height,     4, 1);
    layout->addWidget(QS_height,     4, 2);
    layout->addWidget(QL_width,      5, 1);
    layout->addWidget(QS_width,      5, 2);
    layout->addWidget(QL_length,     6, 1);
    layout->addWidget(QS_length,     6, 2);
    layout->addWidget(QP_generate,   7, 1 );
    layout->addWidget(QP_create_round_tissue,   8, 1 );
    
    setLayout(layout);
    
    // connect the signals and slots
	connect(QC_isPeriodic,SIGNAL(stateChanged(int)),this,SLOT(slot_isPeriodic(int)));
	connect(QS_noRows,SIGNAL(valueChanged(int)),this,SLOT(slot_noRows(int)));
	connect(QS_noCols,SIGNAL(valueChanged(int)),this,SLOT(slot_noCols(int)));
	connect(QS_height,SIGNAL(valueChanged(double)),this,SLOT(slot_height(double)));
	connect(QS_width,SIGNAL(valueChanged(double)),this,SLOT(slot_width(double)));
	connect(QS_length,SIGNAL(valueChanged(double)),this,SLOT(slot_length(double)));
	connect(QP_generate,SIGNAL(clicked()),this,SLOT(slot_generate()));
    connect(QP_create_round_tissue, SIGNAL(clicked()), this, SLOT(slot_createRoundHexagonalTissue()));
    
    //createRoundHexagonalTissue(int radiusInCells, double height, double length, double MaximalRadius)
    //raise();
    //show();
    
    QApplication::processEvents();
}

void HexagonsWindow::slot_generate()
{
    signal_generateTissue(false, QC_isPeriodic->isChecked(),noRows,noCols,height,width,length);
}

// slots for special tissues with special patterning that can be predefined
void HexagonsWindow::slot_generate_cyst()
{
    signal_generateTissue(1);
}

void HexagonsWindow::slot_generate_bigCyst()
{
    signal_generateTissue(2);
}

void HexagonsWindow::slot_show()
{
    if(isVisible()) {hide();}
    else {raise();show();}
}


void HexagonsWindow::slot_createRoundHexagonalTissue()
{
    int radius = QS_radiusInCells->value();
    double height = QS_height->value();
    double length = QS_length->value();
    
    signal_createRoundHexagonalTissue(radius, height, length, 1000);
}

#endif // _USE_QT_