#ifdef _USE_QT_

#include "MainWindow.h"

//detailled constructor of the main window, tissueWindow and the iterator buttons (first box, left hand side) are defined
MainWindow::MainWindow()
: integrator(new Integrator())
{
    // define Geometry
	setWindowTitle(tr("Tissue3D"));
    setGeometry(40, 40, 1250, 760);
	setUnifiedTitleAndToolBarOnMac(true);
	tissueButtons=new QWidget(this); // central widget
	setCentralWidget(tissueButtons);
    //tissueButtons->setMaximumWidth(400);
	int MaxWidth=80;
	CurrentRecordPath=QDir::currentPath()+"frame";
    
    /************************* create the windows ***************************/
	hexagonsWindow = new HexagonsWindow(this);
#ifdef _USE_OPENGL_
    tissue3DWindow = new Tissue3DWindow(integrator);
#endif
    profileWindow = new ProfileWindow(integrator);
    connect(integrator,SIGNAL(signal_set_propertiesOutputFile(QString)),profileWindow,SLOT(slot_set_propertiesOutputFile(QString)));
    
    /************* Define the Push Buttons and other elements of the tissue window **************/
	
	// set push button Integrator
	QP_Integrator=new QPushButton(tr("Integrate"),tissueButtons);
	QP_Integrator->setMaximumWidth(MaxWidth);

    // push button to stop minimization
	QP_Stop=new QPushButton(tr("Stop"),tissueButtons);
	QP_Stop->setMaximumWidth(MaxWidth);
    
    // push buttons to display the cross section
    QPushButton *QP_showProfile= new QPushButton("Cross section",this);
    QP_showProfile->setMaximumWidth(120);
    
    // push button to display the 3D view
    QPushButton *QP_3DView= new QPushButton("3D View",this);
    QP_3DView->setMaximumWidth(120);
    
    connect(QP_showProfile,SIGNAL(clicked()),profileWindow,SLOT(slot_showProfile()));
#ifdef _USE_OPENGL_
    connect(QP_3DView,SIGNAL(pressed()),tissue3DWindow,SLOT(slot_show3DTissue()));
#endif
    /************************ create all the actions and boxes of the main window ****************/
	createActions();
	createMenus();
	createIntegrationBox();
	createGastrulationBox();
	createCystBox();
	createMechanicsBox();
	createTissueBox();
    createCellBox();
    createEdgeBox();
    createVertexBox();
    //createRecordBox();
    
    /************** status bar ***************/
    connect(integrator,SIGNAL(signal_show_status(QString)), this, SLOT(change_status(QString)));
    change_status("Ready");
    
    /************************ Define property tabs shown on left hand side of the tissue window ****************/
	PropertiesTab=new QTabWidget(this);
    PropertiesTab->setUsesScrollButtons(true);

    PropertiesTab->setUsesScrollButtons(true);
    IndexMechanicsBox=PropertiesTab->addTab(MechanicsBox, "Mechanics");
    IndexIntegrationBox=PropertiesTab->addTab(IntegrationBox, "Integration");
	IndexCystBox=PropertiesTab->addTab(CystBox, "Cyst");
	IndexGastrulationBox=PropertiesTab->addTab(GastrulationBox, "Gastrulation");
    IndexTissueBox=PropertiesTab->addTab(TissueBox, "Tissue");
	IndexCellBox=PropertiesTab->addTab(CellBox, "Cell");
	IndexEdgeBox=PropertiesTab->addTab(EdgeBox, "Edge");
	IndexVertexBox=PropertiesTab->addTab(VertexBox, "Vertex");
	connect(PropertiesTab,SIGNAL(currentChanged(int)),integrator->T,SLOT(slot_setMainWindowTab(int)));
	
    
    /************************ create text box to execute scripts on the run ****************/
    QT_scriptTextBox = new QTextEdit();
    QP_executeScriptBox  = new QPushButton("Run Script");
    connect(QP_executeScriptBox,SIGNAL(pressed()),this,SLOT(runScriptFromBox()));
    
    /************************ Define the tabs of the record box under the cell/vertex/... properties ****************/
//    ToolsTab=new QTabWidget(this);
//    ToolsTab->setMaximumWidth(MaxWidth);
//    IndexDisplayBox=ToolsTab->addTab(DisplayBox,"Display");
//    IndexRecordBox=ToolsTab->addTab(RecordBox,"Record");
//    IndexExtractInfosBox=ToolsTab->addTab(ExtractInfosBox,"Extract informations");
//    IndexSelectGroupBox=ToolsTab->addTab(SelectGroupBox,"Select Group");
	
    /************************ connect the defined buttons to the actions that are performed when pressed ****************/
    
	connect(QP_Integrator, SIGNAL(clicked()),integrator,SLOT(integrate()));
	connect(QP_Stop, SIGNAL(clicked()),integrator,SLOT(slot_stop()));
    connect(this,SIGNAL(signal_loadFromTxt(QString)),integrator,SLOT(loadFromTxt(QString)));
    // connect the integrator to the tissue window
    
    /************************ drawing updates ****************/
    
//    connect(this,SIGNAL(update_all()),this,SLOT(update()));
//    connect(this,SIGNAL(update_all()),profileWindow,SLOT(slot_update()));
    connect(integrator, SIGNAL(signal_updateMainWindow()),this,SLOT(slot_update()));
    connect(integrator, SIGNAL(signal_updateCrossSection()),profileWindow,SLOT(slot_update()));
#ifdef _USE_OPENGL_
//    connect(this,SIGNAL(update_all()),tissue3DWindow,SLOT(slot_update()));
    //connect(integrator,SIGNAL(signal_updateWindow()),tissue3DWindow->glWidget,SLOT(slot_TissueChanged()));
    connect(integrator, SIGNAL(signal_update3DWindow()),tissue3DWindow->glWidget,SLOT(slot_TissueChanged()));
    connect(integrator, SIGNAL(signal_update3DWindow()),tissue3DWindow,SLOT(slot_update()));
#endif

    
    
    //connect(PropertiesTab,currentChanged(int),tissueWindow,SLOT(slot_tabChanged(int)));
    /************************ How are the buttons organized? ****************/
    tissueWindow = new TissueWindow(this,integrator);
    
	//define the grid layout of the tissue window
    QGridLayout* gridLayout = new QGridLayout(tissueButtons);

	gridLayout->setHorizontalSpacing(5);
	gridLayout->setVerticalSpacing (5);
	
    gridLayout->addWidget(QP_Integrator,  0,0,1,1);
	gridLayout->addWidget(QP_Stop,        1,0,1,1);
    gridLayout->addWidget(QP_showProfile, 0,1,1,1);
	gridLayout->addWidget(QP_3DView,      1,1,1,1);

    gridLayout->addWidget(PropertiesTab,3,0,5,2);
    
    gridLayout->addWidget(QT_scriptTextBox,8,0,2,2);
    gridLayout->addWidget(QP_executeScriptBox,10,0,1,1);
    
    gridLayout->addWidget(tissueWindow,0,4,12,8);
    
	//gridLayout->addWidget(PropertiesTab,8,0,1,2);
	//gridLayout->addWidget(ToolsTab,9,0,1,2);

    
    setLayout(gridLayout);
	
    
    /************************ set the title of the window and open it, push tension window to front and set the current file name0 ****************/
	CurrentPoint=Point(0,0,0);
    integrator->slot_createHexagonalTissue(0, 1, 4, 4, 100,10,10);
   
    connect(tissueWindow,SIGNAL(signal_addCurrentPoint(double,double)),integrator->T,SLOT(slot_addCurrentPoint(double,double)));
    connect(tissueWindow,SIGNAL(signal_currentPoint(double,double)),integrator->T,SLOT(slot_currentPoint(double,double)));
    connect(tissueWindow,SIGNAL(signal_addCurrentPoint(double,double)),profileWindow->drawProfileWindow,SLOT(slot_currentPoint(double,double)));
    connect(tissueWindow,SIGNAL(signal_currentPoint(double,double)),profileWindow->drawProfileWindow,SLOT(slot_currentPoint(double,double)));
    connect(tissueWindow,SIGNAL(signal_addToCurrentPoint(double,double)),profileWindow->drawProfileWindow,SLOT(slot_addToCurrentPoint(double, double)));
    connect(tissueWindow,SIGNAL(signal_printCellData()),this,SLOT(slot_printCellData()));
    connect(integrator, SIGNAL(signal_writeImage(Point,Point,QString,QString)), profileWindow->drawProfileWindow, SLOT(writeCurrentImage(Point,Point,QString, QString)));
    connect(this, SIGNAL(signal_writeImage(Point,Point,QString,QString)), profileWindow->drawProfileWindow, SLOT(writeCurrentImage(Point,Point,QString, QString)));
    
    update_all();
    
    raise();
	show();
    
    //hexagonsWindow->show();
    //hexagonsWindow->raise();
	CurrentFileName=QDir::currentPath();
    std::cout << "MainWindow created" << std::endl;
    
}

// destructor of the main window
MainWindow::~MainWindow()
{
	delete tissueWindow;
    delete profileWindow;
#ifdef _USE_OPENGL_
    delete tissue3DWindow;
#endif
    delete integrator;
}


// definition of the actions that are carried out after shortcut combinations
void MainWindow::createActions()
{
    /************************ set the actions, that are accessed by short cuts ****************/
	// define action
	evolveAct = new QAction(tr("&Evolve"), this);
	// define short cut
    evolveAct->setShortcut(tr("Ctrl+E"));
	// define status tip
    evolveAct->setStatusTip(tr("Find Minimum Configuration"));
	// connect action to the according function
	connect(evolveAct, SIGNAL(triggered()),integrator,SLOT(integrate()));
	
	//...
	quitAct = new QAction(tr("&Quit"), this);
    quitAct->setShortcut(tr("Ctrl+Q"));
    quitAct->setStatusTip(tr("Quit Application"));
    connect(quitAct, SIGNAL(triggered()), qApp, SLOT(quit()));
	
	openAct = new QAction(tr("&Open..."), this);
    openAct->setShortcut(tr("Ctrl+O"));
    openAct->setStatusTip(tr("Open an existing file"));
    connect(openAct, SIGNAL(triggered()), this, SLOT(slot_open()));
	
	createHexagonalTissueAct = new QAction(tr("&Create hexagonal tissue..."), this);
    createHexagonalTissueAct->setShortcut(tr("Ctrl+P"));
    createHexagonalTissueAct->setStatusTip(tr("Creates a hexagonal tissue"));
    connect(createHexagonalTissueAct, SIGNAL(triggered()), hexagonsWindow, SLOT(slot_show()));
	
	
    saveCrossSectionsAsPDFAct = new QAction(tr("Open all files and save cross sections as PDF..."), this);
    saveCrossSectionsAsPDFAct->setShortcut(tr("Ctrl+I"));
    saveCrossSectionsAsPDFAct->setStatusTip(tr("Saves images as PDFs"));
    connect(saveCrossSectionsAsPDFAct, SIGNAL(triggered()), this, SLOT(slot_saveCrossSectionsAsPDFActFromSubDirectories()));

    
//	  openImageAct = new QAction(tr("&Open Binary Image..."), this);
//    openImageAct->setShortcut(tr("Ctrl+I"));
//    openImageAct->setStatusTip(tr("Open a binary image"));
//    connect(openImageAct, SIGNAL(triggered()), this, SLOT(slot_openImage()));
//	
//	  saveImageAct = new QAction(tr("&Save as Image..."), this);
//    saveImageAct->setShortcut(tr("Ctrl+Shift+I"));
//    saveImageAct->setStatusTip(tr("Save as image"));
//    connect(saveImageAct, SIGNAL(triggered()), this, SLOT(slot_saveImage()));
//	
//	  savePDFImageAct = new QAction(tr("&Save as PDF..."), this);
//    savePDFImageAct->setShortcut(tr("Ctrl+Shift+P"));
//    savePDFImageAct->setStatusTip(tr("Save as PDF"));
//    connect(savePDFImageAct, SIGNAL(triggered()), this, SLOT(slot_savePDFImage()));
	
	saveAct = new QAction(tr("&Save..."), this);
    saveAct->setShortcut(tr("Ctrl+S"));
    saveAct->setStatusTip(tr("Save file"));
    connect(saveAct, SIGNAL(triggered()), this, SLOT(slot_save()));
	connect(this, SIGNAL(signal_writeToTxt(const char*)), integrator->T, SLOT(slot_writeFileToTxt(const char*)));
    
    
    extractDataFromDir = new QAction(tr("&Extract properties from File..."), this);
    extractDataFromDir->setShortcut(tr("Ctrl+D"));
    extractDataFromDir->setStatusTip(tr("Extract Properties from dir"));
    connect(extractDataFromDir, SIGNAL(triggered()), this, SLOT(slot_extractDataFromSubDir()));
	
//	  showTensionWindowAct = new QAction(tr("&Show tensions window"), this);
//	  showTensionWindowAct->setShortcut(tr("Ctrl+T"));
//	  connect(showTensionWindowAct, SIGNAL(triggered()), tensionWindow, SLOT(show()));
//	  connect(showTensionWindowAct, SIGNAL(triggered()), tensionWindow, SLOT(raise()));
	
//	showHexagonsWindowAct = new QAction(tr("&Show hexagons window"), this);
//	showHexagonsWindowAct->setShortcut(tr("Ctrl+H"));
//	connect(showHexagonsWindowAct, SIGNAL(triggered()), hexagonsWindow, SLOT(show()));
//	connect(showHexagonsWindowAct, SIGNAL(triggered()), hexagonsWindow, SLOT(raise()));
//	
//	
//	showGlobalPropertiesAct = new QAction(tr("&Show global properties window"), this);
//	connect(showGlobalPropertiesAct, SIGNAL(triggered()), GlobalPropertiesWindow, SLOT(show()));
//	connect(showGlobalPropertiesAct, SIGNAL(triggered()), GlobalPropertiesWindow, SLOT(raise()));
//	
//	showAddNoiseAct = new QAction(tr("&Show Add Noise window"), this);
//	connect(showAddNoiseAct, SIGNAL(triggered()), AddNoiseWindow, SLOT(show()));
//	connect(showAddNoiseAct, SIGNAL(triggered()), AddNoiseWindow, SLOT(raise()));
//	
	runScriptAct = new QAction(tr("&Run Script"), this);
    runScriptAct->setShortcut(tr("Ctrl+R"));
    runScriptAct->setStatusTip(tr("Run Script"));
    connect(runScriptAct, SIGNAL(triggered()), this, SLOT(runScript()));
}

// toolbar for evolve and quit is created and according actions are defined
void MainWindow::createToolBars()
{
    fileToolBar = addToolBar(tr("evolve"));
    fileToolBar->addAction(evolveAct);
	fileToolBar->addAction(quitAct);
}

// create the menu bar on top of the window (File - Windows - Script) and the commands connected to the buttons
void MainWindow::createMenus()
{
    fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(openAct);
    fileMenu->addAction(createHexagonalTissueAct);
    fileMenu->addAction(saveCrossSectionsAsPDFAct);
    fileMenu->addAction(extractDataFromDir);
    //fileMenu->addAction(openImageAct);
    fileMenu->addAction(saveAct);
    //fileMenu->addAction(writeTissueAct);
    //fileMenu->addAction(saveImageAct);
    //fileMenu->addAction(savePDFImageAct);
    
    
    //WindowsMenu=menuBar()->addMenu(tr("&Windows"));
    //WindowsMenu->addAction(showTensionWindowAct);
    //WindowsMenu->addAction(showHexagonsWindowAct);
    //WindowsMenu->addAction(showGlobalPropertiesAct);
    //WindowsMenu->addAction(showAddNoiseAct);
    
    ScriptMenu=menuBar()->addMenu(tr("&Script"));
    ScriptMenu->addAction(runScriptAct);
    
}

void MainWindow::createTissueBox()
{
    /************** create the box and define the layout object ***************/
	//CellBox = new QGroupBox(tr("Cell Properties"));
	TissueBox = new QGroupBox();

    QGridLayout *layout = new QGridLayout;
	double MaxWidth=60;
    
    /************** add the labels and boxes for number, area, K, preferred area, type and pressure ***************/
	// button to open a window that enables the definition of the geometry of a new tissue
    QPushButton *QP_create_newTissueWindow = new QPushButton("Create Tissue");

	// button to open a window that enables the definition of the mechanical properties of a new window
    //QPushButton *QP_create_mechPropWindow = new QPushButton("Mechanics");
    //QP_create_mechPropWindow->setFont(QFont("Times", 12, QFont::Bold));
    
    // Number of cells in the tissue
	QLabel *QL_NoCells= new QLabel("#cells");
	QL_NoCells->setMaximumWidth(MaxWidth);
	SpecialQSpinBox *QS_NoCells=new SpecialQSpinBox;
	QS_NoCells->setMaximum(10000);
    QS_NoCells->setReadOnly(true);
    
	// System sizes Lx and Ly
	QLabel * QL_Lx=new QLabel("Lx");
    QL_Lx->setMaximumWidth(MaxWidth);
	QDoubleSpinBox *QS_Lx=new QDoubleSpinBox;
	QS_Lx->setMaximum(999999999);
	QS_Lx->setMaximumWidth(MaxWidth);
	QLabel * QL_Ly=new QLabel("Ly");
    QL_Ly->setMaximumWidth(MaxWidth);
	QDoubleSpinBox *QS_Ly=new QDoubleSpinBox;
	QS_Ly->setMaximum(999999999);
	QS_Ly->setMaximumWidth(MaxWidth);
	
   	// force on the system size
	QLabel * QL_FLx=new QLabel("F on Lx");
    QL_FLx->setMaximumWidth(MaxWidth);
	QDoubleSpinBox *QS_FLx=new QDoubleSpinBox;
	QS_FLx->setMinimum(-999999999);
	QS_FLx->setMaximum(999999999);
	QS_FLx->setMaximumWidth(MaxWidth);
	QS_FLx->setReadOnly(true);
    QLabel * QL_FLy=new QLabel("F on Ly");
	QL_FLy->setMaximumWidth(MaxWidth);
	QDoubleSpinBox *QS_FLy=new QDoubleSpinBox;
	QS_FLy->setMinimum(-999999999);
	QS_FLy->setMaximum(999999999);
	QS_FLy->setMaximumWidth(MaxWidth);
	QS_FLy->setReadOnly(true);
   
//    // external force on the system size
//	QLabel * QL_Text=new QLabel("T_ext");
//	QL_Text->setFont(QFont("Times", 10));
//	QDoubleSpinBox *QS_Text=new QDoubleSpinBox;
//	QS_Text->setFont(QFont("Times", 10));
//	QS_Text->setMinimum(-999999999);
//	QS_Text->setMaximum(999999999);
//	QS_Text->setMaximumWidth(MaxWidth);
//    
//    // fix system size
//    QCheckBox *QC_fixLx = new QCheckBox("fix Lx");
//    QC_fixLx->setCheckState(Qt::Unchecked);
//    QC_fixLx->setFont(QFont("Times", 10));
//    QCheckBox *QC_fixLy = new QCheckBox("fix Ly");
//    QC_fixLy->setCheckState(Qt::Unchecked);
//    QC_fixLy->setFont(QFont("Times", 10));
    
    
    // energy
    QLabel * QL_Energy=new QLabel("Energy");
    QL_Energy->setMaximumWidth(MaxWidth);
	QDoubleSpinBox *QS_Energy=new QDoubleSpinBox;
	QS_Energy->setMinimum(-999999999);
	QS_Energy->setMaximum(999999999);
	QS_Energy->setMaximumWidth(MaxWidth);
    QS_Energy->setReadOnly(true);
    /************** define the locations of the labels and spin boxes ***************/
	layout->addWidget(QP_create_newTissueWindow,     1,0);
    layout->addWidget(QL_Lx,     2,0);
	layout->addWidget(QS_Lx,     2,1);
	layout->addWidget(QL_Ly,     2,2);
	layout->addWidget(QS_Ly,     2,3);
	layout->addWidget(QL_FLx,    3,0);
	layout->addWidget(QS_FLx,    3,1);
	layout->addWidget(QL_FLy,    3,2);
	layout->addWidget(QS_FLy,    3,3);
//	layout->addWidget(QL_Text,   4,0);
//	layout->addWidget(QS_Text,   4,1);
//	layout->addWidget(QC_fixLx,  5,0);
//	layout->addWidget(QC_fixLy,  5,2);
    layout->addWidget(QL_NoCells,6,0);
	layout->addWidget(QS_NoCells,6,1);
	layout->addWidget(QL_Energy, 6,2);
	layout->addWidget(QS_Energy, 6,3);
	
	layout->setHorizontalSpacing (0);
	layout->setVerticalSpacing (0);
	
    /************** connect the spin boxes to the functions, i.e. changes made in the window are translated into changes in the objects ***************/
	connect(integrator->T, SIGNAL(signal_currentTissue(int,double,double,double,double,double,bool,bool,double)), this, SLOT(slot_currentTissue(int,double,double,double,double,double,bool,bool,double)) );

    connect(this,SIGNAL(signal_setLx(double)),QS_Lx,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_setLy(double)),QS_Ly,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_setFLx(double)),QS_FLx,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_setFLy(double)),QS_FLy,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_setEnergy(double)),QS_Energy,SLOT(setValue(double)));
//    connect(this,SIGNAL(signal_setText(double)),QS_Text,SLOT(setValue(double)));
//    connect(this,SIGNAL(signal_setFixLx(Qt::CheckState)),QC_fixLx,SLOT(setCheckState(Qt::CheckState)));
//    connect(this,SIGNAL(signal_setFixLy(Qt::CheckState)),QC_fixLy,SLOT(setCheckState(Qt::CheckState)));

    connect(this,SIGNAL(signal_setEnergy(double)),QS_Energy,SLOT(setValue(double)));
    
    // connect outgoing signals
    //connect(QS_Lx,SIGNAL(valueChanged(double)),integrator->T,SLOT(slot_setLx(double)));
    //connect(QS_Ly,SIGNAL(valueChanged(double)),integrator->T,SLOT(slot_setLy(double)));
    //connect(QS_Text,SIGNAL(valueChanged(double)),integrator->T,SLOT(slot_setText(double)));
    connect(QC_fixLx,SIGNAL(stateChanged(int)),integrator->T,SLOT(slot_setFixLx(int)));
    connect(QC_fixLy,SIGNAL(stateChanged(int)),integrator->T,SLOT(slot_setFixLy(int)));
    connect(QP_create_newTissueWindow,SIGNAL(clicked()),hexagonsWindow,SLOT(slot_show()));
    //connect(QP_create_mechPropWindow,SIGNAL(clicked()),mechPropWindow,SLOT(slot_show()));

    
    // connect the hexagonsWindow to the integrator
    connect(hexagonsWindow,SIGNAL(signal_generateTissue(int,bool, int, int, double, double, double)), integrator, SLOT(slot_createHexagonalTissue(int,bool, int, int, double, double, double)));
    connect(hexagonsWindow, SIGNAL(signal_createRoundHexagonalTissue(int, double, double, double)), integrator->T, SLOT(createRoundHexagonalTissue(int, double, double, double)));

    /************** what is done here? ***************/
    layout->setColumnStretch(1, 1);
    //layout->setColumnStretch(2, 1);
    // layout is set
    TissueBox->setLayout(layout);
}


// define the cell box (1st tab, second box on the left hand side)
void MainWindow::createCellBox()
{
    /************** create the box and define the layout object ***************/
	//CellBox = new QGroupBox(tr("Cell Properties"));
	CellBox = new QGroupBox();
	//! [8]
    QGridLayout *layout = new QGridLayout;
	double MaxWidth=60;
    
    /************** add the labels and boxes for number, area, K, preferred area, type and pressure ***************/
	// Number
	QLabel *QL_No= new QLabel("Cell Number");
	QL_No->setMaximumWidth(MaxWidth);
	SpecialQSpinBox *QS_No=new SpecialQSpinBox;
	QS_No->setMaximum(1000000);

    // T_a
	QLabel *QL_Ta= new QLabel("T_a");
	QL_Ta->setMaximumWidth(MaxWidth);
	SpecialQDoubleSpinBox *QS_Ta=new SpecialQDoubleSpinBox;
	QS_Ta->setMaximum(1000000000);
	QS_Ta->setMinimum(-1000000000);
    QS_Ta->setMaximumWidth(MaxWidth);

    // A_a
	QLabel *QL_Aa= new QLabel("apical area");
	QL_Aa->setMaximumWidth(MaxWidth);
	SpecialQDoubleSpinBox *QS_Aa=new SpecialQDoubleSpinBox;
	QS_Aa->setMaximum(1000000000);
	QS_Aa->setMinimum(0);
    QS_Aa->setMaximumWidth(MaxWidth);

    // A_a
	QLabel *QL_AR= new QLabel("aspect ratio");
	QL_AR->setMaximumWidth(MaxWidth);
	SpecialQDoubleSpinBox *QS_AR=new SpecialQDoubleSpinBox;
	QS_AR->setMaximum(1000000000);
	QS_AR->setMinimum(0);
    QS_AR->setMaximumWidth(MaxWidth);
    
    // T_b
    QLabel *QL_Tb= new QLabel("T_b");
    QL_Tb->setMaximumWidth(MaxWidth);
    SpecialQDoubleSpinBox *QS_Tb=new SpecialQDoubleSpinBox;
	QS_Tb->setMaximum(1000000000);
	QS_Tb->setMinimum(-1000000000);
    QS_Tb->setMaximumWidth(MaxWidth);
    
    // Volume
	QLabel * QL_V=new QLabel("Volume");
	QL_V->setMaximumWidth(MaxWidth);
	QDoubleSpinBox *QS_V=new QDoubleSpinBox;
	QS_V->setMaximum(20000000);
	QS_V->setMinimum(-20000000);
    QS_V->setMaximumWidth(MaxWidth);
	QS_V->setReadOnly(true);
	
    // Preferred Volume
	QLabel * QL_V0=new QLabel("V0");
	QL_V0->setMaximumWidth(MaxWidth);
	SpecialQDoubleSpinBox *QS_V0=new SpecialQDoubleSpinBox;
	QS_V0->setMaximum(20000000);
	QS_V0->setMinimum(0);
	QS_V0->setMaximumWidth(MaxWidth);
	QS_V0->setKeyboardTracking(false);
   
    // Preferred apical area
	QLabel * QL_A0=new QLabel("apical A0");
	QL_A0->setMaximumWidth(MaxWidth);
	SpecialQDoubleSpinBox *QS_A0=new SpecialQDoubleSpinBox;
	QS_A0->setMaximum(100000000000);
	QS_A0->setMinimum(0);
	QS_A0->setMaximumWidth(MaxWidth);
	QS_A0->setKeyboardTracking(false);
    
    
	// Type
	QLabel *QL_Type=new QLabel("Type");
	QL_Type->setMaximumWidth(MaxWidth);
	SpecialQSpinBox *QS_Type=new SpecialQSpinBox;
	QS_Type->setMaximum(5);
	QS_Type->setMinimum(0);
	QS_Type->setKeyboardTracking(false);
    QS_Type->setMaximumWidth(MaxWidth);
	// K
	QLabel *QL_K=new QLabel("K_3D");
	QL_K->setMaximumWidth(MaxWidth);
	SpecialQDoubleSpinBox *QS_K=new SpecialQDoubleSpinBox;
	QS_K->setMaximum(1000000);
	QS_K->setMinimum(0);
	QS_K->setValue(1);
	QS_K->setDecimals(6);
	QS_K->setSingleStep(0.01);
	QS_K->setMaximumWidth(MaxWidth);
	QS_K->setKeyboardTracking(false);
	
    // K_2D
    QLabel *QL_K2D=new QLabel("K_2D");
	QL_K2D->setMaximumWidth(MaxWidth);
	SpecialQDoubleSpinBox *QS_K2D=new SpecialQDoubleSpinBox;
	QS_K2D->setMaximum(1000000);
	QS_K2D->setMinimum(0);
	QS_K2D->setValue(1);
	QS_K2D->setDecimals(6);
	QS_K2D->setSingleStep(0.01);
	QS_K2D->setMaximumWidth(MaxWidth);
	QS_K2D->setKeyboardTracking(false);
	
    
    // Pressure
	QLabel *QL_P=new QLabel("Pressure");
	QL_P->setMaximumWidth(MaxWidth);
	QDoubleSpinBox *QS_P=new QDoubleSpinBox;
	QS_P->setMaximum(1000000);
	QS_P->setMinimum(-1000000);
	QS_P->setMaximumWidth(MaxWidth);
	QS_P->setReadOnly(true);
	QPushButton* QP_Apoptosis = new QPushButton(tr("Apoptosis"));
	QP_Apoptosis->setMaximumWidth(MaxWidth);
    
    QPushButton* QP_Division = new QPushButton(tr("Division"));
	QP_Division->setMaximumWidth(MaxWidth);
	
    //QApoptosis->resize(20,10);
	//QApoptosis->setMaximumWidth(MaxWidth);
	
	
	
    /************** define the locations of the labels and spin boxes ***************/
	layout->addWidget(QL_No,1, 0);
	layout->addWidget(QS_No,1, 1);
    layout->addWidget(QL_Type, 1,2);
	layout->addWidget(QS_Type, 1,3);
    layout->addWidget(QL_Ta, 2, 0);
    layout->addWidget(QS_Ta, 2, 1);
    layout->addWidget(QL_Tb, 2, 2);
    layout->addWidget(QS_Tb, 2, 3);
	layout->addWidget(QL_V0, 3, 0);
	layout->addWidget(QS_V0, 3, 1);
	layout->addWidget(QL_A0, 3, 2);
	layout->addWidget(QS_A0, 3, 3);
    layout->addWidget(QL_K,  4, 0);
	layout->addWidget(QS_K,  4, 1);
    layout->addWidget(QL_K2D,  4, 2);
	layout->addWidget(QS_K2D,  4, 3);
	layout->addWidget(QL_V, 5,0);
	layout->addWidget(QS_V, 5,1);
	layout->addWidget(QL_P, 5,2);
	layout->addWidget(QS_P, 5,3);
    layout->addWidget(QL_Aa, 6, 0);
    layout->addWidget(QS_Aa, 6, 1);
    layout->addWidget(QL_AR, 6, 2);
    layout->addWidget(QS_AR, 6, 3);
	layout->addWidget(QP_Apoptosis,7, 1);
	layout->addWidget(QP_Division, 7, 2);
	
    /************** connect the spin boxes to the functions, i.e. changes made in the window are translated into changes in the objects ***************/
	// connect incoming signals
    connect(integrator->T,SIGNAL(signal_currentVertex(int, Point&, Point&, Point&, Point&, double)),this, SLOT(slot_currentVertex(int, Point&, Point&, Point&, Point&, double)));
    
    connect(this,SIGNAL(signal_c_No(int)),QS_No,SLOT(setValue(int)));
    connect(this,SIGNAL(signal_c_Ta(double)),QS_Ta,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_c_Tb(double)),QS_Tb,SLOT(setValue(double)));

    connect(this,SIGNAL(signal_c_Aa(double)),QS_Aa,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_c_AR(double)),QS_AR,SLOT(setValue(double)));
    
    connect(this,SIGNAL(signal_c_V(double)),QS_V,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_c_V0(double)),QS_V0,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_c_A0(double)),QS_A0,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_c_Type(int)),QS_Type,SLOT(setValue(int)));
    connect(this,SIGNAL(signal_c_K(double)),QS_K,SLOT(setValue(double)));
    //connect(this,SIGNAL(signal_c_K2D(double)),QS_K2D,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_c_P(double)),QS_P,SLOT(setValue(double)));
    
    // connect outgoing signals
    connect(QS_No,SIGNAL(valueChanged(int)),integrator->T,SLOT(slot_setFirstCellNumber(int)));
    //connect(QS_Ta,SIGNAL(valueChanged(double)),integrator->T,SLOT(slot_setCellsTa(double)));
    //connect(QS_Tb,SIGNAL(valueChanged(double)),integrator->T,SLOT(slot_setCellsTb(double)));
    connect(QS_V0,SIGNAL(valueChanged(double)),integrator->T,SLOT(slot_setCellsPreferredVolume(double)));
    connect(QS_K,SIGNAL(valueChanged(double)),integrator->T,SLOT(slot_setCellsK(double)));
    connect(QS_Type,SIGNAL(valueChanged(int)), integrator->T, SLOT(slot_setCellsType(int)) );
    connect(QP_Apoptosis,SIGNAL(clicked()),integrator->T,SLOT(slot_firstCellApoptosis()));
    connect(QP_Division,SIGNAL(clicked()),integrator->T,SLOT(slot_firstCellRandomDivision()));
    
    /************** what is done here? ***************/
    layout->setColumnStretch(1, 1);
    //layout->setColumnStretch(2, 1);
    // layout is set
    CellBox->setLayout(layout);
}

// define the edge box (2nd tab, second box on the left hand side)
void MainWindow::createEdgeBox()
{
    EdgeBox = new QGroupBox();
    double MaxWidth=60;
    
    QGridLayout *layout = new QGridLayout;
    
    // Number
    QLabel *QL_No=new QLabel("Number");
    QL_No->setMaximumWidth(MaxWidth);
    QSpinBox *QS_No=new QSpinBox;
    QS_No->setMaximum(1000000);
    QS_No->setMinimum(0);
    QS_No->setFont(QFont("Times", 10));
    QS_No->setMaximumWidth(MaxWidth);
    
    // T_l
    QLabel *QL_Tl= new QLabel("T_l");
    QL_Tl->setMaximumWidth(MaxWidth);
    QDoubleSpinBox *QS_Tl=new QDoubleSpinBox;
    QS_Tl->setMaximum(10000);
    QS_Tl->setMinimum(-10000);
    QS_Tl->setMaximumWidth(MaxWidth);
    
    // G_a
    QLabel *QL_Ga= new QLabel("G_a");
    QL_Ga->setMaximumWidth(MaxWidth);
    QDoubleSpinBox *QS_Ga=new QDoubleSpinBox;
    QS_Ga->setMaximum(10000);
    QS_Ga->setMinimum(-10000);
    QS_Ga->setMaximumWidth(MaxWidth);
    
    // G_a
    QLabel *QL_Gb= new QLabel("G_b");
    QL_Gb->setMaximumWidth(MaxWidth);
    QDoubleSpinBox *QS_Gb=new QDoubleSpinBox;
    QS_Gb->setMaximum(10000);
    QS_Gb->setMinimum(-10000);
    QS_Gb->setMaximumWidth(MaxWidth);
    
    // length_a
    QLabel *QL_length_a= new QLabel("a_b");
    QL_length_a->setMaximumWidth(MaxWidth);
    QDoubleSpinBox *QS_length_a=new QDoubleSpinBox;
    QS_length_a->setMaximum(10000);
    QS_length_a->setMinimum(0);
    QS_length_a->setDecimals(2);
    QS_length_a->setMaximumWidth(MaxWidth);
    QS_length_a->setReadOnly(true);
    
    // length_b
    QLabel *QL_length_b= new QLabel("l_b");
    QL_length_b->setMaximumWidth(MaxWidth);
    QDoubleSpinBox *QS_length_b=new QDoubleSpinBox;
    QS_length_b->setMaximum(10000);
    QS_length_b->setMinimum(0);
    QS_length_b->setDecimals(2);
    QS_length_b->setMaximumWidth(MaxWidth);
    QS_length_b->setReadOnly(true);
    
    // cell 1
    QLabel *QL_c1= new QLabel("c1");
    QL_c1->setMaximumWidth(MaxWidth);
    QSpinBox *QS_c1=new QSpinBox;
    QS_c1->setMaximum(10000);
    QS_c1->setMaximumWidth(MaxWidth);
    QS_c1->setReadOnly(true);
    
    // cell 1
    QLabel *QL_c2= new QLabel("c1");
    QL_c2->setMaximumWidth(MaxWidth);
    QSpinBox *QS_c2=new QSpinBox;
    QS_c2->setMaximum(10000);
    QS_c2->setMaximumWidth(MaxWidth);
    QS_c2->setReadOnly(true);
    
    QPushButton* QP_T1 = new QPushButton(tr("remove Edge"));
    QP_T1->resize(10,10);
    
    layout->addWidget(QL_No, 0, 1);
    layout->addWidget(QS_No, 0, 2);
    layout->addWidget(QL_Tl, 0, 3);
    layout->addWidget(QS_Tl, 0, 4);
    layout->addWidget(QL_Ga, 1, 1);
    layout->addWidget(QS_Ga, 1, 2);
    layout->addWidget(QL_Gb, 1, 3);
    layout->addWidget(QS_Gb, 1, 4);
    layout->addWidget(QL_length_a, 2, 1);
    layout->addWidget(QS_length_a, 2, 2);
    layout->addWidget(QL_length_b, 2, 3);
    layout->addWidget(QS_length_b, 2, 4);
    layout->addWidget(QL_c1, 3, 1);
    layout->addWidget(QS_c1, 3, 2);
    layout->addWidget(QL_c2, 3, 3);
    layout->addWidget(QS_c2, 3, 4);
      
    layout->addWidget(QP_T1, 4, 1);
    
    
    layout->setHorizontalSpacing (0);
    layout->setVerticalSpacing (0);
    
    // outgoing signals
//    connect(QS_No,SIGNAL(valueChanged(int)),integrator->T,SLOT(slot_setFirstEdge_Number(int)));
    connect(QS_Tl,SIGNAL(valueChanged(double)),integrator->T,SLOT(slot_setFirstEdge_Tl(double)));
    connect(QS_Ga,SIGNAL(valueChanged(double)),integrator->T,SLOT(slot_setFirstEdge_Ga(double)));
    connect(QS_Gb,SIGNAL(valueChanged(double)),integrator->T,SLOT(slot_setFirstEdge_Gb(double)));
    connect(QP_T1,SIGNAL(clicked()),integrator->T,SLOT(slot_setFirstEdge_T1()));
    
    // incoming signals
    connect(integrator->T,SIGNAL(signal_currentEdge(int, int, int, double, double, double, double, double)),this, SLOT(slot_currentEdge(int, int, int, double, double, double, double, double)));
    
    connect(this,SIGNAL(signal_e_No(int)),QS_No,SLOT(setValue(int)));
    connect(this,SIGNAL(signal_e_Tl(double)),QS_Tl,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_e_length_a(double)),QS_length_a,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_e_length_b(double)),QS_length_b,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_e_Ga(double)),QS_Ga,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_e_Gb(double)),QS_Gb,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_e_c1(int)),QS_c1,SLOT(setValue(int)));
    connect(this,SIGNAL(signal_e_c2(int)),QS_c2,SLOT(setValue(int)));
        
    EdgeBox->setLayout(layout);
    
}

// define the vertex box (3rd tab, second box on the left hand side)
void MainWindow::createVertexBox()
{
    VertexBox = new QGroupBox();
    double MaxWidth=80;
    double MaxHeight=20;
    QGridLayout *layout = new QGridLayout;
    
    QLabel *LabelVertexNumber= new QLabel("Number of Vertex");
    LabelVertexNumber->setFont(QFont("Times", 12));
    LabelVertexNumber->setMaximumWidth(100000);
    
    QSpinBox *SBVertexNumber=new QSpinBox;
    //SBVertexNumber->setMaximum(tissueWindow->T->ListVertex.size());
    SBVertexNumber->setMinimum(0);
    SBVertexNumber->setMaximum(99999);
    SBVertexNumber->setMaximumWidth(MaxWidth);
    
    
    QLabel *xcolumn= new QLabel("x");
    xcolumn->setMaximumWidth(MaxWidth);
    xcolumn->setMaximumHeight(MaxHeight);
    
    QLabel *ycolumn= new QLabel("y");
    ycolumn->setMaximumWidth(MaxWidth);
    ycolumn->setMaximumHeight(MaxHeight);
    
    QLabel *zcolumn= new QLabel("z");
    zcolumn->setMaximumWidth(MaxWidth);
    zcolumn->setMaximumHeight(MaxHeight);
    
    QLabel *LabelCoordinatesApical= new QLabel("Apical coord:");
    LabelCoordinatesApical->setMaximumWidth(MaxWidth);
    
    
    QLabel *LabelCoordinatesBasal= new QLabel("Basal coord:");
    LabelCoordinatesBasal->setMaximumWidth(MaxWidth);
    
    QDoubleSpinBox *coordxa=new QDoubleSpinBox;
    coordxa->setMaximum(10000);
    coordxa->setMinimum(-10000);
    coordxa->setMaximumWidth(MaxWidth);
    
    QDoubleSpinBox *coordya=new QDoubleSpinBox;
    coordya->setMaximum(10000);
    coordya->setMinimum(-10000);
    coordya->setMaximumWidth(MaxWidth);
    
    QDoubleSpinBox *coordza=new QDoubleSpinBox;
    coordza->setMaximum(10000);
    coordza->setMinimum(-10000);
    coordza->setMaximumWidth(MaxWidth);
    
    QDoubleSpinBox *coordxb=new QDoubleSpinBox;
    coordxb->setMaximum(10000);
    coordxb->setMinimum(-10000);
    coordxb->setMaximumWidth(MaxWidth);
    
    QDoubleSpinBox *coordyb=new QDoubleSpinBox;
    coordyb->setMaximum(10000);
    coordyb->setMinimum(-10000);
    coordyb->setMaximumWidth(MaxWidth);
    
    QDoubleSpinBox *coordzb=new QDoubleSpinBox;
    coordzb->setMaximum(10000);
    coordzb->setMinimum(-10000);
    coordzb->setMaximumWidth(MaxWidth);
    
    QLabel *ApicalVertexFixed= new QLabel("Apical Fixed");
    
    QCheckBox * xaFixed=new QCheckBox(" ",VertexBox);
    xaFixed->setMaximumWidth(MaxWidth);
    xaFixed->setCheckState(Qt::Unchecked);
    xaFixed->setMaximumWidth(MaxWidth);
    
    QCheckBox * yaFixed=new QCheckBox(" ",VertexBox);
    yaFixed->setCheckState(Qt::Unchecked);
    yaFixed->setMaximumWidth(MaxWidth);
    
    QCheckBox * zaFixed=new QCheckBox(" ",VertexBox);
    zaFixed->setCheckState(Qt::Unchecked);
    zaFixed->setMaximumWidth(MaxWidth);
    
    QLabel *BasalVertexFixed= new QLabel("Basal fixed:");
    
    QCheckBox * xbFixed=new QCheckBox(" ",VertexBox);
    xbFixed->setMaximumWidth(MaxWidth);
    xbFixed->setCheckState(Qt::Unchecked);
    xbFixed->setMaximumWidth(MaxWidth);
    
    QCheckBox * ybFixed=new QCheckBox(" ",VertexBox);
    ybFixed->setCheckState(Qt::Unchecked);
    ybFixed->setMaximumWidth(MaxWidth);
    
    QCheckBox * zbFixed=new QCheckBox(" ",VertexBox);
    zbFixed->setCheckState(Qt::Unchecked);
    zbFixed->setMaximumWidth(MaxWidth);
    
    // display FORCES
    
    QLabel *LabelForcesApical= new QLabel("Apical forces:");
    LabelForcesApical->setMaximumWidth(MaxWidth);
    
    QLabel *LabelForcesBasal= new QLabel("Basal forces:");
    LabelForcesBasal->setMaximumWidth(MaxWidth);
    
    QDoubleSpinBox *forcesxa=new QDoubleSpinBox;
    forcesxa->setMaximum(10000);
    forcesxa->setMinimum(-10000);
    forcesxa->setMaximumWidth(MaxWidth);
    
    QDoubleSpinBox *forcesya=new QDoubleSpinBox;
    forcesya->setMaximum(10000);
    forcesya->setMinimum(-10000);
    forcesya->setMaximumWidth(MaxWidth);
    
    QDoubleSpinBox *forcesza=new QDoubleSpinBox;
    forcesza->setMaximum(10000);
    forcesza->setMinimum(-10000);
    forcesza->setMaximumWidth(MaxWidth);
    
    QDoubleSpinBox *forcesxb=new QDoubleSpinBox;
    forcesxb->setMaximum(10000);
    forcesxb->setMinimum(-10000);
    forcesxb->setMaximumWidth(MaxWidth);
    
    QDoubleSpinBox *forcesyb=new QDoubleSpinBox;
    forcesyb->setMaximum(10000);
    forcesyb->setMinimum(-10000);
    forcesyb->setMaximumWidth(MaxWidth);
    
    QDoubleSpinBox *forceszb=new QDoubleSpinBox;
    forceszb->setMaximum(10000);
    forceszb->setMinimum(-10000);
    forceszb->setMaximumWidth(MaxWidth);
    
    // show the tension inbetween
    QLabel *LabelTension= new QLabel("Line Tension");
    LabelTension->setMaximumWidth(100);
    
    QDoubleSpinBox *SBTension=new QDoubleSpinBox;
    SBTension->setMaximum(1000000);
    SBTension->setMinimum(-1000000);
    SBTension->setMaximumWidth(MaxWidth);
    
    QDoubleSpinBox *zForce=new QDoubleSpinBox;
    zForce->setMaximum(1000000);
    zForce->setMinimum(-1000000);
    zForce->setMaximumWidth(MaxWidth);
    
    // ARRANGE the layout
    layout->setVerticalSpacing (0);
    layout->addWidget(LabelVertexNumber,0, 1);
    layout->addWidget(SBVertexNumber,0, 2);
    
    layout->addWidget(xcolumn,1, 2);
    layout->addWidget(ycolumn,1, 3);
    layout->addWidget(zcolumn,1, 4);
    
    layout->addWidget(LabelCoordinatesApical,2, 1);
    layout->addWidget(coordxa,2, 2);
    layout->addWidget(coordya,2, 3);
    layout->addWidget(coordza,2, 4);
    
    
    layout->addWidget(LabelCoordinatesBasal,3, 1);
    layout->addWidget(coordxb,3, 2);
    layout->addWidget(coordyb,3, 3);
    layout->addWidget(coordzb,3, 4);
    
    
    layout->addWidget(ApicalVertexFixed,4, 1);
    layout->addWidget(xaFixed,4, 2);
    layout->addWidget(yaFixed,4, 3);
    layout->addWidget(zaFixed,4, 4);
    
    layout->addWidget(BasalVertexFixed,5, 1);
    layout->addWidget(xbFixed,5, 2);
    layout->addWidget(ybFixed,5, 3);
    layout->addWidget(zbFixed,5, 4);
    
    
    layout->addWidget(LabelForcesApical,6, 1);
    layout->addWidget(forcesxa,6, 2);
    layout->addWidget(forcesya,6, 3);
    layout->addWidget(forcesza,6, 4);
    
    layout->addWidget(LabelForcesBasal,7, 1);
    layout->addWidget(forcesxb,7, 2);
    layout->addWidget(forcesyb,7, 3);
    layout->addWidget(forceszb,7, 4);
    
    // manage incoming signals from the tissue that has changed
	connect(integrator->T,SIGNAL(signal_currentCell(int, double, double, int, double, double, double, double, double, double)),this, SLOT(slot_currentCell(int, double, double, int, double, double, double, double, double, double)));
    
    connect(this,SIGNAL(signal_v_no(int)),SBVertexNumber,SLOT(setValue(int)));
    //connect(this,SIGNAL(signal_v_G_l(double)),VertexNumber,SLOT(setValue(int)));
    connect(this,SIGNAL(signal_v_Ra_x(double)),coordxa,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_v_Ra_y(double)),coordya,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_v_Ra_z(double)),coordza,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_v_Rb_x(double)),coordxb,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_v_Rb_y(double)),coordyb,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_v_Rb_z(double)),coordzb,SLOT(setValue(double)));
    
    connect(this,SIGNAL(signal_v_Fa_x(double)),forcesxa,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_v_Fa_y(double)),forcesya,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_v_Fa_z(double)),forcesza,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_v_Fb_x(double)),forcesxb,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_v_Fb_y(double)),forcesyb,SLOT(setValue(double)));
    connect(this,SIGNAL(signal_v_Fb_z(double)),forceszb,SLOT(setValue(double)));
    
    connect(this,SIGNAL(signal_v_Fb_z(double)),forceszb,SLOT(setValue(double)));
    
    //	connect(coordy,SIGNAL(valueChanged(double)),tissueWindow,SLOT(setyCurrentVertex(double)));
    //	connect(tissueWindow,SIGNAL(xFixedCurrentVertexChanged(bool)),xFixed,SLOT(setChecked(bool)));
    //	connect(tissueWindow,SIGNAL(yFixedCurrentVertexChanged(bool)),yFixed,SLOT(setChecked(bool)));
    //	connect(showForces,SIGNAL(stateChanged(int)),tissueWindow,SLOT(setShowForces(int)));
    VertexBox->setLayout(layout);
}

//void MainWindow::createCystBox()
//{
//    CystBox = new QGroupBox();
//    double MaxWidth=80;
//    QGridLayout *layout = new QGridLayout;
//    
//    // apical surface tension
//    QLabel *QL_Ta= new QLabel("T_a");
//    QL_Ta->setFont(QFont("Times", 12));
//    QDoubleSpinBox *QS_Ta = new QDoubleSpinBox;
//    QS_Ta->setValue(10);
//    QS_Ta->setDecimals(7);
//	QS_Ta->setMinimum(-10000);
//	QS_Ta->setMaximum(10000);
//	QS_Ta->setSingleStep(1);
//    QS_Ta->setMinimum(0);
//    QS_Ta->setFont(QFont("Times", 12));
//    QS_Ta->setMaximumWidth(MaxWidth);
//    
//    // basal surface tension
//    QLabel *QL_Tb= new QLabel("T_b");
//    QL_Tb->setFont(QFont("Times", 12));
//    QDoubleSpinBox *QS_Tb = new QDoubleSpinBox;
//    QS_Tb->setValue(10);
//	QS_Tb->setMaximum(10000);
//	QS_Tb->setSingleStep(1);
//    QS_Tb->setMinimum(0);
//    QS_Tb->setFont(QFont("Times", 12));
//    QS_Tb->setMaximumWidth(MaxWidth);
//    
//    // preferred volume
//    QLabel *QL_V0= new QLabel("V0");
//    QL_V0->setFont(QFont("Times", 12));
//    QDoubleSpinBox *QS_V0 = new QDoubleSpinBox;
//    QS_V0->setValue(10);
//	QS_V0->setSingleStep(1);
//    QS_V0->setMinimum(0);
//    QS_V0->setFont(QFont("Times", 12));
//    QS_V0->setMaximumWidth(MaxWidth);
//    
//    // bulk modulus
//    QLabel *QL_K= new QLabel("K");
//    QL_K->setFont(QFont("Times", 12));
//    QDoubleSpinBox *QS_K = new QDoubleSpinBox;
//    QS_K->setDecimals(7);
//    QS_K->setValue(0.001);
//	QS_K->setSingleStep(0.001);
//    QS_K->setMinimum(0);
//    QS_K->setFont(QFont("Times", 12));
//    QS_K->setMaximumWidth(MaxWidth);
//    
//    // apical line tension
//    QLabel *QL_Ga= new QLabel("G_a");
//    QL_Ga->setFont(QFont("Times", 12));
//    QDoubleSpinBox *QS_Ga = new QDoubleSpinBox;
//    QS_Ga->setValue(10);
//	QS_Ga->setMaximum(10000);
//	QS_Ga->setSingleStep(1);
//    QS_Ga->setMinimum(0);
//    QS_Ga->setFont(QFont("Times", 12));
//    QS_Ga->setMaximumWidth(MaxWidth);
//    
//    // basal line tension
//    QLabel *QL_Gb= new QLabel("G_b");
//    QL_Gb->setFont(QFont("Times", 12));
//    QDoubleSpinBox *QS_Gb = new QDoubleSpinBox;
//    QS_Gb->setValue(10);
//	QS_Gb->setMaximum(10000);
//    QS_Gb->setSingleStep(1);
//    QS_Gb->setMinimum(0);
//    QS_Gb->setFont(QFont("Times", 12));
//    QS_Gb->setMaximumWidth(MaxWidth);
//    
//    // lateral surface tension
//    QLabel *QL_Tl= new QLabel("T_l");
//    QL_Tl->setFont(QFont("Times", 12));
//    QDoubleSpinBox *QS_Tl = new QDoubleSpinBox;
//    QS_Tl->setValue(10);
//    QS_Tl->setMaximum(10000);
//    QS_Tl->setSingleStep(1);
//    QS_Tl->setMinimum(0);
//    QS_Tl->setFont(QFont("Times", 12));
//    QS_Tl->setMaximumWidth(MaxWidth);
//    
//    
//    // prefactor in apical line tension for clone boundaries
//    QLabel *QL_factor_Ga= new QLabel("Factor G_a");
//    QL_factor_Ga->setFont(QFont("Times", 12));
//    QL_factor_Ga->setMaximumWidth(100);
//    QDoubleSpinBox *QS_factor_Ga = new QDoubleSpinBox;
//    QS_factor_Ga->setValue(1.);
//	QS_factor_Ga->setSingleStep(0.1);
//    QS_factor_Ga->setMinimum(0);
//    QS_factor_Ga->setFont(QFont("Times", 12));
//    QS_factor_Ga->setMaximumWidth(MaxWidth);
//    
//    // prefactor in basal line tension for clone boundaries
//    QLabel *QL_factor_Gb= new QLabel("Factor G_b");
//    QL_factor_Gb->setFont(QFont("Times", 12));
//    QL_factor_Gb->setMaximumWidth(100);
//    QDoubleSpinBox *QS_factor_Gb = new QDoubleSpinBox;
//    QS_factor_Gb->setValue(1.);
//	QS_factor_Gb->setSingleStep(0.1);
//    QS_factor_Gb->setMinimum(0);
//    QS_factor_Gb->setFont(QFont("Times", 12));
//    QS_factor_Gb->setMaximumWidth(MaxWidth);
//    
//    // prefactor in lateral surface tension at clone boundaries
//    QLabel *QL_factor_Tl= new QLabel("Factor T_l");
//    QL_factor_Tl->setFont(QFont("Times", 12));
//    QL_factor_Tl->setMaximumWidth(100);
//    QDoubleSpinBox *QS_factor_Tl = new QDoubleSpinBox;
//    QS_factor_Tl->setValue(1.);
//	QS_factor_Tl->setSingleStep(0.1);
//    QS_factor_Tl->setMinimum(0);
//    QS_factor_Tl->setFont(QFont("Times", 12));
//    QS_factor_Tl->setMaximumWidth(MaxWidth);
//    
//    // push button to mutate a round cyst in the middle of the tissue (works only if tissue is periodic)
//    QPushButton *QP_mutate_circle =new QPushButton(tr("circle"),tissueButtons);
//	QP_mutate_circle->setFont(QFont("Times", 12, QFont::Bold));
//	QP_mutate_circle->resize(10,10);
//	QP_mutate_circle->setMaximumWidth(MaxWidth);
//
//    // push button to mutate a stripe cyst in the middle of the tissue (works only if tissue is periodic)
//    QPushButton *QP_mutate_stripe =new QPushButton(tr("stripe"),tissueButtons);
//	QP_mutate_stripe->setFont(QFont("Times", 12, QFont::Bold));
//	QP_mutate_stripe->resize(10,10);
//	QP_mutate_stripe->setMaximumWidth(MaxWidth);
//
//    // push button to mutate two stripe cyst in the middle of the tissue (works only if tissue is periodic)
//    QPushButton *QP_mutate_two_stripes =new QPushButton(tr("2 stripes"),tissueButtons);
//	QP_mutate_two_stripes->setFont(QFont("Times", 12, QFont::Bold));
//	QP_mutate_two_stripes->resize(10,10);
//	QP_mutate_two_stripes->setMaximumWidth(MaxWidth);
//
//    // prefactor in basal line tension for clone boundaries
//    QLabel *QL_clone_size= new QLabel("clone size");
//    QL_clone_size->setFont(QFont("Times", 12));
//    QL_clone_size->setMaximumWidth(100);
//    
//    QS_clone_size = new QDoubleSpinBox;
//    QS_clone_size->setValue(1.);
//	QS_clone_size->setSingleStep(1);
//    QS_clone_size->setMinimum(0);
//    QS_clone_size->setMaximum(100000);
//    QS_clone_size->setFont(QFont("Times", 12));
//    QS_clone_size->setMaximumWidth(MaxWidth);
//    
//    connect(QS_Ta, SIGNAL(valueChanged(double)), integrator->T, SLOT(slot_Ta(double)));
//    connect(QS_Tb, SIGNAL(valueChanged(double)), integrator->T, SLOT(slot_Tb(double)));
//    connect(QS_V0, SIGNAL(valueChanged(double)), integrator->T, SLOT(slot_V0(double)));
//    connect(QS_K, SIGNAL(valueChanged(double)), integrator->T, SLOT(slot_K(double)));
//
//    connect(QS_Tl, SIGNAL(valueChanged(double)), integrator->T, SLOT(slot_Tl(double)));
//    connect(QS_Ga, SIGNAL(valueChanged(double)), integrator->T, SLOT(slot_Ga(double)));
//    connect(QS_Gb, SIGNAL(valueChanged(double)), integrator->T, SLOT(slot_Gb(double)));
//
//    connect(QS_factor_Ga, SIGNAL(valueChanged(double)), integrator->T, SLOT(slot_Ga_factor(double)));
//    connect(QS_factor_Gb, SIGNAL(valueChanged(double)), integrator->T, SLOT(slot_Gb_factor(double)));
//    connect(QS_factor_Tl, SIGNAL(valueChanged(double)), integrator->T, SLOT(slot_Tl_factor(double)));
//    
//    connect(QP_mutate_circle, SIGNAL(clicked()),this,SLOT(slot_mutate_cells_in_middle()));
//    connect(QP_mutate_stripe, SIGNAL(clicked()),this,SLOT(slot_mutate_cells_in_stripe()));
//    connect(QP_mutate_two_stripes, SIGNAL(clicked()),this,SLOT(slot_mutate_cells_in_two_stripes()));
//    
//    connect(this, SIGNAL(signal_mutate_cells_in_middle(int)),integrator->T,SLOT(mutate_cells_in_middle(int)));
//    connect(this, SIGNAL(signal_mutate_cells_in_stripe(double)),integrator->T,SLOT(mutate_cells_in_stripe(double)));
//    connect(this, SIGNAL(signal_mutate_cells_in_two_stripes(double)),integrator->T,SLOT(mutate_cells_in_two_stripes(double)));
//    
//    connect(this, SIGNAL(signal_Ta(double)), QS_Ta, SLOT(setValue(double)));
//    connect(this, SIGNAL(signal_Tb(double)), QS_Tb, SLOT(setValue(double)));
//    connect(this, SIGNAL(signal_Tl(double)), QS_Tl, SLOT(setValue(double)));
//    connect(this, SIGNAL(signal_Ga(double)), QS_Ga, SLOT(setValue(double)));
//    connect(this, SIGNAL(signal_Gb(double)), QS_Gb, SLOT(setValue(double)));
//    connect(this, SIGNAL(signal_factor_Tl(double)), QS_factor_Tl, SLOT(setValue(double)));
//    connect(this, SIGNAL(signal_factor_Ga(double)), QS_factor_Ga, SLOT(setValue(double)));
//    connect(this, SIGNAL(signal_factor_Gb(double)), QS_factor_Gb, SLOT(setValue(double)));
//
//    //connect(integrator, SIGNAL(signal_update_cyst_data(double,double,double,double,double,double,double,double,double)), this, SLOT(slot_update_cyst_data(double,double,double,double,double,double,double,double,double,double)));
//    connect(integrator->T, SIGNAL(signal_update_cyst_data(double,double,double,double,double,double,double,double,double,double)), this, SLOT(slot_update_cyst_data(double,double,double,double,double,double,double,double,double,double,double)));
//    
//    //connect(QP_mutate_circle, SIGNAL(clicked()),integrator->T,SLOT(mutate_cells_in_middle(double)));
//    //connect(QP_mutate_stripe, SIGNAL(clicked()),integrator->T,SLOT(mutate_stripe_in_middle(double)));
//    
//    layout->addWidget(QL_Ta,0, 1);
//    layout->addWidget(QS_Ta,0, 2);
//    layout->addWidget(QL_Tb,0, 3);
//    layout->addWidget(QS_Tb,0, 4);
//    
//    layout->addWidget(QL_V0,1, 1);
//    layout->addWidget(QS_V0,1, 2);
//    layout->addWidget(QL_K, 1, 3);
//    layout->addWidget(QS_K, 1, 4);
//    
//    layout->addWidget(QL_Ga,2, 1);
//    layout->addWidget(QS_Ga,2, 2);
//    layout->addWidget(QL_Gb,2, 3);
//    layout->addWidget(QS_Gb,2, 4);
//    
//    layout->addWidget(QL_Tl,3, 1);
//    layout->addWidget(QS_Tl,3, 2);
//    
//    layout->addWidget(QL_factor_Ga, 4, 1);
//    layout->addWidget(QS_factor_Ga, 4, 2);
//    layout->addWidget(QL_factor_Gb, 4, 3);
//    layout->addWidget(QS_factor_Gb, 4, 4);
//    layout->addWidget(QL_factor_Tl, 5, 1);
//    layout->addWidget(QS_factor_Tl, 5, 2);
//    
//    layout->addWidget(QP_mutate_circle, 6, 1);
//    layout->addWidget(QP_mutate_stripe, 6, 2);
//    layout->addWidget(QP_mutate_two_stripes, 6, 3);
//    layout->addWidget(QS_clone_size,    6, 4);
//    
//    CystBox->setLayout(layout);
//}


void MainWindow::createMechanicsBox()
{
    MechanicsBox = new QGroupBox();
    double MaxWidth=70;
    QGridLayout *layout = new QGridLayout;
    
      // header of the properties table
    QLabel *QL_headerProp = new QLabel("Prop");
    QL_headerProp->setMaximumWidth(MaxWidth/2);
    QLabel *QL_headerType1 = new QLabel("Type 1");
    QL_headerType1->setMaximumWidth(MaxWidth);
    QLabel *QL_headerType2 = new QLabel("Type 2");
    QL_headerType2->setMaximumWidth(MaxWidth);
    
    QLabel *QL_headerEdge11 = new QLabel("1-1");
    QL_headerEdge11->setMaximumWidth(MaxWidth/2);
    QLabel *QL_headerEdge22 = new QLabel("2-2");
    QL_headerEdge22->setMaximumWidth(MaxWidth/2);
    QLabel *QL_headerEdge12 = new QLabel("1-2");
    QL_headerEdge12->setMaximumWidth(MaxWidth/2);
    
    // write the property names underneath each other on the left side of the box
    QLabel *QL_K= new QLabel("K");
    QLabel *QL_V0= new QLabel("V0");
    QLabel *QL_Ta= new QLabel("Ta");
    QLabel *QL_Ta0= new QLabel("T0a");
    QLabel *QL_A0a= new QLabel("A0a");
    QLabel *QL_Tb= new QLabel("Tb");
    QLabel *QL_Tb0= new QLabel("T0b");
    QLabel *QL_A0b= new QLabel("A0b");
    QLabel *QL_Tl= new QLabel("Tl");
    QLabel *QL_Ga= new QLabel("Ga");
    QLabel *QL_Gb= new QLabel("Gb");

    
    // apical surface tensions
    // type 1
    QS_Ta1 = new QDoubleSpinBox;
    QS_Ta1->setValue(1);
    QS_Ta1->setDecimals(3);
	QS_Ta1->setMinimum(-1000000);
	QS_Ta1->setMaximum(1000000);
	QS_Ta1->setSingleStep(1);
    QS_Ta1->setMaximumWidth(MaxWidth);
    // type 2
    QS_Ta2 = new QDoubleSpinBox;
    QS_Ta2->setValue(1);
    QS_Ta2->setDecimals(3);
	QS_Ta2->setMinimum(-1000000);
	QS_Ta2->setMaximum(1000000);
	QS_Ta2->setSingleStep(1);
    QS_Ta2->setMaximumWidth(MaxWidth);

    // apical surface tensions
    // type 1
    QS_Tb1 = new QDoubleSpinBox;
    QS_Tb1->setValue(1);
    QS_Tb1->setDecimals(3);
	QS_Tb1->setMinimum(-1000000);
	QS_Tb1->setMaximum(1000000);
	QS_Tb1->setSingleStep(1);
    QS_Tb1->setMaximumWidth(MaxWidth);
    // type 2
    QS_Tb2 = new QDoubleSpinBox;
    QS_Tb2->setValue(1);
    QS_Tb2->setDecimals(3);
	QS_Tb2->setMinimum(-1000000);
	QS_Tb2->setMaximum(1000000);
	QS_Tb2->setSingleStep(1);
    QS_Tb2->setMaximumWidth(MaxWidth);
    
    // T0_a
    // type 1
    QS_Ta01 = new QDoubleSpinBox;
    QS_Ta01->setValue(1);
    QS_Ta01->setDecimals(3);
	QS_Ta01->setMinimum(-1000000);
	QS_Ta01->setMaximum(1000000);
	QS_Ta01->setSingleStep(1);
    QS_Ta01->setMaximumWidth(MaxWidth);
    // type 2
    QS_Ta02 = new QDoubleSpinBox;
    QS_Ta02->setValue(1);
    QS_Ta02->setDecimals(3);
	QS_Ta02->setMinimum(-1000000);
	QS_Ta02->setMaximum(1000000);
	QS_Ta02->setSingleStep(1);
    QS_Ta02->setMaximumWidth(MaxWidth);
    
    // A0_a
    // type 1
    QS_A0a1 = new QDoubleSpinBox;
	QS_A0a1->setMinimum(0);
	QS_A0a1->setMaximum(100000);
    QS_A0a1->setValue(100);
    QS_A0a1->setDecimals(3);
	QS_A0a1->setSingleStep(1);
    QS_A0a1->setMaximumWidth(MaxWidth);
    // type 2
    QS_A0a2 = new QDoubleSpinBox;
	QS_A0a2->setMinimum(0);
	QS_A0a2->setMaximum(100000);
    QS_A0a2->setValue(100);
    QS_A0a2->setDecimals(3);
	QS_A0a2->setSingleStep(1);
    QS_A0a2->setMaximumWidth(MaxWidth);

    // T0_b
    // type 1
    QS_Tb01 = new QDoubleSpinBox;
    QS_Tb01->setValue(1);
    QS_Tb01->setDecimals(3);
	QS_Tb01->setMinimum(-1000000);
	QS_Tb01->setMaximum(1000000);
	QS_Tb01->setSingleStep(1);
    QS_Tb01->setMaximumWidth(MaxWidth);
    // type 2
    QS_Tb02 = new QDoubleSpinBox;
    QS_Tb02->setValue(1);
    QS_Tb02->setDecimals(3);
	QS_Tb02->setMinimum(-1000000);
	QS_Tb02->setMaximum(1000000);
	QS_Tb02->setSingleStep(1);
    QS_Tb02->setMaximumWidth(MaxWidth);
    
    // A0_b
    // type 1
    QS_A0b1 = new QDoubleSpinBox;
	QS_A0b1->setMinimum(0);
	QS_A0b1->setMaximum(100000);
    QS_A0b1->setValue(100);
    QS_A0b1->setDecimals(3);
	QS_A0b1->setSingleStep(1);
    QS_A0b1->setMaximumWidth(MaxWidth);
    // type 2
    QS_A0b2 = new QDoubleSpinBox;
	QS_A0b2->setMinimum(0);
	QS_A0b2->setMaximum(100000);
    QS_A0b2->setValue(100);
    QS_A0b2->setDecimals(3);
	QS_A0b2->setSingleStep(1);
    QS_A0b2->setMaximumWidth(MaxWidth);

    // preferred volumes
    // type 1
    QS_V01 = new QDoubleSpinBox;
	QS_V01->setMinimum(-1000000);
	QS_V01->setMaximum(1000000);
    QS_V01->setValue(10000);
    QS_V01->setDecimals(3);
	QS_V01->setSingleStep(1);
    QS_V01->setMaximumWidth(MaxWidth);
    // type 2
    QS_V02 = new QDoubleSpinBox;
	QS_V02->setMinimum(-1000000);
	QS_V02->setMaximum(1000000);
    QS_V02->setValue(10000);
    QS_V02->setDecimals(3);
	QS_V02->setSingleStep(1);
    QS_V02->setMaximumWidth(MaxWidth);
    
    // bulk moduli
    // type 1
    QS_K1 = new QDoubleSpinBox;
    QS_K1->setValue(0.001);
    QS_K1->setDecimals(5);
	QS_K1->setMinimum(0);
	QS_K1->setMaximum(100);
	QS_K1->setSingleStep(0.01);
    QS_K1->setMaximumWidth(MaxWidth);
    // type 2
    QS_K2 = new QDoubleSpinBox;
    QS_K2->setValue(0.001);
    QS_K2->setDecimals(5);
	QS_K2->setMinimum(0);
	QS_K2->setMaximum(100);
	QS_K2->setSingleStep(0.01);
    QS_K2->setMaximumWidth(MaxWidth);
    

    // apical line tensions
    // 1-1
    QS_Ga11 = new QDoubleSpinBox;
	QS_Ga11->setMinimum(-100000);
	QS_Ga11->setMaximum(100000);
    QS_Ga11->setValue(30);
    QS_Ga11->setDecimals(1);
	QS_Ga11->setSingleStep(5);
    QS_Ga11->setMaximumWidth(MaxWidth);
    // 2-2
    QS_Ga22 = new QDoubleSpinBox;
	QS_Ga22->setMinimum(-100000);
	QS_Ga22->setMaximum(100000);
    QS_Ga22->setValue(30);
    QS_Ga22->setDecimals(1);
	QS_Ga22->setSingleStep(5);
    QS_Ga22->setMaximumWidth(MaxWidth);
    // 1-2
    QS_Ga12 = new QDoubleSpinBox;
	QS_Ga12->setMinimum(-100000);
	QS_Ga12->setMaximum(100000);
    QS_Ga12->setValue(30);
    QS_Ga12->setDecimals(1);
	QS_Ga12->setSingleStep(5);
    QS_Ga12->setMaximumWidth(MaxWidth);

    // basal line tensions
    // 1-1
    QS_Gb11 = new QDoubleSpinBox;
	QS_Gb11->setMinimum(-100000);
	QS_Gb11->setMaximum(100000);
    QS_Gb11->setValue(30);
    QS_Gb11->setDecimals(1);
	QS_Gb11->setSingleStep(5);
    QS_Gb11->setMaximumWidth(MaxWidth);
    // 2-2
    QS_Gb22 = new QDoubleSpinBox;
	QS_Gb22->setMinimum(-100000);
	QS_Gb22->setMaximum(100000);
    QS_Gb22->setValue(30);
    QS_Gb22->setDecimals(1);
	QS_Gb22->setSingleStep(5);
    QS_Gb22->setMaximumWidth(MaxWidth);
    // 1-2
    QS_Gb12 = new QDoubleSpinBox;
	QS_Gb12->setMinimum(-100000);
	QS_Gb12->setMaximum(100000);
    QS_Gb12->setValue(30);
    QS_Gb12->setDecimals(1);
	QS_Gb12->setSingleStep(5);
    QS_Gb12->setMaximumWidth(MaxWidth);

    // lateral surface tensions
    // 1-1
    QS_Tl11 = new QDoubleSpinBox;
	QS_Tl11->setMinimum(-100000);
	QS_Tl11->setMaximum(100000);
    QS_Tl11->setValue(1);
    QS_Tl11->setDecimals(3);
	QS_Tl11->setSingleStep(1);
    QS_Tl11->setMaximumWidth(MaxWidth);
    // 2-2
    QS_Tl22 = new QDoubleSpinBox;
	QS_Tl22->setMinimum(-1000);
	QS_Tl22->setMaximum(1000);
    QS_Tl22->setValue(1);
    QS_Tl22->setDecimals(3);
	QS_Tl22->setSingleStep(1);
    QS_Tl22->setMaximumWidth(MaxWidth);
    // 1-2
    QS_Tl12 = new QDoubleSpinBox;
	QS_Tl12->setMinimum(-1000);
	QS_Tl12->setMaximum(1000);
    QS_Tl12->setValue(1);
    QS_Tl12->setDecimals(3);
	QS_Tl12->setSingleStep(1);
    QS_Tl12->setMaximumWidth(MaxWidth);

    // external force on the system size
	QL_Text=new QLabel("T_ext");
	QS_Text=new QDoubleSpinBox;
	QS_Text->setMinimum(-999999999);
	QS_Text->setMaximum(999999999);
	QS_Text->setMaximumWidth(MaxWidth);
    
    // fix system size
    QC_fixLx = new QCheckBox("fix Lx");
    QC_fixLx->setCheckState(Qt::Unchecked);
    
    QC_fixLy = new QCheckBox("fix Ly");
    QC_fixLy->setCheckState(Qt::Unchecked);

    
    // set the layout
    layout->addWidget(QL_headerProp, 0, 1);
    layout->addWidget(QL_headerType1,0, 2);
    layout->addWidget(QL_headerType2,0, 3);
    
    layout->addWidget(QL_headerEdge11, 9, 2);
    layout->addWidget(QL_headerEdge22, 9, 3);
    layout->addWidget(QL_headerEdge12, 9, 4);

    layout->addWidget(QL_K,    1, 1); layout->addWidget(QS_K1,     1, 2); layout->addWidget(QS_K2,    1, 3);
    layout->addWidget(QL_V0,   2, 1); layout->addWidget(QS_V01,    2, 2); layout->addWidget(QS_V02,   2, 3);
    layout->addWidget(QL_Ta,   3, 1); layout->addWidget(QS_Ta1,    3, 2); layout->addWidget(QS_Ta2,   3, 3);
    layout->addWidget(QL_Ta0,  4, 1); layout->addWidget(QS_Ta01,   4, 2); layout->addWidget(QS_Ta02,  4, 3);
    layout->addWidget(QL_A0a,  5, 1); layout->addWidget(QS_A0a1,   5, 2); layout->addWidget(QS_A0a2,  5, 3);
    layout->addWidget(QL_Tb,   6, 1); layout->addWidget(QS_Tb1,    6, 2); layout->addWidget(QS_Tb2,   6, 3);
    layout->addWidget(QL_Tb0,  7, 1); layout->addWidget(QS_Tb01,   7, 2); layout->addWidget(QS_Tb02,  7, 3);
    layout->addWidget(QL_A0b,  8, 1); layout->addWidget(QS_A0b1,   8, 2); layout->addWidget(QS_A0b2,  8, 3);
    
    
    layout->addWidget(QL_Tl,   10, 1); layout->addWidget(QS_Tl11, 10, 2); layout->addWidget(QS_Tl22,   10, 3); layout->addWidget(QS_Tl12,  10, 4);
    layout->addWidget(QL_Ga,   11, 1); layout->addWidget(QS_Ga11, 11, 2); layout->addWidget(QS_Ga22,   11, 3); layout->addWidget(QS_Ga12,  11, 4);
    layout->addWidget(QL_Gb,   12, 1); layout->addWidget(QS_Gb11, 12, 2); layout->addWidget(QS_Gb22,   12, 3); layout->addWidget(QS_Gb12,  12, 4);
    
    layout->addWidget(QL_Text,    13, 1); layout->addWidget(QS_Text,  13, 2);
    layout->addWidget(QC_fixLx,   13, 3); layout->addWidget(QC_fixLy, 13, 4);

    
    MechanicsBox->setLayout(layout);

    connect(QS_K1, SIGNAL(valueChanged(double)), this, SLOT(changeK1(double)));
    connect(QS_K2, SIGNAL(valueChanged(double)), this, SLOT(changeK2(double)));
    connect(QS_V01, SIGNAL(valueChanged(double)), this, SLOT(changeV01(double)));
    connect(QS_V02, SIGNAL(valueChanged(double)), this, SLOT(changeV02(double)));
    
    
    connect(QS_Ta1, SIGNAL(valueChanged(double)), this, SLOT(changeTa1(double)));
    connect(QS_Ta2, SIGNAL(valueChanged(double)), this, SLOT(changeTa2(double)));
    connect(QS_Ta01, SIGNAL(valueChanged(double)), this, SLOT(changeTa01(double)));
    connect(QS_Ta02, SIGNAL(valueChanged(double)), this, SLOT(changeTa02(double)));
    connect(QS_A0a1, SIGNAL(valueChanged(double)), this, SLOT(changeAa01(double)));
    connect(QS_A0a2, SIGNAL(valueChanged(double)), this, SLOT(changeAa02(double)));
    
    connect(QS_Tb1, SIGNAL(valueChanged(double)), this, SLOT(changeTb1(double)));
    connect(QS_Tb2, SIGNAL(valueChanged(double)), this, SLOT(changeTb2(double)));
    connect(QS_Tb01, SIGNAL(valueChanged(double)), this, SLOT(changeTb01(double)));
    connect(QS_Tb02, SIGNAL(valueChanged(double)), this, SLOT(changeTb02(double)));
    connect(QS_A0b1, SIGNAL(valueChanged(double)), this, SLOT(changeAb01(double)));
    connect(QS_A0b2, SIGNAL(valueChanged(double)), this, SLOT(changeAb02(double)));

    connect(QS_Tl11, SIGNAL(valueChanged(double)), this, SLOT(changeTl11(double)));
    connect(QS_Tl22, SIGNAL(valueChanged(double)), this, SLOT(changeTl22(double)));
    connect(QS_Tl12, SIGNAL(valueChanged(double)), this, SLOT(changeTl21(double)));
    
    connect(QS_Ga11, SIGNAL(valueChanged(double)), this, SLOT(changeGa11(double)));
    connect(QS_Ga22, SIGNAL(valueChanged(double)), this, SLOT(changeGa22(double)));
    connect(QS_Ga12, SIGNAL(valueChanged(double)), this, SLOT(changeGa21(double)));

    connect(QS_Gb11, SIGNAL(valueChanged(double)), this, SLOT(changeGb11(double)));
    connect(QS_Gb22, SIGNAL(valueChanged(double)), this, SLOT(changeGb22(double)));
    connect(QS_Gb12, SIGNAL(valueChanged(double)), this, SLOT(changeGb21(double)));
   
    connect(QS_Text,SIGNAL(valueChanged(double)),integrator->T,SLOT(slot_setText(double)));
    connect(QC_fixLx,SIGNAL(stateChanged(int)),integrator->T,SLOT(slot_setFixLx(int)));
    connect(QC_fixLy,SIGNAL(stateChanged(int)),integrator->T,SLOT(slot_setFixLy(int)));
    connect(this,SIGNAL(signal_setText(double)),QS_Text,SLOT(setValue(double)));

    connect(integrator->T, SIGNAL(updateMainWindowData()), this, SLOT(updateMechanicalData()));
    connect(integrator, SIGNAL(updateMainWindowData()), this, SLOT(updateMechanicalData()));
}

void MainWindow::createIntegrationBox()
{
    IntegrationBox = new QGroupBox();
    double MaxWidth=80;
    QGridLayout *layout = new QGridLayout;
    
    /* initialize all spin boxes used in the box (will change values in Integrator.minPar) */
    QLabel *QL_minimizationMode = new QLabel("minimizationMode");
    QL_minimizationMode->setMaximumWidth(MaxWidth);
    QS_minimizationMode = new QSpinBox();
    QS_minimizationMode->setValue(2);
	QS_minimizationMode->setSingleStep(1);
    QS_minimizationMode->setMinimum(0);
    QS_minimizationMode->setMaximum(2);
    QS_minimizationMode->setMaximumWidth(MaxWidth);
    
    QLabel *QL_ftol = new QLabel("ftol");
    QL_ftol->setMaximumWidth(MaxWidth);
    QS_ftol = new QDoubleSpinBox();
    QS_ftol->setValue(1e-12);
	QS_ftol->setSingleStep(1e-12);
    QS_ftol->setDecimals(10);
    QS_ftol->setMinimum(0);
    QS_ftol->setMaximum(2);
    QS_ftol->setMaximumWidth(MaxWidth);
    
    QLabel *QL_tol = new QLabel("tol");
    QL_tol->setMaximumWidth(MaxWidth);
    QS_tol = new QDoubleSpinBox();
    QS_tol->setValue(1e-5);
    QS_tol->setDecimals(6);
	QS_tol->setSingleStep(1e-5);
    QS_tol->setMinimum(0);
    QS_tol->setMaximum(2);
    QS_tol->setMaximumWidth(MaxWidth);
    
    QLabel *QL_GTOL = new QLabel("GTOL");
    QL_GTOL->setMaximumWidth(MaxWidth);
    QS_GTOL = new QDoubleSpinBox();
    QS_GTOL->setValue(1e-8);
    QS_GTOL->setDecimals(10);
	QS_GTOL->setSingleStep(1e-8);
    QS_GTOL->setMinimum(0);
    QS_GTOL->setMaximum(2);
    QS_GTOL->setMaximumWidth(MaxWidth);
    
    QLabel *QL_rateT1check = new QLabel("checkForT1");
    QL_rateT1check->setMaximumWidth(MaxWidth);
    QS_rateT1check = new QSpinBox();
    QS_rateT1check->setValue(10);
	QS_rateT1check->setSingleStep(1);
    QS_rateT1check->setMinimum(0);
    QS_rateT1check->setMaximum(1000);
    QS_rateT1check->setMaximumWidth(MaxWidth);
    
    QLabel *QL_ITMAX = new QLabel("maxIterations");
    QL_ITMAX->setMaximumWidth(MaxWidth);
    QS_ITMAX = new QSpinBox();
    QS_ITMAX->setValue(1000);
	QS_ITMAX->setSingleStep(100);
    QS_ITMAX->setMinimum(1);
    QS_ITMAX->setMaximum(1000000);
    QS_ITMAX->setMaximumWidth(MaxWidth);

    QLabel *QL_noise_stdDev = new QLabel("noise_stdDev");
    QL_noise_stdDev->setMaximumWidth(MaxWidth);
    QS_noise_stdDev = new QDoubleSpinBox();
    QS_noise_stdDev->setValue(0.1);
	QS_noise_stdDev->setSingleStep(0.1);
    QS_noise_stdDev->setDecimals(2);
    QS_noise_stdDev->setMinimum(0);
    QS_noise_stdDev->setMaximum(100);
    QS_noise_stdDev->setMaximumWidth(MaxWidth);

    QLabel *QL_numberOfNoiseApplications = new QLabel("noise_noApplications");
    QL_numberOfNoiseApplications->setMaximumWidth(MaxWidth);
    QS_numberOfNoiseApplications = new QSpinBox();
    QS_numberOfNoiseApplications->setValue(10);
	QS_numberOfNoiseApplications->setSingleStep(1);
    QS_numberOfNoiseApplications->setMinimum(0);
    QS_numberOfNoiseApplications->setMaximum(1000000);
    QS_numberOfNoiseApplications->setMaximumWidth(MaxWidth);

    QLabel *QL_minimization_including_SystemSize = new QLabel("minimizeSystemSize");
    QC_minimization_including_SystemSize = new QCheckBox();
    QC_minimization_including_SystemSize->setCheckState(Qt::Unchecked);
    
    // set the layout
    layout->addWidget(QL_minimizationMode, 0, 0);layout->addWidget(QS_minimizationMode, 0, 1);
    layout->addWidget(QL_ftol, 1, 0);layout->addWidget(QS_ftol, 1, 1);
    layout->addWidget(QL_GTOL, 1, 2);layout->addWidget(QS_GTOL, 1, 3);
    layout->addWidget(QL_tol, 2, 0);layout->addWidget(QS_tol, 2, 1);
    layout->addWidget(QL_noise_stdDev, 3, 0);layout->addWidget(QS_noise_stdDev, 3, 1);
    layout->addWidget(QL_numberOfNoiseApplications, 3, 2);layout->addWidget(QS_numberOfNoiseApplications, 3, 3);
    layout->addWidget(QL_ITMAX, 4, 0);layout->addWidget(QS_ITMAX, 4, 1);
    layout->addWidget(QL_rateT1check, 4, 2);layout->addWidget(QS_rateT1check, 4, 3);
    
    IntegrationBox->setLayout(layout);
    
    // outoing signals
    connect(QS_minimizationMode, SIGNAL(valueChanged(int)), integrator, SLOT(set_MinimizationMode(int)));
    connect(QS_rateT1check, SIGNAL(valueChanged(int)), integrator, SLOT(set_rateT1check(int)));
    connect(QS_ITMAX, SIGNAL(valueChanged(int)), integrator, SLOT(set_ITMAX(int)));
    
    connect(QS_ftol, SIGNAL(valueChanged(double)), integrator, SLOT(set_ftol(double)));
    connect(QS_tol, SIGNAL(valueChanged(double)), integrator, SLOT(set_tol(double)));
    connect(QS_GTOL, SIGNAL(valueChanged(double)), integrator, SLOT(set_GTOL(double)));
    connect(QS_noise_stdDev, SIGNAL(valueChanged(double)), integrator, SLOT(set_noise_stdDev(double)));
    connect(QS_numberOfNoiseApplications, SIGNAL(valueChanged(int)), integrator, SLOT(set_noise_numberOfNoiseApplications(int)));
    
    // incoming signal
    connect(integrator, SIGNAL(set_integration_properties(int,double,double, double,int,double,double,int)), this, SLOT(set_integration_properties(int,double,double, double,int,double,double,int)));
    
    // initialize the values
    set_integration_properties(integrator->minPar.MinimizationMode, integrator->minPar.ftol, integrator->minPar.tol, integrator->minPar.GTOL, integrator->minPar.ITMAX, integrator->minPar.noise_numberOfNoiseApplications, integrator->minPar.noise_stdDev, integrator->minPar.rateT1check);
}


void MainWindow::createGastrulationBox()
{
    GastrulationBox = new QGroupBox();
    double MaxWidth=80;
    QGridLayout *layout = new QGridLayout;

    // boxes to define the extensions in the different directions of the ellipse
    QLabel *QL_Ellipse_x = new QLabel("x- & y-extension of ellipsoid");
    QS_Ellipse_x = new QDoubleSpinBox();
    QS_Ellipse_x->setValue(0);
    QS_Ellipse_x->setDecimals(2);
	QS_Ellipse_x->setSingleStep(1e-2);
    QS_Ellipse_x->setMinimum(0);
    QS_Ellipse_x->setMaximum(10000);
    QS_Ellipse_x->setMaximumWidth(MaxWidth);
    QLabel *QL_Ellipse_y = new QLabel("y-extension of ellipsoid");
    QS_Ellipse_y = new QDoubleSpinBox();
    QS_Ellipse_y->setValue(0);
    QS_Ellipse_y->setDecimals(2);
	QS_Ellipse_y->setSingleStep(1e-2);
    QS_Ellipse_y->setMinimum(0);
    QS_Ellipse_y->setMaximum(10000);
    QS_Ellipse_y->setMaximumWidth(MaxWidth);
    QLabel *QL_Ellipse_z = new QLabel("z-extension of ellipsoid");
    QS_Ellipse_z = new QDoubleSpinBox();
    QS_Ellipse_z->setValue(0);
    QS_Ellipse_z->setDecimals(2);
	QS_Ellipse_z->setSingleStep(1e-2);
    QS_Ellipse_z->setMinimum(0);
    QS_Ellipse_z->setMaximum(10000);
    QS_Ellipse_z->setMaximumWidth(MaxWidth);
    QLabel *QL_Ellipse_k = new QLabel("K of ellipsoid");
    QS_Ellipse_k = new QDoubleSpinBox();
    QS_Ellipse_k->setValue(0);
    QS_Ellipse_k->setDecimals(2);
	QS_Ellipse_k->setSingleStep(1e-2);
    QS_Ellipse_k->setMinimum(0);
    QS_Ellipse_k->setMaximum(10000);
    QS_Ellipse_k->setMaximumWidth(MaxWidth);

    // set the lumen properties
    QLabel *QL_Lumen_V0 = new QLabel("V0 of Lumen");
    QS_Lumen_V0 = new QDoubleSpinBox();
    QS_Lumen_V0->setValue(1000000);
    QS_Lumen_V0->setDecimals(0);
	QS_Lumen_V0->setSingleStep(1);
    QS_Lumen_V0->setMinimum(0);
    QS_Lumen_V0->setMaximum(100000000);
    QS_Lumen_V0->setMaximumWidth(MaxWidth);

    QLabel *QL_Lumen_K = new QLabel("K of Lumen");
    QS_Lumen_K = new QDoubleSpinBox();
    QS_Lumen_K->setValue(0);
    QS_Lumen_K->setDecimals(10);
	QS_Lumen_K->setSingleStep(0.00001);
    QS_Lumen_K->setMinimum(0);
    QS_Lumen_K->setMaximum(10);
    QS_Lumen_K->setMaximumWidth(MaxWidth);

    QLabel *QL_Lumen_currentVolume = new QLabel("Current volume of Lumen");
    QS_Lumen_currentVolume = new QDoubleSpinBox();
    QS_Lumen_currentVolume->setValue(0);
    QS_Lumen_currentVolume->setDecimals(1);
    QS_Lumen_currentVolume->setMaximum(100000000);
    QS_Lumen_currentVolume->setMaximumWidth(MaxWidth);
    
    QLabel *QL_K2D = new QLabel("K_2D cells");
    QS_K2D = new QDoubleSpinBox();
    QS_K2D->setValue(0);
    QS_K2D->setDecimals(1);
    QS_K2D->setMaximum(100000000);
    QS_K2D->setMinimum(-100000000);
    QS_K2D->setMaximumWidth(MaxWidth);
    
    QLabel *QL_A0 = new QLabel("A0 cells");
    QS_A0 = new QDoubleSpinBox();
    QS_A0->setValue(0);
    QS_A0->setDecimals(1);
    QS_A0->setMaximum(100000000);
    QS_A0->setMinimum(-100000000);
    QS_A0->setMaximumWidth(MaxWidth);

    QLabel *QL_bulkModulus = new QLabel("bulk modulus");
    QS_bulkModulus = new QDoubleSpinBox();
    QS_bulkModulus->setValue(0);
    QS_bulkModulus->setDecimals(1);
    QS_bulkModulus->setMinimum(-100000000);
    QS_bulkModulus->setMaximum(100000000);
    QS_bulkModulus->setMaximumWidth(MaxWidth);

    QLabel *QL_shearModulus = new QLabel("shear modulus");
    QS_shearModulus = new QDoubleSpinBox();
    QS_shearModulus->setValue(0);
    QS_shearModulus->setDecimals(1);
    QS_shearModulus->setMaximum(100000000);
    QS_shearModulus->setMaximumWidth(MaxWidth);

    QLabel *QL_kappa = new QLabel("kappa");
    QS_kappa = new QDoubleSpinBox();
    QS_kappa->setValue(0);
    QS_kappa->setDecimals(1);
    QS_kappa->setMinimum(-100000000);
    QS_kappa->setMaximum(100000000);
    QS_kappa->setMaximumWidth(MaxWidth);
    
    // create QPushButton which enforces update of gastrulation values
    QPushButton *QP_updateGastrulationProperties = new QPushButton(tr("Update"),tissueButtons);
    QP_updateGastrulationProperties->resize(10,10);
    QP_updateGastrulationProperties->setMaximumWidth(MaxWidth);

    
    layout->addWidget(QL_Lumen_V0,0,0);layout->addWidget(QS_Lumen_V0,0,1);
    layout->addWidget(QL_Lumen_K,1,0); layout->addWidget(QS_Lumen_K,1,1);
    layout->addWidget(QL_Lumen_currentVolume,2,0); layout->addWidget(QS_Lumen_currentVolume,2,1);
    layout->addWidget(QL_Ellipse_x,3,0);layout->addWidget(QS_Ellipse_x,3,1);
    //layout->addWidget(QL_Ellipse_y,4,0);layout->addWidget(QS_Ellipse_y,4,1);
    layout->addWidget(QL_Ellipse_z,4,0);layout->addWidget(QS_Ellipse_z,4,1);
    layout->addWidget(QL_Ellipse_k,5,0);layout->addWidget(QS_Ellipse_k,5,1);
    
    layout->addWidget(QL_Ellipse_z,4,0);layout->addWidget(QS_Ellipse_z,4,1);
    layout->addWidget(QL_Ellipse_k,5,0);layout->addWidget(QS_Ellipse_k,5,1);
    
    layout->addWidget(QL_K2D,6,0);layout->addWidget(QS_K2D,6,1);
    layout->addWidget(QL_A0,6,2);layout->addWidget(QS_A0,6,3);

    layout->addWidget(QL_shearModulus,7,0);layout->addWidget(QS_shearModulus,7,1);
    layout->addWidget(QL_bulkModulus,7,2);layout->addWidget(QS_bulkModulus,7,3);

    layout->addWidget(QL_kappa,8,0);layout->addWidget(QS_kappa,8,1);
    
    layout->addWidget(QP_updateGastrulationProperties,9,0);
    
    GastrulationBox->setLayout(layout);
    
    // connect the boxes to the relevant variables in the tissue
    connect(QS_Ellipse_x, SIGNAL(valueChanged(double)), integrator->T, SLOT(slot_change_ellipsoid_ex(double)));
    connect(QS_Ellipse_y, SIGNAL(valueChanged(double)), integrator->T, SLOT(slot_change_ellipsoid_ey(double)));
    connect(QS_Ellipse_z, SIGNAL(valueChanged(double)), integrator->T, SLOT(slot_change_ellipsoid_ez(double)));
    connect(QS_Ellipse_k, SIGNAL(valueChanged(double)), integrator->T, SLOT(slot_change_ellipsoid_K(double)));
    connect(QS_Lumen_V0, SIGNAL(valueChanged(double)), integrator->T, SLOT(slot_change_lumenV0_a(double)));
    connect(QS_Lumen_K, SIGNAL(valueChanged(double)), integrator->T, SLOT(slot_change_lumenK_a(double)));
    connect(integrator->T, SIGNAL(signal_currentLumen_a(double)), QS_Lumen_currentVolume, SLOT(setValue(double)));
                                                                   
    connect(QP_updateGastrulationProperties, SIGNAL(clicked()), integrator->T, SLOT(slot_updateGastrulationProperties()));
                                                                   
    connect(integrator, SIGNAL(signal_changeLumen_a(double,double,double)), this, SLOT(slot_changeLumen_a(double,double,double)));
    connect(integrator, SIGNAL(signal_changeLumen_b(double,double,double)), this, SLOT(slot_changeLumen_b(double,double,double)));
    connect(integrator->T, SIGNAL(signal_changeEllipsoid(double,double,double,double)), this, SLOT(slot_changeEllipsoid(double,double,double,double)));
    
    connect(integrator->T, SIGNAL(signal_updateGastrulationProperties(double,double,double,double,double)), this, SLOT(slot_updateGastrulationProperties(double,double,double,double,double)));
    
    QS_Ellipse_x->setValue(0);
    QS_Ellipse_y->setValue(0);
    QS_Ellipse_z->setValue(0);
    QS_Ellipse_k->setValue(0);
    
}

void MainWindow::createCystBox()
{
    CystBox = new QGroupBox();
    double MaxWidth=80;
    QGridLayout *layout = new QGridLayout;
    
    // push button to mutate a round cyst in the middle of the tissue (works only if tissue is periodic)
    QPushButton *QP_mutate_circle =new QPushButton(tr("circle"),tissueButtons);
    QP_mutate_circle->resize(10,10);
    QP_mutate_circle->setMaximumWidth(MaxWidth);
    
    // push button to mutate a stripe cyst in the middle of the tissue (works only if tissue is periodic)
    QPushButton *QP_mutate_stripe =new QPushButton(tr("stripe"),tissueButtons);
    QP_mutate_stripe->resize(10,10);
    QP_mutate_stripe->setMaximumWidth(MaxWidth);
    
    // push button to mutate two stripe cyst in the middle of the tissue (works only if tissue is periodic)
    QPushButton *QP_mutate_two_stripes =new QPushButton(tr("2 stripes"),tissueButtons);
    QP_mutate_two_stripes->resize(10,10);
    QP_mutate_two_stripes->setMaximumWidth(MaxWidth);
    
    // set the size of the round clone/width of the stripe
    QLabel *QL_clone_size = new QLabel("Clone size");
    QL_clone_size->setMaximumWidth(MaxWidth);
    QS_clone_size = new QDoubleSpinBox;
    QS_clone_size->setValue(1.);
	QS_clone_size->setSingleStep(1);
    QS_clone_size->setMinimum(0);
    QS_clone_size->setMaximum(100000);
    QS_clone_size->setMaximumWidth(MaxWidth);

    // add the basement membrane
    QLabel *QL_ECM_apicalPosition = new QLabel("apical position of ECM");
    QS_ECM_apicalPosition = new QDoubleSpinBox();
    QS_ECM_apicalPosition->setValue(0.0);
	QS_ECM_apicalPosition->setSingleStep(1);
    QS_ECM_apicalPosition->setMinimum(-10000);
    QS_ECM_apicalPosition->setMaximum(100000);
    QS_ECM_apicalPosition->setMaximumWidth(MaxWidth);
    QLabel *QL_ECM_basalPosition = new QLabel("basal position of ECM");
    QS_ECM_basalPosition = new QDoubleSpinBox();
    QS_ECM_basalPosition->setValue(0.0);
	QS_ECM_basalPosition->setSingleStep(1);
    QS_ECM_basalPosition->setMinimum(-10000);
    QS_ECM_basalPosition->setMaximum(100000);
    QS_ECM_basalPosition->setMaximumWidth(MaxWidth);
    QLabel *QL_ECM_stiffness = new QLabel("stiffness of ECM");
    QS_ECM_stiffness = new QDoubleSpinBox();
    QS_ECM_stiffness->setValue(0.0);
	QS_ECM_stiffness->setSingleStep(1);
    QS_ECM_stiffness->setMinimum(-10000);
    QS_ECM_stiffness->setMaximum(100000);
    QS_ECM_stiffness->setMaximumWidth(MaxWidth);

    
    
    layout->addWidget(QP_mutate_circle,         1, 0);
    layout->addWidget(QP_mutate_stripe,         2, 0);
    layout->addWidget(QP_mutate_two_stripes,    3, 0);
    layout->addWidget(QL_clone_size,            4, 0);layout->addWidget(QS_clone_size,            4, 1);
    layout->addWidget(QL_ECM_apicalPosition,    5, 0);layout->addWidget(QS_ECM_apicalPosition,    5, 1);
    layout->addWidget(QL_ECM_basalPosition,     6, 0);layout->addWidget(QS_ECM_basalPosition,     6, 1);
    layout->addWidget(QL_ECM_stiffness,         7, 0);layout->addWidget(QS_ECM_stiffness,         7, 1);
 
    connect(QP_mutate_circle,       SIGNAL(clicked()), this,    SLOT(slot_mutate_cells_in_middle()));
    connect(QP_mutate_stripe,       SIGNAL(clicked()), this,    SLOT(slot_mutate_cells_in_stripe()));
    connect(QP_mutate_two_stripes,  SIGNAL(clicked()), this,    SLOT(slot_mutate_cells_in_two_stripes()));
    
    connect(this, SIGNAL(signal_mutate_cells_in_middle(int)),integrator->T,SLOT(mutate_cells_in_circle(int)));
    connect(this, SIGNAL(signal_mutate_cells_in_stripe(double)),integrator->T,SLOT(mutate_cells_in_stripe(double)));
    connect(this, SIGNAL(signal_mutate_cells_in_two_stripes(double)),integrator->T,SLOT(mutate_cells_in_two_stripes(double)));
    
    connect(QS_ECM_apicalPosition, SIGNAL(valueChanged(double)),this, SLOT(change_QS_ECM_apicalPosition(double)));
    connect(QS_ECM_basalPosition, SIGNAL(valueChanged(double)),this, SLOT(change_QS_ECM_basalPosition(double)));
    connect(QS_ECM_stiffness, SIGNAL(valueChanged(double)),this, SLOT(change_QS_ECM_stiffness(double)));
    
    connect(integrator->T, SIGNAL(signal_changeFlatECM(double,double,double)), this, SLOT(slot_changeFlatECM(double,double,double)));
    
    CystBox->setLayout(layout);
}

void MainWindow::createRecordBox()
{
	
	double MaxWidth=60;
	RecordBox = new QGroupBox("Record movie");
    QGridLayout *layout = new QGridLayout;
	
	QLabel *frameRatePrefix= new QLabel("Record 1 frame every ");
	frameRatePrefix->setFont(QFont("Times", 10));
	//frameRatePrefix->setMaximumWidth(MaxWidth);
	QSpinBox * frameRate=new QSpinBox(this);
	frameRate->setMinimum(1);
	frameRate->setMaximum(10000);
	frameRate->setMaximumWidth(MaxWidth);
	QLabel *frameRateSuffix= new QLabel(" iterations");
	frameRateSuffix->setFont(QFont("Times", 10));
	frameRateSuffix->setMaximumWidth(MaxWidth);
    
	recordStatus= new QLabel("Recording");
	recordStatus->hide();
	recordStatus->setMaximumWidth(MaxWidth);

	layout->addWidget(frameRatePrefix,0, 1);
	layout->addWidget(frameRate,0, 2);
	layout->addWidget(frameRateSuffix,0, 3);
	layout->addWidget(recordStatus,2, 1);
	RecordBox->setLayout(layout);
//	connect(startRecord,SIGNAL(clicked()),this,SLOT(startRecord()));
//	connect(stopRecord,SIGNAL(clicked()),this,SLOT(stopRecord()));
	connect(frameRate,SIGNAL(valueChanged(int)),tissueWindow,SLOT(setFrameRate(int)));
	frameRate->setValue(10);
    
    //
}

//void MainWindow::createExtractInfosBox()
//{
//	double MaxWidth=60;
//	ExtractInfosBox = new QGroupBox("Extract informations");
//    QGridLayout *layout = new QGridLayout;
//    
//	QPushButton *guessTensions= new QPushButton("Guess Tensions");
//	guessTensions->setFont(QFont("Times", 10));
//    
//	guessTensions->setMaximumWidth(MaxWidth);
//	//QPushButton *guessTensionsRG= new QPushButton("Guess Tensions+RG");
//	//guessTensionsRG->setMaximumWidth(MaxWidth);
//	QPushButton *induceSelection= new QPushButton("Induce Sel.");
//	induceSelection->setMaximumWidth(MaxWidth);
//	induceSelection->setFont(QFont("Times", 10));
//    
//	QPushButton *FindColors= new QPushButton("Find Colors");
//	FindColors->setMaximumWidth(MaxWidth);
//	FindColors->setFont(QFont("Times", 10));
//    
//    
//	//QPushButton *addNoise= new QPushButton("Add Noise");
//	//addNoise->setMaximumWidth(MaxWidth);
//	//QDoubleSpinBox * AmplitudeNoise=new QDoubleSpinBox(this);
//	//AmplitudeNoise->setMaximumWidth(MaxWidth);
//	QSpinBox * nbIterations=new QSpinBox(this);
//	nbIterations->setMaximumWidth(MaxWidth);
//	nbIterations->setMinimum(1);
//	nbIterations->setMaximum(100000);
//	//guessTensions->setFont(QFont("Times", 13));
//	//guessTensionsRG->setFont(QFont("Times", 13));
//	nbIterations->setFont(QFont("Times", 10));
//	//layout->addWidget(guessTensions,0, 0);
//    //	layout->addWidget(guessTensionsRG,1, 0);
//	layout->addWidget(induceSelection,0, 0);
//	layout->addWidget(FindColors,1, 0);
//	layout->addWidget(guessTensions,0, 2);
//	layout->addWidget(nbIterations,1,1);
//	//layout->addWidget(addNoise,2,0);
//	//layout->addWidget(AmplitudeNoise,2,1);
//	layout->setHorizontalSpacing (0);
//	layout->setVerticalSpacing (0);
//	ExtractInfosBox->setLayout(layout);
//	//connect(addNoise,SIGNAL(clicked()),this,SLOT(addNoise()));
//	//connect(AmplitudeNoise,SIGNAL(valueChanged(double)),this,SLOT(setAmplitudeNoise(double)));
//	nbIterations->setValue(1);
//}



// action: File->Open, uses loadFile
void MainWindow::slot_open()
{
    QString fileName = QFileDialog::getOpenFileName(this);
    CurrentFileName=fileName;
    if (!fileName.isEmpty()) {
        signal_loadFromTxt(fileName);
    }
}

// action: File->Open, uses loadFile
void MainWindow::slot_saveCrossSectionsAsPDFAct()
{
    // enable selection of the directory and open it as dir
    QFileDialog getDirDialog(this);
    getDirDialog.setFileMode(QFileDialog::Directory);
    getDirDialog.setOptions(QFileDialog::ShowDirsOnly);
    QString dirName = getDirDialog.getExistingDirectory();
    QDir dir = QDir(dirName);
    
    
    std::cout << dirName.toStdString() << "\n--------------\n";
        
    // list of all files in the selected directory
    QStringList listOfFiles = dir.entryList(QDir::Files);
    
    for(int i = 0; i<listOfFiles.size(); i++) { // iterate through the files and save the snap shots
        // open the files
        QString fileName = listOfFiles.at(i);

        //std::cout << fileName.remove(".txt").append(".pdf").toStdString() << " in folder " << dirName.toStdString() << std::endl;
        
        //if(fileName.endsWith(".txt")) {
        if(true) {
            int flag = integrator->T->loadFileFromTxt(dirName+"/"+fileName);
            signal_updateWindow();
            if(flag==1) { // if the file could be opened properly
                //std::cout << dirName.append("/").remove(".txt").append(".pdf").toStdString() << std::endl;
                signal_writeImage(Point(integrator->T->SystemSize.x/2, integrator->T->SystemSize.y/2, 0),Point(1,0,0),fileName.remove(".txt"),dirName);
            }
        }
    
    }
    
    
    // iterate through all files and save the
    
    //std::cout << folderName << std::endl;
    //if (!fileName.isEmpty()) {
        //signal_writeImage(Point(T->SystemSize.x/2, T->SystemSize.y/2, 0),PI/2,PI/2,"Tlc_"+QString::number(Tl_c)+"_Gc_"+QString::number(G_cyst),QString::fromStdString(topFolderName));
    //}
}

// action: File->Open, uses loadFile
void MainWindow::slot_saveCrossSectionsAsPDFActFromSubDirectories()
{
    // enable selection of the directory and open it as dir
    QFileDialog getDirDialog(this);
    getDirDialog.setFileMode(QFileDialog::Directory);
    getDirDialog.setOptions(QFileDialog::ShowDirsOnly);
    QString dirName = getDirDialog.getExistingDirectory();
    QDir dir = QDir(dirName);
    dir.setFilter(QDir::Dirs);
    
    QStringList dirList = dir.entryList();

    
    for(QStringList::iterator it = dirList.begin(); it != dirList.end(); ++it)
    {
        QString subdirName = dirName+"/"+(*it);
    
        std::cout << subdirName.toStdString() << "\n--------------\n";
        
        QDir subdir(subdirName);
        subdir.setFilter(QDir::Files);
        
        // list of all files in the selected directory
        QStringList listOfFiles = subdir.entryList();
        
        for(int i = 0; i<listOfFiles.size(); i++) { // iterate through the files and save the snap shots
            // open the files
            QString fileName = listOfFiles.at(i);
            
            //std::cout << fileName.remove(".txt").append(".pdf").toStdString() << " in folder " << dirName.toStdString() << std::endl;
            
            //if(fileName.endsWith(".txt")) {
            if(true) {
                int flag = integrator->T->loadFileFromTxt(dirName+"/"+fileName);
                signal_updateWindow();
                if(flag==1) { // if the file could be opened properly
                    //std::cout << dirName.append("/").remove(".txt").append(".pdf").toStdString() << std::endl;
                    signal_writeImage(Point(integrator->T->SystemSize.x/2, integrator->T->SystemSize.y/2, 0),Point(1,0,0),fileName.remove(".txt"),dirName);
                }
            }
            
        }
    }
    
    
    // iterate through all files and save the
    
    //std::cout << folderName << std::endl;
    //if (!fileName.isEmpty()) {
    //signal_writeImage(Point(T->SystemSize.x/2, T->SystemSize.y/2, 0),PI/2,PI/2,"Tlc_"+QString::number(Tl_c)+"_Gc_"+QString::number(G_cyst),QString::fromStdString(topFolderName));
    //}
}

void MainWindow::slot_extractDataFromSubDir() {
    
    // enable selection of the directory and open it as dir
    QFileDialog getDirDialog(this);
    getDirDialog.setFileMode(QFileDialog::Directory);
    getDirDialog.setOptions(QFileDialog::ShowDirsOnly);
    QString dirName = getDirDialog.getExistingDirectory();
    QDir dir = QDir(dirName);
    dir.setFilter(QDir::Dirs);
    
    QStringList dirList = dir.entryList();
    
    for (QStringList::iterator it = dirList.begin(); it != dirList.end(); ++it) {
        QString subdirName = dirName+"/"+(*it);
        std::cout << dirName.toStdString() << std::endl;
        integrator->T->extractDataFromDir(subdirName.toStdString());
    
    }
}

void MainWindow::slot_extractDataFromSubDir2() {
    
    // enable selection of the directory and open it as dir
    QFileDialog getDirDialog(this);
    getDirDialog.setFileMode(QFileDialog::Directory);
    getDirDialog.setOptions(QFileDialog::ShowDirsOnly);
    QString dirName = getDirDialog.getExistingDirectory();
    QDir dir = QDir(dirName);
    dir.setFilter(QDir::Dirs);
    
    // save list of all subdirs
    QStringList dirList = dir.entryList();
    
    // now iterate through all the subdirs
    for (QStringList::iterator it = dirList.begin(); it != dirList.end(); ++it) {
        QString subdirName = dirName+"/"+(*it);
        QDir subDir = QDir(subdirName);
        
        std::cout << "subdir: " << subdirName.toStdString() << std::endl;

        subDir.setFilter(QDir::Dirs);
        QStringList subdirList = subDir.entryList();
        
        for (QStringList::iterator it2 = subdirList.begin(); it2 != subdirList.end(); ++it2) {
            QString subsubdirName = subdirName+"/"+(*it2);
            std::cout << "subsubdir: " <<  subsubdirName.toStdString() << std::endl;
            integrator->T->extractDataFromDir(subsubdirName.toStdString());
        }
    }
}

void MainWindow::slot_extractDataFromDir() {
    // enable selection of the directory and open it as dir
    QFileDialog getDirDialog(this);
    getDirDialog.setFileMode(QFileDialog::Directory);
    getDirDialog.setOptions(QFileDialog::ShowDirsOnly);
    QString dirName = getDirDialog.getExistingDirectory();
    
    if(!dirName.isEmpty()) integrator->T->extractDataFromDir(dirName.toStdString());
}


// action: File->Save, uses saveFile
void MainWindow::slot_save()
{
    QString fileName = QFileDialog::getSaveFileName(this);
    
    int pointPosition = fileName.lastIndexOf(".");
    
    // if there is no file ending, set it to be .3dvm
    if(fileName.lastIndexOf(".") ==-1) fileName = fileName + ".3dvm";
    
    if (!fileName.isEmpty())
    {
        // transform the QString fileName to a const char *
        QByteArray ba = fileName.toLocal8Bit();
        const char *fileName_char = ba.data();
        signal_writeToTxt(fileName_char);
    }
}

// action: File->Write Data, uses T->writeTissueData()
void MainWindow::writeToTxt()
{
    QString fileName = QFileDialog::getSaveFileName(this);
    
    std::cout << "save2\n";

    if (!fileName.isEmpty())
    {
        // transform the QString fileName to a const char *
        QByteArray ba = fileName.toLocal8Bit();
        
        const char *fileName_char = ba.data();
        signal_writeToTxt(fileName_char);
    }
}


// action: Script->Run script
void MainWindow::runScript()
{
    //myEngine.globalObject().setProperty("mainWindow",mainWindowValue);
    QString fileName = QFileDialog::getOpenFileName(this,"Run a script");
    
    integrator->runScript(fileName);
    
//    //start execution...
//    QScriptEngine myEngine;
//    QScriptValue integratorValue = myEngine.newQObject(integrator);
//    myEngine.globalObject().setProperty("Integrator", integratorValue);
//    QScriptValue mainWindowValue = myEngine.newQObject(this);
//
//    if (!fileName.isEmpty())
//    {
//        QFile scriptFile(fileName);
//        CurrentRecordPath=fileName;//pour que les images soient enregistres par defaut la ou est le script
//        //tissueWindow->CurrentRecordPath=CurrentRecordPath;
//        if (!scriptFile.open(QIODevice::ReadOnly))
//            return;
//        QTextStream stream(&scriptFile);
//        QString contents = stream.readAll();
//        scriptFile.close();
//        QScriptValue result=myEngine.evaluate(contents, fileName);//en fait l'evaluation en seconde place ne sert a rien a part pour le debug
//        if (myEngine.hasUncaughtException()) {
//            int line = myEngine.uncaughtExceptionLineNumber();
//            qDebug() << "uncaught exception at line" << line << ":" << result.toString();
//        }
//        
//    }
//    QScriptValue three=myEngine.evaluate("Tissue.CurrentCellType()");
//    //std::cout<<three.toNumber();
    
}

// action: Script->Run script
void MainWindow::runScriptFromBox()
{
    integrator->runScriptFromBox(QT_scriptTextBox->toPlainText());
}




// used in open file - takes a txt file (fullfilling special requirements) and reads out the tissue setup
void MainWindow::loadFile(const QString &fileName)
{
    signal_loadFromTxt(fileName); // returns a signal if the loading failed
}

// slot that receives a signal if the profile window is closed
void MainWindow::slot_showProfileOff()
{
    QPushShowProfile->setChecked(false);
}

void MainWindow::slot_currentCell(int c_no, double c_V0, double c_V, int c_type, double c_K, double c_P, double c_T_a, double c_T_b, double c_Aa, double c_AR)
{
    signal_c_No(c_no);
    signal_c_V0(c_V0);
    signal_c_V(c_V);
    signal_c_Type(c_type);
    signal_c_K(c_K);
    signal_c_P(c_P);
    signal_c_Ta(c_T_a);
    signal_c_Tb(c_T_b);
    signal_c_Aa(c_Aa);
    signal_c_AR(c_AR);
}

void MainWindow::slot_currentTissue(int t_NoCells, double Lx, double Ly, double FLx, double FLy, double T_ext, bool fixLx, bool fixLy, double Energy)
{
    signal_t_NoCells(t_NoCells);
    signal_setLx(Lx);
    signal_setLy(Ly);
    signal_setFixLx(fixLx? Qt::Checked : Qt::Unchecked);
    signal_setFixLx(fixLy? Qt::Checked : Qt::Unchecked);
    signal_setFLx(FLx);
    signal_setFLy(FLy);
    signal_setText(T_ext);
    signal_setEnergy(Energy);
}

void MainWindow::slot_currentVertex(int v_no, Point &v_Ra, Point &v_Rb, Point &v_Fa, Point &v_Fb, double v_G_l)
{
    signal_v_no(v_no);
    signal_v_G_l(v_G_l);
    signal_v_Ra_x(v_Ra.x);
    signal_v_Ra_y(v_Ra.y);
    signal_v_Ra_z(v_Ra.z);
    signal_v_Rb_x(v_Rb.x);
    signal_v_Rb_y(v_Rb.y);
    signal_v_Rb_z(v_Rb.z);
    
    signal_v_Fa_x(v_Fa.x);
    signal_v_Fa_y(v_Fa.y);
    signal_v_Fa_z(v_Fa.z);
    signal_v_Fb_x(v_Fb.x);
    signal_v_Fb_y(v_Fb.y);
    signal_v_Fb_z(v_Fb.z);
    
};

void MainWindow::slot_currentEdge(int e_no, int e_c1, int e_c2,  double e_T_l, double e_G_a, double e_G_b, double e_length_a, double e_length_b)
{
    signal_e_No(e_no);
    signal_e_Tl(e_T_l);
    signal_e_Ga(e_G_a);
    signal_e_Gb(e_G_b);
    signal_e_length_a(e_length_a);
    signal_e_length_b(e_length_b);
    signal_e_c1(e_c1);
    signal_e_c2(e_c2);
}

void MainWindow::slot_printCellData()
{
    if(!integrator->T->currentCells.empty())
    {
        Cell* c = &(*integrator->T->currentCells.front());
        std::cout << "--- Cell " << c->Number << " ---\n";
        std::cout << "number of vertices: " << c->ListEdges.size() << std::endl;
        std::cout << "T_a = " << c->T_a() << ", T_b = " << c->T_b() << std::endl;
        std::cout << "quadrant of basal midpoint:"; c->q_BC_b.Print();
    }
            
}

void MainWindow::slot_update_cyst_data(double Ta, double Tb, double Tl, double V0, double K, double Ga, double Gb, double factor_Tl, double factor_Ga, double factor_Gb)
{
    signal_K(K);
    signal_V0(V0);
    signal_Ta(Ta);
    signal_Tb(Tb);
    signal_Tl(Tl);
    signal_Ga(Ga);
    signal_Gb(Gb);
    signal_factor_Tl(Tl);
    signal_factor_Ga(factor_Ga);
    signal_factor_Gb(factor_Gb);
    signal_K(K);
    signal_V0(V0);
}

void MainWindow::updateMechanicalData()
{
    QS_Ta1->setValue(integrator->T->MechProp.CellPropVector[1].T_a);
    QS_Ta2->setValue(integrator->T->MechProp.CellPropVector[2].T_a);

    QS_Ta01->setValue(integrator->T->MechProp.CellPropVector[1].T0_a);
    QS_Ta02->setValue(integrator->T->MechProp.CellPropVector[2].T0_a);

    QS_A0a1->setValue(integrator->T->MechProp.CellPropVector[1].A0_a);
    QS_A0a2->setValue(integrator->T->MechProp.CellPropVector[2].A0_a);

    QS_Tb1->setValue(integrator->T->MechProp.CellPropVector[1].T_b);
    QS_Tb2->setValue(integrator->T->MechProp.CellPropVector[2].T_b);

    QS_Tb01->setValue(integrator->T->MechProp.CellPropVector[1].T0_b);
    QS_Tb02->setValue(integrator->T->MechProp.CellPropVector[2].T0_b);
    
    QS_A0b1->setValue(integrator->T->MechProp.CellPropVector[1].A0_b);
    QS_A0b2->setValue(integrator->T->MechProp.CellPropVector[2].A0_b);

    
    QS_Tb1->setValue(integrator->T->MechProp.CellPropVector[1].T_b);
    QS_Tb2->setValue(integrator->T->MechProp.CellPropVector[2].T_b);
    
    QS_K1->setValue(integrator->T->MechProp.CellPropVector[1].K);
    QS_K2->setValue(integrator->T->MechProp.CellPropVector[2].K);
    
    QS_V01->setValue(integrator->T->MechProp.CellPropVector[1].V0);
    QS_V02->setValue(integrator->T->MechProp.CellPropVector[2].V0);

    
    QS_Tl11->setValue(integrator->T->MechProp.Intercell_SurfaceTension[1][1]);
    QS_Tl22->setValue(integrator->T->MechProp.Intercell_SurfaceTension[2][2]);
    QS_Tl12->setValue(integrator->T->MechProp.Intercell_SurfaceTension[2][1]);

    QS_Ga11->setValue(integrator->T->MechProp.Intercell_ApicalLineTension[1][1]);
    QS_Ga22->setValue(integrator->T->MechProp.Intercell_ApicalLineTension[2][2]);
    QS_Ga12->setValue(integrator->T->MechProp.Intercell_ApicalLineTension[2][1]);

    QS_Gb11->setValue(integrator->T->MechProp.Intercell_BasalLineTension[1][1]);
    QS_Gb22->setValue(integrator->T->MechProp.Intercell_BasalLineTension[2][2]);
    QS_Gb12->setValue(integrator->T->MechProp.Intercell_BasalLineTension[2][1]);
    
}

void MainWindow::changeK1(double d)
{
    integrator->T->MechProp.CellPropVector[1].K = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeK2(double d)
{
    integrator->T->MechProp.CellPropVector[2].K = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeV01(double d)
{
    integrator->T->MechProp.CellPropVector[1].V0 = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeV02(double d)
{
    integrator->T->MechProp.CellPropVector[2].V0 = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeTa1(double d)
{
    integrator->T->MechProp.CellPropVector[1].T_a = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeTa2(double d)
{
    integrator->T->MechProp.CellPropVector[2].T_a = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeTa01(double d)
{
    integrator->T->MechProp.CellPropVector[1].T0_a = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeTa02(double d)
{
    integrator->T->MechProp.CellPropVector[2].T0_a = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeAa01(double d)
{
    integrator->T->MechProp.CellPropVector[1].A0_a = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeAa02(double d)
{
    integrator->T->MechProp.CellPropVector[2].A0_a = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeTb1(double d)
{
    integrator->T->MechProp.CellPropVector[1].T_b = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeTb2(double d)
{
    integrator->T->MechProp.CellPropVector[2].T_b = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeTb01(double d)
{
    integrator->T->MechProp.CellPropVector[1].T0_b = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeTb02(double d)
{
    integrator->T->MechProp.CellPropVector[2].T0_b = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeAb01(double d)
{
    integrator->T->MechProp.CellPropVector[1].A0_b = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeAb02(double d)
{
    integrator->T->MechProp.CellPropVector[2].A0_b = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeTl11(double d)
{
    integrator->T->MechProp.Intercell_SurfaceTension[1][1] = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeTl21(double d)
{
    integrator->T->MechProp.Intercell_SurfaceTension[2][1] = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeTl22(double d)
{
    integrator->T->MechProp.Intercell_SurfaceTension[2][2] = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeGa11(double d)
{
    integrator->T->MechProp.Intercell_ApicalLineTension[1][1] = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeGa21(double d)
{
    integrator->T->MechProp.Intercell_ApicalLineTension[2][1] = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeGa22(double d)
{
    integrator->T->MechProp.Intercell_ApicalLineTension[2][2] = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeGb11(double d)
{
    integrator->T->MechProp.Intercell_BasalLineTension[1][1] = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeGb21(double d)
{
    integrator->T->MechProp.Intercell_BasalLineTension[2][1] = d;
    integrator->T->setMechanicalProperties();
}
void MainWindow::changeGb22(double d)
{
    integrator->T->MechProp.Intercell_BasalLineTension[2][2] = d;
    integrator->T->setMechanicalProperties();
}

void MainWindow::change_QS_ECM_apicalPosition(double d)
{
    integrator->T->flatECM_apicalPosition = d;
}

void MainWindow::change_QS_ECM_basalPosition(double d)
{
    integrator->T->flatECM_basalPosition = d;
}

void MainWindow::change_QS_ECM_stiffness(double d)
{
    integrator->T->flatECM_stiffness = d;
}

//void MainWindow::RecordImageFile()
//{
//	//convert number
//	std::string s;
//	std::stringstream out;
//	int number=tissueWindow->nbImagesRecorded;
//    
//	//cette partie doit servir a avoir un nombre de zeros constant
//	if(number<10) out<<0<<0<<number;
//	else{
//        if(number<100) out<<0<<number;	else out<<number;
//	}
//
//	s = out.str();
//    QString fileName = CurrentRecordPath+CurrentRecordPath.fromStdString(s)+".png";
//	tissueWindow->CurrentImage.fill(QColor(245,208,123).rgb());//remet l'image a 0
//	tissueWindow->writeCurrentImage();//ecrit l'image en cours dans QImage
//    tissueWindow->CurrentImage.save(fileName,"PNG",-1);
//}
//
//void MainWindow::startRecord()
//{
//    // QString path = QDir::currentPath() + "/Frame";
//	CurrentRecordPath = QFileDialog::getSaveFileName(this, tr("Save movie frames"));
//	tissueWindow->CurrentRecordPath=CurrentRecordPath;
//    
//	tissueWindow->DoRecord=true;
//	tissueWindow->nbIterationsDone=0;
//	tissueWindow->nbImagesRecorded=0;
//	recordStatus->show();
//	//tissueWindow->CurrentImage->fill(QColor(245,208,123).rgb());//remet l'image a 0
//	//tissueWindow->writeCurrentImage();//ecrit l'image en cours dans la QImage CurrentImage
//    //tissueWindow->CurrentImage->save(fileName,"PNG",-1);
//}
//
//void MainWindow::stopRecord()
//{
//    tissueWindow->DoRecord=false;
//    recordStatus->hide();
//}


void MainWindow::set_integration_properties(int MinimizationMode, double ftol, double tol, double GTOL, int ITMAX, double noise_numberOfNoiseApplications, double noise_stdDev, int checkForT1)
{
    QS_ftol->setValue(ftol);
    QS_tol->setValue(tol);
    QS_GTOL->setValue(GTOL);
    QS_rateT1check->setValue(checkForT1);
    QS_ITMAX->setValue(ITMAX);
    QS_minimizationMode->setValue(MinimizationMode);
    QS_numberOfNoiseApplications->setValue(noise_numberOfNoiseApplications);
    QS_noise_stdDev->setValue(noise_stdDev);
}

#endif // _USE_QT_