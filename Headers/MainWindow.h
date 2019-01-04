#ifdef _USE_QT_

#ifndef MainWindow_H
#define MainWindow_H

#include <QMainWindow>
#include <QApplication>
#include <QSpinBox>
#include <QBoxLayout>
#include <list>
#include <QFont>
#include <QGridLayout>
#include <QTextEdit>
#include <QPushButton>
#include <QToolBar>
#include <QAction>
#include <QDialog>
#include <QMouseEvent>
#include <QCheckBox>
#include <QTabWidget>
#include "SpecialQDoubleSpinBox.h"
#include "SpecialQSpinBox.h"
#include "Integrator.h"
#include "Tissue3DWindow.h"

QT_BEGIN_NAMESPACE
class QAction;
class QMenu;
class QTextEdit;
class QGroupBox;
class QLabel;
class QLineEdit;
//class MechPropWindow;
class HexagonsWindow;
class TensionWindow;
class ProfileWindow;
class TissueWindow;
class Statistiques;
QT_END_NAMESPACE


class MainWindow : public QMainWindow
{
    Q_OBJECT
        
public slots:
    //void slot_open();
    void slot_showProfileOff();

    void slot_currentTissue(int t_NoCells, double Lx, double Ly, double FLx, double FLy, double T_ext, bool fixLx, bool fixLy, double Energy);
    void slot_currentCell(int c_no, double c_V0, double c_V, int c_type, double c_K, double c_P, double c_T_a, double c_T_b, double c_Aa, double c_AR);
    void slot_currentVertex(int v_no, Point &v_Ra, Point &v_Rb, Point &v_Fa, Point &v_Fb, double G_l);
    void slot_currentEdge(int e_no, int e_c1, int e_c2,  double e_T_l, double e_G_a, double e_G_b, double e_length_a, double e_length_b);
    void slot_open();
    void slot_save();
    void runScript();
    void runScriptFromBox();
    void slot_saveCrossSectionsAsPDFAct();
    void slot_saveCrossSectionsAsPDFActFromSubDirectories();
    void slot_extractDataFromDir();
    void slot_extractDataFromSubDir();
    void slot_extractDataFromSubDir2();
    void slot_printCellData();
    void slot_mutate_cells_in_middle() {signal_mutate_cells_in_middle((int)QS_clone_size->value());}
    void slot_mutate_selected_cells() {signal_mutate_selected_cells();}
    void slot_mutate_cells_in_stripe() {signal_mutate_cells_in_stripe(QS_clone_size->value());}
    void slot_mutate_cells_in_two_stripes() {signal_mutate_cells_in_two_stripes(QS_clone_size->value());}
    void slot_update_cyst_data(double Ta, double Tb, double Tl, double V0, double K, double Ga, double Gb, double factor_Tl, double factor_Ga, double factor_Gb);
    void slot_update(){update(); QApplication::processEvents();}
    void slot_changeFlatECM(double flatECM_apicalPosition, double flatECM_basalPosition, double flatECM_stiffness)
    {
        QS_ECM_apicalPosition->setValue(flatECM_apicalPosition);
        QS_ECM_basalPosition->setValue(flatECM_basalPosition);
        QS_ECM_stiffness->setValue(flatECM_stiffness);
    };
    void slot_changeEllipsoid(double x, double y, double z, double K)
    {
        QS_Ellipse_x->setValue(x);
        QS_Ellipse_y->setValue(y);
        QS_Ellipse_z->setValue(z);
        QS_Ellipse_k->setValue(K);
    };
    
    void slot_updateGastrulationProperties(double preferredArea_2D,double areaElasticity_2D,double bulkModulus_2D,double shearModulus_2D,double kappa_2D)
    {
        QS_K2D->setValue(areaElasticity_2D);
        QS_A0->setValue(preferredArea_2D);
        QS_bulkModulus->setValue(bulkModulus_2D);
        QS_shearModulus->setValue(shearModulus_2D);
        QS_kappa->setValue(kappa_2D);
    }
    
    void slot_changeLumen_a(double V0, double K, double P0)
    {
        QS_Lumen_V0->setValue(V0);
        QS_Lumen_K->setValue(K);
    };
    
    void slot_changeLumen_b(double V0, double K, double P0)
    {
        //QS_Lumen_V0->setValue(V0);
        //QS_Lumen_K->setValue(K);
    };

    void updateMechanicalData();
    void changeK1(double);
    void changeK2(double);
    void changeV01(double);
    void changeV02(double);
    void changeTa1(double);
    void changeTa2(double);
    void changeTb1(double);
    void changeTb2(double);
    void changeTa01(double);
    void changeTa02(double);
    void changeTb01(double);
    void changeTb02(double);
    void changeAa01(double);
    void changeAa02(double);
    void changeAb01(double);
    void changeAb02(double);
    void changeTl11(double);
    void changeTl21(double);
    void changeTl22(double);
    void changeGa11(double);
    void changeGa21(double);
    void changeGa22(double);
    void changeGb11(double);
    void changeGb21(double);
    void changeGb22(double);
    void change_QS_ECM_apicalPosition(double);
    void change_QS_ECM_basalPosition(double);
    void change_QS_ECM_stiffness(double);
    

    void set_integration_properties(int MinimizationMode, double ftol, double tol, double GTOL, int ITMAX, double noise_numberOfNoiseApplications, double noise_stdDev, int checkForT1); // signal coming from Integrator to change the information displayed

    void change_status(QString newText) {statusBar()->showMessage(newText);}

signals:
    void signal_runIntegrator();
    void signal_stopIntegrator();
    void signal_writeToTxt(const char *fileName_char);
    void signal_loadFromTxt(const char *fileName_char);
    void signal_loadFromTxt(QString fileName_qstring);
   
    void signal_v_no(int);
    void signal_v_G_l(double);
    void signal_v_Ra_x(double);
    void signal_v_Ra_y(double);
    void signal_v_Ra_z(double);
    void signal_v_Rb_x(double);
    void signal_v_Rb_y(double);
    void signal_v_Rb_z(double);
    
    void signal_v_Fa_x(double);
    void signal_v_Fa_y(double);
    void signal_v_Fa_z(double);
    void signal_v_Fb_x(double);
    void signal_v_Fb_y(double);
    void signal_v_Fb_z(double);
    
    void signal_c_No(int);
    void signal_c_V0(double);
    void signal_c_A0(double);
    void signal_c_V(double);
    void signal_c_Type(int);
    void signal_c_K(double);
    void signal_c_P(double);
    void signal_c_Ta(double);
    void signal_c_Tb(double);
    void signal_c_Aa(double);
    void signal_c_AR(double);
    
    void signal_e_No(int);
    void signal_e_Tl(double);
    void signal_e_Ga(double);
    void signal_e_Gb(double);
    void signal_e_c1(int);
    void signal_e_c2(int);
    void signal_e_length_a(double);
    void signal_e_length_b(double);

    
    void signal_t_NoCells(int);
    void signal_setLx(double);
    void signal_setLy(double);
    void signal_setFLx(double);
    void signal_setFLy(double);
    void signal_setText(double);
    void signal_setFixLx(Qt::CheckState);
    void signal_setFixLy(Qt::CheckState);
    
    void signal_setEnergy(double);
    
    void signal_V0(double);
    void signal_A0(double);
    void signal_K(double);
    void signal_Ta(double);
    void signal_Tb(double);
    void signal_Tl(double);
    void signal_Ga(double);
    void signal_Gb(double);
    void signal_factor_Tl(double);
    void signal_factor_Ga(double);
    void signal_factor_Gb(double);

    void signal_mutate_cells_in_stripe(double);
    void signal_mutate_cells_in_two_stripes(double);
    void signal_mutate_cells_in_middle(int);
    void signal_mutate_selected_cells();
    void signal_writeImage(Point,Point,QString,QString);
    
    void signal_updateWindow();
    void update_all();
    
public:
    MainWindow();
    ~MainWindow();
    TissueWindow *tissueWindow;
    HexagonsWindow *hexagonsWindow;
    //MechPropWindow *mechPropWindow;
    ProfileWindow *profileWindow;
#ifdef _USE_OPENGL_
    Tissue3DWindow *tissue3DWindow;
#endif
    Integrator* integrator;
    
    QWidget* tissueButtons;
    QPushButton *QP_Integrator;
    QPushButton *QP_Analysis;
    QPushButton *QP_Stop;
    QPushButton *QP_3DView;
    QSpinBox *QS_numSteps;
    QDoubleSpinBox *QS_stepSize;
    
    QSpinBox *StretchEdit;
    QSpinBox *StretchFrequency;
    QDoubleSpinBox *DistanceEq;
    QLabel *DistanceEqLabel;
    //Tissue3DWindow* tissue3DWindow;
    
    QGroupBox * GlobalPropertiesWindow;//acces aux proprietes globales
    QGroupBox * AddNoiseWindow;//acces aux proprietes globales
    
    //QGroupBox *TensionBox;
    //QLabel* labels[5][2];
    //QDoubleSpinBox* TensionEdits[5][5];
    QGroupBox *TissueBox;
    QGroupBox *CellBox;
    QGroupBox *EdgeBox;
    QGroupBox *VertexBox;
    QGroupBox *MechanicsBox;
    QGroupBox *IntegrationBox;
    QGroupBox *GastrulationBox;
    QGroupBox *CystBox;
    QGroupBox *HorizontalTensionBox;
    QGroupBox *RecordBox;
    QGroupBox *ExtractInfosBox;
    QGroupBox *SelectGroupBox;
    QGroupBox *DisplayBox;
    
    // text box to execute cripts on the run
    QTextEdit *QT_scriptTextBox;
    QPushButton *QP_executeScriptBox;
    
    QDoubleSpinBox *HorizontalExtTension;
    QDoubleSpinBox *VerticalExtTension;
    //QCheckBox * TestCross;//boite qui demande si il faut tester si les cotes se croisent
    QCheckBox * DoBestStep;//boite qui demande si il faut utiliser la recherche du meilleur pas pour l'evolution des cotes
    QCheckBox * DoStretch;//boite qui demande si il faut etirer en fonction de la tension
    QCheckBox *DoxStretch;//boite qui demande si il faut etirer en x
    QCheckBox *DoyStretch;//boite qui demande si il faut etirer en y
    QCheckBox *StretchFromStress;//option pour deformer en fonction du stress
    QCheckBox *Noise;//option pour deformer en fonction du stress
    QDoubleSpinBox *AmpNoise;//amplitude du bruit
    QSpinBox *AddNoiseFrequency;//frequence a laquelle le bruit est ajoute
    QCheckBox *TensionAngleDependent;//option pour faire varier la tension en fonction de l'angle des cotes
    QDoubleSpinBox *TensionAngleDependentFactor;
    QDoubleSpinBox *Volume;
    QDoubleSpinBox *QS_clone_size;
    //QCheckBox *ProbabilisticPop;//option pour faire varier la tension en fonction de l'angle des cotes
    
    QPushButton* RemoveShortEdges;
    QDoubleSpinBox * RemoveShortEdgesEdit;
    QPushButton* Mitosis;
    QDoubleSpinBox* MitosisAngleEdit;
    QCheckBox *MitosisRandom;
    QPushButton* Apoptosis;
    QTabWidget *PropertiesTab;
    QTabWidget *ToolsTab;
    
    /// indices of the used boxes
    int IndexCellBox, IndexTissueBox, IndexEdgeBox, IndexVertexBox,IndexMechanicsBox, IndexIntegrationBox, IndexTensionBox,  IndexRecordBox, IndexDisplayBox, IndexExtractInfosBox, IndexSelectGroupBox, IndexGastrulationBox, IndexCystBox;
    
    QString CurrentRecordPath;
    QString CurrentFileName;
    QLabel *recordStatus;//indique si on est en train d'enregistrer ou pas
    int ExtractTensionsNbIterations;
    double ExtractTensionsTemp;
    double AmplitudeNoise;
    double DistanceCloseLI;
    int StripeColor;
    int FillColor1;
    int FillColor2;
    double StripeLIx1;
    double StripeLIx2;
    
    Point CurrentPoint;
    
    QPushButton *QPushShowProfile;
    QPushButton *QPushShow3DTissue;
    
private:
    /// create actions which can be executed from the MainWindow (especially key commands)
    void createActions();
    /// create the tool bars on the top of the applications
    void createToolBars();
    /// create all the popup menus
    void createMenus();
    /// create the box that contains the tissue information
    void createTissueBox();
    /// create the box that allows to manipulate single cells
    void createCellBox();
    /// create the box that allows to manipulate single edges
    void createEdgeBox();
    /// create the box that allows to manipulate single vertices
    void createVertexBox();
    /// create the box that enables to set global mechanical properties of cells and edges depending on their type and to create cysts
    void createMechanicsBox();
    /// create box where integration parameters can be set
    void createIntegrationBox();
    /// create box where gastrulation parameters can be set
    void createGastrulationBox();
    /// create box where cyst parameters can be set
    void createCystBox();
    /// create box, where records folder and speed can be defined
    void createRecordBox();
    
    
    //void RecordImageFile();
    
    QMenu *fileMenu;
    QMenu *WindowsMenu;
    QMenu *ScriptMenu;
    
    
    QToolBar *fileToolBar;
    QGridLayout *gridLayout;
    QAction *evolveAct;//action de faire evoluer le systeme
    QAction *quitAct;//quitter l'application
    QAction *openAct;//ouvrir un fichier
    QAction *createHexagonalTissueAct;
    QAction *openImageAct;//ouvrir une image
    QAction *saveImageAct;//enregistrer une image
    QAction *savePDFImageAct;//enregistrer un PDF de l'image
    QAction *saveCrossSectionsAsPDFAct; // save the cross sections of all images in one folder as PDFs
    /// save the current tissue as a txt - file
    QAction *saveAct;
    /// write the current properties of the tissue into a txt-file
    QAction *extractDataFromDir;
    QAction *showTensionWindowAct;
    QAction *showHexagonsWindowAct;
    QAction *showGlobalPropertiesAct;
    QAction *showAddNoiseAct;
    
    QAction *runScriptAct;
    void loadImageFile(const QString &fileName);//cette fonction charge une image
    void saveImageFile(const QString &fileName);//cette fonction enregistre l'image en cours
    void writeToTxt();
    void loadFile(const QString &fileName);
    
    /// spin boxes used in MechanicsBox
    QDoubleSpinBox *QS_Ta1, *QS_Ta2, *QS_Ta01, *QS_Ta02, *QS_A0a1, *QS_A0a2, *QS_Tb1, *QS_Tb2, *QS_Tb01, *QS_Tb02, *QS_A0b1, *QS_A0b2, *QS_K1, *QS_K2, *QS_V01, *QS_V02, *QS_Ga11, *QS_Ga12, *QS_Ga22, *QS_Gb11, *QS_Gb12, *QS_Gb22, *QS_Tl11, *QS_Tl12, *QS_Tl22;
    
    /// spin boxes for the CystBox
    QDoubleSpinBox *QS_ECM_apicalPosition, *QS_ECM_basalPosition, *QS_ECM_stiffness;
    
    /// spin boxes used in IntegrationBox
    QDoubleSpinBox *QS_ftol, *QS_tol, *QS_GTOL, *QS_noise_stdDev;
    QSpinBox *QS_rateT1check, *QS_ITMAX, *QS_minimizationMode, *QS_numberOfNoiseApplications;
    QCheckBox *QC_minimization_including_SystemSize;
    
    /// spin boxes used in GastrulationBox
    QDoubleSpinBox *QS_Ellipse_x, *QS_Ellipse_y, *QS_Ellipse_z, *QS_Ellipse_k, *QS_Lumen_V0, *QS_Lumen_K, *QS_Lumen_currentVolume;
    QDoubleSpinBox *QS_K2D, *QS_A0, *QS_bulkModulus, *QS_shearModulus, *QS_kappa;
    
        /// spin boxes used in the MechanicsBox
    QLabel * QL_Text; QDoubleSpinBox *QS_Text; // external tension
    QCheckBox *QC_fixLx; QCheckBox *QC_fixLy; // fixation of the system size
};

#include "HexagonsWindow.h"
//#include "MechPropWindow.h"
#include "ProfileWindow.h"
#include "TissueWindow.h"
#endif

#endif // _USE_QT_