#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <list>
#include <QString>
#include <QtScript>
#include <QScriptString>
#include <QDir>

#include "Tissue.h"

struct MinimizationParameters
{
    int MinimizationMode; // 1 - Fletcher Reese, 2 - Polak Ribiere
    double ftol; // convergence criterion - relative error in successive energy steps
    double GTOL; // convergence criterion - max gradient contribution
    int rateT1check; // number of steps after which the system is checked for T1s
    double tol; // tolerance in line minimization
    int ITMAX; // maximum number of iterations before abort
    double EPS; //
    unsigned noise_distribution; // is the noise distributed normally around 0 with std dev noise_stdDev (1) or uniformly in [-noise_stdDev,+noise_stdDev]
    double noise_stdDev; // strength of noise in the course of the minimization process
    unsigned noise_kind; // how is the noise applied? no noise (0), throughout the minimization (1) or noise_numberOfNoiseApplications times after the minimization process (2)
    unsigned noise_numberOfNoiseApplications; // how often is noise applied in the course of a minization process? for noise_kind==1 during the minimization and for noise_kind==2 only at the end with a successive relaxation after each shot noise
    bool minimization_including_SystemSize;
};


#ifdef _USE_QT_
class Integrator : public QObject
{
    Q_OBJECT
#else
class Integrator
{
#endif //_USE_QT_
public:
    // constructor
    Integrator();
    Integrator(int mode, double maxStepSize, double epsilon);
    // destructor
    ~Integrator(){delete T;};
    
    MinimizationParameters minPar;

    // integrate the system N steps
    int evolve_N_steps(unsigned N);
    int evolve_N_steps(unsigned N, double maxStep);
    void evolve();
    
    // current state of the tissue (public such that the other objects can access it easily)
    Tissue *T;
    
    void createHexagonalTissue(int kind, bool isPeriodic, int noRows, int noCols, double height, double width, double length);
    void createHexagonalTissue(int kind, bool isPeriodic, int noRows, int noCols, double height, double width, double length, double K, double T_a, double T_b, double T_l, double G_a, double G_b, double G_l);

    int test(unsigned int iterations);
    
    int test_evolve();
    
    double inline SIGN(double a,double b){return (b>0 ? std::abs(a):-1.0*std::abs(a));}
    inline void shft3(double &a, double &b, double &c, const double d)
    {
		a=b;
		b=c;
		c=d;
    }
    inline void SWAP(double &x, double &y) { double temp=x; x=y; y=temp; }
    inline void TIMES(double d, std::vector<double> &v_in, std::vector<double> &v_res)
    {
        if(v_res.size()!=v_in.size()) throw;
        for(unsigned i=0;i<v_in.size();i++) v_res[i]=d*v_in[i];
    }
    inline void RESCALE(std::vector<double> &dir,double factor_x, double factor_y ){
        if(dir.size() != 6*T->ListVertex.size()+2) throw("in RESCALE: vectors incompatible");
        
        int N_times_3 = T->ListVertex.size()*3;
        
        if(factor_x==1 && factor_y==1) return; // nothing to be done
        
        for(unsigned i=0; i<T->ListVertex.size(); i++){ // adjust the gradient since the system size has changed
            dir[N_times_3+3*i+2] *= factor_x;
            dir[N_times_3+3*i+3] *= factor_y;
            dir[N_times_3+3*i+2] *= factor_x;
            dir[N_times_3+3*i+3] *= factor_y;
        }
    }

//    void analysis_varying_Gab();
//    void analysis_varying_Tl();
//    void analysis_vary_both();
//    void analysis_symmetric();
//    void analysis_general(unsigned rows, unsigned cols, double radius_cyst, double Ta, double Tb, double Ga, double Gb, double Tl_c);
//    void analysis_general(unsigned rows, unsigned cols, double radius_cyst, double Ta, double Tb, double Ga, double Gb, double Tl_c, QString folderName);
//    void analysis_general(double Tl_c, QString tissueName, QString folderName);
//    void analysis_general(QString tissueName, QString fileName_save, double Ta_c, double Tb_c, double Tl_c, double Ga_c, double Gb_c, MinimizationParameters minParameters, bool ECM);
//    void analysis_asymmetric();
private:
    // optimization parameters
    double maxStepSize; // maximal allowed minimization step size for vertices
    double maxStepSize_systemSize; // maximal allowed minimization step size for the system size
    unsigned numberOfSteps;
    int mode; // kind of minimization (steepest descent vs. gradient methods)
    double time;
    double epsilon;
    bool stopNow;
    
    QString currentRecordPath;
    
    QString folderName_script;
    
    // recording routine
    QString record_path;
    QString record_name; // what is the name of the record
    bool record_3dvm, record_top, record_crossSection, record_3D;
    bool record_nextIntegration;
    unsigned record_rate;
    
    // tissue properties
    std::list<double> EnergyOverTime;

    /// how often shall the TissueWindow be updated
    unsigned updateRate;
    
    // the vizualization can be turned off in the process of integration
    bool useGUI;

    int minimize();

    void bracket(std::vector<double> &dir, double &ax, double &bx, double &cx, double &fb, double &current_u);
    double line_minimization_Brent(std::vector<double> &dir);

#ifdef _USE_QT_
public slots:
    void slot_numberOfSteps(int N) {numberOfSteps=N;};
    void slot_maxStepSize(double s) {maxStepSize = s; maxStepSize_systemSize=s/10.;}
    void slot_maxStepSize_systemSize(double s) {maxStepSize_systemSize = s;}
    //void slot_analysis(){analysis_symmetric();};
    int integrate();
    void slot_stop(){stopNow=true;std::cout << "Stop button has been pressed.\n";}
    void loadFromTxt(QString fileName){
        
        int flag = T->loadFileFromTxt(fileName);
        if(flag!=1) std::cout << "Reading Tissue " << fileName.toStdString() << " failed with flag " << flag << "!\n";
        //resetTime();
        updateWindows();
    }
    void writeToTxt(QString fileName) {
        
        std::cout << "Saving tissue to " << fileName.toStdString() << "!\n";
                
        if(fileName.lastIndexOf(".")==-1) fileName = fileName+".3dvm"; // if there is no file type set it to be .3dvm
        
        if (!fileName.isEmpty())
        {
            // transform the QString fileName to a const char *
            QByteArray ba = fileName.toLocal8Bit();
            const char *fileName_char = ba.data();
            T->writeFileToTxt(fileName_char);
        }
    }
    
    
    // functions to be used in scripts and in MainWindow
    void set_updateRate(int updateRate_) {updateRate=updateRate_; set_integration_properties(minPar.MinimizationMode, minPar.ftol, minPar.tol, minPar.GTOL, minPar.ITMAX, minPar.noise_numberOfNoiseApplications, minPar.noise_stdDev, minPar.rateT1check); }
    void set_MinimizationMode(int MinimizationMode) {minPar.MinimizationMode = MinimizationMode;    set_integration_properties(minPar.MinimizationMode, minPar.ftol, minPar.tol, minPar.GTOL, minPar.ITMAX, minPar.noise_numberOfNoiseApplications, minPar.noise_stdDev, minPar.rateT1check);}
    void set_ftol(double ftol) {minPar.ftol = ftol; set_integration_properties(minPar.MinimizationMode, minPar.ftol, minPar.tol, minPar.GTOL, minPar.ITMAX, minPar.noise_numberOfNoiseApplications, minPar.noise_stdDev, minPar.rateT1check);}
    void set_GTOL(double GTOL) {minPar.GTOL = GTOL;set_integration_properties(minPar.MinimizationMode, minPar.ftol, minPar.tol, minPar.GTOL, minPar.ITMAX, minPar.noise_numberOfNoiseApplications, minPar.noise_stdDev, minPar.rateT1check);}
    void set_tol(double tol) {minPar.tol = tol;set_integration_properties(minPar.MinimizationMode, minPar.ftol, minPar.tol, minPar.GTOL, minPar.ITMAX, minPar.noise_numberOfNoiseApplications, minPar.noise_stdDev, minPar.rateT1check);}
    void set_rateT1check(int checkForT1) {minPar.rateT1check = checkForT1;set_integration_properties(minPar.MinimizationMode, minPar.ftol, minPar.tol, minPar.GTOL, minPar.ITMAX, minPar.noise_numberOfNoiseApplications, minPar.noise_stdDev, minPar.rateT1check);}
    void set_ITMAX(int ITMAX) {minPar.ITMAX = ITMAX;set_integration_properties(minPar.MinimizationMode, minPar.ftol, minPar.tol, minPar.GTOL, minPar.ITMAX, minPar.noise_numberOfNoiseApplications, minPar.noise_stdDev, minPar.rateT1check);}
    void set_EPS(double EPS) {minPar.EPS = EPS;set_integration_properties(minPar.MinimizationMode, minPar.ftol, minPar.tol, minPar.GTOL, minPar.ITMAX, minPar.noise_numberOfNoiseApplications, minPar.noise_stdDev, minPar.rateT1check);}
    void set_noise_stdDev(double noise_stdDev) {minPar.noise_stdDev = noise_stdDev;set_integration_properties(minPar.MinimizationMode, minPar.ftol, minPar.tol, minPar.GTOL, minPar.ITMAX, minPar.noise_numberOfNoiseApplications, minPar.noise_stdDev, minPar.rateT1check);}
    void set_noise_numberOfNoiseApplications(int noise_numberOfNoiseApplications) {minPar.noise_numberOfNoiseApplications = noise_numberOfNoiseApplications;set_integration_properties(minPar.MinimizationMode, minPar.ftol, minPar.tol, minPar.GTOL, minPar.ITMAX, minPar.noise_numberOfNoiseApplications, minPar.noise_stdDev, minPar.rateT1check);}
    
    void set_cellProperties(int type, double Ta, double T0_a, double A0_a, double Tb, double T0_b, double A0_b, double K, double V0);
    void set_edgeProperties(int type1, int type2, double Tl, double Ga, double Gb);
    void set_vertexProperties(double Gl);
    void set_periodicityProperties(double T_ext, bool Lx_fix, bool Ly_fix);
    void setECM(bool exists){T->STIFF_BASAL_MEMBRANE = exists; if(exists) T->basal_wall = false; };
    
    void slot_setCurrentRecordPath(QString crp) {currentRecordPath = crp;};

    void set_ECM_Properties(double pos_a, double pos_b, double K, bool basallyBothWays);
    
    /// run a script saved in the file filename
    void runScript(QString fileName);
    
    /// run the script that is given by boxContent
    void runScriptFromBox(QString boxContent);
    
    /// add "Integrator." at the beginning and after every semicolon of the script to make it easier readible
    QString modifyScript(QString s);

    void mutate_cells_in_middle(int number_of_closest_Cells){T->mutate_cells_in_middle(number_of_closest_Cells);}
    void mutate_cells_in_circle(int number_of_closest_Cells){T->mutate_cells_in_circle(number_of_closest_Cells);}
    void mutate_cells_in_annulus(double r1, double r2){T->mutate_cells_in_annulus(r1,r2);}
    void mutate_cells_in_stripe(double stripeWidth){T->mutate_cells_in_stripe(stripeWidth);}

    /// signal to mutate randomly N cells in the tissue
    void randomCellDivisions(unsigned N){T->randomCellDivisions(N);}
    
    /// signal to divide every cell in the tissue once
    void allCellsDivide(){T->allCellsDivide();};
    void divideCells(int cellType, int maxDivisions=0){T->divideCells(cellType,maxDivisions);};
    
    /// signal to mutate all cells in a box
    void mutate_cells_in_box(int cell_type, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max) {T->mutate_cells_in_box(cell_type,x_min,x_max,y_min,y_max,z_min,z_max);}
    
    void recordNextIntegration(QString folderName, QString style, unsigned rate, QString fileName = QString("step"));
    
    /// move all vertices in the tissue, both apical and basal sides, by D
    void moveAllVertices(double dx, double dy, double dz);
    
    void GUI(bool g) {useGUI = g;}
    
    void updateWindows();
    
    void setPopLengths(double minLength, double initLength){T->pop_minLength = minLength; T->pop_initLength = initLength;}

    void change_lumen_a(double K, double V0, double P0){
        T->V0_Lumen_a = V0;
        T->K_Lumen_a = K;
        T->P0_Lumen_a = P0;
        
        // signal for main window
        signal_changeLumen_a(K,V0,P0);
    }
    
    void change_lumen_b(double K, double V0, double P0){
        T->V0_Lumen_b = V0;
        T->K_Lumen_b = K;
        T->P0_Lumen_b = P0;
        
        // signal for main window
        signal_changeLumen_b(K,V0,P0);

    }
    
    /// fix all cells which are in the box [x_min,x_max] x [y_min,y_max] x [z_min,zs_max]
    void fix_cells_in_box(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max);

    /// the type 2 cells with the biggest and the smallest z-value get fixed in their positions
    void fix_outmost_Type2_cells();
    
    /// the type 2 cells which are in the box get fixed in their positions
    void fix_type2_cells_in_box(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max);

    /// unfix all cells which are in the box [x_min,x_max] x [y_min,y_max] x [z_min,zs_max]
    void unfix_cells_in_box(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max);

    /// unfix all cells in the tissue
    void unfix_all_cells();
    
    /// define an additional anisotropic line tension between cells of Type type1 and type2
    void set_anisotropic_lineTension(int type1, int type2, double value)
    {
        T->MechProp.Intercell_ApicalLineTension_anisotropic[std::max(type1,type2)][std::min(type1,type2)] = value;
        T->setMechanicalProperties();
    };
    
    void set_propertiesOutputFile(QString fileName)
    {
        signal_set_propertiesOutputFile(fileName);
    }
    
    /// a function to set the apical tissue properties such that the cells have the right preferred area and cell elasticity
    void set_for_A0_and_K2D(double K2D, double A0);
    void set_for_A0_and_K2D(int cellType, double K2D, double A0);
    
    void saveCrossSectionsAsPDFActFromSubDirectories(QString dirName, QString extension, double px, double py, double pz, double nx, double ny, double nz);
    void saveCrossSectionsForGastrulation(QString dirName);

#endif // _USE_QT_
    void slot_createHexagonalTissue(int,bool, int, int, double, double, double);
    
#ifdef _USE_QT_    
public: signals:
    void signal_writeImage(Point,Point,QString,QString);
    //void signal_update_cyst_data(double,double,double,double,double,double,double,double,double,double);
    //void updateMainWindowData(MechanicalTissueProperties);
    void signal_updateWindow();
    void set_integration_properties(int MinimizationMode, double ftol, double tol, double GTOL, int ITMAX, double noise_numberOfNoiseApplications, double noise_stdDev, int rateT1check);
    void signal_show_status(QString someString);
    void signal_updateMainWindow();
    void signal_update3DWindow();
    void signal_updateCrossSection();
    void signal_currentPointChanged(double x,double y);
    void signal_changeLumen_a(double K,double V0,double P0);
    void signal_changeLumen_b(double K,double V0,double P0);
    void signal_set_propertiesOutputFile(QString fileName);
    void updateMainWindowData();
#else
    //void signal_update_cyst_data(double,double,double,double,double,double,double,double,double,double,double){};
    void signal_currentPointChanged(double x,double y){};
#endif // _USE_QT_
    
    
    
};

#endif