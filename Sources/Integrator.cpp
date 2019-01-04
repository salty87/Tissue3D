#include "Integrator.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <iomanip>
#include <limits> // numeric limits
#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include <direct.h> // needed to create folders using mkdi
#define S_IRWXU 0000700 // chmod mask to create full access for owner
#define PI 3.141592653589793238462643383279

Integrator::Integrator()
{
    // create the empty tissue
    T = new Tissue();
    updateRate=10;
    
    minPar = MinimizationParameters();
    minPar.MinimizationMode = 2;
    minPar.ftol = 1e-10;
    minPar.GTOL = 1e-7;
    minPar.rateT1check = 0;
    minPar.tol = 1e-6;
    minPar.ITMAX = 10000;
    minPar.EPS = 1e-18;
    minPar.noise_stdDev = 0;
    minPar.noise_numberOfNoiseApplications = 0; // at how many steps shall noise be applied at the beginning of the simulations
    minPar.minimization_including_SystemSize = true;
    record_nextIntegration = false;
    useGUI = true;
};

Integrator::Integrator(int mode, double maxStepSize, double epsilon):
mode(mode), maxStepSize(maxStepSize), epsilon(epsilon), maxStepSize_systemSize(maxStepSize)
{
    T = new Tissue();
    updateRate=10;
    
    minPar = MinimizationParameters();
    minPar.MinimizationMode = 2;
    minPar.ftol = 1e-12;
    minPar.GTOL = 1e-8;
    minPar.rateT1check = 0;
    minPar.tol = 1e-5;
    minPar.ITMAX = 10000;
    minPar.EPS = 1e-18;
    minPar.noise_stdDev = 0.1;
    minPar.noise_numberOfNoiseApplications = 1; // at how many steps shall noise be applied at the beginning of the simulations
    minPar.minimization_including_SystemSize = true;
    record_nextIntegration = false;
};

int Integrator::integrate()
{
    signal_show_status(QString("Integrating..."));

    int maxT1s = 500;
    
    if(minPar.noise_numberOfNoiseApplications==0)
    {
        for(int T1s_accomplished=0; T1s_accomplished<maxT1s; T1s_accomplished++)
        {
            int flag = minimize();
            if(flag==1 || flag==0) break; // minimization found minimium configuration or stopped after too many steps
            else continue; // topological transition happened and minimization needs to be restarted
        }
    }
    
    for(unsigned its=0; its<minPar.noise_numberOfNoiseApplications; its++)
    {
        T->applyPositionalNoise(minPar.noise_stdDev,1);
    
        for(int T1s_accomplished=0; T1s_accomplished<=maxT1s; T1s_accomplished++)
        {
            bool flag = minimize();
            if(flag==1 || flag==0) break; // minimization found minimium configuration or stopped after too many steps
            else continue; // topological transition happened and minimization needs to be restarted
        }
    }
    
    signal_show_status(QString("Ready!"));
    
    // make sure that recording is only done for one integration process
    record_nextIntegration = false;
}

int Integrator::evolve_N_steps(unsigned N)
{
    
    stopNow = false;
    
    for(unsigned int step=1;step<=N; ++step)
    {
        if(stopNow)
        {
            stopNow=false;
            updateWindows();
            return 1;
        }
        
        // update time
        time += maxStepSize;
        
        //if the tissue has been changed since the last force update, update it again
        if(T->Changed)
        {
            T->update();
        }
        
        if(mode==1)
        {
            // non-periodic BCs
            if(!T->isPeriodic)
            {
                // now that the gradient has been calculated, the system is iterated accordingly
                for(std::list<Vertex>::iterator it_v=T->ListVertex.begin();it_v!=T->ListVertex.end();it_v++)
                {
                    (*it_v).movePoints(maxStepSize);
                }
                
                T->Changed = true;
            }
            // periodic BCs
            else
            {
                maxStepSize_systemSize = maxStepSize;
                
                // new system sizes
                double Lx_new = T->SystemSize.x;
                double Ly_new = T->SystemSize.y;
             
                // append the current energy of the tissue to the vector of energies
                EnergyOverTime.push_back(T->Energy);
                
                double maxMove = 0;
                double move_it;
                // after all the resulting forces have been calculated, the vertices are moved using the steepest descent method
                for(std::list<Vertex>::iterator it=T->ListVertex.begin();it!=T->ListVertex.end();it++)
                {
                    move_it = (*it).movePoints(maxStepSize,Lx_new,Ly_new);
                    if(move_it>maxMove) maxMove = move_it;
                }
                
                maxMove = std::max(std::max(std::abs(T->dW_dL.x),std::abs(T->dW_dL.y)),maxMove);
                
                // set the new size
                T->SystemSize.x = Lx_new;
                T->SystemSize.y = Ly_new;
                
                T->Changed = true;
            }
            
            T->PopAndUnpopAllEdges();
#if 1
            
            // extract parameters
            double a = T->ListEdge.begin()->l_a();
            double h = T->ListVertex.begin()->l_l;
            double G_l = T->ListVertex.begin()->G_l;
            double P = T->ListCell.begin()->Pressure;
            double T_l = T->ListEdge.begin()->T_l;
            double G_a = T->ListEdge.begin()->G_a;
            double G_b = T->ListEdge.begin()->G_b;
            double T_a = T->ListCell.begin()->Ta;
            double T_b = T->ListCell.begin()->Tb;
            double V = T->ListCell.begin()->Volume;
            double V0 = T->ListCell.begin()->V0;
            double K = T->ListCell.begin()->K;
            
            double dW_da = -P*std::sqrt(3.)*3*a*h+3*T_l*h+3*std::sqrt(3.)*(T_a+T_b)+3*(G_a+G_b);
            double dW_dh = -P*std::sqrt(3.)*3/2*a*a+3*T_l*a+3*G_l;
            
            double V_calculated = h*3*std::sqrt(3)/2*a*a;
            double E_calculated = K/2*std::pow(V0-3*std::sqrt(3)/2*h*a*a,2)+3*T_l*a*h+(T_a+T_b)*3*std::sqrt(3)/2*a*a+3*a*(G_a+G_b)+3*h*G_l;
            
            std::cout << "dW_da= "<< dW_da << " ,dW_dh =" << dW_dh << ", V-h*3*sqrt(3)/2*a**2 =" << V-V_calculated << ", Energy difference = " << T->Energy/4 - E_calculated << std::endl;
#endif
            
        }
        
        T->updateMainWindow();
        
        if(step%updateRate==0) signal_updateWindow();
    }
    
    signal_updateWindow();
    
    return 1;
    
}

// evolve until system has converged, that is until max(dR/dt <= maxStep)
int Integrator::evolve_N_steps(unsigned N, double maxStep)
{
    maxStep = 0.001;
    
    stopNow = false;
    
    for(unsigned int step=1;step<=N; ++step)
    {
        if(stopNow)
        {
            stopNow=false;
            signal_updateWindow();
            return 1;
        }
        
        // update time
        time += maxStepSize;
        
        //if the tissue has been changed since the last force update, update it again
        if(T->Changed)
        {
            T->update();
        }
        
        if(mode==1)
        {
            // non-periodic BCs
            if(!T->isPeriodic)
            {
                // now that the gradient has been calculated, the system is iterated accordingly
                for(std::list<Vertex>::iterator it_v=T->ListVertex.begin();it_v!=T->ListVertex.end();it_v++)
                {
                    (*it_v).movePoints(maxStepSize);
                }
                
                T->Changed = true;
            }
            // periodic BCs
            else
            {
                maxStepSize_systemSize = maxStepSize;
                
                // new system sizes
                double Lx_new = T->SystemSize.x;
                double Ly_new = T->SystemSize.y;
                
                //if(!T->Lx_fix) Lx_new += maxStepSize_systemSize * (T->dW_dL.x+T->ext_dW_dL.x*T->SystemSize.y);
                //if(!T->Ly_fix) Ly_new += maxStepSize_systemSize * (T->dW_dL.y+T->ext_dW_dL.y*T->SystemSize.x);

                // append the current energy of the tissue to the vector of energies
                EnergyOverTime.push_back(T->Energy);
                
                double maxMove = 0;
                double move_it;
                // after all the resulting forces have been calculated, the vertices are moved using the steepest descent method
                for(std::list<Vertex>::iterator it=T->ListVertex.begin();it!=T->ListVertex.end();it++)
                {
                    move_it = (*it).movePoints(maxStepSize,Lx_new,Ly_new);
                    if(move_it>maxMove) maxMove = move_it;
                }
                
                // set the new size
                T->SystemSize.x = Lx_new;
                T->SystemSize.y = Ly_new;
                
                T->Changed = true;

                maxMove = std::max(std::max(std::abs(T->dW_dL.x),std::abs(T->dW_dL.y)),maxMove);

                if(maxMove<maxStep) break;
            }
            T->PopAndUnpopAllEdges();
#if 0
            
            // extract parameters
            double a = T->ListEdge.begin()->l_a();
            double h = T->ListVertex.begin()->l_l;
            double G_l = T->ListVertex.begin()->G_l;
            double P = T->ListCell.begin()->Pressure;
            double T_l = T->ListEdge.begin()->T_l;
            double G_a = T->ListEdge.begin()->G_a;
            double G_b = T->ListEdge.begin()->G_b;
            double T_a = T->ListCell.begin()->T_a;
            double T_b = T->ListCell.begin()->T_b;
            double V = T->ListCell.begin()->Volume;
            double V0 = T->ListCell.begin()->V0;
            double K = T->ListCell.begin()->K;
            
            double dW_da = -P*std::sqrt(3.)*3*a*h+3*T_l*h+3*std::sqrt(3.)*(T_a+T_b)+3*(G_a+G_b);
            double dW_dh = -P*std::sqrt(3.)*3/2*a*a+3*T_l*a+3*G_l;
            
            double V_calculated = h*3*std::sqrt(3)/2*a*a;
            double E_calculated = K/2*std::pow(V0-3*std::sqrt(3)/2*h*a*a,2)+3*T_l*a*h+(T_a+T_b)*3*std::sqrt(3)/2*a*a+3*a*(G_a+G_b)+3*h*G_l;
            
            std::cout << "dW_da= "<< dW_da << " ,dW_dh =" << dW_dh << ", V-h*3*sqrt(3)/2*a**2 =" << V-V_calculated << ", Energy difference = " << T->Energy/4 - E_calculated << std::endl;
#endif
            
        }
        
        
        T->updateMainWindow();
        
        if(step%updateRate==0) signal_updateWindow();
    }
    
    signal_updateWindow();
    
}


void Integrator::slot_createHexagonalTissue(int kind, bool isPeriodic, int noRows, int noCols, double height, double width, double length)
{
    if(kind == 0)
    {
        T->clear();
        
        //resetTime();
        
        MechanicalTissueProperties M;
        
        // calculate the length scale of the system
//        double V0 = width*length*height;
//        double l0 = std::pow(V0,1./3.);
//        double l0_pow4 = std::pow(l0,4);
//        double l0_pow5 = std::pow(l0,5);
//        
//        double T_a = 1;
//        double T_b = 1;
//        double T_l = 1;
//        double G_a = 1;
//        double G_b = 1;
//        double G_l = 1;
//        
//        // now rescale all the mechanical parameters
//        double K = 0.0001;
//
//        T_a /= K*l0_pow4;
//        T_b /= K*l0_pow4;
//        T_l /= K*l0_pow4;
//        
//        G_a /= K*l0_pow5;
//        G_b /= K*l0_pow5;
//        G_l /= K*l0_pow5;
        
        
//        double K = 0.001;
//        
//        double T_a = 0.5;
//        double T_b = 2;
//        
//        double T_l = 1;
//        
//        double G_a = 50;
//        double G_b = 0;
//        double G_l = 0;
//        
//        double A0_a = 0;
//        double T0_a = T_a;
//        double A0_b = 0;
//        double T0_b = T_b;
//        
//        M.CellPropVector[1].V0=0;//V0=0 -> is calculated to equal the initial volume of the cells
//        M.CellPropVector[1].K=K;
//        M.CellPropVector[1].T_a = T_a;
//        M.CellPropVector[1].T_b = T_b;
//        M.CellPropVector[1].A0_a = A0_a;
//        M.CellPropVector[1].A0_b = A0_b;
//        M.CellPropVector[1].T0_a = T0_a;
//        M.CellPropVector[1].T0_b = T0_b;
//        
//        M.CellPropVector[2].V0=0;//V0=0 -> is calculated to equal the initial volume of the cells
//        M.CellPropVector[2].K=K;
//        M.CellPropVector[2].T_a = T_a;
//        M.CellPropVector[2].T_b = T_b;
//        M.CellPropVector[2].A0_a = A0_a;
//        M.CellPropVector[2].A0_b = A0_b;
//        M.CellPropVector[2].T0_a = T0_a;
//        M.CellPropVector[2].T0_b = T0_b;
//
//        // define the surface tensions on the lateral sides inbetween the different cell types
//        M.Intercell_SurfaceTension[1][0] = T_l; // surface if border edge
//        M.Intercell_SurfaceTension[1][1] = T_l; // surface to other cell 1
//        M.Intercell_SurfaceTension[2][0] = T_l; // surface if 2 has outline edge
//        M.Intercell_SurfaceTension[2][1] = T_l; // surface between cell 1 and cell 2
//        M.Intercell_SurfaceTension[2][2] = T_l; // surface to another cell 2
//        
//        // define the tensions on the apical lines inbetwen the different cell types
//        M.Intercell_ApicalLineTension[1][0] = G_a; // surface if border edge
//        M.Intercell_ApicalLineTension[1][1] = G_a; // surface to other cell 1
//        M.Intercell_ApicalLineTension[2][0] = G_a; // surface if border edge
//        M.Intercell_ApicalLineTension[2][1] = G_a; // apical line tension between cell 1 and cell 2
//        M.Intercell_ApicalLineTension[2][2] = G_a; // surface to another cell 2
//        
//        // define the tensions on the basal lines inbetwen the different cell types
//        M.Intercell_BasalLineTension[1][0] = G_b; // surface if border edge
//        M.Intercell_BasalLineTension[1][1] = G_b; // surface to other cell 1
//        M.Intercell_BasalLineTension[2][0] = G_b; // surface if border edge
//        M.Intercell_BasalLineTension[2][1] = G_b; // basal line tension between cell 1 and cell 2
//        M.Intercell_BasalLineTension[2][2] = G_b; // surface to another cell 2
//        M.Vertex_LineTension = G_l;
//        
//        // set the tissue properties
//        T->MechProp=M;
//            
//        T->T_ext = 0;
//        
        // ... and create the tissue
        if(isPeriodic)
        {
            T->createPeriodicHexagonalTissue(noRows,noCols,height,length,width);
        }
        else
        {
            T->createHexagonalTissue(noRows,noCols,length,width,height);
        }
        
        
        // now write the V0 into the cell properties of types 1 and 2
        T->MechProp.CellPropVector[1].V0 = T->ListCell.begin()->Volume;
        T->MechProp.CellPropVector[2].V0 = T->ListCell.begin()->Volume;
                
        T->setMechanicalProperties();
        
        T->updateMainWindow();
        
        updateWindows();
    }
    else if(kind == 1)
    {
        // ... do something else
    }
    
};

void Integrator::createHexagonalTissue(int kind, bool isPeriodic, int noRows, int noCols, double height, double width, double length, double K, double T_a, double T_b, double T_l, double G_a, double G_b, double G_l)
{
    if(kind == 0)
    {
        T->clear();
                
        MechanicalTissueProperties M;
        
        M.CellPropVector[1].V0=0;//V0=0 -> is calculated to equal the initial volume of the cells
        M.CellPropVector[1].K=K;
        M.CellPropVector[1].T_a = T_a;
        M.CellPropVector[1].T_b = T_b;
        M.CellPropVector[2].V0=0;//V0=0 -> is calculated to equal the initial volume of the cells
        M.CellPropVector[2].K=K;
        M.CellPropVector[2].T_a = T_a;
        M.CellPropVector[2].T_b = T_b;
        
        // define the surface tensions on the lateral sides inbetween the different cell types
        M.Intercell_SurfaceTension[1][0] = T_l; // surface if border edge
        M.Intercell_SurfaceTension[1][1] = T_l; // surface to other cell 1
        M.Intercell_SurfaceTension[2][0] = T_l; // surface if 2 has outline edge
        M.Intercell_SurfaceTension[2][1] = T_l; // surface between cell 1 and cell 2
        M.Intercell_SurfaceTension[2][2] = T_l; // surface to another cell 2
        
        // define the tensions on the apical lines inbetwen the different cell types
        M.Intercell_ApicalLineTension[1][0] = G_a; // surface if border edge
        M.Intercell_ApicalLineTension[1][1] = G_a; // surface to other cell 1
        M.Intercell_ApicalLineTension[2][0] = G_a; // surface if border edge
        M.Intercell_ApicalLineTension[2][1] = G_a; // apical line tension between cell 1 and cell 2
        M.Intercell_ApicalLineTension[2][2] = G_a; // surface to another cell 2
        
        // define the tensions on the basal lines inbetwen the different cell types
        M.Intercell_BasalLineTension[1][0] = G_b; // surface if border edge
        M.Intercell_BasalLineTension[1][1] = G_b; // surface to other cell 1
        M.Intercell_BasalLineTension[2][0] = G_b; // surface if border edge
        M.Intercell_BasalLineTension[2][1] = G_b; // basal line tension between cell 1 and cell 2
        M.Intercell_BasalLineTension[2][2] = G_b; // surface to another cell 2
        M.Vertex_LineTension = G_l;
        
        // set the tissue properties
        T->MechProp=M;
        T->T_ext = 0;
                
        // ... and create the tissue
        if(isPeriodic)
        {
            T->createPeriodicHexagonalTissue(noRows,noCols,height,length,width);
        }
        else
        {
            T->createHexagonalTissue(noRows,noCols,length,width,height);
        }
        
        // now write the V0 into the cell properties of types 1 and 2
        T->MechProp.CellPropVector[1].V0 = T->ListCell.begin()->Volume;
        T->MechProp.CellPropVector[2].V0 = T->ListCell.begin()->Volume;
        
        T->updateMainWindow();
        
        updateWindows();
    }
    else if(kind == 1)
    {
        // ... do something else
    }
    
};

// function to bracket the minimum (minimize Energy(T->pos + u*dir) for u)
void Integrator::bracket(std::vector<double> &dir, double &ax, double &bx, double &cx, double &fb, double &current_u)
{
    const double GOLD = 1.618034; // default ratio by which successive intervals are magnified
    const double GLIMIT = 100.0; //maximum magnification allowed for a parabolic-fit step
    const double TINY = 1.0e-20;
    
    double fa,fc,fu,u;
    
    // evaluate energy at T->pos+ax*dir and at T->pos+bx*dir
    for(unsigned long i=0;i<T->pos.size();i++) T->pos[i] += (ax-current_u)*dir[i];
    //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
    T->set_pos(minPar.minimization_including_SystemSize);
    T->W();
    fa=T->Energy;
    current_u=ax;
    
    for(unsigned long i=0;i<T->pos.size();i++) T->pos[i] += (bx-current_u)*dir[i];
    //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
    T->set_pos(minPar.minimization_including_SystemSize);
    
    fb=T->W();;
    current_u=bx;
    
    if(fb>fa) { // switch roles of a and b so that we can go downhill in the direction from a to b
        SWAP(ax,bx);
        SWAP(fa,fb);
    }
    
    cx = bx + GOLD*(bx-ax); // first guess for c - now bx is golden section between ax and cx
    // evaluate energy at T->pos+cx*dir
    for(unsigned long i=0;i<T->pos.size();i++) T->pos[i] += (cx-current_u)*dir[i];
    //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
    T->set_pos(minPar.minimization_including_SystemSize);
    fc=T->W();;
    current_u=cx;
    
    while(fb>fc) { // keep returning here until we bracket
        double r=(bx-ax)*(fb-fc);
        double q=(bx-cx)*(fb-fa);
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0*SIGN(std::max(std::abs(q-r),TINY),q-r));
        double ulim=bx+GLIMIT*(cx-bx);
        
        // now test various possibilities
        if((bx-u)*(u-cx)>0.0) { // parabolic u is between b and c: try it
            // fu=f(u):
            for(unsigned long i=0;i<T->pos.size();i++) T->pos[i] += (u-current_u)*dir[i];
            //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
            T->set_pos(minPar.minimization_including_SystemSize);
            fu=T->W();;
            current_u=u;
            if(fu<fc){ // got a minumum between b and c
                ax=bx;
                bx=u;
                fa=fb;
                fb=fu;
                return;
            } else if(fu>fb){ // got a minimum between a and u
                cx=u;
                fc=fu;
                return;
            }
            
            u=cx+GOLD*(cx-bx); // parabolic fit was no use; use default magnification
            for(unsigned long i=0;i<T->pos.size();i++) T->pos[i] += (u-current_u)*dir[i];
            //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
            T->set_pos(minPar.minimization_including_SystemSize);
            fu=T->W();
            current_u=u;
        } else if((cx-u)*(u-ulim)>0.0) { // parabolic fit is between c and its allowed limit
            // fu=f(u):
            for(unsigned long i=0;i<T->pos.size();i++) T->pos[i] += (u-current_u)*dir[i];
            //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
            T->set_pos(minPar.minimization_including_SystemSize);
            fu=T->W();
            current_u=u;
            if(fu<fc) {
                shft3(bx,cx,u,u+GOLD*(u-cx));
                fb=fc;
                fc=fu;
                for(unsigned long i=0;i<T->pos.size();i++) T->pos[i] += (u-current_u)*dir[i];
                //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
                T->set_pos(minPar.minimization_including_SystemSize);
                fu=T->W();
                current_u=u;
            }
            
        } else if((u-ulim)*(ulim-cx) >= 0.0) { // limit parabolic u to maximum allowed value
            u=ulim;
            for(unsigned long i=0;i<T->pos.size();i++) T->pos[i] += (u-current_u)*dir[i];
            //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
            T->set_pos(minPar.minimization_including_SystemSize);
            fu=T->W();
            current_u=u;
            
        } else { // reject parabolic u, use default magnification
            u=cx+GOLD*(cx-bx);
            for(unsigned long i=0;i<T->pos.size();i++) T->pos[i] += (u-current_u)*dir[i];
            //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
            T->set_pos(minPar.minimization_including_SystemSize);
            fu=T->W();
            current_u=u;
        }
        
        shft3(ax,bx,cx,u); // eliminate oldest point and continue
        shft3(fa,fb,fc,fu);
        
    }
    
}


// minimize the energy along the line p+lambda*dir, the result is written into T->pos
double Integrator::line_minimization_Brent(std::vector<double> &dir)
{
    const int ITMAX = minPar.ITMAX; // max number of iterations before the algorithm stops
    const double CGOLD = 0.3819660; // golden ratio
    const double ZEPS = std::numeric_limits<double>::epsilon()*10e-3; // machine limit
    double a,b,d(0),etemp,fu,fv,fw,fx;
    double p,q,r,tol1,tol2,u,v,w,x,xm;
    double e=0; // this will be the distance on the step moved before last
    double tol = minPar.tol;
    T->update_pos(minPar.minimization_including_SystemSize);
    std::vector<double> pos_init = T->pos;
    
    std::vector<double> dir_vec = dir; // take copy of the direction vector to be able to modify it
    
    // call Brents bracket method to enclose the minumum along the gradient
    double ax(-0.0000001),bx(0.0000001), cx;
    bracket(dir,ax,bx,cx,fx,u);
    
    a = (ax<cx ? ax:cx); // a and b must be in ascending order
    b = (ax>cx ? ax:cx);
    x=w=v=bx; // initializations...
    
    double u_previous=u;
    
    fw = fv = fx; // energy at middle point of interval (bx)
    
    for(int iter = 0; iter<ITMAX;iter++){ // main program loop
        
        xm=0.5*(a+b);
        tol2=2.0*(tol1=tol*std::abs(x)+ZEPS);
        if(std::abs(x-xm)<=(tol2-0.5*(b-a))){ // test for done here
            //std::cout << "line minimization finished" << std::endl;
            return 1;
        }
        if(std::abs(e)>tol1){ // construct a trial parabolic fit
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
            p=(x-v)*q - (x-w)*r;
            q=2.0*(q-r);
            if(q>0) p=-p;
            q=std::abs(q);
            etemp=e;
            e=d;
            if(std::abs(p)>=std::abs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x)) d=CGOLD*(e=(x>=xm ? a-x:b-x));
            // the above conditions determine the acceptability of the parabolic fit. Here we take the golden section step into the larger of the two segments
            else{
                d=p/q;
                u=x+d;
                if(u-a<tol2 || b-u<tol2) d = SIGN(tol1,xm-x);
            }
        } else {
            d=CGOLD*(e=(x>=xm ? a-x:b-x));
        }
        
        u= (std::abs(d)>=tol1 ? x+d:x+SIGN(tol1,d));
        
        // now evaluate W for the only time in the iteration
        for(unsigned long i=0;i<T->pos.size();i++) { // move the tissue along dir
            T->pos[i] += (u-u_previous)*dir[i];
        }
        
        u_previous = u;
        
        //RESCALE(dir,T->pos[T->pos.size()-2]/T->SystemSize.x, T->pos[T->pos.size()-1]/T->SystemSize.y );
        
        T->set_pos(T->isPeriodic && minPar.minimization_including_SystemSize);
        fu = T->W();
        
        //std::cout << "a= " << a << ", u= " << u  << ", b=" << b << ", fu= " << fu << std::endl;
        
        if(fu<=fx){
            if(u>=x) a=x; else b=x;
            shft3(v,w,x,u);
            shft3(fv,fw,fx,fu);
        } else {
            if(u<x) a=u; else b=u;
            if(fu<=fw || w==x) {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            } else if(fu<=fv || v==x || v==w){
                v=u;
                fv=fu;
            }
        }
    }
    std::cout << "Too many iterations in brent" << std::endl;
    throw("Too many iterations in brent");
    
}

////
//// function to bracket the minimum (minimize Energy(T->pos + u*dir) for u)
//void Integrator::bracket(std::vector<double> &dir, double &ax, double &bx, double &cx, double &fb, double &current_u)
//{
//    const double GOLD = 1.618034; // default ratio by which successive intervals are magnified
//    const double GLIMIT = 100.0; //maximum magnification allowed for a parabolic-fit step
//    const double TINY = 1.0e-20;
//    
//    double fa,fc,fu,u;
//    
//    // evaluate energy at T->pos+ax*dir and at T->pos+bx*dir
//    for(int i=0;i<T->pos.size();i++) T->pos[i] += (ax-current_u)*dir[i];
//    //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
//    T->set_pos(minPar.minimization_including_SystemSize);
//    T->update_energy();
//    fa=T->Energy;
//    current_u=ax;
//        
//    for(int i=0;i<T->pos.size();i++) T->pos[i] += (bx-current_u)*dir[i];
//    //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
//    T->set_pos(minPar.minimization_including_SystemSize);
//    T->update_energy();
//    fb=T->Energy;
//    current_u=bx;
//    
//    if(fb>fa) { // switch roles of a and b so that we can go downhill in the direction from a to b
//        SWAP(ax,bx);
//        SWAP(fa,fb);
//    }
//    
//    cx = bx + GOLD*(bx-ax); // first guess for c - now bx is golden section between ax and cx
//    // evaluate energy at T->pos+cx*dir
//    for(int i=0;i<T->pos.size();i++) T->pos[i] += (cx-current_u)*dir[i];
//    //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
//    T->set_pos(minPar.minimization_including_SystemSize);
//    T->update_energy();
//    fc=T->Energy;
//    current_u=cx;
//    
//    while(fb>fc) { // keep returning here until we bracket
//        double r=(bx-ax)*(fb-fc);
//        double q=(bx-cx)*(fb-fa);
//        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0*SIGN(std::max(std::abs(q-r),TINY),q-r));
//        double ulim=bx+GLIMIT*(cx-bx);
//        
//        // now test various possibilities
//        if((bx-u)*(u-cx)>0.0) { // parabolic u is between b and c: try it
//            // fu=f(u):
//            for(int i=0;i<T->pos.size();i++) T->pos[i] += (u-current_u)*dir[i];
//            //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
//            T->set_pos(minPar.minimization_including_SystemSize);
//            T->update_energy();
//            fu=T->Energy;
//            current_u=u;            
//            if(fu<fc){ // got a minumum between b and c
//                ax=bx;
//                bx=u;
//                fa=fb;
//                fb=fu;
//                return;
//            } else if(fu>fb){ // got a minimum between a and u
//                cx=u;
//                fc=fu;
//                return;
//            }
//            
//            u=cx+GOLD*(cx-bx); // parabolic fit was no use; use default magnification
//            for(int i=0;i<T->pos.size();i++) T->pos[i] += (u-current_u)*dir[i];
//            //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
//            T->set_pos(minPar.minimization_including_SystemSize);
//            T->update_energy();
//            fu=T->Energy;
//            current_u=u;
//        } else if((cx-u)*(u-ulim)>0.0) { // parabolic fit is between c and its allowed limit
//            // fu=f(u):
//            for(int i=0;i<T->pos.size();i++) T->pos[i] += (u-current_u)*dir[i];
//            //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
//            T->set_pos(minPar.minimization_including_SystemSize);
//            T->update_energy();
//            fu=T->Energy;
//            current_u=u;
//            if(fu<fc) {
//                shft3(bx,cx,u,u+GOLD*(u-cx));
//                fb=fc;
//                fc=fu;
//                for(int i=0;i<T->pos.size();i++) T->pos[i] += (u-current_u)*dir[i];
//                //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
//                T->set_pos(minPar.minimization_including_SystemSize);
//                T->update_energy();
//                fu=T->Energy;
//                current_u=u;
//            }
//            
//        } else if((u-ulim)*(ulim-cx) >= 0.0) { // limit parabolic u to maximum allowed value
//            u=ulim;
//            for(int i=0;i<T->pos.size();i++) T->pos[i] += (u-current_u)*dir[i];
//            //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
//            T->set_pos(minPar.minimization_including_SystemSize);
//            T->update_energy();
//            fu=T->Energy;
//            current_u=u;
//            
//        } else { // reject parabolic u, use default magnification
//            u=cx+GOLD*(cx-bx);
//            for(int i=0;i<T->pos.size();i++) T->pos[i] += (u-current_u)*dir[i];
//            //RESCALE(dir,T->pos[0]/T->SystemSize.x, T->pos[1]/T->SystemSize.y );
//            T->set_pos(minPar.minimization_including_SystemSize);
//            T->update_energy();
//            fu=T->Energy;
//            current_u=u;
//        }
//        
//        shft3(ax,bx,cx,u); // eliminate oldest point and continue
//        shft3(fa,fb,fc,fu);
//        
//    }
//    
//}
//
//
//// minimize the energy along the line p+lambda*dir, the result is written into T->pos
//double Integrator::line_minimization_Brent(std::vector<double> &dir)
//{
//    const int ITMAX = minPar.ITMAX; // max number of iterations before the algorithm stops
//    const double CGOLD = 0.3819660; // golden ratio
//    const double ZEPS = std::numeric_limits<double>::epsilon()*10e-3; // machine limit
//    double a,b,d(0),etemp,fu,fv,fw,fx;
//    double p,q,r,tol1,tol2,u,v,w,x,xm;
//    double e=0; // this will be the distance on the step moved before last
//    double tol = minPar.tol;
//    T->update_pos(minPar.minimization_including_SystemSize);
//    std::vector<double> pos_init = T->pos;
//    
//    std::vector<double> dir_vec = dir; // take copy of the direction vector to be able to modify it
//    
//    // call Brents bracket method to enclose the minumum along the gradient
//    double ax(-0.0000001),bx(0.0000001), cx;
//    bracket(dir,ax,bx,cx,fx,u);
//    
//    a = (ax<cx ? ax:cx); // a and b must be in ascending order
//    b = (ax>cx ? ax:cx);
//    x=w=v=bx; // initializations...
//
//    double u_previous=u;
//    
//    fw = fv = fx; // energy at middle point of interval (bx)
//
//    for(int iter = 0; iter<ITMAX;iter++){ // main program loop
//        
//        xm=0.5*(a+b);
//        tol2=2.0*(tol1=tol*std::abs(x)+ZEPS);
//        if(std::abs(x-xm)<=(tol2-0.5*(b-a))){ // test for done here
//            //std::cout << "line minimization finished" << std::endl;
//            return 1;
//        }
//        if(std::abs(e)>tol1){ // construct a trial parabolic fit
//            r=(x-w)*(fx-fv);
//            q=(x-v)*(fx-fw);
//            p=(x-v)*q - (x-w)*r;
//            q=2.0*(q-r);
//            if(q>0) p=-p;
//            q=std::abs(q);
//            etemp=e;
//            e=d;
//            if(std::abs(p)>=std::abs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x)) d=CGOLD*(e=(x>=xm ? a-x:b-x));
//            // the above conditions determine the acceptability of the parabolic fit. Here we take the golden section step into the larger of the two segments
//            else{
//                d=p/q;
//                u=x+d;
//                if(u-a<tol2 || b-u<tol2) d = SIGN(tol1,xm-x);
//            }
//        } else {
//            d=CGOLD*(e=(x>=xm ? a-x:b-x));
//        }
//
//        u= (std::abs(d)>=tol1 ? x+d:x+SIGN(tol1,d));
//
//        // now evaluate W for the only time in the iteration
//        for(int i=0;i<T->pos.size();i++) { // move the tissue along dir
//            T->pos[i] += (u-u_previous)*dir[i];
//        }
//                        
//        u_previous = u;
//        
//        //RESCALE(dir,T->pos[T->pos.size()-2]/T->SystemSize.x, T->pos[T->pos.size()-1]/T->SystemSize.y );
//        
//        T->set_pos(minPar.minimization_including_SystemSize);
//        T->update_energy();
//        fu = T->Energy;
//        
//        //std::cout << "a= " << a << ", u= " << u  << ", b=" << b << ", fu= " << fu << std::endl;
//        
//        if(fu<=fx){
//            if(u>=x) a=x; else b=x;
//            shft3(v,w,x,u);
//            shft3(fv,fw,fx,fu);
//        } else {
//            if(u<x) a=u; else b=u;
//            if(fu<=fw || w==x) {
//                v=w;
//                w=u;
//                fv=fw;
//                fw=fu;
//            } else if(fu<=fv || v==x || v==w){
//                v=u;
//                fv=fu;
//            }
//        }
//    }
//    std::cout << "Too many iterations in brent" << std::endl;
//    throw("Too many iterations in brent");
//    
//}

// multidimensional minimization by the Fletcher-Reese-Polak-Ribiere method; algorithm following "Numerical Recipes"
int Integrator::minimize()
{
    unsigned step = 0;
    
    const int ITMAX = minPar.ITMAX;
    const double EPS = minPar.EPS;
    const double GTOL = minPar.GTOL; // maximal gradient size as stop criterion
    const double ftol = minPar.ftol;//1.0e-8; // relative max change in energy as stop criterion
    
    double gg, dgg;
    int iter;
    double fret;

    Point recordCs_normalPlane = Point(1,0,0);
    
    int rateT1check= minPar.rateT1check; // number of steps after each which the system is checked for possible T1s
    
    T->update();
    T->update_pos_and_der(T->isPeriodic && minPar.minimization_including_SystemSize);
    int n = T->pos.size();
    std::vector<double> g(n), h(n);
    
    double fp = T->Energy; // energy of starting point
    std::vector<double> xi = T->der; // derivative at starting point of minimization
    std::vector<double> pp = T->pos; // position of the starting point
    
    
    for(int j=0; j<n; j++){ // initial values for gradient calculation
        g[j]=-xi[j];
        xi[j]=h[j]=g[j];
    }
    
    stopNow = false;
    
    // create a directory if it doesn't exist
#ifdef _USE_QT_
    //QString folderName = "/Users/silvanus/MPI/EpithelialMechanics/Tissue3D/Movies/recording";
    QByteArray ba = record_path.toLocal8Bit();
    const char *folderName_char = ba.data();
    mkdir(folderName_char, S_IRWXU);
#endif // _USE_QT_
    
    
    for(int its=0;its<ITMAX;its++){// loop over iterations
        
        if(stopNow) // signal to stop the calculations
        {
            stopNow=false;
            updateWindows();
            return 1;
        }
        
        iter=its; // current iteration
        line_minimization_Brent(xi); // do a line minimization along vector xi and the Tissue will be in this optimal position afterwards
        T->update(); // update the derivatives
        T->update_pos_and_der(minPar.minimization_including_SystemSize); // save the derivatives in der
                
        fret = T->Energy;

        if(2.0*std::abs(fret-fp) < ftol){ // one possible return - check absolute change
        //if(2.0*std::abs(fret-f p) < ftol*(std::abs(fret)+std::abs(fp)+EPS)){ // one possible return - check for relative change
            std::cout << "Minimization successful - difference between successive energies sufficiently low, error = " << 2.0*std::abs(fret-fp) << std::endl;
            
            // write the final result
            if(record_nextIntegration && record_3dvm) {
                QString recName = record_name+QString::number(its)+"_final";
                writeToTxt(record_path+"/"+recName);
            }
            
#ifdef _USE_QT_
            if(record_nextIntegration && record_crossSection) {
                QString recName = record_name+QString::number(its)+"_final";
                signal_writeImage(Point(0, 0, 0),recordCs_normalPlane,recName,record_path);
            }
            updateWindows();
#endif // _USE_QT_
            
            return 1;
        }
        
        fp=fret; // save previous energy fp

        xi = T->der;
        double test = 0.0; // test for convergence on zero gradient
        double den = std::max(fp, 1.0);
        for(int j=0; j<n; j++){
            double temp = std::abs(xi[j])*std::max(std::abs(T->pos[j]),1.0)/den; // Standard way
          //  double temp = std::abs(xi[j])*(n/6.)/den; // here: seems more useful to set something like: change/(energy per vertex) < GTOL   as a criterion since the tissue behaves translational invariant
            if(temp>test) test=temp;
        }
            
        if(test<GTOL){ // the other possible return
            std::cout << "Minimization successful - gradient short enough" << std::endl;

            // write the final result as txt file
            if(record_nextIntegration && record_3dvm) {
                QString recName = record_name+QString::number(its)+"_final";
                writeToTxt(record_path+"/"+recName);
            }

#ifdef _USE_QT_
            if(record_nextIntegration && record_crossSection) {
                QString recName = record_name+QString::number(its)+"_final";
                signal_writeImage(Point(0, 0, 0),recordCs_normalPlane,recName,record_path);
            }
            updateWindows();
#endif //_USE_QT_
            
            return 1;
        }
        
        dgg=gg=0;
        for(int j=0; j<n; j++){ // calculate the conjugate gradient according to Fletcher-Reese or Polak-Ribiere
            gg += g[j]*g[j];
            if(minPar.MinimizationMode == 1) {
                dgg += xi[j]*xi[j]; // statement for Fletcher-Reese
            } else if(minPar.MinimizationMode == 2) {
                dgg += (xi[j]+g[j])*xi[j]; // statement for Polak-Ribiere
            }
        }
    
        if(gg==0){ // unlikely: if gradient is exactly zero, then we are already done
        std::cout << "Minimization successful - gradient zero" << std::endl;

        // write the final result as txt file
        if(record_nextIntegration && record_3dvm) {
            QString recName = record_name+QString::number(its)+"_final";
            writeToTxt(record_path+"/"+recName);
        }

#ifdef _USE_QT_
            updateWindows();
            if(record_nextIntegration && record_crossSection) {
                QString recName = record_name+QString::number(its)+"_final";
                signal_writeImage(Point(0, 0, 0),recordCs_normalPlane,recName,record_path);
            }
#endif //_USE_QT_
            
            return 1;
        }
        
        double gam = dgg/gg;
        
        for(int j=0; j<n; j++) {
            g[j]=-xi[j];
            xi[j]=h[j]=g[j]+gam*h[j];
        }
        
        step++;
        
        // write the current state at a rate record_rate
        if(record_nextIntegration && record_3dvm && its%record_rate==0) {
            QString recName = QString("step%1").arg(QString::number(its),5,QLatin1Char('0'));
            writeToTxt(record_path+"/"+recName);
        }

        
#ifdef _USE_QT_
        if(record_nextIntegration && record_crossSection && its%record_rate==0) { // save another image
            QString recName = QString(record_name+"%1").arg(QString::number(its),5,QLatin1Char('0'));
            signal_writeImage(Point(T->SystemSize.x/2, T->SystemSize.y/2, 0),recordCs_normalPlane,recName,record_path);
        }
        
        updateWindows();
#endif // _USE_QT_
        
        if(rateT1check>0 && its>0 && its%rateT1check==0 && T->PopAndUnpopAllEdges()){
            return 2;
        }
    }
    
    T->update();
    
    std::cout << "Minimization did not succeed in predefined maximal number of iterations!" << std::endl;
    
    return 0; // too many iterations in frpr - method
    
}


void Integrator::evolve()
{
    int number_divisions = 0; // number of random cell divisions in the tissue
    
    //minimize();
    
    for(int i=0;i<number_divisions;i++)  {
        (T->get_random_cell())->RandomDivision();
        minimize(); // relax the system
    }
    
    int flag;
    
    int max_T1s = 500;
    //for(int i=0; i<number_divisions; i++){ // main loop (divide, relax, divide, ...)
        
        int iter=0;
        do { // relax the system
            
            flag = minimize();
            //T->PopAndUnpopAllEdges(0.5,2);
            if(flag==2) std::cout << "T1 took place!" << std::endl;
            if(flag==0) {std::cout << "Minimization did not succeed in predefined maximal number of iterations!" << std::endl; return;}
            iter++;
        } while (iter<max_T1s && flag==2);
    
    if(iter==max_T1s ){std::cout << "More than " << max_T1s << " T1 transitions" << std::endl;}
}

// only for gastrulation: set apical tissue properties of all cell types to account for elastic behavior of apical surfaces: W = K2D/2*(A-A0)^2
void Integrator::set_for_A0_and_K2D(double K2D, double A0)
{
    // calculate the right values for Ta and T0...
    double A0_star = 9999999;
    double Ta = K2D*(A0_star-A0);
    double T0 = -K2D*A0;
    
    for(int i=1; i<=4; ++i)
    {
        // and change the tissue propertieso
        T->MechProp.CellPropVector[i].T_a = Ta;
        T->MechProp.CellPropVector[i].T0_a = T0;
        T->MechProp.CellPropVector[i].A0_a = A0_star;
    }
    
    updateMainWindowData();

}

// only for gastrulation: set apical tissue properties for one cell type to account for elastic behavior of apical surfaces: W = K2D/2*(A-A0)^2
void Integrator::set_for_A0_and_K2D(int cellType, double K2D, double A0)
{
    // only allowed if allowed cell type
    if(cellType<=0 || cellType>=5) return;
    
    // calculate the right values for Ta and T0...
    double A0_star = 9999999;
    double Ta = K2D*(A0_star-A0);
    double T0 = -K2D*A0;
    
    // and change the tissue propertieso
    T->MechProp.CellPropVector[cellType].T_a = Ta;
    T->MechProp.CellPropVector[cellType].T0_a = T0;
    T->MechProp.CellPropVector[cellType].A0_a = A0_star;
    
    updateMainWindowData();
    
}

// action: File->Open, uses loadFile
void Integrator::saveCrossSectionsAsPDFActFromSubDirectories(QString dirName, QString extension, double px, double py, double pz, double nx, double ny, double nz)
{
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
            //if(true){
            if(fileName.endsWith(".3dvm") | fileName.endsWith(".3dmv")) {
                int flag = T->loadFileFromTxt(dirName+"/"+fileName);
                signal_updateWindow();
                if(flag==1) { // if the file could be opened properly
                    //std::cout << dirName.append("/").remove(".txt").append(".pdf").toStdString() << std::endl;
                    signal_writeImage(Point(px,py,pz),Point(nx,ny,nz),extension+fileName.left(fileName.size()-5),dirName);
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


void Integrator::saveCrossSectionsForGastrulation(QString dirName)
{
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
            //if(true){
            if(fileName.endsWith(".3dvm") | fileName.endsWith(".3dmv")) {
                int flag = T->loadFileFromTxt(dirName+"/"+fileName);
                signal_updateWindow();
                if(flag==1) { // if the file could be opened properly
                    //std::cout << dirName.append("/").remove(".txt").append(".pdf").toStdString() << std::endl;
                    // save sagittal cross section
                    signal_writeImage(Point(0,0,0),Point(1,0,0),"sagCS_1"+fileName.left(fileName.size()-5),dirName);
                    // save DV sections
                    signal_writeImage(Point(0,0,0),Point(0,0,1),"dvCS_1_0"+fileName.left(fileName.size()-5),dirName);
                    signal_writeImage(Point(0,0,-50),Point(0,0,1),"dvCS_2_m50_"+fileName.left(fileName.size()-5),dirName);
                    signal_writeImage(Point(0,0,50),Point(0,0,1),"dvCS_3_p50_"+fileName.left(fileName.size()-5),dirName);
                    signal_writeImage(Point(0,0,-100),Point(0,0,1),"dvCS_4_m100_"+fileName.left(fileName.size()-5),dirName);
                    signal_writeImage(Point(0,0,100),Point(0,0,1),"dvCS_5_p100_"+fileName.left(fileName.size()-5),dirName);
                    signal_writeImage(Point(0,0,-150),Point(0,0,1),"dvCS_6_m150_"+fileName.left(fileName.size()-5),dirName);
                    signal_writeImage(Point(0,0,150),Point(0,0,1),"dvCS_7_p150_"+fileName.left(fileName.size()-5),dirName);

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

//void Integrator::analysis_varying_Tl()
//{
//    // folder to save the results
//    QString folderName = "/Users/silvanus/Programming/Tissue/results/report/vary_Tl_c/";
//    
//    // parameters of the analysis
//    int rows = 20;
//    int cols = 20;
//    double width = 10;
//    double length = 10;
//    double radius = rows*width/8. + 20;
//    double height = 100;
//    double K=0.001;
//    double T_l = 1;
//    double T_a = 0.5;
//    double T_b = 0.5;
//    double G_a= 0;
//    double G_b= 0;
//    double G_l = 0;
//    
//    bool STIFF_BASAL_MEMBRANE = 0; // dont want a stiff membrane
//    
//    // calculate the normalized parameter
//    double V0 = width*length*height; // volume of the cells
//    double l0 = std::pow(V0,1./3.); // length scale in the system
//    T_a *= T_l;
//    T_b *= T_l;
//    G_a *= T_l*l0;
//    G_b *= T_l*l0;
//    
//    // vary the interface surface tension around the cyst
//    double T_l_cyst_init = 1;
//    double T_l_cyst_last = 4;
//    double T_l_cyst_steps = 0.2;
//
//    // write output to file
//    std::ofstream output;
//    QString fileName = folderName;
//    fileName.append(QString("result.txt"));
//    
//    QByteArray ba = fileName.toLocal8Bit();
//    const char *fileName_char = ba.data();
//	output.open(fileName_char);// initial position is set to be at the end of the file
//	
//    // Write into file
//	if(output.is_open())
//	{
//        output << "Results for cyst formation in hexagonal tissue \n Parameters:\n";
////        output << "T_a = " << T_a << std::endl;
////        output << "T_b = " << T_b << std::endl;
////        output << "K = " << K << std::endl;
////        output << "G_a = " << G_a << std::endl;
////        output << "G_b = " << G_b << std::endl;
////        output << "rows = " << rows << std::endl;
////        output << "columns = " << cols << std::endl;
////        output << "width = " << width << std::endl;
////        output << "length = " << length << std::endl;
////        output << "ftol = " << minPar.ftol << std::endl;
////        output << "basalMembrane is straight? " << T->STIFF_BASAL_MEMBRANE << std::endl;
//        output << "T = " << T_a << std::endl;
//        output << "G = " << G_a << std::endl;
//        output << "Simulation results for the cyst for varying interface surface tensions\n(T_l^c apicalCircumference basalCircumference interfaceArea height apicalIndentation basalIndentation apicalDisplacement basalDisplacement)" << std::endl;
//        
//    }
//    
//    createHexagonalTissue(0, true, rows, cols, height, width, length, K, T_a, T_b, T_l, G_a, G_b, G_l);
//    T->mutate_cells_in_middle(radius);
//
//    std::vector<double> T_l_cyst_vector, cyst_height, apical_circumference, basal_circumference;
//        
//    for(double T_l_cyst = T_l_cyst_init; T_l_cyst <= T_l_cyst_last+0.5*T_l_cyst_steps; T_l_cyst += T_l_cyst_steps) { // main analysis loop
//        
//        // start by creating a hexagonally packed tissue
//        
//        T->STIFF_BASAL_MEMBRANE = STIFF_BASAL_MEMBRANE;
//        
//        // change interface properties
//        T->MechProp.Intercell_SurfaceTension[2][1] = T_l_cyst;
//        T->setMechanicalProperties();
//        
//        // relax the tissue
//        //minimize();
//        
//        // change the cells in the middle to type 2
//        //T->mutate_cells_in_middle(radius);
//        
//        // relax the tissue
//        minimize();
//        
//        // save current tissue state
//        QString fileName = folderName;
//        fileName.append(QString::number(T_l_cyst));
//        QByteArray ba = fileName.toLocal8Bit();
//        const char *fileName_char = ba.data();
//        T->writeFileToTxt(fileName_char);
//        
//        // get readout of cyst properties
//        T->updateData();
//        if(output.is_open()) { // write current cyst shape
//            output << T_l_cyst << " " << T->tissueData.apicalCircumference << " " << T->tissueData.basalCircumference << " " << T->tissueData.interfaceArea << " " << T->tissueData.max_height_2 << " " << T->tissueData.apicalIndentation << " " << T->tissueData.basalIndentation << " " << T->tissueData.apicalMaxZDistance << " " << T->tissueData.basalMaxZDistance <<  std::endl;
//        }
//
//        std::cout<< "T_l_cyst = " << T_l_cyst << "; h = " << T->tissueData.max_height_2 << "; AC = " << T->tissueData.apicalCircumference << "; BC = " << T->tissueData.basalCircumference << std::endl;
//    }
//    
//        // extract the statistics from the tissue data
//        //double mean_height_1, var_height_1, mean_height_2, var_height_2;
//        //Tissue::vector_statistics(T->tissueData.height_type1, mean_height_1, var_height_1);
//        //Tissue::vector_statistics(T->tissueData.height_type2, mean_height_2, var_height_2);
//        
//        //height_ratio.push_back(mean_height_1/mean_height_2);
//        //std::cout << "factor = " << G_a*Ga_factor << " , height ratio = " << mean_height_1/mean_height_2 << std::endl;
//       
////       QString fileName = ;
////       
////       //std::cout << "save\n";
////       
////       if (!fileName.isEmpty())
////       {
////           // transform the QString fileName to a const char *
////           QByteArray ba = fileName.toLocal8Bit();
////           const char *fileName_char = ba.data();
////           signal_writeToTxt(fileName_char);
////       }
////    }
//
//    if(output.is_open())  output.close();
//
//   
//}

//void Integrator::analysis_symmetric()
//{
//    // parameters of the analysis that are not varied
//    int rows = 30;
//    int cols = 30;
//    double width = 10;
//    double length = 10;
//    double radius = 55;//rows*width/8. + 20;
//    double height = 100;
//    double K=0.001;
//    double T_l = 1;
//    double G_l = 0;
//    double V0 = width*length*height; // volume of the cells
//    double l0 = std::pow(V0,1./3.); // length scale in the system
//    bool STIFF_BASAL_MEMBRANE = 0; // dont want a stiff membrane
//
//    // arrays of surface and line tensions over which the parameter space shall be explored
//    double T_ab[] = {0,0.5,1,2,3,5,10};
//    double G_ab[] = {0,0.5,1,2,3,5,10};
//    
//    // ----------------- first part vary lateral surface tensions for different parameters -------------------- //
//    // folder to save the results
//    std::string topFolderName("/Users/silvanus/Programming/Tissue/results/symmetryBreaking/cyst38_new/");
//
//    // range to vary the interface surface tension around the cyst
//    double T_l_cyst_init = 1;
//    double T_l_cyst_last = 10;
//    double T_l_cyst_steps = 0.3;
//    
//    // now iterate over those lists
//    for(int i = 0; i<7; i++) { // iterate over surface tensions
//        
//        double Tab = T_ab[i]*T_l; // denormalized surface tension
//        
//        for(int j = 0; j<7; j++) {  // iterate over line tensions
//        
//            double Gab = G_ab[j]*l0; // denormalized line tension
//
//            if(Gab==0 && Tab ==0) continue; // not useful when no apical forces whatsoever
//            
//            // define name of current folder where all the simulation results shall be placed
//            std::ostringstream strs;
//            strs << topFolderName << "T=" << T_ab[i] << ",G=" << G_ab[j] << "/";
//            std::string folderName  = strs.str();
//            mkdir(folderName.c_str(), S_IRWXU);
//            std::string fileName = folderName;
//            fileName.append("result_Tl_c.txt");
//            
//            std::ofstream output;
//            output.open(fileName.c_str());// initial position is set to be at the end of the file
//            
//            // Write into file
//            if(output.is_open())  {
//                output << "Results for cyst formation in hexagonal tissue with apical-basal symmetry \n Parameters:\n";
//                output << "T = " << T_ab[i] << std::endl;
//                output << "G = " << G_ab[j] << std::endl;
//                output << "Simulation results for the cyst for varying interface surface tensions\n(T_l^c apicalCircumference basalCircumference interfaceArea height apicalIndentation basalIndentation apicalDisplacement basalDisplacement)" << std::endl;
//            }
//            
//            // create tissue for the first time, relax it and create cyst
//            createHexagonalTissue(0, true, rows, cols, height, width, length, K, Tab, Tab, T_l, Gab, Gab, G_l);
//            T->STIFF_BASAL_MEMBRANE = STIFF_BASAL_MEMBRANE;
//            minimize();
//            T->mutate_cells_in_middle(radius);
//            
//            // iterate through interface tensions
//            for(double T_l_cyst = T_l_cyst_init; T_l_cyst <= T_l_cyst_last+0.6*T_l_cyst_steps; T_l_cyst += T_l_cyst_steps) { // main analysis loop
//                
//                // change interface properties
//                T->MechProp.Intercell_SurfaceTension[2][1] = T_l_cyst*T_l;
//                T->setMechanicalProperties();
//                
//                // relax the tissue
//                minimize();
//                
//                // save current tissue state
//                strs.str("");
//                strs << folderName << "Tl_c=" << T_l_cyst << ".txt";
//                fileName = strs.str();
//                T->writeFileToTxt(fileName.c_str());
//                
//#ifdef _USE_QT_
//                signal_writeImage(Point(T->SystemSize.x/2, T->SystemSize.y/2, 0),PI/2,PI/2,"Tl_c="+QString::number(T_l_cyst),QString::fromStdString(folderName));
//#endif //8_USE_QT_
//                // get readout of cyst properties
//                T->updateData();
//                if(output.is_open()) { // write current cyst shape
//                    //output << T_l_cyst << " " << T->tissueData.apicalCircumference << " " << T->tissueData.basalCircumference << " " << T->tissueData.interfaceArea << " " << T->tissueData.max_height_2 << " " << T->tissueData.apicalIndentation << " " << T->tissueData.basalIndentation << " " << T->tissueData.apicalMaxZDistance << " " << T->tissueData.basalMaxZDistance << std::endl;
//                }
//            }
//            
//            if(output.is_open())  output.close();
//        }
//        
//    }
//    
//    // ----------------- second part: vary apical and basal line tensions for different parameters -------------------- //
//    // folder to save the results
//    topFolderName = "/Users/silvanus/Programming/Tissue/results/symmetric/big/";
//
//    // range to vary the interface surface tension around the cyst
//    double G_cyst_init = 1;
//    double G_cyst_last = 10;
//    double G_cyst_steps = 0.3;
//    
//    
//    // list of surface and line tensions over which the parameter space shall be explored
//    
//    // now iterate over those lists
//    for(int i = 0; i<7; i++) { // iterate over surface tensions
//        
//        double Tab = T_ab[i]*T_l; // denormalized surface tension
//        
//        for(int j = 0; j<7; j++) {  // iterate over line tensions
//            
//            double Gab = G_ab[j]*l0; // denormalized line tension
//            
//            if(T_ab[i] == 0 & G_ab[j]<10) continue;
//            
//            // define name of current folder where all the simulation results shall be placed
//            std::ostringstream strs;
//            strs << topFolderName << "T=" << T_ab[i] << "G=" << G_ab[j] << "/";
//            std::string folderName = strs.str();
//            // create directory to put the data there
//            mkdir(folderName.c_str(), S_IRWXU);
//            
//            std::string fileName = folderName;
//            fileName.append("result_G_c.txt");
//            // define the output name and write output to file
//            std::ofstream output;
//            output.open(fileName.c_str());// initial position is set to be at the end of the file
//            
//            // Write into file
//            if(output.is_open())  {
//                output << "Results for cyst formation in hexagonal tissue with apical-basal symmetry \n Parameters:\n";
//                output << "T = " << T_ab[i] << std::endl;
//                output << "G = " << G_ab[j] << std::endl;
//                output << "Simulation results for the cyst for varying apical and basal line tensions around the cysts\n(G^c apicalCircumference basalCircumference interfaceArea height apicalIndentation basalIndentation apicalDisplacement basalDisplacement)" << std::endl;
//            }
//            
//            // create tissue for the first time, relax it and create cyst
//            createHexagonalTissue(0, true, rows, cols, height, width, length, K, Tab, Tab, T_l, Gab, Gab, G_l);
//            T->STIFF_BASAL_MEMBRANE = STIFF_BASAL_MEMBRANE;
//            minimize();
//            T->mutate_cells_in_middle(radius);
//            
//            // iterate through interface tensions
//            for(double G_cyst = G_cyst_init; G_cyst <= G_cyst_last+0.6*G_cyst_steps; G_cyst += G_cyst_steps) { // main analysis loop
//                
//                // change interface properties
//                T->MechProp.Intercell_ApicalLineTension[2][1] = G_cyst*Gab;
//                T->MechProp.Intercell_BasalLineTension[2][1] = G_cyst*Gab;
//                T->setMechanicalProperties();
//                
//                // relax the tissue
//                minimize();
//                
//                // save current tissue state
//                strs.str("");
//                strs << folderName << "G_c=" << G_cyst<<".txt";
//                fileName = strs.str();
//                T->writeFileToTxt(fileName.c_str());
//                
//#ifdef _USE_QT_
//                signal_writeImage(Point(T->SystemSize.x/2, T->SystemSize.y/2, 0),PI/2,PI/2,"G_c="+QString::number(G_cyst),QString::fromStdString(folderName));
//#endif //_USE_QT_
//                // get readout of cyst properties
//                T->updateData();
//                if(output.is_open()) { // write current cyst shape
//                    //output << G_cyst << " " << T->tissueData.apicalCircumference << " " << T->tissueData.basalCircumference << " " << T->tissueData.interfaceArea << " " << T->tissueData.max_height_2 << " " << T->tissueData.apicalIndentation << " " << T->tissueData.basalIndentation << " " << T->tissueData.apicalMaxZDistance << " " << T->tissueData.basalMaxZDistance << std::endl;
//                }
//            }
//            
//            if(output.is_open())  output.close();
//        }
//        
//    }
//    
//    
//    
//    
//    
//}
//
//void Integrator::analysis_asymmetric()
//{
//    // parameters of the analysis that are not varied
//    int rows = 10;
//    int cols = 10;
//    double width = 10;
//    double length = 10;
//    double radius = 55;// rows*width/8. + 20;
//    double height = 100;
//    double K=0.001;
//    double T_l = 1;
//    double G_l = 0;
//    double V0 = width*length*height; // volume of the cells
//    double l0 = std::pow(V0,1./3.); // length scale in the system
//    bool STIFF_BASAL_MEMBRANE = 0; // dont want a stiff membrane
//    
//    // arrays of surface and line tensions over which the parameter space shall be explored
//    double T_ab[] = {0.5,1,3,5};
//    double G_a[] = {0.5,1,3,5};
//    
//    // ----------------- first part vary lateral surface tensions for different parameters -------------------- //
//    // folder to save the results
//    std::string topFolderName("/Users/silvanus/Programming/Tissue/results/asymmetric/Gb=0,Ta=Tb,big/");
//    // range to vary the interface surface tension around the cyst
//    double T_l_cyst[] = {1, 1.2, 1.5, 2, 3};
//    double G_a_cyst[] = {1, 1.2, 1.5, 2, 3};
//    
//    // now iterate over those lists
//    for(unsigned i = 0; i<4; i++) { // iterate over surface tensions
//        
//        double Tab = T_ab[i]*T_l; // denormalized apical surface tension
//        //double Tb = Ta; // basal surface tension with the same magnitude as the apical surface tension
//        
//        for(unsigned j = 0; j<4; j++) {  // iterate over line tensions
//            
//            double Ga = G_a[j]*l0*T_l; // denormalized line tension
//            
//            // define name of current folder where all the simulation results shall be placed
//            std::ostringstream strs;
//            strs << topFolderName << "Tab=" << T_ab[i] << ",Ga=" << G_a[j] << "/";
//            
//            std::string folderName  = strs.str();
//            
//            mkdir(folderName.c_str(), S_IRWXU);
//            
//            
//            for(unsigned k=0; k<5; k++) {  // now loop over the varying lateral mut-wt surface tension
//                
//                double Tl_c = T_l_cyst[k]*T_l;
//                
//                strs.str("");
//                
//                // create a folder for this configuration: Tab, Ga, T_l_cyst
//                strs << folderName << "Tl_c=" << T_l_cyst[k] << "/";
//                
//                std::string folderName_Tlc = strs.str();
//                
//                // create directory to put the data there
//                mkdir(folderName_Tlc.c_str(), S_IRWXU);
//
//                 // define the output name and write output to file
//                
//                std::string fileName = folderName_Tlc;
//                fileName.append("results.txt");
//                std::ofstream output;
//                output.open(fileName.c_str());// initial position is set to be at the end of the file
//                
//                // Write into file
//                if(output.is_open())  {
//                    output << "Results for cyst formation in hexagonal tissue without basal line tension, and same apical and basal surface tension for varying boundary effects \n Parameters:\n";
//                    output << "Ta = Tb = " << T_ab[i] << std::endl;
//                    output << "Tl_cyst = " << T_l_cyst[k] << std::endl;
//                    output << "Ga = " << G_a[j] << std::endl;
//                    output << "Gb = 0" << std::endl;
//                    output << "Simulation results for the cyst for varying apical line tensions Ga_c around the cyst\n(Ga_c apicalCircumference basalCircumference interfaceArea height apicalIndentation basalIndentation apicalDisplacement basalDisplacement apicalSurface basalSurface)" << std::endl;
//                }
//                
//                // create tissue for the first time, relax it and create cyst
//                createHexagonalTissue(0, true, rows, cols, height, width, length, K, Tab, Tab, T_l, Ga, 0, G_l);
//                
//                T->MechProp.Intercell_SurfaceTension[2][1] = Tl_c;
//                T->setMechanicalProperties();
//                
//                T->STIFF_BASAL_MEMBRANE = false;
//                minimize();
//                T->mutate_cells_in_middle(radius);
//                
//                // iterate through interface tensions
//                for(unsigned l = 0; l<5; l++) { // main analysis loop
//                    
//                    // change interface properties
//                    T->MechProp.Intercell_ApicalLineTension[2][1] = G_a_cyst[l]*Ga;
//                    T->setMechanicalProperties();
//                    
//                    // relax the tissue
//                    minimize();
//                    
//                    // save current tissue state
//                    strs.str("");
//                    strs << folderName_Tlc << "Ga_c=" << G_a_cyst[l] << ".txt";
//                    fileName = strs.str();
//                    T->writeFileToTxt(fileName.c_str());
//
//#ifdef _USE_QT_
//                    signal_writeImage(Point(T->SystemSize.x/2, T->SystemSize.y/2, 0),PI/2,PI/2,"Ga_c="+QString::number(G_a_cyst[l]),QString::fromStdString(folderName_Tlc));
//#endif // _USE_QT_
//                    // get readout of cyst properties
//                    T->updateData();
//                    if(output.is_open()) { // write current cyst shape
//                        //output << G_a_cyst[l] << " " << T->tissueData.apicalCircumference << " " << T->tissueData.basalCircumference << " " << T->tissueData.interfaceArea << " " << T->tissueData.max_height_2 << " " << T->tissueData.apicalIndentation << " " << T->tissueData.basalIndentation << " " << T->tissueData.apicalMaxZDistance << " " << T->tissueData.basalMaxZDistance << " " << T->tissueData.apicalArea_2 << " " << T->tissueData.basalArea_2 << std::endl;
//                    }
//                    
//                    // stop criterion: if the apical side of the cyst is completely vanished the energy does not change with changing intercell apical line tension but the simulations there have some instabilities since the cells are 'one point'
//                    //if(T->tissueData.apicalCircumference*T->MechProp.Intercell_ApicalLineTension[2][1]<1) break;
//                }
//                
//                if(output.is_open())  output.close();
//            }
//        }
//    }
//}
//
//void Integrator::analysis_general(unsigned rows, unsigned cols, double radius_cyst, double Ta, double Tb, double Ga, double Gb, double Tl_c)
//{
//    // parameters of the analysis that are not varied
//    double width = 10;
//    double length = 10;
//    double height = 100;
//    double K=0.001;
//    double T_l = 1;
//    double G_l = 0;
//    double V0 = width*length*height; // volume of the cells
//    double l0 = std::pow(V0,1./3.); // length scale in the system
//    
//    Ta  *= T_l; // denormalized surface tension
//    Tb  *= T_l;
//    Ga *= l0*T_l; // denormalized line tension
//    Gb *= l0*T_l;
//    Tl_c *= T_l;
//    
//    // ----------------- first part vary lateral surface tensions for different parameters -------------------- //
//    // folder to save the results
//    std::string topFolderName("/Users/silvanus/Programming/Tissue/results/Gac_vs_Tlc/");
//    
//    // range to vary the interface surface tension around the cyst
//    double G_cyst_init = 1;
//    double G_cyst_last = 4;
//    double G_cyst_steps = 0.2;
//    
//    // define name of current folder where all the simulation results shall be placed
//    std::ostringstream strs;
//    strs << topFolderName << "Tlc_" << std::setprecision(2) << Tl_c;
//    std::string fileName = strs.str();
//    // define the output name and write output to file
//    std::ofstream output;
//    output.open(fileName.c_str());// initial position is set to be at the end of the file
//    
//    // Write into file
//    if(output.is_open())  {
//        output << "Results for cyst formation in hexagonal tissue with apical-basal symmetry \n Parameters:\n";
//        output << "Ta = " << Ta << std::endl;
//        output << "Tb = " << Tb << std::endl;
//        output << "Ga = " << Ga << std::endl;
//        output << "Gb = " << Gb << std::endl;
//        output << "Tlc =" << Tl_c << std::endl;
//        output << "Simulation results for the cyst for varying line and surface tensions around the cysts\n(G^c apicalCircumference basalCircumference interfaceArea height apicalIndentation basalIndentation apicalDisplacement basalDisplacement)" << std::endl;
//    }
//    
//    // create tissue for the first time, relax it and create cyst
//    createHexagonalTissue(0, true, rows, cols, height, width, length, K, Ta, Tb, T_l, Ga, Gb, G_l);
//    T->STIFF_BASAL_MEMBRANE = 0;
//    minimize();
//    T->mutate_cells_in_middle(radius_cyst);
//    
//    // iterate through interface tensions
//    for(double G_cyst = G_cyst_init; G_cyst <= G_cyst_last+0.6*G_cyst_steps; G_cyst += G_cyst_steps) { // main analysis loop
//		
//		// change interface properties
//		T->MechProp.Intercell_ApicalLineTension[2][1] = G_cyst*Ga;
//		T->MechProp.Intercell_SurfaceTension[2][1] = Tl_c*T_l;
//        //T->MechProp.Intercell_BasalLineTension[2][1] = G_cyst*Gb;
//		T->setMechanicalProperties();
//		
//		// relax the tissue
//		minimize();
//		
//		// save current tissue state
//		strs.str("");
//		strs << topFolderName << "Tlc_" << Tl_c << "_Gc_" << G_cyst << ".txt";
//		fileName = strs.str();
//		T->writeFileToTxt(fileName.c_str());
//		
//#ifdef _USE_QT_
//		signal_writeImage(Point(T->SystemSize.x/2, T->SystemSize.y/2, 0),PI/2,PI/2,"Tlc_"+QString::number(Tl_c)+"_Gc_"+QString::number(G_cyst),QString::fromStdString(topFolderName));
//#endif //_USE_QT_
//		// readout of cyst properties
//		T->updateData();
//		if(output.is_open()) { // write current cyst shape
//		    output << G_cyst << " " << T->tissueData.apicalCircumference << " " << T->tissueData.basalCircumference << " " << T->tissueData.interfaceArea << " " << T->tissueData.max_height_2 << " " << T->tissueData.apicalIndentation << " " << T->tissueData.basalIndentation << " " << T->tissueData.apicalMaxZDistance << " " << T->tissueData.basalMaxZDistance << std::endl;
//		}
//    }
//    
//    if(output.is_open())  output.close();
//}
//

//void Integrator::analysis_general(unsigned rows, unsigned cols, double radius_cyst, double Ta, double Tb, double Ga, double Gb, double Tl_c, QString folderName)
//{
//    // parameters of the analysis that are not varied
//    double width = 10;
//    double length = 10;
//    double height = 100;
//    double K=0.001;
//    double T_l = 1;
//    double G_l = 0;
//    double V0 = width*length*height; // volume of the cells
//    double l0 = std::pow(V0,1./3.); // length scale in the system
//    
//    Ta  *= T_l; // denormalized surface tension
//    Tb  *= T_l;
//    Ga *= l0*T_l; // denormalized line tension
//    Gb *= l0*T_l;
//    Tl_c *= T_l;
//    
//    // ----------------- first part vary lateral surface tensions for different parameters -------------------- //
//    // range to vary the interface surface tension around the cyst
//    double G_cyst_init = 1;
//    double G_cyst_last = 2;
//    double G_cyst_steps = 0.1;
//    
//    // define name of current folder where all the simulation results shall be placed
//    // and create the directory if it doesnt exist
//    QDir dir;
//    dir.mkdir(folderName);
//    
//    std::string fileName = folderName.toStdString()+"/results.txt";
//    // define the output name and write output to file
//    std::ofstream output;
//    output.open(fileName.c_str());// initial position is set to be at the end of the file
//    
//    // Write into file
//    if(output.is_open())  {
//        output << "Results for cyst formation in hexagonal tissue with apical-basal symmetry \n Parameters:\n";
//        output << "Ta = " << Ta << std::endl;
//        output << "Tb = " << Tb << std::endl;
//        output << "Ga = " << Ga << std::endl;
//        output << "Gb = " << Gb << std::endl;
//        output << "Tlc =" << Tl_c << std::endl;
//        output << "Simulation results for the cyst for varying line and surface tensions around the cysts\n(G^c apicalCircumference basalCircumference interfaceArea height apicalIndentation basalIndentation apicalDisplacement basalDisplacement)" << std::endl;
//    }
//    
//    // create tissue for the first time, relax it and create cyst
//    createHexagonalTissue(0, true, rows, cols, height, width, length, K, Ta, Tb, T_l, Ga, Gb, G_l);
//    T->STIFF_BASAL_MEMBRANE = 0;
//    T->MechProp.Intercell_SurfaceTension[2][1] = Tl_c*T_l;
//    minimize();
//    T->mutate_cells_in_middle(radius_cyst);
//    
//    // iterate through interface tensions
//    for(double G_cyst = G_cyst_init; G_cyst <= G_cyst_last+0.6*G_cyst_steps; G_cyst += G_cyst_steps) { // main analysis loop
//		
//		// change interface properties
//		T->MechProp.Intercell_ApicalLineTension[2][1] = G_cyst*Ga;
//		//T->MechProp.Intercell_BasalLineTension[2][1] = G_cyst*Gb;
//		T->MechProp.Intercell_BasalLineTension[2][1] = Gb;
//		T->setMechanicalProperties();
//		
//		// relax the tissue
//		minimize();
//		
//		// save current tissue state
//        std::stringstream strs;
//		strs <<folderName.toStdString() << std::setprecision(2) << "Tlc_" << Tl_c << "_Gc_" << G_cyst << ".txt";
//		fileName = strs.str();
//		T->writeFileToTxt(fileName.c_str());
//		
//#ifdef _USE_QT_
//		//signal_writeImage(Point(T->SystemSize.x/2, T->SystemSize.y/2, 0),PI/2,PI/2,"Tlc_"+QString::number(Tl_c)+"_Gc_"+QString::number(G_cyst),QString::fromStdString(folderName));
//#endif //_USE_QT_
//		// readout of cyst properties
//		T->updateData();
//		if(output.is_open()) { // write current cyst shape
//		    //output << G_cyst << " " << T->tissueData.apicalCircumference << " " << T->tissueData.basalCircumference << " " << T->tissueData.interfaceArea << " " << T->tissueData.max_height_2 << " " << T->tissueData.apicalIndentation << " " << T->tissueData.basalIndentation << " " << T->tissueData.apicalMaxZDistance << " " << T->tissueData.basalMaxZDistance << std::endl;
//		}
//    }
//    
//    if(output.is_open())  output.close();
//}
//

//
//
//void Integrator::analysis_general(double Tl_c, QString tissueName, QString folderName)
//{
//    // range in which to vary the interface surface tension around the cyst
//    double G_cyst_init = 1;
//    double G_cyst_last = 4;
//    double G_cyst_steps = 0.25;
//    
//    int flag = T->loadFileFromTxt(tissueName); // load the tissue from file tissueName
//    if(flag!=1) return; // make sure the tissue has been loaded
//    
//
//    // in the tissue set the surface tension around the clone to Tl_c
//    T->MechProp.Intercell_SurfaceTension[2][1] = Tl_c*T->MechProp.Intercell_SurfaceTension[1][1];
//    
//    
//    // define name of current folder where all the simulation results shall be placed
//    // and create the directory if it doesnt exist
//    QDir dir;
//    dir.mkdir(folderName);
//    std::string fileName = folderName.toStdString()+"/results.txt";
//    
//    // define the output name and write output to file
//    std::ofstream output;
//    output.open(fileName.c_str());// initial position is set to be at the end of the file
//    // Write into file
//    if(output.is_open())  {
//        output << "Results for cyst formation in hexagonal tissue with apical-basal symmetry. Parameters are:\n";
//        output << "Ta = " << T->MechProp.CellPropVector[1].T_a << std::endl;
//        output << "Tb = " << T->MechProp.CellPropVector[1].T_b << std::endl;
//        output << "Ga = " << T->MechProp.Intercell_ApicalLineTension[1][1] << std::endl;
//        output << "Gb = " << T->MechProp.Intercell_BasalLineTension[1][1] << std::endl;
//        output << "Text = " << T->T_ext << std::endl;
//        output << "Tlc = " << Tl_c << std::endl;
//        output << "Simulation results for the cyst for varying line and surface tensions around the cysts\n(G^c apicalCircumference basalCircumference interfaceArea height apicalIndentation basalIndentation apicalDisplacement basalDisplacement)" << std::endl;
//    }
//    
//    // iterate through interface tensions
//    for(double G_cyst = G_cyst_init; G_cyst <= G_cyst_last+0.6*G_cyst_steps; G_cyst += G_cyst_steps) { // main analysis loop
//		
//		// change interface properties
//		T->MechProp.Intercell_ApicalLineTension[2][1] = G_cyst*T->MechProp.Intercell_ApicalLineTension[1][1];
//		//T->MechProp.Intercell_BasalLineTension[2][1] = G_cyst*Gb;
//		T->MechProp.Intercell_BasalLineTension[2][1] = G_cyst*T->MechProp.Intercell_BasalLineTension[1][1];
//		T->setMechanicalProperties();
//		
//		// relax the tissue
//		minimize();
//		
//		// save current tissue state
//        std::stringstream strs;
//        strs << std::setprecision(2) << std::fixed; // set the precision such that the output always has the form 1.00, 1.25, ...
//		strs <<folderName.toStdString() << "/Gc=" << std::setprecision(2) << std::fixed << G_cyst << ".txt";
//		fileName = strs.str();
//		T->writeFileToTxt(fileName.c_str());
//		
//#ifdef _USE_QT_
//		//signal_writeImage(Point(T->SystemSize.x/2, T->SystemSize.y/2, 0),PI/2,PI/2,"Tlc_"+QString::number(Tl_c)+"_Gc_"+QString::number(G_cyst),QString::fromStdString(folderName));
//#endif //_USE_QT_
//		// readout of cyst properties
//		T->updateData();
//		if(output.is_open()) { // write current cyst shape
//		    //output << G_cyst << " " << T->tissueData.apicalCircumference << " " << T->tissueData.basalCircumference << " " << T->tissueData.interfaceArea << " " << T->tissueData.max_height_2 << " " << T->tissueData.apicalIndentation << " " << T->tissueData.basalIndentation << " " << T->tissueData.apicalMaxZDistance << " " << T->tissueData.basalMaxZDistance << std::endl;
//		}
//    }
//    
//    if(output.is_open())  output.close();
//}
//
//
//void Integrator::analysis_general(QString tissueName, QString fileName_save, double Ta_c, double Tb_c, double Tl_c, double Ga_c, double Gb_c, MinimizationParameters minParameters, bool ECM)
//{
//    int flag = T->loadFileFromTxt(tissueName); // load the tissue from file tissueName
//    if(flag!=1) {std::cout << "Loading Tissue failed!!" << std::endl; return;} // make sure the tissue has been loaded
//    
//    T->STIFF_BASAL_MEMBRANE = ECM;
//    
//    // set the cyst properties
//    T->MechProp.Intercell_SurfaceTension[2][1] = Tl_c;
//    T->MechProp.Intercell_ApicalLineTension[2][1] = Ga_c;
//    T->MechProp.Intercell_BasalLineTension[2][1] = Gb_c;
//    T->MechProp.CellPropVector[2].T_a = Ta_c;
//    T->MechProp.CellPropVector[2].T_b = Tb_c;
//    
//    T->setMechanicalProperties();
//    
//    minPar = minParameters;
//    
//	// relax the tissue
//	minimize();
//	
//    // save current tissue state
//    std::stringstream strs;
//    strs <<fileName_save.toStdString();
//    std::string fileName = strs.str();
//    T->writeFileToTxt(fileName.c_str());
//
//    T->updateData();
//
//    // define the output name and write output to file
//    std::ofstream output;
//    strs.str("");
//    strs << fileName_save.append("_result.txt").toStdString();
//    fileName = strs.str();
//    output.open(fileName.c_str());
//    // Write into file
//    if(output.is_open())  {
//        output << "Ta = " << T->MechProp.CellPropVector[1].T_a << ", Tb = " << T->MechProp.CellPropVector[1].T_b << "Ta_c = " << T->MechProp.CellPropVector[2].T_a << ", Tb_c = " << T->MechProp.CellPropVector[2].T_b << ", Text = " << T->T_ext << ", Tlc = " << Tl_c << std::endl;
//        output << "Ga = " << T->MechProp.Intercell_ApicalLineTension[1][1] << ", Gb = " << T->MechProp.Intercell_BasalLineTension[1][1] << "Ga_c = " << T->MechProp.Intercell_ApicalLineTension[2][1] << ", Gb_c = " << T->MechProp.Intercell_BasalLineTension[2][1] << std::endl;
//        output << "K = " << T->MechProp.CellPropVector[1].K << ", K_c = " << T->MechProp.CellPropVector[2].K << ", V0 = " << T->MechProp.CellPropVector[1].V0 << ", V0_c = " << T->MechProp.CellPropVector[2].V0 << std::endl;
//        output << "apicalCircumference basalCircumference interfaceArea height apicalIndentation basalIndentation apicalDisplacement basalDisplacement)" << std::endl;
//       // output << T->tissueData.apicalCircumference << " " << T->tissueData.basalCircumference << " " << T->tissueData.interfaceArea << " " << T->tissueData.max_height_2 << " " << T->tissueData.apicalIndentation << " " << T->tissueData.basalIndentation << " " << T->tissueData.apicalMaxZDistance << " " << T->tissueData.basalMaxZDistance << std::endl;
//        output.close();
//    }
//}

void Integrator::set_cellProperties(int type, double Ta, double T0_a, double A0_a, double Tb, double T0_b, double A0_b, double K, double V0)
{
    if(!std::isnan(Ta))      T->MechProp.CellPropVector[type].T_a = Ta;
    if(!std::isnan(T0_a)) T->MechProp.CellPropVector[type].T0_a = T0_a;
    else T->MechProp.CellPropVector[type].T0_a = Ta; // make sure that no area elasticity is added involuntarily
    if(!std::isnan(A0_a))    T->MechProp.CellPropVector[type].A0_a = A0_a;

    if(!std::isnan(Tb))      T->MechProp.CellPropVector[type].T_b = Tb;
    if(!std::isnan(T0_b))    T->MechProp.CellPropVector[type].T0_b = T0_b;
    else T->MechProp.CellPropVector[type].T0_b = Tb; // make sure that no area elasticity is added involuntarily
    if(!std::isnan(A0_b))    T->MechProp.CellPropVector[type].A0_b = A0_b;
    
    if(K>0)  T->MechProp.CellPropVector[type].K = K;
    if(V0>0) T->MechProp.CellPropVector[type].V0 = V0;
    
    T->setMechanicalProperties();
}

void Integrator::set_edgeProperties(int type1, int type2, double Tl, double Ga, double Gb)
{
    if(type1<1 || type1>NO_CELL_TYPES || type2<1 || type2>NO_CELL_TYPES) return;
    
    if(!std::isnan(Tl)) T->MechProp.Intercell_SurfaceTension[std::max(type1, type2)][std::min(type1, type2)] = Tl;
    if(!std::isnan(Ga)) T->MechProp.Intercell_ApicalLineTension[std::max(type1, type2)][std::min(type1, type2)] = Ga;
    if(!std::isnan(Gb)) T->MechProp.Intercell_BasalLineTension[std::max(type1, type2)][std::min(type1, type2)] = Gb;
    
    T->setMechanicalProperties();
}

//void Integrator::set_multiplyEdgeProperties(int type1, int type2, double Tl, double Ga, double Gb)
//{
//    if(Tl>=0) T->MechProp.Intercell_SurfaceTension[std::max(type1, type2)][std::min(type1, type2)] = Tl;
//    if(Ga>=0) T->MechProp.Intercell_ApicalLineTens§ion[std::max(type1, type2)][std::min(type1, type2)] = Ga;
//    if(Gb>=0) T->MechProp.Intercell_BasalLineTension[std::max(type1, type2)][std::min(type1, type2)] = Gb;
//    T->setMechanicalProperties();
//}


void Integrator::set_vertexProperties(double Gl)
{
    T->MechProp.Vertex_LineTension = Gl;
    T->setMechanicalProperties();
}

void Integrator::set_periodicityProperties(double T_ext, bool Lx_fix, bool Ly_fix)
{
    if(!std::isnan(T_ext)) T->T_ext = T_ext;
    if(!std::isnan(Lx_fix)) T->Lx_fix = Lx_fix;
    if(!std::isnan(Ly_fix)) T->Ly_fix = Ly_fix;
}

void Integrator::runScript(QString fileName)
{
    //start execution...
    QScriptEngine myEngine;
    
    // obtain the folder name of the script which is being run
    folderName_script = fileName.left(fileName.lastIndexOf('/'));
    folderName_script += '/';
    QByteArray ba = folderName_script.toLocal8Bit();
    const char *c_str2 = ba.data();

    myEngine.globalObject().setProperty("folderName", c_str2);
    
    // one argument has to be the current instance of the Integrator class
    QScriptValue integratorValue = myEngine.newQObject(this);
    std::cout << folderName_script.toStdString() << std::endl;
    
    myEngine.globalObject().setProperty("Integrator", integratorValue);
    
    QScriptValue mainWindowValue = myEngine.newQObject(this);
    
    if (!fileName.isEmpty())
    {
        QFile scriptFile(fileName);
        //CurrentRecordPath=fileName;//pour que les images soient enregistres par defaut la ou est le script
        //tissueWindow->CurrentRecordPath=CurrentRecordPath;
        if (!scriptFile.open(QIODevice::ReadOnly))
            return;
        QTextStream stream(&scriptFile);
        QString contents = stream.readAll();
        scriptFile.close();
        
        //QString contents_mod = modifyScript(contents);
        
        QScriptValue result=myEngine.evaluate(contents, fileName);
        if (myEngine.hasUncaughtException()) {
            int line = myEngine.uncaughtExceptionLineNumber();
            qDebug() << "uncaught exception at line" << line << ":" << result.toString();
        }
        
    }
}

void Integrator::runScriptFromBox(QString boxContent)
{
    if(boxContent.isEmpty()) return;
    QScriptEngine myEngine;
    QScriptValue integratorValue = myEngine.newQObject(this);
    myEngine.globalObject().setProperty("Integrator", integratorValue);
    QScriptValue mainWindowValue = myEngine.newQObject(this);
    
    QScriptValue result=myEngine.evaluate(modifyScript(boxContent));
    
    if (myEngine.hasUncaughtException()) {
        int line = myEngine.uncaughtExceptionLineNumber();
        qDebug() << "uncaught exception at line" << line << ":" << result.toString();
    }
}

QString Integrator::modifyScript(QString s)
{
    // split the string in its single commands
    QStringList semicolon_list = s.split(";", QString::SkipEmptyParts);
    
    QString s_new, currentSection;
    
    // now iterate through the string list, and add Integrator. at the  beginning if necessary
    for(int i=0; i<semicolon_list.size();i++)
    {
        if(semicolon_list.at(i).startsWith("Integrator.")) s_new.append(semicolon_list.at(i)+";");
        else s_new.append("Integrator."+semicolon_list.at(i)+";");
    }
    
    // modified script
    std::cout << s_new.toStdString() << std::endl;
    
    return s_new;
}

void Integrator::recordNextIntegration(QString folderName, QString style, unsigned rate, QString fileName)
{
    record_nextIntegration = true;
    record_path = folderName;
    record_rate = rate;

    record_name = fileName;
    
    // which information have to be saved
    record_3D = style.contains("3D");
    record_crossSection = style.contains("cross") || style.contains("section") || style.contains("cs");
    record_3dvm = style.contains("3dvm") || style.contains("txt");

}

void Integrator::updateWindows()
{
    // update all windows
    signal_updateMainWindow();
    
    if(useGUI){
        signal_update3DWindow();
        signal_updateCrossSection();
    }
}

void Integrator::moveAllVertices(double dx, double dy, double dz)
{
    Point D = Point(dx,dy,dz);
    for(std::list<Vertex>::iterator it_v = T->ListVertex.begin(); it_v!=T->ListVertex.end(); it_v++) it_v->moveVertex(D,D);
    T->update();
}

void Integrator::fix_outmost_Type2_cells()
{
    // create list of type 2 cells
    std::list<Cell*> Type2_cells;
    for(std::list<Cell>::iterator it_c = T->ListCell.begin();it_c!=T->ListCell.end();it_c++)
    {
        it_c->positionFixed = false;
        if(it_c->Type==2) Type2_cells.push_back(&(*it_c));
    }
    
    if(Type2_cells.empty()) return;
    
    // first find the outmost type 2 cells and then fix them
    double min_z_position = (*Type2_cells.begin())->BC_a.z;
    double max_z_position = (*Type2_cells.begin())->BC_a.z;
    Cell *min_z_cell = (*Type2_cells.begin());
    Cell *max_z_cell = (*Type2_cells.begin());
    
    for(std::list<Cell*>::iterator it_c = Type2_cells.begin();it_c!=Type2_cells.end();it_c++)
    {
        if((*it_c)->BC_a.z>max_z_position)
        {
            max_z_position = (*it_c)->BC_a.z;
            max_z_cell = (*it_c);
        }
        
        if((*it_c)->BC_a.z<min_z_position)
        {
            min_z_position = (*it_c)->BC_a.z;
            min_z_cell = (*it_c);
        }
    }
    
    // now fix the two outmost cells in z direction
    max_z_cell->positionFixed = true;
    min_z_cell->positionFixed = true;
}

void Integrator::fix_type2_cells_in_box(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max)
{
    // iterate through all cells, check if they are in the box and if so, fix them
    for(std::list<Cell>::iterator it_c = T->ListCell.begin();it_c!=T->ListCell.end();it_c++)
    {
        if(it_c->Type==2
           && it_c->BC_a.x > x_min && it_c->BC_a.x < x_max
           && it_c->BC_a.y > y_min && it_c->BC_a.y < y_max
           && it_c->BC_a.z > z_min && it_c->BC_a.z < z_max)
        {
            it_c->positionFixed = true;
        }
    }
}

void Integrator::fix_cells_in_box(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max)
{
    // iterate through all cells, check if they are in the box and if so, fix them
    for(std::list<Cell>::iterator it_c = T->ListCell.begin();it_c!=T->ListCell.end();it_c++)
    {
        
        if(   it_c->BC_a.x > x_min && it_c->BC_a.x < x_max
           && it_c->BC_a.y > y_min && it_c->BC_a.y < y_max
           && it_c->BC_a.z > z_min && it_c->BC_a.z < z_max)
        {
            it_c->positionFixed = true;
        }
    }
}

void Integrator::unfix_all_cells()
{
    // iterate through all cells and unfix them
    for(std::list<Cell>::iterator it_c = T->ListCell.begin();it_c!=T->ListCell.end();it_c++)
    {
        it_c->positionFixed = false;
    }
}

void Integrator::unfix_cells_in_box(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max)
{
    // iterate through all cells, check if they are in the box and if so, unfix them
    for(std::list<Cell>::iterator it_c = T->ListCell.begin();it_c!=T->ListCell.end();it_c++)
    {
        
        if(   it_c->BC_a.x > x_min && it_c->BC_a.x < x_max
           && it_c->BC_a.y > y_min && it_c->BC_a.y < y_max
           && it_c->BC_a.z > z_min && it_c->BC_a.z < z_max)
        {
            it_c->positionFixed = false;
        }
    }
}


void Integrator::set_ECM_Properties(double pos_a, double pos_b, double K, bool basallyBothWays)
{
    T->flatECM_apicalPosition = pos_a;
    T->flatECM_basalPosition = pos_b;

    T->flatECM_stiffness = K;
    T->flatECM_basalPosition_bothways = basallyBothWays;
}