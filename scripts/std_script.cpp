//
//  std_script.c
//  
//
//  Created by Silvanus Alt on 27.02.14.
//
//

folderName = '/Users/silvanus/MPI/EpithelialMechanics/CystFormation/Simulations/';

/* --------------- mechanical parameters ------------------ */
Ta =   1.5
T0a =  1.5
A0a = 100
Tb =   2.5
T0b =  2.5
A0b = 100
Text = 3

Tl = 1
Tl_cyst = 1.5

Ga = 30
Ga_cyst = 45

Gb = 20
Gb_cyst = 20

Lxfix = 0
Lyfix = 0

// set K and V0 to zero such that they are not changed
K = 0
V0 = 0


/* --------------- set integration properties ------------------ */
Integrator.set_updateRate(1);                          // rate at which the widgets are updated
Integrator.set_MinimizationMode(2);                    // conjugate gradient == 2
Integrator.set_ftol(1e-14);                            // minimization convergence criterion - difference between succinct function evaluations
Integrator.set_GTOL(1e-12);                            // minimization convergence criterion - steepness of the gradient
Integrator.set_tol(1e-5);                              // tolerance in Brent line minimization
Integrator.set_check_for_T1(0);                        // after how many steps shall the algorithm check for T1s
Integrator.set_ITMAX(10000);                           // maximum number of steps before the algorithm stops
Integrator.set_EPS(1e-18);                             // small value
Integrator.set_noise_stdDev(0.0);                    // std deviation of the positional noise applied on the vertices
Integrator.set_noise_numberOfNoiseApplications(0);    // number of steps in the beginning of the minimization before which the noise is applied




/* --------------- load or initialize tissue ------------------ */
// load tissue
//Integrator.slot_loadFromTxt('/Users/silvanus/MPI/EpithelialMechanics/CystFormation/Simulations/SizeDependency/startingTissue');
Integrator.slot_createHexagonalTissue(0, 1, 20, 20, 100,10,10);

// wt
Integrator.set_cellProperties( 1,
                              
                              /*Ta*/      Ta,
                              /*T0_a*/    T0a,
                              /*A0_a*/    A0a,
                              
                              /*Tb*/      Tb,
                              /*T0_b*/    T0b,
                              /*A0_b*/    A0b,
                              
                              /*K*/       K,
                              /*V0*/      V0 );

// mut
Integrator.set_cellProperties( 2,
                              
                              /*Ta*/      Ta,
                              /*T0_a*/    T0a,
                              /*A0_a*/    A0a,
                              
                              /*Tb*/      Tb,
                              /*T0_b*/    T0b,
                              /*A0_b*/    A0b,
                              
                              /*K*/       K,
                              /*V0*/      V0 );


// set wt-wt interface properties
Integrator.set_edgeProperties(1,1,
                              /*Tl*/        Tl,
                              /*Ga*/        Ga,
                              /*Gb*/        Gb);

// set wt-mut interface properties
Integrator.set_edgeProperties(2,1,
                              /*Tl*/        Tl_cyst,
                              /*Ga*/        Ga_cyst,
                              /*Gb*/        Gb_cyst);

// set mut-mut interface properties
Integrator.set_edgeProperties(2,2,
                              /*Tl*/        Tl,
                              /*Ga*/        Ga,
                              /*Gb*/        Gb);


// mechanical properties of the tissue
Integrator.set_periodicityProperties(/*T_ext*/  Text,
                                     /*fix Lx*/ Lxfix,
                                     /*fix Ly*/ Lyfix );


/* --------------- relax tissue and save result ------------------ */
Integrator.minimize();
initFileName = folderName+'init';
Integrator.writeToTxt(initFileName);


/* --------------- change tissue, relax it and save the result ------------------ */
resultFileName = folderName+'res';
for(i=20;i<21;i=i+1)
{
    // mutate i cells in the middle
    Integrator.mutate_cells_in_middle(i);
    
    // ... and relax
    Integrator.minimize();
    
    // save file
    Integrator.writeToTxt(resultFileName+i);
    
    // load initial tissue
    Integrator.loadFromTxt(initFileName);
}



