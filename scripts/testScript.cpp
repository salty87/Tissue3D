//
//  testScript.cpp
//  
//
//  Created by Silvanus Alt on 25.02.14.
//
//

// set integration properties
Integrator.set_updateRate(1);                          // rate at which the widgets are updated
Integrator.set_MinimizationMode(2);                    // conjugate gradient == 2
Integrator.set_ftol(1e-14);                            // minimization convergence criterion - difference between succinct function evaluations
Integrator.set_GTOL(1e-12);                            // minimization convergence criterion - steepness of the gradient
Integrator.set_tol(1e-5);                              // tolerance in Brent line minimization    
Integrator.set_check_for_T1(0);                        // after how many steps shall the algorithm check for T1s
Integrator.set_ITMAX(10000);                           // maximum number of steps before the algorithm stops
Integrator.set_EPS(1e-18);                             // small value
Integrator.set_noise_stdDev(0.001);                    // std deviation of the positional noise applied on the vertices
Integrator.set_noise_numberOfNoiseApplications(10);    // number of steps in the beginning of the minimization before which the noise is applied


// create/load tissue
//Integrator.slot_createHexagonalTissue(0, 1, 6, 30, 100,10,10);
// alternatively: Integrator.slot_loadFromTxt(QString)




// mechanical properties of the tissue
// wt
Integrator.set_cellProperties( 1,
    
    /*Ta*/      1.5,
    /*T0_a*/    1.5,
    /*A0_a*/    0,
    
    /*Tb*/      1.5,
    /*T0_b*/    1.5,
    /*A0_b*/    0,
    
    /*K*/       0,
    /*V0*/      0 );

// mut
Integrator.set_cellProperties( 2,
                              
                              /*Ta*/      1.5,
                              /*T0_a*/    1.5,
                              /*A0_a*/    0,
                              
                              /*Tb*/      1.5,
                              /*T0_b*/    1.5,
                              /*A0_b*/    0,
                              
                              /*K*/       0,
                              /*V0*/      0 );

// wt-wt
Integrator.set_edgeProperties(1,1,
                              /*Tl*/        4,
                              /*Ga*/        20,
                              /*Gb*/        20);

// wt-mut
Integrator.set_edgeProperties(2,1,
                              /*Tl*/        5,
                              /*Ga*/        20,
                              /*Gb*/        20);

// mut-mut
Integrator.set_edgeProperties(2,2,
                              /*Tl*/        6,
                              /*Ga*/        20,
                              /*Gb*/        20);


Integrator.set_periodicityProperties(/*T_ext*/ 3,
                                     /*fix Lx*/ false,
                                     /*fix Ly*/ false );
Integrator.slot_evolve_N_steps();


// mutate patch in the middle
