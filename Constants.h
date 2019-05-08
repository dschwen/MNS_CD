#pragma once

#include <sundials/sundials_types.h> // definition of realtype

const realtype kb = 8.617E-5; // Boltzmann
const realtype pi = 3.141592; // pi

const int numPhase = 2;     // Number of precipitating phases
const int numComp = 3;      // Number of precipitating components
const int numClass = 50000; // Number of cluster classes/maximum cluster size considered
const int runs = 70;        // Number of loops to run
const int CutoffSize = 65;  // Cutoff size used for output
const realtype T0 = 0.0;    // Initial time

// Method to calculate mean radius of precipitate, radM1 or radM2
#define RadiusCalc radM2

// number of calculating phases, including both homo- and heterogeneously nucleated phases
const int numCalcPhase = (numPhase * 2);
// number of ODEs
const int neq = (numCalcPhase * numClass + numComp);

const realtype RTOL = 1.0E-6; // rel tolerance
const realtype ATOL = 0.0;    // abs tolerance
