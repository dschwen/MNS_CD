/**
 * Code developed by Huibin Ke from UW-Madison for the evolution of Mn-Ni-Si
 * precipitates in Reactor Pressure Vessel steels.
 */

#include "Constants.h" // Constants header file
#include "Input.h"     // Input dd header file

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include <time.h>

// Sundials includes (requires version <= 2.7.0)
#include <cvode/cvode.h>             // main integrator header file
#include <cvode/cvode_band.h>        // band solver header
#include <nvector/nvector_serial.h>  // serial N_Vector types, fct. and macros
#include <sundials/sundials_math.h>  // contains the macros ABS, SQR, and EXP
#include <sundials/sundials_types.h> // definition of realtype

InputCondition ICond;    // Defined in Input.h, including irradiation conditions
InputMaterial IMaterial; // Defined in Input.h, including material information
InputProperty IProp;     // Defined in Input.h, including all other parameters used in model

double sqr(double a)
{
  return a * a;
}

// function defs

typedef struct
{
  realtype size[numClass], radClust[numPhase][numClass], beta[numPhase][numClass], delG[numPhase][numClass];
} * UserData;

void loadData(UserData &data);
void getInitVals(realtype y0[neq]);
void initParams();
void GetRED(realtype D[numComp], realtype Flux);
void printYVector(N_Vector y);
void getOutput(N_Vector y, realtype radM1[numCalcPhase], realtype radM2[numCalcPhase], realtype rhoC[numCalcPhase]);
void getbeta(realtype size[numClass], realtype beta[numPhase][numClass]);
void getSize(realtype size[numClass]);
void getRadClust(realtype size[numClass], realtype radClust[numPhase][numClass]);
void getDelG(realtype size[numClass], realtype radClust[numPhase][numClass], realtype delG[numPhase][numClass]);
void getFlux(UserData data, N_Vector y, realtype J[numCalcPhase][numClass + 1]);
int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
std::string int_to_string(int i);

// global problem parameters

// Diffusion coefficient,effective precipitate lattice constant, square of precipitate composition
realtype D[numComp], aP[numPhase], nu[numPhase][numComp];
// precipitate number density, two kinds of mean radius (see readme for more detail)
realtype rhoC[numCalcPhase], radM1[numCalcPhase], radM2[numCalcPhase];
// Irradiation flux, solute product
realtype Flux, solProd[numPhase];

int main()
{
  realtype y0data[neq], t;
  realtype tout = 1E0; // This is the solution time of each run.  The solver will solve for y0(t) on the range (0,tout)
  double ts = 0.0;     // This is the start time for each run, tf+tout will be the final time ater each run

  // The following block setup the information for the solver and initial all the parameters
  N_Vector y0 = NULL;
  ICond = NULL;
  IMaterial = NULL;
  IProp = NULL;
  UserData data = NULL;
  void *cvode_mem = NULL;
  int flag, itask, mxsteps;
  mxsteps = 2000000;
  long int mu = 3;
  long int ml = 3;
  realtype *yd;
  ICond = (InputCondition)malloc(sizeof *ICond);        // Allocate irradiation condition memory
  IMaterial = (InputMaterial)malloc(sizeof *IMaterial); // Allocate material properties memory
  IProp = (InputProperty)malloc(sizeof *IProp);         // Allocate other properties memory
  LoadInput(ICond, IMaterial, IProp);                   // Load all the input parameters
  initParams();
  getInitVals(y0data);
  y0 = N_VMake_Serial(neq, y0data);

  // Write the output file
  std::ofstream O_file;
  O_file.open("Output");
  O_file << "Run\tCalcTime(s)\tTime(s)\tFluence(n/m2s)\t";

  for (int p = 0; p < numPhase; ++p)
  {
    std::cout << p << " " << (p + 1) << " " << (p + 1) << std::endl;
    O_file << "Mean_Radius_of_Phase_" << (p + 1) << "_Homo(m)\tNumber_Density_of_Phase_" << (p + 1) << "_Homo(1/m3)\t";
  }

  for (int p = 0; p < numPhase; ++p)
    O_file << "Mean_Radius_of_Phase_" << (p + 1) << "_Heter(m)\tNumber_Density_of_Phase_" << (p + 1)
           << "_Heter(1/m3)\t";

  O_file << "Mn\tNi\tSi" << std::endl;

  // Solving all the equations, each run gives output at time ts+tout
  for (int i = 0; i < runs; ++i)
  {
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    data = (UserData)malloc(sizeof *data);
    loadData(data);
    flag = CVodeSetUserData(cvode_mem, data);
    flag = CVodeInit(cvode_mem, f, T0, y0);
    flag = CVodeSStolerances(cvode_mem, RTOL, ATOL);
    flag = CVBand(cvode_mem, neq, mu, ml);
    flag = CVodeSetMaxNumSteps(cvode_mem, mxsteps);
    time_t tik, tok; /*tik is the time before solving equations, tok is the time after solving equations, tok-tik will
                        be the time used in each run*/
    time(&tik);
    flag = CVode(cvode_mem, tout, y0, &t, CV_NORMAL);
    time(&tok);
    getOutput(y0, radM1, radM2, rhoC); // Obtain the output after each run
    yd = NV_DATA_S(y0);

    // The following few lines write the output into file Output
    O_file << i << "\t" << difftime(tok, tik) << "\t" << ts + tout << "\t" << (ts + tout) * Flux << "\t";
    for (int p = 0; p < numCalcPhase; ++p)
      O_file << RadiusCalc[p] << "\t" << rhoC[p] << "\t";

    O_file << yd[neq - 3] << "\t" << yd[neq - 2] << "\t" << yd[neq - 1];
    O_file << std::endl;

    // Write the particle size distribution of each phase into a seperate file
    printYVector(y0);

    CVodeFree(&cvode_mem);
    ts = ts + tout;             // Calculate the start time of next run
    tout = std::pow(10, i / 9); // Define the solution time for next run
  }
  N_VDestroy_Serial(y0);

  return 0;
}

/*************************************************

This function initialize parameters used in model

*************************************************/
void initParams()
{
  // Irradiatio flux
  Flux = ICond->Flux;

  for (int i = 0; i < numComp; ++i)
    // Thermal diffusion coefficients
    D[i] = IMaterial->D[i];

  // Calculate radiation enhanced diffusion coefficients
  GetRED(D, Flux);

  for (int p = 0; p < numPhase; ++p)
    // Effective atomic radius in precipitate
    aP[p] = std::cbrt((3 * IMaterial->cVol[p]) / (4 * pi));

  for (int p = 0; p < numPhase; ++p)
    for (int c = 0; c < numComp; c++)
      // Square of precipitate composition
      nu[p][c] = sqr(IMaterial->X[p][c]);
}

/*****************************************************************

This function gives the radiation enhanced diffusion coefficients
Described in SI Sec. C

*****************************************************************/
void GetRED(realtype D[numComp], realtype Flux)
{
  realtype Eta, Gs, Cvr;

  if (Flux > IProp->Rflux)
  {
    // When flux is higher than reference flux, use p-scaling (Eq. SI-18) to calculate gs
    Eta = 16 * pi * IProp->rv * IProp->DCB * IProp->SigmaDpa * IProp->Rflux / IMaterial->aVol / IProp->DV /
          sqr(IProp->DDP);
    Gs = 2.0 / Eta * (std::sqrt(1 + Eta) - 1.0) * std::pow((IProp->Rflux / Flux), IProp->p_factor);
  }
  else
  {
    // When flux is lower than reference flux, use Eq. SI-19 and SI-20 to calculate gs
    Eta = 16 * pi * IProp->rv * IProp->DCB * IProp->SigmaDpa * Flux / IMaterial->aVol / IProp->DV / sqr(IProp->DDP);
    Gs = 2.0 / Eta * (std::sqrt(1 + Eta) - 1.0);
  }

  // Calculate vacancy concentration under irradiation, Eq. SI-17
  Cvr = IProp->DCB * Flux * IProp->SigmaDpa * Gs / (IProp->DV * IProp->DDP);

  for (int i = 0; i < numComp; ++i)
    // Radiation enhanced diffusion coefficients, Eq. SI-16
    D[i] = D[i] + IProp->DV * Cvr * D[i] / IMaterial->DFe;
}

/*****************************************************************

This function loads data that is as a function of cluster size
All values are calcuate via next few functions.

*****************************************************************/
void loadData(UserData &data)
{
  realtype size[numClass], radClust[numPhase][numClass], beta[numPhase][numClass], delG[numPhase][numClass];
  getSize(size);
  getRadClust(size, radClust);
  getbeta(size, beta);
  getDelG(size, radClust, delG);

  for (int i = 0; i < numClass; ++i)
    data->size[i] = size[i];

  for (int p = 0; p < numPhase; ++p)
    for (int i = 0; i < numClass; ++i)
    {
      data->radClust[p][i] = radClust[p][i];
      data->beta[p][i] = beta[p][i];
      data->delG[p][i] = delG[p][i];
    }
}

/*****************************************************************

This function calculates part of the absorption rate

*****************************************************************/
void getbeta(realtype size[numClass], realtype beta[numPhase][numClass])
{
  for (int p = 0; p < numPhase; ++p)
    for (int i = 0; i < numClass; ++i)
      beta[p][i] = (4 * pi * aP[p] * std::cbrt(size[i])) / IMaterial->aVol;
}

/*****************************************************************

This function calculates number of atoms (size) in each cluster

*****************************************************************/

void getSize(realtype size[numClass])
{
  for (int i = 0; i < numClass; ++i)
    size[i] = double(i + 1);
}

/*****************************************************************

This function calculates radius of each cluster

*****************************************************************/

void getRadClust(realtype size[numClass], realtype radClust[numPhase][numClass])
{

  for (int p = 0; p < numPhase; ++p)
    for (int i = 0; i < numClass; ++i)
      radClust[p][i] = std::cbrt(3.0 * IMaterial->cVol[p] * size[i] / (4.0 * pi));
}

/*********************************************************************************

This function calculates difference of interfacial energy between adjacent clusters

**********************************************************************************/

void getDelG(realtype size[numClass], realtype radClust[numPhase][numClass], realtype delG[numPhase][numClass])
{
  realtype delta;

  for (int p = 0; p < numPhase; ++p)
    for (int i = 1; i < numClass; ++i)
    {
      delta = -((IMaterial->sig[p] * 4.0 * pi * sqr(radClust[p][i - 1])) -
                (IMaterial->sig[p] * 4.0 * pi * sqr(radClust[p][i])));
      delG[p][i] = std::exp(delta / (kb * ICond->Temp));
    }
}

/*****************************************************************

This function gives the inital values of number density of each cluster
1E-30 for clusters >= 2 atoms
monomer concentration is calculated based on Eq. SI-15

*****************************************************************/

void getInitVals(realtype y0[neq])
{
  realtype base = 1E-30;
  for (int i = 0; i < neq; ++i)
    y0[i] = base;

  for (int c = 0; c < numComp; c++)
    // Concentration of solute in matrix
    y0[neq - numComp + c] = IMaterial->C0[c];

  for (int p = 0; p < numPhase; ++p)
  {
    solProd[p] = 1;
    for (int c = 0; c < numComp; c++)
      solProd[p] *= std::pow(IMaterial->C0[c], IMaterial->X[p][c]); // SI-14

    // Effective monomer concentration, based on Eq. SI-15
    y0[p * numClass] = solProd[p];
  }
}

/**************************************************************************

This function calculates flux (Jn->n+1) between adjacent clusters, Eq. SI-2

****************************************************************************/

void getFlux(UserData data, N_Vector y, realtype J[numCalcPhase][numClass + 1])
{
  realtype *yd;
  yd = NV_DATA_S(y);
  realtype solP[numCalcPhase], wp[numComp], sumwp, wpEff;

  // pref is used to refer to the real phase for both homo and heter nucleated phases,
  // pref=0 is T3 phase, pref=1 is T6 phase
  int pref;
  for (int p = 0; p < numPhase; ++p)
  {
    solP[p] = 1;
    for (int c = 0; c < numComp; c++)
      solP[p] *= std::pow(yd[neq - numComp + c], IMaterial->X[p][c]); // Calculate solute product SI-14
  }

  // solute product for heterogeneous nucleation phase, used to calcuate effective monomer concentration
  for (int p = numPhase; p < numCalcPhase; ++p)
    solP[p] = 1E-30;

  for (int p = 0; p < numCalcPhase; ++p)
  {
    pref = p % numPhase;
    J[p][0] = 0.0;
    // Effective monomer concentration
    yd[p * numClass] = solP[p];
    sumwp = 0;
    for (int c = 0; c < numComp; c++)
    {
      wp[c] = yd[neq - numComp + c] * D[c] * data->beta[pref][0];
      sumwp += nu[pref][c] / wp[c];
    }

    // Absorption rate for momoner to dimer
    wpEff = 1.0 / sumwp;

    // flux from monomer to dimer
    J[p][1] = wpEff * (yd[p * numClass] / numPhase -
                       (IMaterial->solPBar[pref] / solP[pref]) * data->delG[pref][1] * yd[p * numClass + 1]);

    for (int i = 2; i < numClass; ++i)
    {
      sumwp = 0;
      for (int c = 0; c < numComp; c++)
      {
        wp[c] = yd[neq - numComp + c] * D[c] * data->beta[pref][i - 1];
        sumwp += (nu[pref][c] / wp[c]);
      }

      wpEff = 1.0 / sumwp;

      // flux from size i to size i+1
      J[p][i] = wpEff * (yd[p * numClass + i - 1] -
                         (IMaterial->solPBar[pref] / solP[pref]) * data->delG[pref][i] * yd[p * numClass + i]);
    }
    J[p][numClass] = 0.0;
  }
}

/**********************************************************************************

This function calculates Eq. SI-1, the change of concentration of each cluster size

***********************************************************************************/

int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype *yd, *ydotd;
  yd = NV_DATA_S(y);
  ydotd = NV_DATA_S(ydot);
  UserData data;

  realtype size[numClass], radClust[numPhase][numClass], beta[numPhase][numClass], delG[numPhase][numClass],
      J[numCalcPhase][numClass + 1];
  int pref;
  realtype solP[numPhase];
  realtype GR;

  data = (UserData)user_data;

  getFlux(data, y, J);

  realtype sumNdot[numCalcPhase], Cdot;
  for (int p = 0; p < numCalcPhase; ++p)
  {
    sumNdot[p] = 0.0;
    ydotd[p * numClass] = 0.0;
    for (int i = 1; i < numClass; ++i)
    {
      ydotd[p * numClass + i] = J[p][i] - J[p][i + 1];       // Eq. 1 in Sec. 2.1 without Rhet term
      sumNdot[p] += ydotd[p * numClass + i] * data->size[i]; // Monomer consumed
    }
  }

  for (int p = numPhase; p < numCalcPhase; ++p)
  {
    pref = p % numPhase;
    solP[pref] = 1;
    for (int c = 0; c < numComp; c++)
      solP[pref] *= std::pow(yd[neq - numComp + c], IMaterial->X[pref][c]); // solute product
  }

  // The next three lines is the generation of clusters in cascade damage

  // Generation rate of clusters, Eq. 5 in Sec. 2.2
  GR = IProp->Alpha * Flux * IProp->ccs * solP[IProp->HGPhase] / IProp->RsolP;

  // Add Rhet term to Eq. (1) in Sec. 2.1
  ydotd[(IProp->HGPhase + numPhase) * numClass + IProp->HGSize - 1] += GR;

  // Calculate the monomers consumed in heterogeneous nucleation
  sumNdot[(IProp->HGPhase + numPhase)] += GR * data->size[IProp->HGSize - 1];

  for (int c = 0; c < numComp; c++)
  {
    Cdot = 0.0;
    for (int p = 0; p < numCalcPhase; ++p)
    {
      pref = p % numPhase;
      Cdot += IMaterial->X[pref][c] * sumNdot[p];
    }
    ydotd[neq - numComp + c] = -Cdot; // change of solute/monomer in matrix
  }

  return 0;
}

/*****************************************************************

This function calculates mean cluster radius and cluster density

*****************************************************************/

void getOutput(N_Vector y, realtype radM1[numCalcPhase], realtype radM2[numCalcPhase], realtype rhoC[numCalcPhase])
{
  realtype *yd, size[numClass], radClust[numPhase][numClass], numC, numCxSize1, numCxSize2;
  yd = NV_DATA_S(y);
  int pref;
  getSize(size);
  getRadClust(size, radClust);
  for (int p = 0; p < numCalcPhase; ++p)
  {
    pref = p % numPhase;

    // Concertration of clusters in unit per lattice site
    numC = 0;
    numCxSize1 = 0;
    numCxSize2 = 0;

    // Mean radius of clusters calculated via sum(r_i*N_i)/sum(N_i)
    // refer to readme for details regarding these two different ways of calculating mean radius
    radM1[p] = 0;

    // Mean radius of clusters calculated via equilivant radius of a cluster with size sum(n_i*N_i)/sum(N_i)
    radM2[p] = 0;
    rhoC[p] = 0;

    for (int i = CutoffSize; i < numClass; ++i)
    {
      // Count clusters larger than CutoffSize for output
      numC += yd[p * numClass + i];
      numCxSize1 += yd[p * numClass + i] * radClust[pref][i];
      numCxSize2 += yd[p * numClass + i] * size[i];
    }
    radM1[p] = numCxSize1 / numC;
    radM2[p] = std::cbrt((numCxSize2 / numC * IMaterial->cVol[pref]) / ((4. / 3.) * pi));

    // Calculate number density of clusters in unit of per m^3
    rhoC[p] = numC / IMaterial->aVol;
  }
}

/**********************************************************************************************

This function prints cluster size distribution in the file Profile for the final solution time.

***********************************************************************************************/
void printYVector(N_Vector y)
{
  realtype *yd, size[numClass], radClust[numPhase][numClass];
  std::ofstream P_file;
  yd = NV_DATA_S(y);
  int pref;

  for (int p = 0; p < numCalcPhase; ++p)
  {
    P_file.open(("Profile_" + int_to_string(p)).c_str());
    P_file << "cluster size (# atoms)\tcluster radius (m)\tcluster density (1/m3)" << std::endl;

    getSize(size);
    getRadClust(size, radClust);

    pref = p % numPhase;
    for (int i = 0; i < numClass; ++i)
      P_file << size[i] << "\t" << radClust[pref][i] << "\t" << yd[p * numClass + i] / IMaterial->aVol << std::endl;

    P_file.close();
  }
}

std::string int_to_string(int i)
{
  std::stringstream ss;
  ss << i;
  return ss.str();
}
