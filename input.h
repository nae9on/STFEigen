/*
 * input.h
 *
 *  Created on: 20 Jul 2018
 *      Author: msshah
 */

#ifndef INPUT_H_
#define INPUT_H_

const double L_flat = 15;        // length of the flat portion of the film
const unsigned long int N_init = pow(2, 8);  // number of grid points on the flat film
const double deltaX = L_flat / N_init;  // grid/mesh size

// Physical properties
const double Tmp = 0.001;     // dimensionless noise strength (= 0, for deterministic)

// Curvature
const double kappa = 0.0;      // dimensionless curvature (= 0 for flat films)
const double L_curv = 0;
const double L = L_curv + L_flat;   // total length of the film (curved+flat)
const unsigned long int N = floor(L / deltaX);   // adjusted number of grid points;
const unsigned long int h_size = N + 5;

// Determine time step
const double c = 2.75;           // exponent used in deciding deltaT = deltaX^c;
const double deltaT = pow(deltaX, c);     // time step size
const double endTime = 122;          // end time of a realization
const unsigned long int seN = 20;              // save every these many time steps
const unsigned int N_Reals = 50;          // number of realizations

//Discretization parameters
const double p1 = deltaT/(deltaX*deltaX);
const double p2 = deltaT/std::pow(deltaX,4);
const double p3 = 1/(deltaX)*sqrt(2*deltaT*Tmp);

// output
static unsigned int numFiles = 0;
const std::string folderName = "realization";
const std::string baseFileName = "./output/data";

#endif /* INPUT_H_ */
