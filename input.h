/*
 * input.h
 *
 *  Created on: 20 Jul 2018
 *      Author: msshah
 */

#ifndef INPUT_H_
#define INPUT_H_

# include <iostream>
# include <fstream>
# include <math.h>

double L_flat = 15;        // length of the flat portion of the film
unsigned long int N_init = pow(2, 8);  // number of grid points on the flat film
double deltaX = L_flat / N_init;  // grid/mesh size

// Physical properties
double Tmp = 0.001;     // dimensionless noise strength (= 0, for deterministic)

// Curvature
double kappa = 0.0;      // dimensionless curvature (= 0 for flat films)
double L_curv = 0;
double L = L_curv + L_flat;   // total length of the film (curved+flat)
unsigned long int N = floor(L / deltaX);   // adjusted number of grid points;
unsigned long int h_size = N + 5;

// Determine time step
double c = 2.75;           // exponent used in deciding deltaT = deltaX^c;
double deltaT = pow(deltaX, c);     // time step size
double endTime = 122;          // end time of a realization
unsigned long int seN = 20;              // save every these many time steps
unsigned int N_Reals = 50;          // number of realizations

//Discretization parameters
double p1 = deltaT/(deltaX*deltaX);
double p2 = deltaT/std::pow(deltaX,4);

// output
static unsigned int numFiles = 0;
static std::string folderName = "realization";
static std::string baseFileName = "./output/data";



#endif /* INPUT_H_ */
