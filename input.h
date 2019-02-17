/*
 * input.h
 *
 *  Created on: 20 Jul 2018
 *      Author: msshah
 */

#ifndef INPUT_H_
#define INPUT_H_

const double global_L_flat = 15;        // length of the flat portion of the film
const unsigned long int global_N_init = pow(2, 8);  // number of grid points on the flat film
const double global_deltaX = global_L_flat / global_N_init;  // grid/mesh size

// Physical properties
const double global_Tmp = 0.001;     // dimensionless noise strength (= 0, for deterministic)

// Curvature
const double global_kappa = 0.0;      // dimensionless curvature (= 0 for flat films)
const double global_L_curv = 0;
const double global_L = global_L_curv + global_L_flat;   // total length of the film (curved+flat)
const unsigned long int global_N = floor(global_L / global_deltaX);   // adjusted number of grid points;
const unsigned long int global_h_size = global_N + 5;

// Determine time step
const double global_c = 2.75;           // exponent used in deciding deltaT = deltaX^c;
const double global_deltaT = pow(global_deltaX, global_c);     // time step size
const double global_endTime = 122;          // end time of a realization
const unsigned long int global_seN = 100;              // save every these many time steps
const unsigned int global_N_Reals = 2;          // number of realizations

//Discretization parameters
const double global_p1 = global_deltaT/(global_deltaX*global_deltaX);
const double global_p2 = global_deltaT/std::pow(global_deltaX,4);
const double global_p3 = 1/(global_deltaX)*sqrt(2*global_deltaT*global_Tmp);

// output
static unsigned int global_numFiles = 0;
const std::string global_folderName = "realization";
const std::string global_baseFileName = "./output/data";

#endif /* INPUT_H_ */
