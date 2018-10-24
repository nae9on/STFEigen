/*
 * initialize_gx.h
 *
 *  Created on: Oct 24, 2018
 *      Author: akadar
 */

#ifndef INITIALIZE_GX_H_
#define INITIALIZE_GX_H_

#include <iostream>

void gx_generator(double * gx)
{
	// Initialize x = 0:deltaX:L;
	double * x = new double[global_N+1] ();
	for (unsigned long int i=0; i<global_N+1; i++) x[i] = i*global_deltaX;

	delete[] x;
	std::cout<<"\n\nFinished initializing gx\n\n";
}



#endif /* INITIALIZE_GX_H_ */
