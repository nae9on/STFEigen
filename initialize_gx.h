/*
 * initialize_gx.h
 *
 *  Created on: Oct 24, 2018
 *      Author: akadar
 */

#ifndef INITIALIZE_GX_H_
#define INITIALIZE_GX_H_

#include <iostream>
#include <math.h>

void gx_generator(double * gx)
{
	// Initialize x = 0:deltaX:L;
	double * x = new double[global_N+1] ();
	for (unsigned long int i=0; i<global_N+1; i++)
		{
		x[i] = i*global_deltaX;
		for (unsigned long int k=0; k<global_N+1; k++)
		{
			int m = k - (global_N/2 + 1);
					if(m < 0)
					{
						gx[i+k*(global_N+1)] = sqrt(2/global_L)*sin(2*M_PI*m*x[i]/global_L);
					}
					else if(m > 0)
					{
						gx[i+k*(global_N+1)] = sqrt(2/global_L)*cos(2*M_PI*m*x[i]/global_L);
					}
					else
					{
						gx[i+k*(global_N+1)] = sqrt(2/global_L);
					}

		}
		}

	delete[] x;
	std::cout<<"\n\nFinished initializing gx\n\n";
}



#endif /* INITIALIZE_GX_H_ */
