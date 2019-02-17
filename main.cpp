/*
 * main.cpp
 *
 *  Created on: 20 Jul 2018
 *      Author: msshah
 */

#include <Eigen/Sparse>
#include <iostream>
#include <math.h>
#include <stack>
#include <ctime>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

// My headers
#include "input.h"
#include "output.h"
#include "update.h"
#include "initialize_gx.h"

std::stack<clock_t> tictoc_stack;

void tic() {
	tictoc_stack.push(clock());
}

void toc() {
	std::cout << "Time elapsed: "
			<< ((double) (clock() - tictoc_stack.top())) / CLOCKS_PER_SEC
			<< " sec" << std::endl;
	tictoc_stack.pop();
}

typedef Eigen::Triplet<double> tripleData;

int main(int argc, char** argv) {
	assert(argc == 1);

	tic();

	// Declaration and initialization of h
	Eigen::VectorXd hLU(global_h_size);
	hLU.setOnes();
	// displayVector(hLU);

	// Declaration and initialization of A
	// Declares a column-major sparse matrix type of double
	Eigen::SparseMatrix<double, Eigen::ColMajor> ALU(global_h_size, global_h_size);
	// Reserve NNZ for A
	ALU.reserve(Eigen::VectorXi::Constant(global_h_size, 5));
	std::vector<tripleData> coefficients;
	updateA(coefficients, hLU, global_h_size);
	ALU.setFromTriplets(coefficients.begin(), coefficients.end());
	// displayFullMatrix(ALU);

	/*// the following is to teach how to read the triple data coefficients
	 for(std::vector<tripleData>::iterator it = coefficients.begin(); it != coefficients.end(); ++it) {
	 tripleData t = *it;
	 std::cout<<t.col()<<" "<<t.row()<<" "<<t.value()<<"\n";
	 }
	 */

	// Declaration and initialization of b
	Eigen::VectorXd bLU(global_h_size);
	bLU.setZero();
	//displayVector(bLU);

	// Initialize solver
	Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor> > solverLU;

	// Analyze pattern
	solverLU.analyzePattern(ALU);

	// Initialize gx
	double * gx = new double[(global_N+1)*(global_N+1)] ();
	double * gx_f = new double[(global_N+1)*(global_N+1)] ();
	double * g_noise = new double[(global_N+1)] ();
	gx_generator(gx);

	unsigned long int counter01 = 0;
	bool flag = 0;

	// Initialize random number generator
	typedef boost::mt19937 RNGType; //boost will use mersenne twister generator.
	RNGType rng(time(0)); // Defining the generator and intializing with time as seed.
	boost::normal_distribution<> nd(0.0, 1.0); // Defining the distribution
	// variate_generator combines a generator with a distribution.
	boost::variate_generator<RNGType, boost::normal_distribution<> > randNum(rng, nd);

	// Random number can be obtained by calling randNum()
	std::cout<<"Random no = "<< randNum();

	for (double time = global_deltaT; time <= global_endTime; time = time + global_deltaT) {

		// Update A using h
		updateA(coefficients, hLU, global_h_size);
		ALU.setFromTriplets(coefficients.begin(), coefficients.end());
		// displayFullMatrix(ALU);

		// Update the gx_f matrix
		for (unsigned long int k = 0; k < (global_N + 1); k++) {
			double R_noise = randNum();
			for (unsigned long int i = 0; i < (global_N + 1); i++) {
				gx_f[i + (global_N + 1) * k] = gx[i + (global_N + 1) * k] * R_noise;

			}
		}
		// Update the g_noise vector
		for (unsigned long int i = 0; i < (global_N + 1); i++) {
			g_noise[i] = 0;
			for (unsigned long int k = 0; k < (global_N + 1); k++) {
				g_noise[i] += gx_f[i + (global_N + 1) * k];

			}
		}

		// Update b using h
		updateRHS(bLU, hLU, g_noise);
		// displayVector(bLU);

		// Factorize A
		solverLU.factorize(ALU);
		//solverLU.compute(ALU);

		// Update h
		hLU = solverLU.solve(bLU);
		// displayVector(hLU);

		// Check if film height is below the solid substrate
		for (unsigned int i = 0; i <= sizeof(hLU); i++) {
			if (hLU[i] <= 0.1) {
				flag = 1; // set termination flag to true
				break;
			}
		}

		counter01 = counter01 + 1;
		if (flag==1 || counter01%global_seN==0) {
			write_h_toFile(hLU, time);
			std::cout << "\nData written to the file at time = " << time
					<< "\n";
		}

		if (flag == 1) {
			std::cout << "\nTerminating time loop at time = "<<time;
			std::cout << "\nTotal no of time iterations = "<<counter01;
			std::cout << "\nTotal no of files written to disk = "<<floor(counter01/global_seN)+1<<"\n";
			break;
		}

	}

	toc();

	delete[] gx;
	delete[] gx_f;
	delete[] g_noise;
	return 0;
}
