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
#include "input.h"
#include "output.h"
#include "update.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

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
	Eigen::VectorXd hLU(h_size);
	hLU.setOnes();
	// displayVector(hLU);

	// Declaration and initialization of A
	// Declares a column-major sparse matrix type of double
	Eigen::SparseMatrix<double, Eigen::ColMajor> ALU(h_size, h_size);
	// Reserve NNZ for A
	ALU.reserve(Eigen::VectorXi::Constant(h_size, 5));
	std::vector<tripleData> coefficients;
	updateA(coefficients, hLU, h_size);
	ALU.setFromTriplets(coefficients.begin(), coefficients.end());
	// displayFullMatrix(ALU);

	/*// the following is to teach how to read the triple data coefficients
	 for(std::vector<tripleData>::iterator it = coefficients.begin(); it != coefficients.end(); ++it) {
	 tripleData t = *it;
	 std::cout<<t.col()<<" "<<t.row()<<" "<<t.value()<<"\n";
	 }
	 */

	// Declaration and initialization of b
	Eigen::VectorXd bLU(h_size);
	bLU.setZero();
	//displayVector(bLU);

	// Initialize solver
	Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor> > solverLU;

	// Analyze pattern
	solverLU.analyzePattern(ALU);

	unsigned long int counter01 = 0;
	bool flag = 1;

	// Initialize random number generator
	typedef boost::mt19937 RNGType; //boost will use mersenne twister generator.
	RNGType rng(time(0)); // Defining the generator and intializing with time as seed.
	boost::normal_distribution<> nd(0.0, 1.0); // Defining the distribution
	// variate_generator combines a generator with a distribution.
	boost::variate_generator<RNGType, boost::normal_distribution<> > randNum(rng, nd);

	// Random number can be obtained by calling randNum()
	std::cout<<"Random no = "<< randNum();

	for (double time = deltaT; time <= endTime; time = time + deltaT) {

		// Update A using h
		updateA(coefficients, hLU, h_size);
		ALU.setFromTriplets(coefficients.begin(), coefficients.end());
		// displayFullMatrix(ALU);

		// Update b using h
		updateRHS(bLU, hLU);
		// displayVector(bLU);

		// Factorize A
		solverLU.factorize(ALU);
		//solverLU.compute(ALU);

		// Update h
		hLU = solverLU.solve(bLU);
		// displayVector(hLU);

		counter01 = counter01 + 1;
		if (counter01 == seN) {
			write_h_toFile(hLU, time);
			std::cout << "\nData written to the file at time = " << time
					<< "\n";
			counter01 = 0;

			// Check if film height is below the solid substrate
			for (unsigned int i = 0; i <= sizeof(hLU); i++) {
				if (hLU[i] <= 0.1) {
					flag = 0;
					break;
				}
			}
		}

		if (flag == 0)
			break;
	}

	toc();

	return 0;
}
