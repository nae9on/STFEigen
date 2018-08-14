/*
 * main.cpp
 *
 *  Created on: Aug 14, 2018
 *      Author: nae9on
 */

#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <math.h>
#include "input.h"
#include "output.h"
#include "update.h"

typedef Eigen::Triplet<double> tripleData;

int main(int argc, char** argv) {
	assert(argc == 1);

	// Declaration and initialization of h
	Eigen::VectorXd hLU(h_size);
	hLU.setOnes();

	// Declaration and initialization of A
	// Declares a column-major sparse matrix type of double
	Eigen::SparseMatrix<double, Eigen::ColMajor> ALU(h_size, h_size);
	// Reserve NNZ for A
	ALU.reserve(Eigen::VectorXi::Constant(h_size, 5));
	std::vector<tripleData> coefficients;
	updateA(coefficients, hLU, h_size);
	ALU.setFromTriplets(coefficients.begin(), coefficients.end());
	//displayFullMatrix(ALU);

	// Declaration and initialization of b
	Eigen::VectorXd bLU(h_size);
	bLU.setZero();
	//displayVector(bLU);

	// Initialize solver
	Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor> > solverLU;

	// Analyze pattern
	solverLU.analyzePattern(ALU);

	std::cout << "\n\nDone initialization \n\n";

	unsigned long int itr = 0;

	for (double time = 0; itr <= 2; time = time + deltaT, itr++) {

		// Update A
		updateA(coefficients, hLU, h_size);
		ALU.setFromTriplets(coefficients.begin(), coefficients.end());
		//displayFullMatrix(ALU);

		// Update b
		updateRHS(bLU, hLU);
		//displayVector(bLU);

		// Factorize A
		solverLU.factorize(ALU);

		// Update h
		hLU = solverLU.solve(bLU);

		if (itr % 10 == 0)
		{
			write_h_toFile(hLU, time);
			std::cout << "Done time = " << time << "\n";
		}

	}
}
