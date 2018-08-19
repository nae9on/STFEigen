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

std::stack<clock_t> tictoc_stack;

void tic() {
    tictoc_stack.push(clock());
}

void toc() {
    std::cout << "Time elapsed: "
              << ((double)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC
              << std::endl;
    tictoc_stack.pop();
}

typedef Eigen::Triplet<double> tripleData;

int main(int argc, char** argv) {
	assert(argc == 1);

	// Declaration and initialization of h
	Eigen::VectorXd hLU(h_size);
	hLU.setOnes();
	// hLU.initialize_h(deltaX);
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



	int counter01 = 0;
	for (double time = deltaT; time <= endTime; time = time + deltaT) {

		// Update A
		updateA(coefficients, hLU, h_size);
		ALU.setFromTriplets(coefficients.begin(), coefficients.end());
		// displayFullMatrix(ALU);

		// Update b
		updateRHS(bLU, hLU);
		// displayVector(bLU);

		// Factorize A
		solverLU.factorize(ALU);

		// Update h
		hLU = solverLU.solve(bLU);
		// displayVector(hLU);
		counter01 = counter01 + 1;
		if(counter01 == seN)
		{
			// Update A
			tic();
			updateA(coefficients, hLU, h_size);
			ALU.setFromTriplets(coefficients.begin(), coefficients.end());
			toc();
			// displayFullMatrix(ALU);

			// Update b
			tic();
			updateRHS(bLU, hLU);
			toc();
			// displayVector(bLU);

			// Factorize A
			tic();
			solverLU.factorize(ALU);
			toc();

			// Update h
			tic();
			hLU = solverLU.solve(bLU);
			toc();
			write_h_toFile(hLU, time);

			std::cout << "\nDone time = " << time << "\n";
			counter01 = 0;
			displayVector(hLU);
			// displayVector(bLU);
			// displayFullMatrix(ALU);
			tic();
			for(unsigned int i = 0; i <= sizeof(hLU); i++)
			{
				while (hLU[i] <= 0.1)
				{
					write_h_toFile(hLU, time);
					displayVector(hLU);
					return 0;
				}
			}
			toc();
		}



	}

	return 0;
}


