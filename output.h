/*
 * output.h
 *
 *  Created on: 20 Jul 2018
 *      Author: msshah
 */

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <iostream>
#include <string>
#include <fstream>

void displayFullMatrix(Eigen::SparseMatrix<double, Eigen::ColMajor> A) {
	std::cout << "\n\nMatrix :";
	for (int i = 0; i < A.rows(); i++) {
		std::cout << "\n";
		for (int j = 0; j < A.cols(); j++) {
			std::cout << A.coeff(i, j) << " ";
		}
	}

}

void displayVector(Eigen::VectorXd b) {
	std::cout << "\n\nVector :\n";
	for (int i = 0; i < b.rows(); i++) {
		// std::cout << b(i) << "\n";
		printf("%16.16f\n",b(i));

	}

}

void write_h_toFile(Eigen::VectorXd hLU, double time) {

	std::string fileName;
	char numstr[21];
	sprintf(numstr, "%d", ++numFiles);
	fileName.append(baseFileName);
	fileName.append(numstr);
	fileName.append(".txt");

	// Convert string to char *
	char *name = new char[fileName.length()+1];
	std::strcpy(name,fileName.c_str());

	std::ofstream myfile(name, std::ofstream::out | std::ofstream::app);
	if (myfile.is_open()) {
		for (unsigned int i = 0; i < hLU.rows(); i++) {
			myfile << hLU(i) << " ";
		}
		myfile.close();
	} else
		std::cout << "Unable to open file";

	delete[] name;
}



#endif /* OUTPUT_H_ */
