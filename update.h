/*
 * update.h
 *
 *  Created on: 20 Jul 2018
 *      Author: msshah
 */

#ifndef UPDATE_H_
#define UPDATE_H_

#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

typedef Eigen::Triplet<double> tripleData;

void insertCoefficient(int id, int i, int j, double w,
		std::vector<tripleData>& coeffs, Eigen::VectorXd& b,
		const Eigen::VectorXd& boundary) {
	int n = int(boundary.size());
	int id1 = i + j * n;

	if (i == -1 || i == n)
		b(id) -= w * boundary(j); // constrained coefficient
	else if (j == -1 || j == n)
		b(id) -= w * boundary(i); // constrained coefficient
	else
		coeffs.push_back(tripleData(id, id1, w));         // unknown coefficient
}

// Fill right hand side
void updateRHS(Eigen::VectorXd& b, Eigen::VectorXd h, boost::variate_generator<boost::mt19937, boost::normal_distribution<> > randNum) {
	for (unsigned int i = 2; i < b.size() - 2; i++) {
		double h1 = 2*h(i+1)*h(i+1)*h(i)*h(i)/(h(i+1)+h(i));
		double h2 = 2*h(i-1)*h(i-1)*h(i)*h(i)/(h(i-1)+h(i));
		b(i) = h(i) -  global_p1 * 0.5 * (( 1/h(i + 1) + 1/ h(i) ) * ( h(i+1) - h(i) ) - ( 1/h(i) + 1/h(i-1) ) * ( h(i) - h(i-1) ))
				+ global_p3*(sqrt(h1)*randNum() - sqrt(h2)*randNum());
	}
}

void initialize_h(double deltaX, Eigen::VectorXd h) {
	for (unsigned int i = 0; i < h.size() ; i++) {
		h(i) = 1 + 0.001*sin(6*(deltaX*i - 7)*(deltaX*i - 7));;
	}

}

// Inefficient
void fillPentaDiagonal_Ineff(Eigen::SparseMatrix<double, Eigen::ColMajor> *A) {
	for (int i = 0; i < A->rows(); i++) {
		for (int j = 0; j < A->cols(); j++) {
			A->insert(i, j) = i + j;
		}
	}

}

// Efficient
void fillPentaDiagonal(std::vector<tripleData>& coefficients, Eigen::VectorXd h, unsigned long int N) {

	// Fill left boundary
	coefficients.push_back(tripleData(0, 0, 1));
	coefficients.push_back(tripleData(0, N - 4, -1));
	coefficients.push_back(tripleData(1, 1, 1));
	coefficients.push_back(tripleData(1, N - 3, -1));

	// Fill the core entries
	for (unsigned long int i = 2; i <= N-3; i++) {
		double h1 = 2*h(i+1)*h(i+1)*h(i)*h(i)/(h(i+1)+h(i));
		double h2 = 2*h(i-1)*h(i-1)*h(i)*h(i)/(h(i-1)+h(i));

		coefficients.push_back(tripleData(i, i - 2, h2*global_p2));
		coefficients.push_back(tripleData(i, i - 1, -global_p2*(h1+3*h2)));
		coefficients.push_back(tripleData(i, i, 3*global_p2*(h1+h2)+1));
		coefficients.push_back(tripleData(i, i + 1, -global_p2*(3*h1+h2)));
		coefficients.push_back(tripleData(i, i + 2, h1*global_p2));
	}

	// Fill right boundary
	coefficients.push_back(tripleData(N - 2, 2, -1));
	coefficients.push_back(tripleData(N - 2, N - 2, 1));
	coefficients.push_back(tripleData(N - 1, 3, -1));
	coefficients.push_back(tripleData(N - 1, N - 1, 1));
}

// Efficient
void updateA(std::vector<tripleData>& coefficients, Eigen::VectorXd h, unsigned long int N) {

	coefficients.clear();

	// Fill left boundary
	coefficients.push_back(tripleData(0, 0, 1));
	coefficients.push_back(tripleData(0, N - 4, -1));
	coefficients.push_back(tripleData(1, 1, 1));
	coefficients.push_back(tripleData(1, N - 3, -1));

	// Fill the core entries
	for (unsigned long int i = 2; i <= N-3; i++) {
		double h1 = 2*h(i+1)*h(i+1)*h(i)*h(i)/(h(i+1)+h(i));
		double h2 = 2*h(i-1)*h(i-1)*h(i)*h(i)/(h(i-1)+h(i));

		coefficients.push_back(tripleData(i, i - 2, h2*global_p2));
		coefficients.push_back(tripleData(i, i - 1, -global_p2*(h1+3*h2)));
		coefficients.push_back(tripleData(i, i, 3*global_p2*(h1+h2)+1));
		coefficients.push_back(tripleData(i, i + 1, -global_p2*(3*h1+h2)));
		coefficients.push_back(tripleData(i, i + 2, h1*global_p2));
	}

	// Fill right boundary
	coefficients.push_back(tripleData(N - 2, 2, -1));
	coefficients.push_back(tripleData(N - 2, N - 2, 1));
	coefficients.push_back(tripleData(N - 1, 3, -1));
	coefficients.push_back(tripleData(N - 1, N - 1, 1));
}


#endif /* UPDATE_H_ */
