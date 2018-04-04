/*
 * schroedinger.cpp
 *
 *  Created on: 02.06.2017
 *      Author: marcel
 */

#define _USE_MATH_DEFINES

#include <iostream>
#include <stdio.h>
#include <time.h>
#include <random>
#include <complex>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <sstream>
#include <utility>
#include <math.h>
#include <random>
#include <string>
#include <fstream>
#include <Eigen/Dense>

using namespace std;

const double dxi = 0.1;
const double dtau = 0.01;
const double sigma = 1.;
const double xi0 = -5.;
const double k0 = 5.5;
const int dim = 201;
const double b = 1.;
const complex<double> i(0, 1);
const double normierung = pow(2 * M_PI * sigma, -0.25);

double theta(double x) {
	if (x >= 0)
		return 1.;
	else
		return 0;
}

Eigen::MatrixXcd hamilton(double V) {
	Eigen::MatrixXcd H(dim, dim);
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			if (i == j) {
				H(i, j) = 2. / (dxi * dxi)
						+ V * theta((b / 2) - abs(-10 + j * dxi));
			}
			if (i == j + 1 || i == j - 1) {
				H(i, j) = -1. / (dxi * dxi);
			}
		}
	}
	return H;
}

Eigen::MatrixXcd SH(double V) {
	Eigen::MatrixXcd M(dim, dim);
	Eigen::MatrixXcd H = hamilton(V);

	for (int n = 0; n < dim; n++) {
		for (int m = 0; m < dim; m++) {
			M(n, m) = H(n, m) * 0.5 * i + dtau;
		}
	}

	return (Eigen::MatrixXcd::Identity(dim,dim) + M).inverse() * (Eigen::MatrixXcd::Identity(dim,dim) - M);
}

Eigen::VectorXcd PSI() {

	Eigen::VectorXcd Phi(dim);
	for (int m = 0; m < dim; m++) {
		Phi(m) = normierung
				* exp(
						-(-10 + m * dxi - xi0) * (-10 + m * dxi - xi0)
								/ (4 * sigma)) * exp(i * k0 * (-10 + m * dxi));
	}
	return Phi;
}

double T(Eigen::VectorXcd psi) {
	double sum = 0;
	for (int j = 100; j < dim; j++) {
		sum += psi.real()(j);
	}
	return sum * dxi;
}

void Crank_Nicolson(double V, string dateiname) {
	// ===============================
	ofstream crank;
	crank.open("psi_" + dateiname + ".txt");
	crank.precision(10);
	// ===============================
	// ===============================
	ofstream data;
	data.open("transmission" + dateiname + ".txt");
	data.precision(10);
	// ===============================
	Eigen::VectorXcd psi = PSI();
	Eigen::MatrixXcd U = SH(V);
	Eigen::VectorXcd psiquad(dim);

	for (int j = 0; j < dim; j++) {
		psiquad(j) = psi.conjugate()(j) * psi(j);
		crank << -10 + j * dxi << "\t" << psiquad.real()(j) << "\n";
	}
	crank << "\n";

	for (double t = 0; t <= 1; t += dtau) {
		psi = U * psi;
		for (int j = 0; j < dim; j++) {
			psiquad(j) = psi.conjugate()(j) * psi(j);
			crank << -10 + j * dxi << "\t" << psiquad.real()(j) << "\n";
		}
		crank << "\n";
		data << t << "\t" << T(psiquad) << "\n";

	}

	crank.close();
	data.close();
}

void Potential_Plot(double V_0) {

	ofstream Potential;
	Potential.open("Potential_Plot.txt");
	Potential.precision(10);

	for(double xi = -10.; xi <= 10; xi += dxi) {
		Potential << xi << "\t" << V_0 * theta((b / 2) - abs(xi)) << "\n";
	}

	Potential.close();
}

int main() {
	Crank_Nicolson(0, "V=0");
	Crank_Nicolson(10, "V=10");
	Crank_Nicolson(30, "V=30");
	Crank_Nicolson(50, "V=50");

	Potential_Plot(10);
	return 0;
}
