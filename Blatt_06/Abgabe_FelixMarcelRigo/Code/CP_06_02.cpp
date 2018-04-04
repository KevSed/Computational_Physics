//============================================================================
// Name        : Poisson-Gleichung.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

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
#include <string>
#include <fstream>
#include <Eigen/Dense>

using namespace std;

const double Delta = 0.05;
const double L = 1;
const int N = L / Delta;
const double eps = 10e-5;

double diff_phi(double phi_1, double phi_2) {
	return (phi_1 - phi_2) / 2;
}

void plot_E(Eigen::MatrixXd phi) {

	ofstream E_Vektorfeld;
	E_Vektorfeld.open("Vektorfeld.txt");
	E_Vektorfeld.precision(10);

	ofstream abs_E_Vektorfeld;
	abs_E_Vektorfeld.open("Vektorfeld_Betrag.txt");
	abs_E_Vektorfeld.precision(10);

	Eigen::MatrixXd E_x(N + 1, N + 1);
	Eigen::MatrixXd E_y(N + 1, N + 1);
	Eigen::MatrixXd E_ges(N + 1, N + 1);

	for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++) {

			E_x(i, j) = diff_phi(phi(i + 1, j), phi(i - 1, j));
			E_y(i, j) = diff_phi(phi(i, j + 1), phi(i, j - 1));
			E_ges(i, j) = sqrt(E_x(i, j) * E_x(i, j) + E_y(i, j) * E_y(i, j));

		}
	}

	for (int i = 0; i <= N; i++) {
		for (int j = 0; j <= N; j++) {
			if (i == 0 || j == 0 || i == N || j == N)
				abs_E_Vektorfeld << i * Delta << "\t" << j * Delta << "\t" << 0
						<< "\n";
			else {
				E_Vektorfeld << i * Delta << "\t" << j * Delta << "\t"
						<< -E_x(i, j) << "\t" << -E_y(i, j) << "\n";
				abs_E_Vektorfeld << i * Delta << "\t" << j * Delta << "\t"
						<< E_ges(i, j) << "\n";
			}
		}
	}

	E_Vektorfeld.close();
}

void Phi_ana(int n_sum) {

	ofstream Poissongleichung_analytisch;
	Poissongleichung_analytisch.open("Poissongleichung_Analytisch.txt");
	Poissongleichung_analytisch.precision(10);

	Eigen::MatrixXd Phi(N + 1, N + 1);

	for (int i = 0; i <= N; i++) {
		for (int j = 0; j <= N; j++) {
			Phi(i, j) = 0;
			for (int n = 1; n <= n_sum; n++) {
				Phi(i, j) += (2 * (1 - cos(n * M_PI)))
						/ (n * M_PI * sinh(n * M_PI))
						* sin(n * M_PI * i * Delta)
						* sinh(n * M_PI * j * Delta);
			}
		}
	}

	for (int i = 0; i <= N; i++) {
		for (int j = 0; j <= N; j++) {
			Poissongleichung_analytisch << i * Delta << "\t" << j * Delta
					<< "\t" << Phi(i, j) << "\n";
		}
	}
	Poissongleichung_analytisch.close();
}

Eigen::MatrixXd Init_Rho() {
	Eigen::MatrixXd rho(N, N);

	rho(N/2,N/2) = 1.;

	//cout << rho << "\n";

	return rho;
}

Eigen::MatrixXd Init_Phi() {

	Eigen::MatrixXd Phi(N + 1, N + 1);

	for (int i = 0; i <= N; i++) {
		for (int j = 0; j <= N; j++) {
			if (j == N)
				Phi(i, j) = 0;
			else
				Phi(i, j) = 0;
		}
	}

	return Phi;
}

Eigen::MatrixXd Gauss_Seidel(Eigen::MatrixXd Phi, Eigen::MatrixXd rho) {

	for (int j = 1; j < N; j++) {
		for (int l = 1; l < N; l++) {
			Phi(j, l) = 0.25
					* ((Phi(j + 1, l) + Phi(j - 1, l) + Phi(j, l + 1)
							+ Phi(j, l - 1)) + rho(j, l));
		}
	}
	return Phi;
}

void Poisson() {

	ofstream Poissongleichung;
	Poissongleichung.open("Poissongleichung.txt");
	Poissongleichung.precision(10);

	Eigen::MatrixXd Phi = Init_Phi();
	Eigen::MatrixXd Rho = Init_Rho();

	for (int i = 0; i < 5000; i++) {

		Phi = Gauss_Seidel(Phi, Rho);

	}

	for (int i = 0; i <= N; i++) {
		for (int j = 0; j <= N; j++) {
			Poissongleichung << i * Delta << "\t" << j * Delta << "\t"
					<< Phi(i, j) << "\n";
		}
	}

	plot_E(Phi);

	Poissongleichung.close();
}


void Poisson_2() {

	ofstream Poissongleichung;
	Poissongleichung.open("Poissongleichung.txt");
	Poissongleichung.precision(10);

	Eigen::MatrixXd Phi1 = Init_Phi();
	Eigen::MatrixXd Phi2(N+1,N+1);
	Eigen::MatrixXd Rho = Init_Rho();

	int cnt = 0;

	while (cnt < (N+1)*(N+1)) {

		cnt = 0;

		Phi2 = Phi1;
		Phi1 = Gauss_Seidel(Phi1, Rho);

		Phi2 = Phi2 - Phi1;
		for(int i = 0; i <= N; i++) {
			for(int j = 0; j <= N; j++) {
				if(abs(Phi2(i,j)) < eps) cnt += 1;
			}
		}
		cout << cnt << "\n";
	}

	for (int i = 0; i <= N; i++) {
		for (int j = 0; j <= N; j++) {
			Poissongleichung << i * Delta << "\t" << j * Delta << "\t"
					<< Phi1(i, j) << "\n";
		}
	}

	plot_E(Phi1);

	Poissongleichung.close();

}

int main() {

	Poisson_2();
	//Phi_ana(200);

	return 0;
}
