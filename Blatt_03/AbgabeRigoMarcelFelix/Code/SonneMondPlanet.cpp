//============================================================================
// Name        : CP_02_02.cpp
// Author      : Rigo Bause,
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#define _USE_MATH_DEFINES

#include <iostream>
#include <stdio.h>
#include <time.h>
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

using namespace std;

const double h = 0.01;

const double G = 1;
const double m_0 = 1;

vector<double> f(vector<double> y_1, vector<double> y_2, double m_1,
		double m_2) {
	vector<double> y_neu(y_1.size());

	double r1 = sqrt(y_1[0] * y_1[0] + y_1[1] * y_1[1] + y_1[2] * y_1[2]);
	double r2 = sqrt(
			(y_1[0] - y_2[0]) * (y_1[0] - y_2[0])
					+ (y_1[1] - y_2[1]) * (y_1[1] - y_2[1])
					+ (y_1[2] - y_2[2]) * (y_1[2] - y_2[2]));
	for (unsigned int i = 0; i < y_1.size() / 2; i++) {
		y_neu[i] = h * y_1[i + y_1.size() / 2];
		y_neu[i + y_1.size() / 2] = -h * G
				* (m_0 * y_1[i] / (r1 * r1 * r1)
						+ m_2 * (y_1[i] - y_2[i]) / (r2 * r2 * r2));
	}
	return y_neu;
}

void runge_kutta(vector<double> y_planet, vector<double> y_mond, double m_1,
		double m_2) {

	// ===============================
	ofstream data_planet;
	data_planet.open("Planet.txt");
	data_planet.precision(10);

	ofstream data_mond;
	data_mond.open("Mond.txt");
	data_mond.precision(10);

	ofstream data_mond_rel;
	data_mond_rel.open("Mond_Rel.txt");
	data_mond_rel.precision(10);
	// ===============================

	vector<double> k1_planet(y_planet.size()), k2_planet(y_planet.size()),
			k3_planet(y_planet.size()), k4_planet(y_planet.size()),
			y_step_planet(y_planet.size());
	vector<double> k1_mond(y_mond.size()), k2_mond(y_mond.size()), k3_mond(
			y_mond.size()), k4_mond(y_mond.size()), y_step_mond(y_mond.size());

	for (int i = 0; i <= 4800; i++) {

		k1_planet = f(y_planet, y_mond, m_1, m_2);
		k1_mond = f(y_mond, y_planet, m_2, m_1);

		for (unsigned int i = 0; i < y_planet.size(); i++) {
			y_step_planet[i] = y_planet[i] + 0.5 * k1_planet[i];
			y_step_mond[i] = y_mond[i] + 0.5 * k1_mond[i];
		}

		k2_planet = f(y_step_planet, y_step_mond, m_1, m_2);
		k2_mond = f(y_step_mond, y_step_planet, m_2, m_1);

		for (unsigned int i = 0; i < y_planet.size(); i++) {
			y_step_planet[i] = y_planet[i] + 0.5 * k2_planet[i];
			y_step_mond[i] = y_mond[i] + 0.5 * k2_mond[i];
		}

		k3_planet = f(y_step_planet, y_step_mond, m_1, m_2);
		k3_mond = f(y_step_mond, y_step_planet, m_2, m_1);

		for (unsigned int i = 0; i < y_planet.size(); i++) {
			y_step_planet[i] = y_planet[i] + k3_planet[i];
			y_step_mond[i] = y_mond[i] + k3_mond[i];
		}

		k4_planet = f(y_step_planet, y_step_mond, m_1, m_2);
		k4_mond = f(y_step_mond, y_step_planet, m_2, m_1);

		for(unsigned int i = 0; i < y_planet.size(); i++) {
			y_planet[i] += (k1_planet[i] + 2*k2_planet[i] + 2*k3_planet[i] + k4_planet[i])/6;
			y_mond[i] += (k1_mond[i] + 2*k2_mond[i] + 2*k3_mond[i] + k4_mond[i])/6;
		}

		data_planet << y_planet[0] << "\t" << y_planet[1] << "\t" << y_planet[2] << "\n";
		data_mond << y_mond[0] << "\t" << y_mond[1] << "\t" << y_mond[2] << "\n";
		data_mond_rel << y_mond[0] - y_planet[0] << "\t" << y_mond[1] - y_planet[1] << "\t" << y_mond[2] - y_mond[2] << "\n";
	}

	data_planet.close();
	data_mond.close();
	data_mond_rel.close();

}

int main() {

	// Laufzeit starten
	clock_t check_clock;
	check_clock = clock();

	vector<double> y_planet = { 1, 0, 0, 0.0, sqrt(1.5), 0 }; // Startwerte, Planet
	vector<double> y_mond = { 1.1, 0, 0, 0, 4.7, 0 };

	runge_kutta(y_planet, y_mond, 1, 0.01);

	//Laufzeit berechnen
	check_clock = clock() - check_clock;
	cout << check_clock;

	return 0;
}
