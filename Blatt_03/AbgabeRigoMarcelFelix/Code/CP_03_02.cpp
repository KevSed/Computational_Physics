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
const double m_1 = 1;

vector<double> kepler_1(vector<double> y, double alpha) {
	vector<double> y_neu(y.size());

	double r = 0;
	for (unsigned int i = 0; i < y.size() / 2; i++) {
		r += y[i] * y[i];
	}
	r = sqrt(r);

	for (unsigned int i = 0; i < y.size() / 2; i++) {
		y_neu[i] = h * y[i + y.size() / 2];
		y_neu[i + y.size() / 2] = -h * G * alpha * y[i] / (pow(r, alpha + 2));
	}

	return y_neu;
}

double E(vector<double> y, double alpha) {

	double r = 0;
	for (unsigned int i = 0; i < y.size() / 2; i++) {
		r += y[i] * y[i];
	}
	r = sqrt(r);

	double E = 0;

	for (unsigned int i = 0; i < y.size() / 2; i++) {

		E += 0.5 * y[i + y.size() / 2] * y[i + y.size() / 2];

	}

	return E - G / pow(r, alpha);

}

// Nur für drei Raumdimensionen wegen Kreuzprodukt
double L(vector<double> y) {
	return sqrt(
			(y[1] * y[5] - y[2] * y[4]) * (y[1] * y[5] - y[2] * y[4])
					+ (y[2] * y[3] - y[0] * y[5]) * (y[2] * y[3] - y[0] * y[5])
					+ (y[0] * y[4] - y[1] * y[3]) * (y[0] * y[4] - y[1] * y[3]));
}

void runge_kutta(vector<double> y, vector<double> (*f)(vector<double>, double),
		double alpha) {

	// ===============================
	ofstream data;
	data.open("runge_kutta.txt");
	data.precision(10);
	// ===============================

	// ===============================
	ofstream erhaltung;
	erhaltung.open("Erhaltung.txt");
	erhaltung.precision(10);
	// ===============================

	vector<double> k1(y.size()), k2(y.size()), k3(y.size()), k4(y.size()),
			y_step(y.size());

	for (int i = 0; i < 100000; i++) {

		k1 = f(y, alpha);

		for (unsigned int i = 0; i < y.size(); i++) {
			y_step[i] = y[i] + 0.5 * k1[i];
		}

		k2 = f(y_step, alpha);

		for (unsigned int i = 0; i < y.size(); i++) {
			y_step[i] = y[i] + 0.5 * k2[i];
		}

		k3 = f(y_step, alpha);

		for (unsigned int i = 0; i < y.size(); i++) {
			y_step[i] = y[i] + k3[i];
		}

		k4 = f(y_step, alpha);

		for (unsigned int i = 0; i < y.size(); i++) {
			y[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
		}

		data << h * i << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << "\n";
		erhaltung << h * i << "\t" << E(y, alpha) << "\t" << L(y) << "\n";
	}

	data.close();
	erhaltung.close();
}


int main() {

	// Laufzeit starten
	clock_t check_clock;
	check_clock = clock();

	vector<double> y = { 1, 0, 0, 0, sqrt(1.5), 0 }; // Startwerte

	runge_kutta(y, kepler_1, 1);

	//Laufzeit berechnen
	check_clock = clock() - check_clock;
	cout << check_clock;

	return 0;
}
