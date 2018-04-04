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

vector<double> a = { 1, 0, 0 };

vector<double> f_a(vector<double> y, double h) {
	vector<double> y_neu(y.size());

	for (unsigned int i = 0; i < y.size() / 2; i++) {
		y_neu[i] = h * y[i + y.size() / 2];
		y_neu[i + y.size() / 2] = -h * y[i];
	}

	return y_neu;
}

vector<double> give_abs(vector<double> y) {
	vector<double> abs_vec(2);
	for (unsigned int i = 0; i < y.size()/2; i++) {
		abs_vec[0] += y[i] * y[i];
		abs_vec[1] += y[i+y.size()/2]*y[i+y.size()/2];
	}
	abs_vec[0] = sqrt(abs_vec[0]);
	abs_vec[1] = sqrt(abs_vec[1]);
	return abs_vec;
}

vector<double> f_b(vector<double> y, double h) {

	vector<double> y_neu(y.size());

	for (unsigned int i = 0; i < y.size() / 2; i++) {
		y_neu[i] = h * y[i + y.size() / 2];
		y_neu[i + y.size() / 2] = -h * (2 * y[i] - a[i]);
	}

	return y_neu;
}

vector<double> E(vector<double> y) {

	vector<double> energy(3);

	for (unsigned int i = 0; i < y.size() / 2; i++) {

		energy[0] += y[i + y.size() / 2] * y[i + y.size() / 2];
		energy[1] += y[i] * y[i];
	}

	energy[0] *= 0.5;
	energy[1] *= 0.5;
	energy[2] = energy[0] + energy[1];

	return energy;

}

void runge_kutta(vector<double> y, vector<double> (*f)(vector<double>, double),
		double h) {

	// ===============================
	ofstream data;
	data.open("runge_kutta.txt");
	data.precision(10);

	ofstream energie;
	energie.open("energie.txt");
	energie.precision(10);

	ofstream r_v;
	r_v.open("r_v_plot.txt");
	r_v.precision(10);
	// ===============================

	vector<double> k1(y.size()), k2(y.size()), k3(y.size()), k4(y.size()),
			y_step(y.size());
	vector<double> energy(3);
	vector<double> r_v_vec(2);

	for (int i = 0; i < 10000; i++) {

		k1 = f(y, h);

		for (unsigned int i = 0; i < y.size(); i++) {
			y_step[i] = y[i] + 0.5 * k1[i];
		}

		k2 = f(y_step, h);

		for (unsigned int i = 0; i < y.size(); i++) {
			y_step[i] = y[i] + 0.5 * k2[i];
		}

		k3 = f(y_step, h);

		for (unsigned int i = 0; i < y.size(); i++) {
			y_step[i] = y[i] + k3[i];
		}

		k4 = f(y_step, h);

		for (unsigned int i = 0; i < y.size(); i++) {
			y[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
		}

		energy = E(y);

		data << h * i << "\t";
		for (unsigned int i = 0; i < y.size(); i++) {
			data << y[i] << "\t";
		}
		data << "\n";

		r_v_vec = give_abs(y);

		energie << h * i << "\t" << energy[0] << "\t" << energy[1] << "\t"
				<< energy[2] << "\n";

		r_v << r_v_vec[0] << "\t" << r_v_vec[1] << "\n";
	}

	data.close();
	energie.close();
	r_v.close();
}

void h_fehler(vector<double> (*f)(vector<double>,double)) {

	// ===============================
	ofstream data;
	data.open("fehler_h.txt");
	data.precision(10);

	ofstream data2;
	data2.open("fehler_energie.txt");
	data2.precision(10);
	// ===============================

	vector<double> y = {1,0};

	double E_ges = E(y)[2];

	vector<double> k1(y.size()), k2(y.size()), k3(y.size()), k4(y.size()),
			y_step(y.size());

	double h;

	for (int k = 10; k < 10000; k+=10) {

		y = {1,0};

		h = M_PI / k;

		for (int i = 0; i*h <= 10*M_PI; i++) {

			k1 = f(y, h);

			for (unsigned int i = 0; i < y.size(); i++) {
				y_step[i] = y[i] + 0.5 * k1[i];
			}

			k2 = f(y_step, h);

			for (unsigned int i = 0; i < y.size(); i++) {
				y_step[i] = y[i] + 0.5 * k2[i];
			}

			k3 = f(y_step, h);

			for (unsigned int i = 0; i < y.size(); i++) {
				y_step[i] = y[i] + k3[i];
			}

			k4 = f(y_step, h);

			for (unsigned int i = 0; i < y.size(); i++) {
				y[i] += (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
			}

		}

		data << h << "\t" << abs(1 - y[0]) << "\n";
		data2 << h << "\t" << abs((E_ges - E(y)[2])/E_ges) << "\n";

	}
	data.close();
}

int main() {

	// Laufzeit starten
	clock_t check_clock;
	check_clock = clock();

	vector<double> y = { 2, 2, 2, 1, 0, 0 }; // Startwerte

	runge_kutta(y, f_b, 0.01);

	//h_fehler(f_a);

	//Laufzeit berechnen
	check_clock = clock() - check_clock;
	cout << check_clock;

	return 0;
}
