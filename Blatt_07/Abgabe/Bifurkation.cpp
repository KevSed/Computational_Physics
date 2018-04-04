/*
 * Bifurkation.cpp
 *
 *  Created on: 12.06.2017
 *      Author: mona
 */

#include <iostream>
#include <string>
#include <cstdlib>
#include <complex>
#include <vector>
#include <cmath>
#include <sstream>
#include <utility>
#include <math.h>
#include <fstream>
#include <functional>
#include "Eigen/Dense"
#include "Eigen/Core"

using namespace std;
using namespace Eigen;


double Step_log_Abb(double x1, double r){
	double x2;

	if(x1<0 or x1>1){
		cout << "x liegt ausserhalb des Intervalls [0,1]!" << endl;
		return 0.0;
	}

	x2=r*x1*(1-x1);
	return x2;
}

double Step_kub_Abb(double x1, double r){
	double x2;

	if(x1<(-sqrt(1+r)) or x1>(sqrt(1+r))){
		cout << "x liegt ausserhalb des Intervalls [-sqrt(1+r),sqrt(1+r)]!" << endl;
		cout << "x1= " << x1 << "    r= " << r << endl;
		return 0.0;
	}


	x2=r*x1-x1*x1*x1;
	return x2;
}

void iteration(double (*F)(double, double), double rmax , double delta_r, string Name, bool kubisch){
	int kalibrierung=300;
	int auswertungen=60;

	// ===============================
	ofstream Data;
	Data.open(Name.c_str());
	Data.precision(10);//
	// ===============================

	for (double r=0; r<rmax; r=r+delta_r){
		for (int j=0; j<auswertungen; j++){
			double x0=0.0;
			if(kubisch==true){
				x0=rand()* (1./RAND_MAX)*sqrt(1+r)*2 ;
				x0=x0-sqrt(1+r);
			}
			else{
				x0=rand()* (1./RAND_MAX);
			}

			for(int i=0; i<kalibrierung; i++){
			x0=F(x0, r);}

			Data << r << "\t" << x0 << endl;

		}
	}


}


int main(){
	//Fuer die logistische Abbildung divergiert der wert der iteration fuer r>4
	iteration(Step_log_Abb, 4.0 , 0.001, "logistische_Abb", false);
	cout << "Ab hier kubische Abbildung" << endl;
	//Fuer die kubische Abbildung divergiert der wert der iteration fuer r>3
	iteration(Step_kub_Abb, 3.00 , 0.001, "kubische_Abb", true);

	cout << "Ende der main-Funktion" << endl;
	return 0;
}