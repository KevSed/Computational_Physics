/*
 * Ising_Magnetisierung.cpp
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
//#include "Eigen/Dense"
//#include "Eigen/Core"

using namespace std;
//using namespace Eigen;



double Funktion(double h, double t, double m){
	double F=m-tanh((h+m)/t);
	return F;
}

double Startwert_select(double h, double t){
	//cout << "startwert wird gesucht" << endl;
	double prod=1;
	double x0;
	double y0;
	while (prod>=0){
		x0=rand()*1./RAND_MAX*2;
		y0=-rand()*1./RAND_MAX*2;
		//cout << "x0= " << x0 << '\t' << "y0= " << y0 << endl;
		prod=Funktion(h,t,x0)*Funktion(h,t,y0);
	}
	double x=x0;
	double y=y0;
	double krit=0.0001;

	while(Funktion(h,t,x)>krit){
		double z=(x+y)/2;
		if(Funktion(h,t,z)*Funktion(h,t,x)<0){y=z;}
		else{x=z;}
	}
	return x;
}

double Newton_Raphson_Step(double m, double h, double t){

	double f=m-tanh((h+m)/t);
	double sech_quad=cosh(2*(h+m)/t)+1;
	sech_quad=2/sech_quad;
	double f_strich=1-(1/t)*sech_quad;

	double m_new=m-(f/f_strich);

	return m_new;
}

double Newton_Raphson_Iterate(double m0, double h, double t){

	double m_krit=0.0001;
	double m_new=0;
	double diff=1;
	double m=m0;
	int counter=0;
	int counter_max=50;

	while(diff>m_krit && counter < counter_max){
		m_new=Newton_Raphson_Step(m, h, t);
		diff=abs(m_new-m);
		m=m_new;
		counter ++;
		if (counter==counter_max){cout << "Achtung, wahrscheinlich liegen mehrere Fixpunkte vor!" << endl;}
	}
	//cout << "counter" << counter << endl;
	return m;
}

void Newton_Raphson_tvari( double h, string Name ){

	double t_start=0.01;
	double t_max=3;
	double tstep=0.001;

	// ===============================
	ofstream Data;
	Data.open(Name.c_str());
	Data.precision(10);//
	// ===============================

	for(double t=t_start; t<t_max; t=t+tstep){
		double m0=Startwert_select(h, t);
		double m=Newton_Raphson_Iterate(m0, h, t);
		Data << t << "\t" << m << endl;
	}

}

void Newton_Raphson_hvari(double t, string Name ){

	double h_start=-3.01;
	double h_max=3;
	double hstep=0.00001;

	// ===============================
	ofstream Data;
	Data.open(Name.c_str());
	Data.precision(10);//
	// ===============================

	for(double h=h_start; h<h_max; h=h+hstep){
		//cout << "beginn der for-Schleife" << endl;
		double m0=Startwert_select(h, t);
		double m=Newton_Raphson_Iterate(m0, h, t);
		Data << h << "\t" << m << endl;
	}

}




int main(){

	double h;
	double t;

	cout << "h=0: " << endl;
	h=0.0;
	Newton_Raphson_tvari(h, "Data_h0" );
	cout << "h=01: " << endl;
	h=0.1;
	Newton_Raphson_tvari( h, "Data_h01" );
	cout << "h=05: " << endl;
	h=0.5;
	Newton_Raphson_tvari(h, "Data_h05" );

	cout << "t=10: " << endl;
	t=1.0;
	Newton_Raphson_hvari(t, "Data_t10" );
	cout << "t=15: " << endl;
	t=1.5;
	Newton_Raphson_hvari(t, "Data_t15" );
	cout << "t=05: " << endl;
	t=0.5;
	Newton_Raphson_hvari( t, "Data_t05" );

	cout << "Ende der main-Funktion" << endl;
	return 0;
}
