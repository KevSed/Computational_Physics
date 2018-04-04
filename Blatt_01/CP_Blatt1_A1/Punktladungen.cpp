/*
 * Punktladungen.cpp
 *
 *  Created on: 24.04.2017
 *      Author: mona
 */



#include <iostream>
#include <complex>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <sstream>
#include <utility>
#include <math.h>
#include <fstream>
#include <functional>
using namespace std;

double Energie(double r, double a){
	double r1=sqrt(2*a*a-2*a*r+r*r);
	double r2=sqrt(2*a*a+2*a*r+r*r);
	double E=-2*((1/r1)+(1/r2));
	return E;
}

double zPunkt(double h, double r, double a){
	double f=Energie(r+h,a)-Energie(r-h,a);
	f=f/(2*h);
	return f;
}

double Kraft(double r, double a){
	double t1=sqrt(2*a*a-2*a*r+r*r)*(2*a*a-2*a*r+r*r); //entspricht ^3/2
	double t2=sqrt(2*a*a+2*a*r+r*r)*(2*a*a+2*a*r+r*r);
	t1=(-2*a+2*r)/t1;
	t2=(2*a+2*r)/t2;
	double F=t1+t2;
	return -F;
}

int main(){
	double a=1;
	double relFehler;
	//Abspeichern der Daten für die Energie
	ofstream b;
	b.open("Energie.txt");
	b.precision(10);
	//

	for(double n=-3*a; n<=3*a; n=n+0.001){
		b << n << "\t" << Energie(n,a) <<  "\n";
	}
	b.close();

	//Abspeichern der Daten
		ofstream c;
		c.open("Kraft_03a.txt");
		c.precision(10);
		//

		for(double n=-3*a; n<=3*a; n=n+0.001){
			relFehler=(-zPunkt(0.3*a,n,a)-Kraft(n,a))/Kraft(n,a);
			c << n << "\t" << -zPunkt(0.3*a,n,a) << "\t" << Kraft(n,a) << "\t" << relFehler <<"\n";
		}
		c.close();

		//Abspeichern der Daten
		ofstream d;
		d.open("Kraft_E4a.txt");
		d.precision(10);
		//
		double h=pow(10,-4);
		for(double n=-3*a; n<=3*a; n=n+0.001){
			relFehler=(-zPunkt(h*a,n,a)-Kraft(n,a))/Kraft(n,a);
			d << n << "\t" << -zPunkt(h*a,n,a) << "\t" << Kraft(n,a)<< "\t" << relFehler  << "\n";
		}
		d.close();

		//Abspeichern der Daten
		ofstream g;
		g.open("Kraft_E15a.txt");
		g.precision(10);
				//
		h=pow(10,-15);
		for(double n=-3*a; n<=3*a; n=n+0.001){
			relFehler=(-zPunkt(h*a,n,a)-Kraft(n,a))/Kraft(n,a);
			g << n << "\t" << -zPunkt(h*a,n,a) << "\t" << Kraft(n,a) << "\t" << relFehler <<"\n";
				}
		g.close();


	return 0;
}

