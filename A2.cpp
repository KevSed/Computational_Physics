/*
 * A2.cpp
 *
 *  Created on: 04.05.2017
 *      Author: mona & kevin
 */


#include <iostream>
#include <iomanip>
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

//****************************************
// Funktionsvorschriften
//****************************************

double F1(double x){
        return exp(x)/x;
}

double F1_sing(double x, double lim){
			return lim+0.5*lim*lim*x;
}

double F2(double x){
        return 2 * exp(-x * x);
}

double F3(double x){
        if(x==0){
                return 1;
        }
        return sin(x)/x;
}

//****************************************
// Trapezregel
//****************************************

double trapez(double a, double b, double N, double(*f)(double)){
        double h=(b-a)/N;
        double summe=(h/2)*(f(a)+f(b));
        for (int n=1; n<N; n++){
                summe=summe+h*f(a+n*h);
        };
        return summe;
}

//****************************************************************
// Trapezregel für den singulären Bereich um 0 für Funktion 2
// benötigt ein zusätzliches Argument in f: lim.
//****************************************************************

double trapezlim(double a, double b, double N, double(*f)(double,double), double lim){
        double h=(b-a)/N;
        double summe=(h/2)*(f(a,lim)+f(b,lim));
        for (int n=1; n<N; n++){
                summe=summe+h*f(a+n*h,lim);
        };
        return summe;
}

//****************************************************************
// GENAUIGKEIT
//****************************************************************

double eps = 1e-5;

//****************************************************************
// Funktion für Aufgabenteil a)
//****************************************************************

double Trapez(double a, double b){
	double lim = 0.00001;
	        double Delta = 1;
	        double N=2;
	        //double links;
	        double integral1;
	        int i = 0;

	    	while(abs(Delta)>eps){
	    				++i;
	    	      integral1 = trapez(a,-lim,N,F1)+trapez(lim,b,N,F1)+trapezlim(-1,1,N,F1_sing,lim);
	    	      Delta = abs(integral1-2.1145018)/2.1145018; // relativer Fehler zum analytischen Wert
	    	      N = 2*N;
	    	      lim=lim/2;}
	    	cout << "Anzahl der Iterationen: " << i << endl;
        cout << "Wert des Integrals 2a): " << integral1 << endl;
	    	return integral1;

	return integral1;
}

//****************************************************************
// Funktion für Aufgabenteil b)
//****************************************************************

double F2Trapez(){
	double up = 2;
	        double Delta = 1;
	        double N=2;
	        double integral2;
	        int i = 0;

	    	while(abs(Delta)>eps){
	    		   ++i;
	    	     integral2 = trapez(0,up,N,F2);
	    	     Delta = abs(1.77245385-integral2)/1.77245385; // relativer Fehler zum "wahren" Wert
	    	     N = 2*N;
	    	     up *= 1.5;}
	    	cout << "\nAnzahl Iterationen: " << i << endl;
        cout << "Wert des Integrals 2b): " << integral2 << endl;
	    	return integral2;

	return integral2;
}

//****************************************************************
// Funktion für Aufgabenteil c)
//****************************************************************

double F3Trapez(){
	double ob = 1;
	        double Delta = 1;
	        double N=2;
	        double integral3;
	        int i = 0;

	    	while(abs(Delta)>eps){
	    				++i;
	    	      integral3 = 2 * trapez(0,ob,N,F3);
	    	      Delta = abs(M_PI-integral3)/M_PI; // relativer Fehler zum "wahren" Wert
	    	      N = 2*N;
	    	      ob *= 1.5;}
              cout << "\nAnzahl Iterationen: " << i << endl;
              cout << "Wert des Integrals 2c): " << integral3 << endl;
	    	return integral3;

	return integral3;
}

//****************************************
// main
//****************************************

int main(){
        int a1 = -1;
        int b1 = 1;
//Aufgabenteil a)
        cout << "****************\nAufgabe 2a)\n****************\n" << endl;
        cout << scientific << setprecision(8) << Trapez(a1,b1) << endl;
//Aufgabenteil b)
        cout << "\n****************\nAufgabe 2b)\n****************" << endl;
        cout << scientific << setprecision(8) << F2Trapez() << endl;
//Aufgabenteil c)
        cout << "\n****************\nAufgabe 2c)\n****************" << endl;
        cout << scientific << setprecision(8) << F3Trapez() << endl;

        return 0;
}
