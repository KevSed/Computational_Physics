/*
 * harmOszillator.cpp
 *
 *  Created on: 08.05.2017
 *      Author: mona
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
// Ein Schritt des Runge-Kutta Verfahrens 4. Ordnung
//****************************************


vector<double> Schritt(vector<double> (*F)(vector<double>), vector<double> y, double h){

	double dimension= (y.size())/2;
	vector<double> f(y.size());
	vector<double> k2(y.size());
	vector<double> k3(y.size());
	vector<double> k4(y.size());
	vector<double> r(y.size()/2);
	vector<double> v(y.size()/2);

	for(int i=0; i<dimension; i++){
		r[i]=y[i];
		v[i]=y[i+dimension];
	}
		//cout << "r= " << r[0]  << "\t" << r[1] << "\t" << r[2] << endl;
		//cout << "v= " << v[0]  << "\t" << v[1] << "\t" << v[2] << endl;

	for(int i=0; i<dimension; i++){
		f[i]=v[i];
		f[i+dimension]=F(r)[i];
		f[i]=h*f[i];
		f[i+dimension]=h*f[i+dimension];
	}

	vector<double> rtemp(y.size()/2);
	vector<double> vtemp(y.size()/2);

	// Berechnung der Hilfsgroessen k
	for(int n=0;n<dimension;n++){
			rtemp[n]=r[n]+0.5*f[n];
			vtemp[n]=v[n]+0.5*f[n+dimension];
		}
	for(int i=0; i<dimension; i++){

			k2[i]=vtemp[i];
			k2[i+dimension]=F(rtemp)[i];
			k2[i]=h*k2[i];
			k2[i+dimension]=h*k2[i+dimension];
	}

	for(int n=0;n<dimension;n++){
			rtemp[n]=r[n]+0.5*k2[n];
			vtemp[n]=v[n]+0.5*k2[n+dimension];
		}
	for(int i=0; i<dimension; i++){

			k3[i]=vtemp[i];
			k3[i+dimension]=F(rtemp)[i];
			k3[i]=h*k3[i];
			k3[i+dimension]=h*k3[i+dimension];
	}

	for(int n=0;n<dimension;n++){
			rtemp[n]=r[n]+k3[n];
			vtemp[n]=v[n]+k3[n+dimension];
		}
	for(int i=0; i<dimension; i++){

			k4[i]=vtemp[i];
			k4[i+dimension]=F(rtemp)[i];
			k4[i]=h*k4[i];
			k4[i+dimension]=h*k4[i+dimension];
	}

	// Berechnung des neuen y
	for(int i=0; i<2*dimension; i++){
			y[i]=y[i]+(1./6.)*(f[i]+2*k2[i]+2*k3[i]+k4[i]);
	}

	return y;
}

//****************************************
// Ausführen des Runge-Kutta Verfahrens 4. Ordnung
//****************************************


void RungeKutta(vector<double> (*F)(vector<double>) , vector<double> ystart, double N, double h){
	//cout << "Start"<< endl;
	double dimension= (ystart.size())/2;
	//cout << "Dimension= " << dimension << endl;
	vector<double> y(ystart.size());

	for(int i=0; i<2*dimension;i++){
		y[i]=ystart[i];
	}

	// ===============================
	ofstream a;
	//a.open("RungeKutta.txt");
	//a.open("TrajektorieV0.txt");
	//a.open("TrajektorieV1.txt");
	//a.open("TrajektorieF2V0.txt");
	a.open("TrajektorieF2V1.txt");
	a.precision(10);//
	// ===============================
	// ===============================
	ofstream Energiedata;
	Energiedata.open("Energie.txt");
	Energiedata.precision(10);//
	// ===============================

	double Energie;


	for(int i=0; i<N; i++){
		//Berechnung der Energie
		double kinE=0;
		double potE=0;
		for(int n=0; n<dimension;n++){
			potE=potE+0.5*y[n]*y[n];
			kinE=kinE+0.5*y[n+dimension]*y[n+dimension];
		}
		Energie=potE+kinE;
		Energiedata << i*h  << "\t" << Energie <<  "\t" << potE << "\t" << kinE <<endl;
		//Werte fur den Zeitpunkt in die ersten 3 Spalten der Datei schreiben:
		a << i << "\t";
		a << h << "\t"<< i*h << "\t";
		//Werte fuer die Koordinaten in Spalte 4 bis 9 schreiben
		//for(int n=0; n<dimension; n++){
		for(int n=0; n<2*dimension; n++){
			a << y[n] << "\t";
		}
		a << endl;
		y=Schritt(F,y,h);
	}

	a.close();
	return;
}


//****************************************
// Runge Kutta mit angepasster Schrittweite
//****************************************

void RungeKuttaSchrittweite(vector<double> (*F)(vector<double>) , vector<double> ystart, double toleranz){
	//cout << "Start"<< endl;
	double dimension= (ystart.size())/2;
	//cout << "Dimension= " << dimension << endl;
	vector<double> y(ystart.size());

	for(int i=0; i<2*dimension;i++){
		y[i]=ystart[i];
	}

	double h=0.1;
	bool erreicht = false;
	double Fehler;
	while(erreicht==false){
	for(int i=0; (i*h)<7; i++){
		y=Schritt(F,y,h);
		int count=0;
		//cout << "count wird genullt" << endl;
			if((i*h)>1){

//hier war der Fehler: wir haben N so klein gewhählt dass nach weniger als einer Periode aufgehoert wurde.
				for (int n=0; n<2*dimension;n++){
					Fehler = abs(y[n]-ystart[n])/abs(ystart[n]);
					//cout << "Fehler von "<< n << " betraegt " << Fehler  << endl;
					if (Fehler < toleranz){
						count=count+1;
					//cout << "count wird auf" << count << " hochgesetzt" << endl;
					}
				}
		if(count==2*dimension){
		erreicht=true;
		cout << "h= " << h << "   bei einer Toleranz von = " << toleranz << endl;
		return;
		}
		}
	}
	h=0.9*h;
	//cout << "h wird auf " << h << " verkleinert" << endl;
	}
	return;
}



//****************************************
// Berechnung eines Kraftfeldes
//****************************************


vector<double> Feld_harmOsz(vector<double> y){
	vector<double> F(y.size()/2);
	double dimension= (y.size());
	//cout << "Dimension Schleife= " << dimension << endl;
	for(int i=0; i<dimension; i++){
			F[i]=-y[i];
	}
	//cout << "F= " << F[0]  << "\t" << F[1] << "\t" << F[2] << endl;
	return F;
}

//****************************************
// Zusätzliches harmonisches Kraftfeld
//****************************************


vector<double> Feld_neu(vector<double> y){
	vector<double> F(y.size()/2);
	vector<double> a={1,0,0};

	double dimension= (y.size());
	//cout << "Dimension Schleife= " << dimension << endl;
	for(int i=0; i<dimension; i++){
			F[i]=-2*y[i]+a[i];
	}
	//cout << "F= " << F[0]  << "\t" << F[1] << "\t" << F[2] << endl;
	return F;
}

int main(){
	vector<double> ystart={3,2,1,0,0,0};
	//v=0
	//RungeKutta(Feld_harmOsz,ystart,10000,0.01);
	vector<double> y2={3,2,1,1,2,2};
	//v ungleich null
	//RungeKutta(Feld_harmOsz,y2,10000,0.01);
	vector<double> y_neq_null={2,1,1,1,1,1};

	//Berechnung der Daten fuer die Schrittweite
	/*
	double x=0.01;
	for(int i=1; i<20; i++){
		RungeKuttaSchrittweite(Feld_harmOsz , y_neq_null, x);
		x=0.5*x;
	}
	*/
	//Berechnung der Groeßen mit zusaetzlichem Kraftfeld
	//RungeKutta(Feld_neu,ystart,10000,0.01);
	RungeKutta(Feld_neu,y2,10000,0.01);

	//Berechnung der Trajektorien

return 0;
}
