/*
 * Doppelpendel.cpp
 *
 *  Created on: 15.05.2017
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

double mu= 0.5;
double g1=9.81;
double g2=9.81;
double lambda=1;

//****************************************
// Funktionen fuer 2. Ableitungen
//****************************************

double phi_1(vector<double> y){
	double phi1;

	double Vorfaktor = mu*cos(y[1]-y[0])*cos(y[1]-y[0]);
	Vorfaktor=1-Vorfaktor;
	Vorfaktor=1/Vorfaktor;

	double Faktor1=mu*g1*sin(y[1])*cos(y[1]-y[0]);
	double Faktor2=mu*y[2]*y[2]*sin(y[1]-y[0])*cos(y[1]-y[0]);
	double Faktor3=-g1*sin(y[0]);
	double Faktor4=(mu/lambda)*y[3]*y[3]*sin(y[1]-y[0]);

	phi1=Vorfaktor*(Faktor1+Faktor2+Faktor3+Faktor4);
	return phi1;
}

double phi_2(vector<double> y){
	double phi2;

	double Vorfaktor = mu*cos(y[1]-y[0])*cos(y[1]-y[0]);
	Vorfaktor=1-Vorfaktor;
	Vorfaktor=1/Vorfaktor;

	double Faktor1=g2*sin(y[0])*cos(y[1]-y[0]);
	double Faktor2=-mu*y[3]*y[3]*sin(y[1]-y[0])*cos(y[1]-y[0]);
	double Faktor3=-g2*sin(y[1]);
	double Faktor4=-lambda*y[2]*y[2]*sin(y[1]-y[0]);

	phi2=Vorfaktor*(Faktor1+Faktor2+Faktor3+Faktor4);
	return phi2;
}

vector<double> Schritt(double (*phi1)(vector<double>), double (*phi2)(vector<double>), vector<double> y, double h){

	vector<double> k1(4);
	vector<double> k2(4);
	vector<double> k3(4);
	vector<double> k4(4);
	vector<double> ytemp(4);

	k1[0]=h*y[2];
	k1[1]=h*y[3];
	k1[2]=h*phi1(y);
	k1[3]=h*phi2(y);

	for(int i=0; i<4; i++){
		ytemp[i]=y[i]+0.5*k1[i];
	}

	k2[0]=h*ytemp[2];
	k2[1]=h*ytemp[3];
	k2[2]=h*phi1(ytemp);
	k2[3]=h*phi2(ytemp);

	for(int i=0; i<4; i++){
		ytemp[i]=y[i]+0.5*k2[i];
	}

	k3[0]=h*ytemp[2];
	k3[1]=h*ytemp[3];
	k3[2]=h*phi1(ytemp);
	k3[3]=h*phi2(ytemp);

	for(int i=0; i<4; i++){
		ytemp[i]=y[i]+k3[i];
	}

	k4[0]=h*ytemp[2];
	k4[1]=h*ytemp[3];
	k4[2]=h*phi1(ytemp);
	k4[3]=h*phi2(ytemp);


	// Berechnung des neuen y
	for(int i=0; i<4; i++){
			y[i]=y[i]+(1./6.)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
	}

	return y;
}

//****************************************
// Funktionen zur Berechnung der Energie
//****************************************

double Epot(vector<double> y){
	double Epot=0;
	double m1=1;
	double m2=1;
	double L=1;

	Epot=(m1+m2)*L*(1-cos(y[0]))+m2*L*(1-cos(y[1]));
	Epot=g1*Epot;
	return Epot;
}

double Ekin(vector<double> y){
	double Ekin=0;
	double m1=1;
	double m2=1;
	double L=1;

	double V1x=(L*cos(y[0])*y[2]);
	double V1y=(L*sin(y[0])*y[2]);
	double V2x=(L*cos(y[0])*y[2])+(L*cos(y[1])*y[3]);
	double V2y=(L*sin(y[0])*y[2])+(L*sin(y[1])*y[3]);

	double Ekin1=0.5*m1*(V1x*V1x+V1y*V1y);
	double Ekin2=0.5*m2*(V2x*V2x+V2y*V2y);

	Ekin=Ekin1+Ekin2;

	return Ekin;
}


int main(){

	double N=10000;
	double h=0.01;
	vector<double> phi(4);

	vector<double> phi_period={0.1,0.1*sqrt(2),0,0};
	vector<double> phi_quasi={0,0,0,4.472};
	vector<double> phi_chaos={0,0,0,11.832};

	// ===============================
	ofstream RKp;
	RKp.open("RungeKutta_period.txt");
	RKp.precision(10);
	ofstream RKpEnergie;
	RKpEnergie.open("RungeKutta_period_E.txt");
	RKpEnergie.precision(10);
	// ===============================

		for(int i=0; i<4;i++){
			phi[i]=phi_period[i];
		}

		for(int i=0; i<N; i++){
			RKp << i*h << "\t";
			RKpEnergie << i*h << "\t" << Epot(phi) << "\t" << Ekin(phi) << "\t" << Epot(phi)+Ekin(phi) << endl;
			//Koordinaten in Datei Schreiben

			for(int n=0; n<4; n++){
				RKp << phi[n] << "\t";
			}

			RKp << endl;
			phi=Schritt(phi_1,phi_2,phi,h);
		}
	 RKp.close();

		// ===============================
		ofstream RKq;
		RKq.open("RungeKutta_quasi.txt");
		RKq.precision(10);//
		ofstream RKqEnergie;
		RKqEnergie.open("RungeKutta_quasi_E.txt");
		RKqEnergie.precision(10);
		// ===============================

		for(int i=0; i<4;i++){
			phi[i]=phi_quasi[i];
		}

		for(int i=0; i<N; i++){
			RKq << i*h << "\t";
			RKqEnergie << i*h << "\t" << Epot(phi) << "\t" << Ekin(phi) << "\t" << Epot(phi)+Ekin(phi)<< endl;
			//Koordinaten in Datei Schreiben

			for(int n=0; n<4; n++){
				RKq << phi[n] << "\t";
			}
			RKq << endl;
			phi=Schritt(phi_1,phi_2,phi,h);
			}
		 RKq.close();

		// ===============================
		ofstream RKc;
		RKc.open("RungeKutta_chaos.txt");
		RKc.precision(10);//
		ofstream RKcEnergie;
		RKcEnergie.open("RungeKutta_chaos_E.txt");
		RKcEnergie.precision(10);
		// ===============================
		double pi=M_PI;
		for(int i=0; i<4;i++){
			phi[i]=phi_chaos[i];
		}

		for(int i=0; i<N; i++){
			RKc << i*h << "\t";
			//Koordinaten in Datei Schreiben
			RKcEnergie << i*h << "\t" << Epot(phi) << "\t" << Ekin(phi) << "\t" << Epot(phi)+Ekin(phi)<< endl;
			//Hier muss darauf geachtet werden, dass sich das Pendel bei der gegebenen Anfangsgeschrindigkeit ueberschlaegt
			for(int n=0; n<2; n++){
				while(phi[n]>pi){phi[n]=phi[n]-2*pi;}
				while(phi[n]<-pi){phi[n]=phi[n]+2*pi;}
				RKc << phi[n] << "\t";
			}
			for(int n=2; n<4; n++){
				RKc << phi[n] << "\t";
						}
			RKc << endl;
			phi=Schritt(phi_1,phi_2,phi,h);
				}
		RKc.close();

	cout << "Epot_period= " << Epot(phi_period) << "\t" << "Ekin_period= " << Ekin(phi_period);
	cout << "\t" << "E_period= " << Ekin(phi_period)+Epot(phi_period) << endl;

	cout << "Epot_quasi= " << Epot(phi_quasi) << "\t" << "Ekin_quasi= " << Ekin(phi_quasi);
	cout << "\t" << "E_quasi= " << Ekin(phi_quasi)+Epot(phi_quasi) << endl;

	cout << "Epot_chaos= " << Epot(phi_chaos) << "\t" << "Ekin_chaos= " << Ekin(phi_chaos);
	cout << "\t" << "E_chaos= " << Ekin(phi_chaos)+Epot(phi_chaos) << endl;

	cout << endl << "Ende der main-Funktion" << endl;
	return 0;
}
