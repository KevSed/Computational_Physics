/*
 * MagnetischesPendel.cpp
 *
 *  Created on: 12.05.2017
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


vector<double> rm1={0.,1.};
double a=sqrt(3)/2.;
vector<double> rm2={-a,-0.5};
vector<double> rm3={a,-0.5};

//****************************************
// Ein Schritt des Runge-Kutta Verfahrens 4. Ordnung
//****************************************


vector<double> Schritt(vector<double> (*F)(vector<double>, vector<double>), vector<double> y, double m, double h){

	double dimension= (y.size())/2;
	double m_1=m;
	// m_1 ist die Masse dehren Bahn verfolgt werden soll,
	//daher muss vor berechnung der Kraft durch m_1 dividiert werden

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
		f[i+dimension]=(1/m_1)*F(r,v)[i];
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
			k2[i+dimension]=(1/m_1)*F(rtemp,vtemp)[i];
			k2[i]=h*k2[i];
			k2[i+dimension]=h*k2[i+dimension];
	}

	for(int n=0;n<dimension;n++){
			rtemp[n]=r[n]+0.5*k2[n];
			vtemp[n]=v[n]+0.5*k2[n+dimension];
		}
	for(int i=0; i<dimension; i++){

			k3[i]=vtemp[i];
			k3[i+dimension]=(1/m_1)*F(rtemp,vtemp)[i];
			k3[i]=h*k3[i];
			k3[i+dimension]=h*k3[i+dimension];
	}

	for(int n=0;n<dimension;n++){
			rtemp[n]=r[n]+k3[n];
			vtemp[n]=v[n]+k3[n+dimension];
		}
	for(int i=0; i<dimension; i++){

			k4[i]=vtemp[i];
			k4[i+dimension]=(1/m_1)*F(rtemp,vtemp)[i];
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
// Potentialfunktion
//****************************************

double Potenzial(vector<double> r){

	double V=0;
	double alpha=1;
	double gamma=1;
	double z=0.1;

	for(int i=0;i<2;i++){
		V=V+(alpha/2)*r[i]*r[i];
	}
	double nenner1;
	nenner1 = (r[0]-rm1[0])*(r[0]-rm1[0])+(r[1]-rm1[1])*(r[1]-rm1[1])+z*z;
	nenner1 =sqrt(nenner1);

	double nenner2 = (r[0]-rm2[0])*(r[0]-rm2[0])+(r[1]-rm2[1])*(r[1]-rm2[1])+z*z;
	nenner2 =sqrt(nenner2);

	double nenner3 = (r[0]-rm3[0])*(r[0]-rm3[0])+(r[1]-rm3[1])*(r[1]-rm3[1])+z*z;
	nenner3 =sqrt(nenner3);

	V=V-gamma*(1/nenner1)-gamma*(1/nenner2)-gamma*(1/nenner3);

	return V;
}

//****************************************
// Kraftfeld
//****************************************

vector<double> Kraft(vector<double> r, vector<double> v){

	vector<double> F(2);
	double alpha=1;
	double gamma=1;
	//double beta=0.25;
	double beta=0.1;
	double z=0.1;


	double nenner1;
	nenner1 = (r[0]-rm1[0])*(r[0]-rm1[0])+(r[1]-rm1[1])*(r[1]-rm1[1])+z*z;
	nenner1 =sqrt(nenner1);
	nenner1=nenner1*nenner1*nenner1;

	double nenner2 = (r[0]-rm2[0])*(r[0]-rm2[0])+(r[1]-rm2[1])*(r[1]-rm2[1])+z*z;
	nenner2 =sqrt(nenner2);
	nenner2=nenner2*nenner2*nenner2;

	double nenner3 = (r[0]-rm3[0])*(r[0]-rm3[0])+(r[1]-rm3[1])*(r[1]-rm3[1])+z*z;
	nenner3 =sqrt(nenner3);
	nenner3=nenner3*nenner3*nenner3;

	for(int i=0;i<2;i++){
		F[i]=-alpha*r[i]-beta*v[i];
		F[i]=F[i]-gamma*(r[i]-rm1[i])/nenner1;
		F[i]=F[i]-gamma*(r[i]-rm2[i])/nenner2;
		F[i]=F[i]-gamma*(r[i]-rm3[i])/nenner3;
	}


	return F;
}


int main(){

	double toleranz=1e-4;

// ==============================================================Daten fuer den Plot des Potentials

	// ===============================
	ofstream Potentialdata;
	Potentialdata.open("Potential_V.txt");
	Potentialdata.precision(10);//
	// ===============================

	double genau=0.01;
	vector<double> r(2);
	vector<double> v(2);

	for(double n=-2.0;n<2.0;n=n+genau){
		for(double j=-2.0;j<2.0;j=j+genau){
			r[0]=n;
			r[1]=j;
			double V=Potenzial(r);
			Potentialdata << r[0] << "\t" << r[1] << "\t" << V << endl;
		}
	}

	Potentialdata.close();
	cout << "ende Potential" << endl;

// ==================================================Daten fuer die Integration der Bewegungsgleichungen

	// ===============================
	ofstream RK;
	RK.open("RungeKutta1.txt");
	RK.precision(10);//
	// ===============================

	double N=10000;
	double h=0.01;
	double m=1;

	vector<double> y(4);

	//Startwerte initialisieren
	r={2,0.99};
	v={0,0};

	for(int i=0; i<2;i++){
		y[i]=r[i];
		y[i+2]=v[i];
	}

	for(int i=0; i<N; i++){
		RK << i*h << "\t";
		//Koordinaten in Datei Schreiben
		for(int n=0; n<2; n++){
			RK << y[n] << "\t";
		}
		RK << endl;
		y=Schritt(Kraft,y,m,h);
	}
 RK.close();

//2.Startwert
	// ===============================
		ofstream RK2;
		RK2.open("RungeKutta2.txt");
		RK2.precision(10);//
		// ===============================

		//Startwerte initialisieren
		r={2,1.01};
		v={0,0};

		for(int i=0; i<2;i++){
			y[i]=r[i];
			y[i+2]=v[i];
		}

		for(int i=0; i<N; i++){
			RK2 << i*h << "\t";
			//Koordinaten in Datei Schreiben
			for(int n=0; n<2; n++){
				RK2 << y[n] << "\t";
			}
			RK2 << endl;
			y=Schritt(Kraft,y,m,h);
		}
//3. Startwert
	// ===============================
		ofstream RK3;
		RK3.open("RungeKutta3.txt");
		RK3.precision(10);//
	// ===============================
	// ===============================
		ofstream Energied;
		Energied.open("Energie.txt");
		Energied.precision(10);//
	// ===============================


		//Startwerte initialisieren
		r={2,1};
		v={0,0};
		double Epot;
		double Ekin;

		for(int i=0; i<2;i++){
			y[i]=r[i];
			y[i+2]=v[i];
		}

		for(int i=0; i<N; i++){
			for(int j=0; j<2;j++){
				r[j]=y[j];
				v[j]=y[j+2];
			}
			Epot=Potenzial(r);
			Ekin=0.5*m*(v[0]*v[0]+v[1]*v[1]);
			Energied << i*h << "\t" << Epot << "\t" << Ekin << "\t" << Epot+Ekin << endl;
			RK3 << i*h << "\t";
			//Koordinaten in Datei Schreiben
			for(int n=0; n<2; n++){
				RK3 << y[n] << "\t";
			}
			RK3 << endl;
			y=Schritt(Kraft,y,m,h);
		}
	cout << "Hier beginnen die Berechnungen für den optionalen Heatplot" << endl;


	// ===============================
	ofstream Heat;
	Heat.open("Heat.txt");
	//Heat.open("Heat_b25.txt");
	Heat.precision(10);//
	// ===============================
	double E_alt;
	double auflosung=1./100.;
	double abstand1;
	double abstand2;
	double abstand3;
	int color=0;
	double E;

	for(double a=-2; a<2; a=a+auflosung){
		for(double b=-2; b<2; b=b+auflosung){
			y[0]=a;
			r[0]=a;
			y[1]=b;
			r[1]=b;
			y[2]=0;
			y[3]=0;
			E=Potenzial(r)+0.5*m*(y[2]*y[2]+y[3]*y[3]);
			do{
				E_alt=E;
				//cout << a << "\t" << b << "\t" << E_alt << "\t" ;
				y=Schritt(Kraft,y,m,h);
				r[1]=y[0];
				r[2]=y[1];
				E=Potenzial(r)+0.5*m*(y[2]*y[2]+y[3]*y[3]);
				//cout << E_alt << endl ;
			}while(abs(E-E_alt)>toleranz);
			abstand1=0;
			abstand2=0;
			abstand3=0;
			for(int i=0; i<2;i++){
				abstand1=abstand1+(y[i]-rm1[i])*(y[i]-rm1[i]);
				abstand2=abstand2+(y[i]-rm2[i])*(y[i]-rm2[i]);
				abstand3=abstand3+(y[i]-rm3[i])*(y[i]-rm3[i]);
			}
			if(abstand1<abstand2 && abstand1<abstand3){color=1;}
			if(abstand2<abstand3 && abstand2<abstand1){color=2;}
			if(abstand3<abstand2 && abstand3<abstand1){color=3;}
			Heat << color <<"\t";
			cout << "Farbe" <<"\t" << a << "\t" << b << "\t" << color << endl;
		}
		Heat << endl;
	}

cout << "Ende der Main-Funktion" << endl;
return 0;
}