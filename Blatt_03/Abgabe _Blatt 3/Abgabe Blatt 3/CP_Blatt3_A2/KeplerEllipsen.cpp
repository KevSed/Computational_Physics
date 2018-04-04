/*
 * KeplerEllipsen.cpp
 *
 *  Created on: 10.05.2017
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
	double m_1=1;
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
		f[i+dimension]=(1/m_1)*F(r)[i];
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
			k2[i+dimension]=(1/m_1)*F(rtemp)[i];
			k2[i]=h*k2[i];
			k2[i+dimension]=h*k2[i+dimension];
	}

	for(int n=0;n<dimension;n++){
			rtemp[n]=r[n]+0.5*k2[n];
			vtemp[n]=v[n]+0.5*k2[n+dimension];
		}
	for(int i=0; i<dimension; i++){

			k3[i]=vtemp[i];
			k3[i+dimension]=(1/m_1)*F(rtemp)[i];
			k3[i]=h*k3[i];
			k3[i+dimension]=h*k3[i+dimension];
	}

	for(int n=0;n<dimension;n++){
			rtemp[n]=r[n]+k3[n];
			vtemp[n]=v[n]+k3[n+dimension];
		}
	for(int i=0; i<dimension; i++){

			k4[i]=vtemp[i];
			k4[i+dimension]=(1/m_1)*F(rtemp)[i];
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
// Ein Schritt des Runge-Kutta Verfahrens 4. Ordnung, wobei die Masse uebergeben wird
//****************************************

vector<double> Schritt2(vector<double> (*F)(vector<double>, vector<double> , double, double), vector<double> y1, vector<double> y2 , double h, double m1, double m2){

	double dimension= (y1.size())/2;
	double m=m1;
	// m_1 ist die Masse dehren Bahn verfolgt werden soll,
	//daher muss vor berechnung der Kraft durch m_1 dividiert werden

	vector<double> f(y1.size());
	vector<double> k2(y1.size());
	vector<double> k3(y1.size());
	vector<double> k4(y1.size());
	vector<double> r(y1.size()/2);
	vector<double> v(y1.size()/2);
	vector<double> Abstand(y1.size()/2);

	for(int i=0; i<dimension; i++){
		r[i]=y1[i];
		v[i]=y1[i+dimension];
		Abstand[i]=y1[i]-y2[i];

	}
		//cout << "r= " << r[0]  << "\t" << r[1] << "\t" << r[2] << endl;
		//cout << "v= " << v[0]  << "\t" << v[1] << "\t" << v[2] << endl;

	for(int i=0; i<dimension; i++){
		f[i]=v[i];
		f[i+dimension]=(1/m)*F(r,Abstand,m1,m2)[i];
		f[i]=h*f[i];
		f[i+dimension]=h*f[i+dimension];
	}

	vector<double> rtemp(y1.size()/2);
	vector<double> vtemp(y1.size()/2);

	// Berechnung der Hilfsgroessen k
	for(int n=0;n<dimension;n++){
			rtemp[n]=r[n]+0.5*f[n];
			vtemp[n]=v[n]+0.5*f[n+dimension];
		}
	for(int i=0; i<dimension; i++){

			k2[i]=vtemp[i];
			k2[i+dimension]=(1/m)*F(rtemp,Abstand,m1,m2)[i];
			k2[i]=h*k2[i];
			k2[i+dimension]=h*k2[i+dimension];
	}

	for(int n=0;n<dimension;n++){
			rtemp[n]=r[n]+0.5*k2[n];
			vtemp[n]=v[n]+0.5*k2[n+dimension];
		}
	for(int i=0; i<dimension; i++){

			k3[i]=vtemp[i];
			k3[i+dimension]=(1/m)*F(rtemp,Abstand,m1,m2)[i];
			k3[i]=h*k3[i];
			k3[i+dimension]=h*k3[i+dimension];
	}

	for(int n=0;n<dimension;n++){
			rtemp[n]=r[n]+k3[n];
			vtemp[n]=v[n]+k3[n+dimension];
		}
	for(int i=0; i<dimension; i++){

			k4[i]=vtemp[i];
			k4[i+dimension]=(1/m)*F(rtemp,Abstand,m1,m2)[i];
			k4[i]=h*k4[i];
			k4[i+dimension]=h*k4[i+dimension];
	}

	// Berechnung des neuen y
	for(int i=0; i<2*dimension; i++){
			y1[i]=y1[i]+(1./6.)*(f[i]+2*k2[i]+2*k3[i]+k4[i]);
	}

	return y1;
}

//****************************************
// Ausführen des Runge-Kutta Verfahrens 4. Ordnung
//****************************************


void RungeKutta(vector<double> (*F)(vector<double>) , vector<double> ystart, double N, double h){
	//cout << "Start"<< endl;
	double dimension= (ystart.size())/2;
	cout << "Dimension= " << dimension << endl;
	vector<double> y(ystart.size());

	for(int i=0; i<2*dimension;i++){
		y[i]=ystart[i];
	}

	// ===============================
	ofstream a;
	a.open("Keplerohnemond.txt");
	a.precision(10);//
	// ===============================
	// ===============================
	ofstream Energiedata;
	Energiedata.open("Energie.txt");
	Energiedata.precision(10);//
	// ===============================

	double Energie;

	double G=1;
	double m_0=1;
	double m_1=1;

	for(int i=0; i<N; i++){


		//Berechnung der Energie
		double kinE=0;
		double potE=0;

		// Betrag von r ausrechnen
		double rbetrag=0;
		for(int i=0; i<dimension; i++){
				rbetrag=rbetrag+y[i]*y[i];
		}
		rbetrag=sqrt(rbetrag);
		//cout << "rbetrag:  " << rbetrag;

		//Berechnung der potentiellen Energie
		potE=potE-G*m_0*m_1/rbetrag;
		//cout << "   Epot:  " << potE;

		//Berechnung der kinetischen Energie
		for(int n=0; n<dimension;n++){
			kinE=kinE+0.5*m_1*y[n+dimension]*y[n+dimension];
		}
		//cout << "   Ekin:  " << kinE << endl;
		Energie=potE+kinE;

		// Berechnung des Dehimpulses
		vector<double> L(3);
		vector<double> r(y.size()/2);
		vector<double> v(y.size()/2);

		for(int i=0; i<dimension; i++){
			r[i]=y[i];
			v[i]=y[i+dimension];
		}
		//Berechnung des Kreuzproduktes
		L[0] = r[1]*v[2]-r[2]*v[1];
		L[1] = r[2]*v[0]-r[0]*v[2];
		L[2] = r[0]*v[1]-r[1]*v[0];

		double Lbetrag=0;

		for(int i=0; i<dimension; i++){
			Lbetrag=Lbetrag+L[i]*L[i];
		}
		Lbetrag=sqrt(Lbetrag);

		//Energiedatei mit daten fuellen:
		Energiedata << i  << "\t" << Energie << "\t" << m_1*Lbetrag << endl;

		//Koordinatendatei mit daten fuellen
		a << i << "\t";
		//Pruefen wo eine Periode zu ende ist:
		a << h << "\t"<< i*h << "\t";
		//Koordinaten in Datei Schreiben
		for(int n=0; n<dimension; n++){
			a << y[n] << "\t";
		}
		a << endl;
		y=Schritt(F,y,h);
	}

	a.close();
	Energiedata.close();
	return;
}

//****************************************
// Ausführen des Runge-Kutta Verfahrens 4. Ordnung fuer 2 Planeten
//****************************************


void RungeKuttaMond(vector<double> (*F)(vector<double>, vector<double> , double, double), vector<double> Mstart, vector<double> Pstart, double N, double h){
	//cout << "Start"<< endl;
	double dimension= 3;
	cout << "Dimension= " << dimension << endl;
	vector<double> Mond(6);
	vector<double> Planet(6);

	for(int i=0; i<2*dimension;i++){
		Planet[i]=Pstart[i];
		Mond[i]=Mstart[i];
	}

	// ===============================
	ofstream a;
	a.open("Kepler_2.txt");
	a.precision(10);//
	// ===============================

	double m_Planet=1;
	double m_Mond=0.01*m_Planet;

	for(int i=0; i<N; i++){

		a << i << "\t";
		a << h << "\t"<< i*h << "\t";

		for(int n=0; n<dimension; n++){
			a << Mond[n] << "\t" << Planet[n] << "\t";
		}
		a << endl;
		Mond=Schritt2(F,Mond,Planet,h,m_Mond,m_Planet);
		Planet=Schritt2(F,Planet,Mond,h,m_Planet,m_Mond);
	}

	a.close();
	return;
}

//****************************************
// Berechnung des Gravitationsfeldes
//****************************************


vector<double> Feld_grav(vector<double> y){
	vector<double> F(y.size()/2);
	double G=1;
	double m_0=1;
	double m_1=1;
	double dimension= (y.size());
	double rbetrag=0;
	//cout << "Dimension Schleife= " << dimension << endl;
	// Betrag von r ausrechnen
	for(int i=0; i<dimension; i++){
			rbetrag=rbetrag+y[i]*y[i];
	}

	rbetrag=sqrt(rbetrag);
	//cout << "betrag von r:  " << rbetrag << endl;
	rbetrag=rbetrag*rbetrag*rbetrag;

	for(int i=0; i<dimension; i++){
			F[i]=-G*m_0*m_1*y[i]/(rbetrag);
	}
	//cout << "F= " << F[0]  << "\t" << F[1] << "\t" << F[2] << endl;
	return F;
}



//****************************************
// Berechnung des Gravitationsfeldes für 2 Massen
//****************************************


vector<double> Feld_2Massen(vector<double> y, vector<double> Abstand, double m1, double m2){
	vector<double> F(y.size()/2);
	double G=1;
	double dimension= (y.size());

	double AbstandBetrag=0;
	double rbetrag=0;
	double mSonne=1;

	//cout << "Dimension Schleife= " << dimension << endl;
	// Betrag von r ausrechnen
	for(int i=0; i<dimension; i++){
			rbetrag=rbetrag+y[i]*y[i];
			AbstandBetrag=AbstandBetrag+Abstand[i]*Abstand[i];
	}

	rbetrag=sqrt(rbetrag);
	AbstandBetrag=sqrt(AbstandBetrag);
	//cout << "betrag von r:  " << rbetrag << endl;
	rbetrag=rbetrag*rbetrag*rbetrag;
	AbstandBetrag=AbstandBetrag*AbstandBetrag*AbstandBetrag;

	for(int i=0; i<dimension; i++){
			F[i]=-(G*mSonne*m1*y[i]/(rbetrag))-(G*m1*m2*Abstand[i]/(AbstandBetrag));
	}
	//cout << "F= " << F[0]  << "\t" << F[1] << "\t" << F[2] << endl;
	return F;
}

//****************************************
// Berechnung des abgeaenderten potentials
//****************************************
vector<double> Feld_alpha(vector<double> y){

	double alpha=1.1;
	vector<double> F(y.size()/2);
	double G=1;
	double m_0=1;
	double m_1=1;
	double dimension= (y.size());
	double rbetrag=0;
	//cout << "Dimension Schleife= " << dimension << endl;
	// Betrag von r ausrechnen
	for(int i=0; i<dimension; i++){
			rbetrag=rbetrag+y[i]*y[i];
	}

	rbetrag=sqrt(rbetrag);
	//cout << "betrag von r:  " << rbetrag << endl;

	rbetrag=exp(log(rbetrag)*(alpha+2));

	for(int i=0; i<dimension; i++){
			F[i]=-G*m_0*m_1*alpha*y[i]/(rbetrag);
	}
	//cout << "F= " << F[0]  << "\t" << F[1] << "\t" << F[2] << endl;
	return F;
}

int main(){

	//v_0 zu klein:
	vector<double> ystart_klein={1,0,0,0,0.35,0};

	//v_0 zu gross:
	vector<double> ystart_gros={1,0,0,0,1,2,0};
	//v_0 geschlossene Ellipse
	vector<double> ystart_geschlossen={1,0,0,0,0.8,0};
	vector<double> ystart_alpha={1,0,0,0,0.8,0};
	//RungeKutta(Feld_grav,ystart_geschlossen,10000,0.01);
	//RungeKutta(Feld_grav,ystart_klein,6459,0.01);
	//RungeKutta(Feld_grav,ystart_gros,3000,0.01);
	//RungeKutta(Feld_alpha,ystart_alpha,10000,0.01);
	//RungeKutta(Feld_grav,ystart_gros,10000,0.01);
	//Kepler mit 2 Massen
	double start = sqrt(3./2.);
	vector<double> Mstart = {1.1,0,0,0,4.7,0};
	vector<double> Pstart = {1,0,0,0,start,0};
	double N=700000;
	double h=24/N;
	//RungeKuttaMond(Feld_2Massen ,  Mstart, Pstart , N, h);
	RungeKutta(Feld_grav,Pstart , N, h);


return 0;
}

