/*
 * Wettervorhersage.cpp
 *
 *  Created on: 02.05.2017
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
// Definition des Startvektors und des DGL-Systems
//****************************************

vector<double> DGL_start={-10,15,19}; //Vektor der Startwerte

vector<double> f(vector<double> DGL, double sigma, double r, double b){

	double fx=sigma*(-DGL[0]+DGL[1]);
	double fy=-DGL[0]*DGL[2]+r*DGL[0]-DGL[1];
	double fz=DGL[0]*DGL[1]-b*DGL[2];
	DGL={fx,fy,fz};
	return DGL;
}


//****************************************
// Aufgabenteil a) Integration des gleichungssystems mit dem Euler verfahren und Ausgabe von 3 Textdateien,
// welche die werte für N={10,100,1000} ausgibt.
//****************************************

void euler(vector<double> DGL, double h1, double h2, double h3, double N, double sigma, double r, double b){

	vector<double> DGL_1=DGL;
	vector<double> DGL_2=DGL;
	vector<double> DGL_3=DGL;
	// ===============================
	ofstream a;
	a.open("Wetter_h1.txt");
	a.precision(10);//
	// ===============================

	for(double n=1; n<=N; n++){
		vector<double> func1=f(DGL_1,sigma,r,b);
		double x=DGL_1[0]+h1*func1[0];
		double y=DGL_1[1]+h1*func1[1];
		double z=DGL_1[2]+h1*func1[2];
		if(n==10 || n==100 || n==1000){
			a << n << "\t" << h1*n << "\t" << x << "\t" << y << "\t" <<  z << endl;
		}
		DGL_1={x,y,z};
	}
	a.close();
	// ===============================
		ofstream c;
		c.open("Wetter_h2.txt");
		c.precision(10);//
	// ===============================

	for(double n=1; n<=N; n++){
		vector<double> func2=f(DGL_2,sigma,r,b);
		double x=DGL_2[0]+h2*func2[0];
		double y=DGL_2[1]+h2*func2[1];
		double z=DGL_2[2]+h2*func2[2];
		if(n==10 || n==100 || n==1000){
			c << n << "\t" << h2*n << "\t" << x << "\t" << y << "\t" <<  z << endl;
		}
		DGL_2={x,y,z};
	}
	c.close();
	// ===============================
		ofstream d;
		d.open("Wetter_h3.txt");
		d.precision(10);//
	// ===============================

	for(double n=1; n<=N; n++){
		vector<double> func3=f(DGL_3,sigma,r,b);
		double x=DGL_3[0]+h3*func3[0];
		double y=DGL_3[1]+h3*func3[1];
		double z=DGL_3[2]+h3*func3[2];
		if(n==10 || n==100 || n==1000){
			d << n << "\t" << h3*n << "\t" << x << "\t" << y << "\t" <<  z << endl;
		}
		DGL_3={x,y,z};
	}
	d.close();
}

//****************************************
// Aufgabenteil b) und c) Ausgabe von Datensätzen für x und y für verschiedene r und Startwerte
//****************************************


void euler_plot(vector<double> DGL, double h, double N, double sigma, double r1, double r2, double b){

	vector<double> DGL_1=DGL;
	vector<double> DGL_2=DGL;

	// ===============================
	ofstream a;
	a.open("Wetter_plotr20.txt");
	a.precision(10);//
	// ===============================

	for(double n=1; n<=N; n++){
		vector<double> func1=f(DGL_1,sigma,r1,b);
		double x=DGL_1[0]+h*func1[0];
		double y=DGL_1[1]+h*func1[1];
		double z=DGL_1[2]+h*func1[2];
		a << n << "\t" << h*n << "\t" << x << "\t" << y  << endl;
		DGL_1={x,y,z};
	}
	a.close();

	// ===============================
	ofstream c;
	c.open("Wetter_plotr28.txt");
	c.precision(10);//
	// ===============================

	for(double n=1; n<=N; n++){
		vector<double> func2=f(DGL_2,sigma,r2,b);
		double x=DGL_2[0]+h*func2[0];
		double y=DGL_2[1]+h*func2[1];
		double z=DGL_2[2]+h*func2[2];
		c << n << "\t" << h*n << "\t" << x << "\t" << y << endl;
		DGL_2={x,y,z};
	}
	c.close();

	// Aufgabenteil c)
	vector<double> DGL_3={-10.01,15,19};

	// ===============================
	ofstream d;
	d.open("Wetter_plotr20_c.txt");
	d.precision(10);//
	// ===============================

	for(double n=1; n<=N; n++){
		vector<double> func3=f(DGL_3,sigma,r1,b);
		double x=DGL_3[0]+h*func3[0];
		double y=DGL_3[1]+h*func3[1];
		double z=DGL_3[2]+h*func3[2];
		d << n << "\t" << h*n << "\t" << x << "\t" << y << endl;
		DGL_3={x,y,z};
	}
	d.close();


}

//****************************************
// Aufgabenteil d) Numerische Ermittlung des Fixpunktes als letzten Wert
//****************************************

vector<double> euler_Fixpunkt(vector<double> DGL, double h, double N, double sigma, double r, double b){

	vector<double> DGL_1=DGL;

	double x_alt=DGL_1[0];
	double y_alt=DGL_1[1];
	double z_alt=DGL_1[2];
	double x;
	double y;
	double z;

	double variation=0.00001;

	for(double n=1; n<=N; n++){
		vector<double> func1=f(DGL_1,sigma,r,b);
		x=DGL_1[0]+h*func1[0];
		y=DGL_1[1]+h*func1[1];
		z=DGL_1[2]+h*func1[2];
		if(abs(x_alt-x)<variation && abs(y_alt-y)<variation && abs(z_alt-z)<variation){
			cout << "fixpunkt bei x = " << x << " , y = " << y << " , z = " << z << endl;
			break;
		}
		x_alt=x;
		y_alt=y;
		z_alt=z;
		DGL_1={x,y,z};
	}
	vector<double> letzterWert= {x,y,z};
	return letzterWert;
}



int main(){
	double r=20;
	double sigma=10;
	double b= 8./3.;   //Wird sonst behandelt wie ein Integer
	euler(DGL_start,0.05,0.005,0.0005,1000,sigma,r,b);
	euler_plot(DGL_start,0.01,1e5,sigma,20,28,b);

	//Berechnung des analytischen Fixpunktes
	double Fix_20= sqrt(b*(20-1));
	double Fix_28= sqrt(b*(28-1));
	cout << scientific << setprecision(8)<< "Analytisch berechneter Fixpunkt bei r=20: " << "\t" << Fix_20  << endl;
	cout << scientific << setprecision(8)<< "Analytisch berechneter Fixpunkt bei r=28: " << "\t" << Fix_28  << endl << endl << endl;

	// Numerische berechnung der Fixpunkte
	cout << scientific << setprecision(8)<< "Fixpunkt bei r=20, Start=(-10,15,19): " << "\t" ;
	vector<double> letzter=euler_Fixpunkt(DGL_start,0.001,1e5,sigma,20,b);

	cout << "Letzter erreichter Wert: " << "\t" << letzter[0] <<"\t" << letzter[1]<<"\t" << letzter[2] << endl ;


	//Fehler für x
	double x_Fehler=abs(-letzter[0]-Fix_20)/Fix_20;
	double y_Fehler=abs(-letzter[1]-Fix_20)/Fix_20;
	cout << "Fehler fuer x: " << "\t" << x_Fehler <<"\t" << "Fehler fuer y:"<<"\t" << y_Fehler << endl << endl;



	cout << "Fixpunkt bei r=20, Start=(10.01,15,19): " << "\t" ;
	vector<double> DGL_start_neu={10.01,15,19}; //Neuer Vektor der Startwerte
	letzter=euler_Fixpunkt(DGL_start_neu,0.001,1e5,sigma,20,b);
	x_Fehler=abs(letzter[0]-Fix_20)/Fix_20;
	y_Fehler=abs(letzter[1]-Fix_20)/Fix_20;
	cout << "Letzter erreichter Wert: " << "\t" << letzter[0] <<"\t" << letzter[1]<<"\t" << letzter[2] << endl << endl;
	cout << "Fehler fuer x: " << "\t" << x_Fehler <<"\t" << "Fehler fuer y:"<<"\t" << y_Fehler << endl << endl;

	return 0;
}