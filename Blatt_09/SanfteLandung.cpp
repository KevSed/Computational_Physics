/*
 * SanfteLandung.cpp
 *
 *  Created on: 29.06.2017
 *      Author: mona
 */

#include <iostream>
#include <iomanip>
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

/*
Globale Konstanten.
Achtung, Laengen muessen in Metern, Zeiten in Sekunden und Massen in kg angegeben werden!
*/

const double A = 10.;
const double m_R = 10000.;
const double m_E = 5.97e24;
const double G = 6.67e-11;
const double R_ERDE=6370000;
const double h_0 = R_ERDE + 130000;


//Luftdichte

double rho(double h){
	double p;
	if(h<0){
		//cout << "Achtung, rho liefert nicht mehr die richtigen Ergebnisse!" << endl << endl;
		return 1.2;
	}
	p=1.2*exp((-1.17e-4)*h);
	return p;
}

Vector2d kraft(Vector4d y, double C){

	//Achtung: Hier wird in Wahrheit nicht F sondern F/m_R berechnet

	Vector2d F_L;
	Vector2d F_G;

	Vector2d r;
	Vector2d v;
	r << y[0],y[1];
	v << y[2],y[3];

	F_G=-r;
	F_G=F_G*(G*m_E);
	F_G=F_G/(r.norm()*r.norm()*r.norm());

	F_L=-v;
	F_L=F_L*(C*A*v.norm());
	F_L=F_L/(2*m_R);

	double Hoehe=r.norm()-R_ERDE;

	F_L=F_L*rho(Hoehe);

	return F_L+F_G;
}

//****************************************
// Ein Schritt des Runge-Kutta Verfahrens 4. Ordnung
//****************************************


Vector4d Schritt(Vector2d (*F)(Vector4d, double), Vector4d y0, double h, double C){


	Vector4d f;
	Vector4d k2;
	Vector4d k3;
	Vector4d k4;
	Vector2d Kraft;


	for(int i=0; i<2; i++){
		f[i]=h*y0[i+2];
		f[i+2]=h*F(y0,C)[i];
	}

	Vector4d ytemp;
	Vector4d y;

	// Berechnung der Hilfsgroessen k

	ytemp=y0+0.5*f;
	for(int i=0; i<2; i++){
			k2[i]=h*ytemp[i+2];
			k2[i+2]=h*F(ytemp,C)[i];
	}

	ytemp=y0+0.5*k2;
	for(int i=0; i<2; i++){
			k3[i]=h*ytemp[i+2];
			k3[i+2]=h*F(ytemp,C)[i];
	}

	ytemp=y0+k3;
	for(int i=0; i<2; i++){
			k4[i]=h*ytemp[i+2];
			k4[i+2]=h*F(ytemp,C)[i];
	}

	// Berechnung des neuen y
	y=y0+(1./6.)*(f+2*k2+2*k3+k4);

	return y;
}

//****************************************
// Ausfuehren des Runge-Kutta Verfahrens 4. Ordnung
// Anders als zuvor gibt die Funktion den Endvektor zurueck
//****************************************


Vector4d RungeKutta(Vector2d (*F)(Vector4d, double), Vector4d ystart, double C, double h){

	Vector4d y;
	Vector2d r;
	Vector2d v;

	y=ystart;

	// ===============================
	ofstream Data;
	Data.open("RK.txt");
	Data.precision(10);//
	// ===============================

	int i=0;
	do{
		Data << i*h << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << y[3] << "\n";
		y=Schritt(F,y,h,C);
		r << y[0],y[1];
		v << y[2],y[3];
		i++;
	}while(r.norm()>R_ERDE);

	Data.close();
	return y;
}

double v_vert(double v_0, double C, double h){
	Vector4d ystart;
	Vector4d yend;
	Vector2d rend;
	Vector2d vend;
	ystart << h_0,0,0.,v_0;

	yend=RungeKutta(kraft,ystart,C,h);
	rend << yend[0],yend[1];
	vend << yend[2],yend[3];

	//Ortogonalprojektion von v auf r, wobei r normiert wird
	double vvert=(rend.dot(vend))/rend.norm();

	// Zum Suchen der Minima:
	return abs(vvert);
	// Zum Suchen der Maxima:
	//return -abs(vvert);
}

double v_norm(double v_0, double C, double h){
	Vector4d ystart;
	Vector4d yend;
	Vector2d rend;
	Vector2d vend;
	ystart << h_0,0,0.,v_0;

	yend=RungeKutta(kraft,ystart,C,h);
	vend << yend[2],yend[3];

	// Zum Suchen der Minima:
	return vend.norm();
	// Zum Suchen der Maxima:
	//return -vend.norm();
}

Vector3d Intervallhalbierung_step(double x, double y, double z, double (*F)(double, double, double), double C, double h){
	Vector3d v;
	double u;
	double diff_xy=abs(y-x);
	double diff_yz=abs(z-y);

	if(diff_xy >= diff_yz){
		u=x+0.5*diff_xy;
		if(F(u,C,h)>=F(y,C,h)){
			x=u;
		}
		else{
			z=y;
			y=u;
		}
	}
	else{
		u=y+0.5*diff_yz;
		if(F(u,C,h)>=F(y,C,h)){
			z=u;
		}
		else{
			x=y;
			y=u;
		}

	}
	v[0]=x;
	v[1]=y;
	v[2]=z;
	//cout << v << endl;
	return v;
}

double Intervallhalbierung(double (*Funk)(double, double, double), double h, double C){

	cout << "Start: " << endl << endl;
	double Minimum=0;
	/*
	 Da wir wissen wo das Intervall seinen Rand hat
	 koennen wir diese Werte gezielt als Randpunkte des Dreieks nehmen.
	 Wenn die Funktion ein Minimum hat, welches nicht in einem der Randpunkte liegt,
	 so muss es einen Punkt y zwischen ihnen geben, dessen Funktionswert
	 kleiner ist als der der Randwerte. Es muss also nur ein dritter Startpunkt gewucht werden
	 */

	double x_0=0;
	double z_0=7800;
	double y_0;

	int counter=10;
	bool found=false;

	do{
		y_0=z_0-counter;
		if(y_0<0){
			cout << "mit dieser Schrittweite konnte kein Minimum gefunden werden" << endl;
			cout << "es wird daher der kleinere Randwert zurueckgegeben" << endl << endl;
			if(Funk(x_0,C,h)<Funk(z_0,C,h)){Minimum = x_0;}
			else{Minimum= z_0;}
			cout << endl << "Das Minimum liegt bei x= " << endl << Minimum;
			cout << endl << "und betraegt F(x)= " << Funk(Minimum,C,h) << endl;
			cout << endl << endl;
			return Minimum;
		}
		if(Funk(y_0,C,h)<Funk(x_0,C,h) && Funk(y_0,C,h)<Funk(z_0,C,h)){found=true;}
		counter=counter+100;
		//cout << "y= " << "\t" << y_0 << endl;
	}while(found==false);


	//---------------------------Die Startwerte sind nun gefunden, hier beginnt die suche nach dem Minimum

	double x=x_0;
	double y=y_0;
	double z=z_0;

	//cout << "Startvektor:" << "\t" << x << "\t" << y << "\t" << z << endl << endl;

	double eps=0.00000001;
	Vector3d Step;
	while(abs(Funk(x,C,h)-Funk(y,C,h))>eps){
		Step=Intervallhalbierung_step(x, y, z, Funk, C, h);
		x=Step[0];
		y=Step[1];
		z=Step[2];
	}
	Minimum=y;

	cout << endl << "Das Minimum liegt bei x= " << endl << y << endl << "und betraegt F(x)= " << Funk(y,C,h) << endl;
	cout << endl << endl;

	return Minimum;
}


int main(){
	double C1=0.2;
	double C2=0.1;
	//double h=0.01;
	double h=0.001;

	// ===============================
		ofstream Datav;
		Datav.open("vvert.txt");
		Datav.precision(10);//
	// ===============================
/*
 	 	for(double v=7800; v>=0; v=v-10){
		Datav << v << "\t" << v_vert(v,C1,h) << "\t" << v_norm(v,C1,h) << "\t" << v_vert(v,C2,h) << "\t" << v_norm(v,C2,h) << endl;
	}
*/
	//Suchen der Minima
	Intervallhalbierung(v_vert, h, C1);
	Intervallhalbierung(v_norm, h, C1);
	Intervallhalbierung(v_vert, h, C2);
	Intervallhalbierung(v_norm, h, C2);



	cout << "Ende der main-Funktion" << endl;
	return 0;
}



