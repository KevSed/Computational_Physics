/*
 * Minimierung_2D.cpp
 *
 *  Created on: 27.06.2017
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


double f(Vector2d x_vec){
	double f;
	double x=x_vec[0];
	double y=x_vec[1];

	f=exp(-10*(x*y-3)*(x*y-3));
	f=f/(x*x+y*y);
	f=1+f;
	f=1/f;
	return f;
}

double g(Vector2d x_vec){
	double G;
	double x=x_vec[0];
	double y=x_vec[1];
	G=x*x+y*y+x*x*y*y;
	return G;
}

MatrixXd Intervallhalbierung_step(Vector2d x, Vector2d y, Vector2d z, double (*F)(Vector2d)){
	MatrixXd Spaltenvektoren(2,3);
	Spaltenvektoren.setZero();
	Vector2d u;
	Vector2d diff_xy=y-x;
	double norm_diff_xy=diff_xy.norm();
	Vector2d diff_yz=z-y;
	double norm_diff_yz=diff_yz.norm();

	if(norm_diff_xy >= norm_diff_yz){
		u=x+0.5*diff_xy;
		if(F(u)>=F(y)){
			x=u;
		}
		else{
			z=y;
			y=u;
		}
	}
	else{
		u=y+0.5*diff_yz;
		if(F(u)>=F(y)){
			z=u;
		}
		else{
			x=y;
			y=u;
		}

	}
	Spaltenvektoren.col(0)=x;
	Spaltenvektoren.col(1)=y;
	Spaltenvektoren.col(2)=z;

	return Spaltenvektoren;
}

Vector2d Intervallhalbierung(Vector2d x, Vector2d p, double (*F)(Vector2d)){
	Vector2d Minimum;
	Minimum.setZero();
	Vector2d y_0;
	Vector2d z_0;
	Vector2d x_0=x;
	Vector2d temp;
	Vector2d y_found;
	Vector2d z_found;
	Vector2d x_found;

	int counter=1;
	bool asymptotik_rechts=false;
	bool asymptotik_links=false;
	bool found=false;
	double norm=p.norm();
	p=p/norm;
	do{
		y_0=x_0+counter*p;
		counter=counter+1;
	}while(F(y_0)==F(x_0));
	counter=1;

	//-------------------------------------------nach rechts gehts bergab
	label_asymptotik:
	if(F(y_0)<F(x_0)){
		do{
			z_0=y_0+p;
			if(F(z_0)>F(y_0)){
				y_found=y_0;
				z_found=z_0;
				x_found=x_0;
				found=true;
				break;
			}
			if(F(z_0)<F(y_0)){
				x_0=y_0;
				y_0=z_0;
				counter++;
				if(counter==1000){
					asymptotik_rechts=true;
					cout << "In positive Richtung geht die Funktion von oben wahrscheinlich";
					cout << " asymptotisch gegen minus unendlich oder eine Konstante" << endl;
					x_0=x;
				}
			}
		}while(counter<1000);
	}

	counter=0;
	if(asymptotik_rechts==true){
		x_0=x;
		do{
			y_0=x_0-p;
			if(F(y_0)<F(x_0)){
				temp=x_0;
				x_0=y_0;
				y_0=temp;
				break;
			}
			if(F(y_0)>F(x_0)){
				x_0=y_0;
				counter++;
			}

			if(counter==1000){
				cout << "Die Funktion hat auf dem intervall (x+1000p , x-1000p)";
				cout << " kein Minimum das breiter ist als p" << endl;
				Minimum.setZero();
				return Minimum;
			}
		}while(counter<1000);
	}

	//-------------------------------------------nach links gehts bergab
	if(F(y_0)>F(x_0)){
		temp=x_0;
		x_0=y_0;
		y_0=temp;
		do{
			z_0=y_0-p;
			if(F(z_0)>F(y_0)){
				x_found=z_0;
				y_found=y_0;
				z_found=x_0;
				found=true;
				break;
			}
			if(F(z_0)<F(y_0)){
				x_0=y_0;
				y_0=z_0;
				counter++;
				if(counter==1000){
					asymptotik_links=true;
					cout << "In negative Richtung geht die Funktion von oben wahrscheinlich";
					cout << " asymptotisch gegen minus unendlich oder eine Konstante" << endl;
					x_0=x;
				}
			}

		}while(counter<1000);
	}

	counter=0;
	if(asymptotik_links==true){
		x_0=x;
		do{
			y_0=x_0+p;
			if(F(y_0)<F(x_0)){
				cout << "goto-Fall tritt ein" << endl;
				goto label_asymptotik;
			}
			if(F(y_0)>F(x_0)){
				x_0=y_0;
				counter++;
			}

			if(counter==1000){
				cout << "Die Funktion hat auf dem intervall (x+1000p , x-1000p)";
				cout << " kein Minimum das breiter ist als p" << endl;
				Minimum.setZero();
				return Minimum;
			}
		}while(counter<1000);
	}

	if(found==true){
		/*
		cout << "Startwerte konnten gefunden werden!" << endl << endl;
		cout << "x ist: " << endl << x_found << endl;
		cout << "F(x) ist: " << F(x_found) << endl << endl;
		cout << "y ist: " << endl << y_found << endl;
		cout << "F(y) ist: " << F(y_found) << endl << endl;
		cout << "z ist: " << endl << z_found << endl;
		cout << "F(z) ist: " << F(z_found) << endl << endl;
		*/
	}
	else{
		cout << "es konnten keine Startwerte gefunden werden, wahrscheinlich komische Funktionseigenschaften" << endl;
		cout << "x ist: " << endl << x_0 << endl;
		cout << "p ist: " << endl << p << endl;
		Minimum.setZero();
		return Minimum;
	}
	//---------------------------Die Startwerte sind nun gefunden, hier beginnt die suche nach dem Minimum

	Vector2d x1=x_found;
	Vector2d y1=y_found;
	Vector2d z1=z_found;

	double eps=0.000000000001;
	MatrixXd Step(2,3);
	while(abs(F(x1)-F(y1))>eps){
		//cout << "while Schleife" << endl;
		Step=Intervallhalbierung_step(x1, y1, z1, F);
		x1=Step.col(0);
		y1=Step.col(1);
		z1=Step.col(2);
	}
	Minimum=y1;
	/*
	cout << endl << "Das Minimum liegt bei x= " << endl << y1 << endl << "und betraegt F(x)= " << F(y1) << endl;
	cout << endl << endl;

	cout << "test" << endl;
	cout << "x ist: " << endl << x1 << endl;
	cout << "F(x) ist: " << F(x1) << endl << endl;
	cout << "y ist: " << endl << y1 << endl;
	cout << "F(y) ist: " << F(y1) << endl << endl;
	cout << "z ist: " << endl << z1 << endl;
	cout << "F(z) ist: " << F(z1) << endl << endl;
	*/

	return Minimum;
}

Vector2d Steepest_Descent (Vector2d x0, double (*F)(Vector2d)){
	//------------------Gradienten berechnen
	int counter=1;
	Vector2d x=x0;
	double min;
	double min_old;
	Vector2d p;
	p.setZero();
	double h=0.000001;
	double eps=0.000000000001;
	Vector2d ex;
	Vector2d ey;
	ex << 1,0;
	ey << 0,1;

	do{
	counter++;
	p[0]=-F(x+h*ex)+F(x-h*ex);
	p[0]=p[0]/(2*h);
	p[1]=-F(x+h*ey)+F(x-h*ey);
	p[1]=p[1]/(2*h);
	//cout << "p ist: " << endl << p << endl;
	min_old=F(x);
	x=Intervallhalbierung(x, p, F);
	min=F(x);
	}while(abs(min_old-min)>eps);

	cout << endl << "Das Minimum liegt bei x= " << endl << x;
	cout << endl << "und betraegt F(x)= " << F(x) << endl;
	cout << endl << "Konvergenz wurde in  " << counter << " Schritten erreicht" << endl;

	return x;

}

Vector2d KonjugierteGradienten (Vector2d x0, double (*F)(Vector2d)){

	int counter=1;
	Vector2d x=x0;
	Vector2d p;
	p.setZero();
	Vector2d p_old;
	p_old.setZero();
	double h=0.000001;
	double mu;
	double eps=0.000000000001;
	double min;
	double min_old;
	Vector2d ex;
	Vector2d ey;
	ex << 1,0;
	ey << 0,1;

	p[0]=-F(x+h*ex)+F(x-h*ex);
	p[0]=p[0]/(2*h);
	p[1]=-F(x+h*ey)+F(x-h*ey);
	p[1]=p[1]/(2*h);
	//cout << "p ist " << endl << p << endl;

	do{
	counter++;
	p_old=p;
	min_old=F(x);
	x=Intervallhalbierung(x, p, F);
	min=F(x);

	p[0]=-F(x+h*ex)+F(x-h*ex);
	p[0]=p[0]/(2*h);
	p[1]=-F(x+h*ey)+F(x-h*ey);
	p[1]=p[1]/(2*h);

	mu=p.dot(p);
	mu=mu/(p_old.norm()*p_old.norm());
	//cout << endl << "mu= " << mu << endl;

	p=p+mu*p_old;

	}while(abs(min_old-min)>eps);

	cout << endl << "Das Minimum liegt bei x= " << endl << x;
	cout << endl << "und betraegt F(x)= " << F(x) << endl;
	cout << endl << "Konvergenz wurde in  " << counter << " Schritten erreicht" << endl;

	return x;

}

int main(){
	Vector2d x;
	x << 1,0;
	Vector2d p;
	p << -2,0;
	//cout << x_vec << "\t" << g(x_vec) << endl;

	cout << "Aufgabenteil a)" << endl << endl << "Erster Startwert x=(1,0): " << endl << endl;Intervallhalbierung(x, p, g);
	cout << endl << endl << "Zweiter Startwert x=(1,1): " << endl << endl;
	x << 1,1;
	p << -2,0;

	Intervallhalbierung(x, p, g);

	cout << "Aufgabenteil b) (SD)" << endl << endl << "Erster Startwert x=(1.5 , 1): " << endl << endl;
	x << 1.5,1;
	Steepest_Descent (x, g);
	cout  << endl << endl << "Zweiter Startwert x=(100 , 80): " << endl << endl;
	x << 100,80;
	Steepest_Descent (x, g);
	cout  << endl << endl << "Dritter Startwert x=(100 , 100): " << endl << endl;
	x << 100,100;
	Steepest_Descent (x, g);

	cout << "Aufgabenteil c) (CG)" << endl << endl << "Erster Startwert x=(1.5 , 1): " << endl << endl;
	x << 1.5,1;
	KonjugierteGradienten (x, g);
	cout  << endl << endl << "Zweiter Startwert x=(100 , 80): " << endl << endl;
	x << 100,80;
	KonjugierteGradienten (x, g);
	cout  << endl << endl << "Zweiter Startwert x=(100 , 100): " << endl << endl;
	x << 100,100;
	KonjugierteGradienten (x, g);

	cout << "Aufgabenteil d) (SD)" << endl << endl << "Erster Startwert x=(1.5 , 2.3): " << endl << endl;
	x << 1.5,2.3;
	Steepest_Descent (x, f);
	cout  << endl << endl << "Zweiter Startwert x=(-1.7 , 1.9): " << endl << endl;
	x << -1.7,1.9;
	Steepest_Descent (x, f);
	cout  << endl << endl << "Dritter Startwert x=(0.5 , 0.6): " << endl << endl;
	x << 0.5,0.6;
	Steepest_Descent (x, f);

	cout << "Aufgabenteil d) (cD)" << endl << endl << "Erster Startwert x=(1.5 , 2.3): " << endl << endl;
	x << 1.5,2.3;
	KonjugierteGradienten (x, f);
	cout  << endl << endl << "Zweiter Startwert x=(-1.7 , 1.9): " << endl << endl;
	x << -1.7,1.9;
	KonjugierteGradienten (x, f);
	cout  << endl << endl << "Dritter Startwert x=(0.5 , 0.6): " << endl << endl;
	x << 0.5,0.6;
	KonjugierteGradienten (x, f);


	cout << "Ende der main-Funktion" << endl;
	return 0;
}
