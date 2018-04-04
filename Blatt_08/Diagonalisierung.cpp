/*
 * Diagonalisierung.cpp
 *
 *  Created on: 19.06.2017
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

VectorXd Potenzmethode(Ref<MatrixXd> A_ref, double eps){

	MatrixXd A(A_ref.rows(),A_ref.cols());
	MatrixXd B(A_ref.rows(),A_ref.cols());
	A=A_ref;
	int Anzahl_EW=A.rows();
	VectorXd EW(Anzahl_EW);
	EW.setZero();
	if(A.rows()!=A.cols()){
		cout << "Matrix ist nicht quadratisch und somit nicht diagonalisierbar"<< endl;
		return EW;
	}

	VectorXd v0(Anzahl_EW);
	v0.setZero();
	VectorXd diff(Anzahl_EW);
	double norm_w;
	double norm_diff;
	double lambda;

	//Erstellen eines zufaelligen, normierten Startvektors
	for(int i=0; i<Anzahl_EW;i++){
		v0[i]=rand()* (1./RAND_MAX);
	}
	v0=v0/(v0.norm());
	VectorXd v(Anzahl_EW);
	VectorXd w(Anzahl_EW);

	//Berechnung der Eigenwerte
	for(int i=0; i<Anzahl_EW;i++){
		//cout << "A:" << endl << A << endl;
		v=v0;

		do{
			w=A*v;
			norm_w=w.norm();
			w=w/norm_w;

			//Es kann passieren dass v und w im Vorzeichen alternieren
			//wenn der Eigenwert negantiv ist.
			for(int i=0; i<Anzahl_EW;i++){
				diff[i]=abs(v[i])-abs(w[i]);
			}

			norm_diff=diff.norm();
			v=w;
			//cout << "norm " << endl << norm_diff << endl;
		}while(norm_diff>eps);

		//Berechnung von v^T * A * V
		w=A*v;
		lambda=v.dot(w);
		EW[i]=lambda;

		B=v*v.transpose();
		//cout << "B:" << endl << B << endl;
		A=A-lambda*B;
		//cout << "A:" << endl << A << endl;
	}

	return EW;
}


int main(){

	//Aufgabenteil a)-----------------------------------------------------------------

	Matrix4d A; //Erstellen einer 4x4 Matrix A
	A << 	 1,-2, -3, 4,
			-2, 2, -1, 7,
			-3,-1,  3, 6,
			 4, 7,  6, 4;

	//Berechnung der Eigenwerte mit vorimplementierter Funktion:

	Vector4d A_EW; //Erstellen des Vektors fuer die Eigenwerte

	A_EW=A.eigenvalues().real();

	cout << "Die Eigenwerte der Matrix A berechnet mit Eigen sind" << endl;
	cout << setprecision(10) << A_EW << endl;

	//Aufgabenteil b))-----------------------------------------------------------------

	VectorXd EW_pot(4);
	EW_pot=Potenzmethode(A, 0.0001);
	cout << endl << "Die Eigenwerte der Matrix A berechnet mit der Potenzmethode sind" << endl;
	cout << EW_pot << endl;
	// Die Eigenwerte stimmen bis auf die 6. Nachkommastelle ueberein.
	cout << "Ende der main-Funktion" << endl;
	return 0;
}

