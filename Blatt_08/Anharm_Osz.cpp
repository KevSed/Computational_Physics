/*
 * Anharm_Osz.cpp
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

double E_n(int N, double L, int n){
	double delta_E = L/N;
	double E_n=-L+n*delta_E;
	return E_n;
}

MatrixXd Ham(double lambda, int N, double L){
	MatrixXd H(2*N+1,2*N+1);
	H.setZero();
	double delta_E = L/N;
	for(int n=0; n<2*N+1;n++){
		for(int m=0; m<2*N+1;m++){
				if(n==m-1 or n==m+1){
					H(n,m)=-1/(delta_E*delta_E);
				}
				if(n==m){
					H(n,m)=2/(delta_E*delta_E);
					double En=E_n(N,L,n);
					H(n,m)=H(n,m)+En*En+lambda*En*En*En*En;
				}
			}
	}
	return H;
}

VectorXd compute_EV(Ref<MatrixXd> H_ref, int Anzahl_kleine_Ew){
	int Anzahl_EW=H_ref.rows();
	VectorXd EW(Anzahl_EW);
	EW.setZero();
	VectorXd Ausgabe(Anzahl_kleine_Ew);

	if(H_ref.rows()!=H_ref.cols()){
		cout << "Matrix ist nicht quadratisch und somit nicht diagonalisierbar"<< endl;
		return EW;
	}
	MatrixXd H(Anzahl_EW,Anzahl_EW);
	H=H_ref;

	EW=H.eigenvalues().real();
	// Vorimplementiert: Sortiert Vektorelemente in aufsteigender Groesse
	sort(EW.data(),EW.data()+EW.size());

	for (int i=0;i<Anzahl_kleine_Ew;i++){
		Ausgabe[i]=EW[i];
	}
	return Ausgabe;
}

MatrixXd Ep(double lambda, int N, double L){
	MatrixXd H(2*N+1,2*N+1);
	H.setZero();
	for(int n=0; n<2*N+1;n++){
		for(int m=0; m<2*N+1;m++){
			if(n==m-4){
				H(n,m)=lambda*sqrt(m*(m-1)*(m-2)*(m-3));
			}
			if(n==m+4){
				H(n,m)=lambda*sqrt((m+1)*(m+2)*(m+3)*(m+4));
			}
			if(n==m-2){
				H(n,m)=lambda*sqrt(m*(m-1))*(4*m-2);
			}
			if(n==m+2){
				H(n,m)=lambda*sqrt((m+1)*(m+2))*(4*m+6);
			}
			if(n==m){
				H(n,m)=lambda*(6*m*m+6*m+3);
				H(n,m)=H(n,m)+8*m+4;
			}
			}
	}
	return H/4;
}


int main(){
	int L=10;
	double delta_E=0.1;
	double lambda=0.2;
	int N= L/delta_E;
	cout << "N= " << N << endl << endl;

	cout << "Eigenwerte fuer lamda= " << setprecision(15)<< lambda << endl << endl;
	MatrixXd H(2*N+1,2*N+1);
	H=Ham(lambda, N, L);
	cout << compute_EV(H,10) << endl;

	lambda=0;
	cout << "Eigenwerte fuer lamda= " << lambda << " und N= " << 2*N+1 << endl << endl;
	H=Ham(lambda, N, L);
	cout << compute_EV(H,10) << endl;

	lambda=0.2;
	cout << "Eigenwerte fuer lamda= " << lambda << " und N= " << 2*N+1 << endl << endl;
	MatrixXd H_besetzung(2*N+1,2*N+1);
	H_besetzung=Ep(lambda, N, L);
	VectorXd EigenW(10);
	EigenW=compute_EV(H_besetzung,10);
	cout << EigenW << endl;

	cout << "Differenzen:" << endl << endl;

	for(int i=0; i<9; i++){
	cout << EigenW[i+1]-EigenW[i] << endl;
	}
	cout << endl << endl;

	lambda=0;
	cout << "Eigenwerte fuer lamda= " << lambda << " und N= " << 2*N+1 << endl << endl;
	H_besetzung=Ep(lambda, N, L);
	EigenW=compute_EV(H_besetzung,10);
	cout << EigenW << endl;

	cout << "Differenzen:" << endl << endl;

	for(int i=0; i<9; i++){
	cout << EigenW[i+1]-EigenW[i] << endl;
	}
	cout << endl << endl;

	lambda=0.2;
	N=25; //entspricht dem N=50 aus der Aufgabenstellung
	cout << "Eigenwerte fuer lamda= " << lambda << " und N= " << 2*N+1 << endl << endl;
	H_besetzung=Ep(lambda, N, L);
	cout << compute_EV(H_besetzung,10) << endl;

	lambda=0;
	N=25; //entspricht dem N=50 aus der Aufgabenstellung
	cout << "Eigenwerte fuer lamda= " << lambda << " und N= " << 2*N+1 << endl << endl;
	H_besetzung=Ep(lambda, N, L);
	cout << compute_EV(H_besetzung,10) << endl;


	cout << "Ende der main-Funktion" << endl;
	return 0;
}

