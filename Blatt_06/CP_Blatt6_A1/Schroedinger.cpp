/*
 * Schroedinger.cpp
 *
 *  Created on: 06.06.2017
 *      Author: mona
 */

#include <iostream>
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

const complex<double> i(0, 1);
const complex<double> eins(1, 0);

double Epsillon(int j, double dE, double epsmin){
	double epsillon=j*dE;
	epsillon=epsillon+epsmin;
	return epsillon;
}

MatrixXcd Hamiltonian(double V, int dim, double b, double dE, double epsmin){
	MatrixXcd Ham(dim,dim);
	Ham.setZero();
	double frac=1./dE;
	frac=(-1)*frac*frac;

	for(int n=0; n<201;n++){
		for(int m=0; m<201;m++){
			if(n==m-1){Ham(n,m)=Ham(n,m)+frac;}
			if(n==m+1){Ham(n,m)=Ham(n,m)+frac;}
			if(n==m){
				Ham(n,m)=Ham(n,m)-2*frac;
				if((0.5*b-abs(Epsillon(n, dE, epsmin)))>0){Ham(n,m)=Ham(n,m)+V;}
			}
		}
	}


	return Ham;
}

MatrixXcd Zeitentwicklungsop(double V, int dim, MatrixXcd H, double dtau){
	MatrixXcd M(dim, dim);
	MatrixXcd Einheitsmatrix(dim, dim);
	Einheitsmatrix=MatrixXcd::Identity(dim,dim);
	M.setZero();

	for(int n=0; n<201;n++){
		for(int m=0; m<201;m++){
			M(n,m) = 0.5 * i * H(n, m) * dtau;
		}
	}

	M=(Einheitsmatrix + M).inverse()*(Einheitsmatrix - M);

	return M;
}

VectorXcd Psi(double sigma, double epsnull, double knull, int dim, double dE, double epsmin){
	VectorXcd psi(dim);
	psi.setZero();
	double faktor=2*M_PI*sigma;
	faktor=pow(faktor, 0.25);
	faktor=1./faktor;

	for(int j=0; j<dim;j++){
		psi(j)=exp(-(Epsillon(j,dE,epsmin)-epsnull)*(Epsillon(j,dE,epsmin)-epsnull)/(4*sigma));
		psi(j)=psi(j)*exp(i*Epsillon(j,dE,epsmin)*knull);
	}
	return faktor*psi;
}

VectorXcd Zeitschritt(VectorXcd psi, MatrixXcd Sh){
	VectorXcd psi_nplus1=Sh*psi;

	return psi_nplus1;
}

void Zeitentwicklung(VectorXcd psi, MatrixXcd Sh, double Tau, double dTau, string Name, string trans_name, int dim, double dE, double epsmin){
	VectorXcd psi_entw=psi;
	int counter=0;
	// ===============================
	ofstream Data;
	Data.open(Name.c_str());
	Data.precision(10);//
	// ===============================
	// ===============================
	ofstream Datatrans;
	Datatrans.open(trans_name.c_str());
	Datatrans.precision(10);//
	// ===============================
	double q=Tau/dTau;

	for(int f=0; f<q; f++){
		double Transkoeff=0;
		psi_entw=Zeitschritt(psi_entw, Sh);
		counter++;
		for(int j=0;j<dim;j++){
			if(Epsillon(j,dE,epsmin)>0){
				Transkoeff=Transkoeff+dE*abs(psi_entw(j))*abs(psi_entw(j));
			}
		}
		Datatrans << f*dTau << "\t" << Transkoeff << endl;
	}
	Datatrans.close();

	for(int j=0;j<dim;j++){
		Data << Epsillon(j,dE,epsmin) << "\t" << abs(psi_entw(j))*abs(psi_entw(j)) << endl;
	}
	cout << Name << "  counter:  " << counter << endl;
	Data.close();
}


int main(){
	int dim=201;
	double epsmin=-10;
	double b=1;
	double V=10;
	double dTau=0.01;
	double dE=0.1;
	double sigma=1;
	double epsnull=-5.;
	double knull=5.5;

	MatrixXcd H=Hamiltonian(V, dim, b, dE, epsmin);
	MatrixXcd SH=Zeitentwicklungsop(V, dim, H, dTau);

	// ===============================
		ofstream psidata;
		psidata.open("Anfangszustand.txt");
		psidata.precision(10);//
	// ===============================

	VectorXcd psi=Psi(sigma, epsnull, knull, dim, dE, epsmin);

	for(int j=0;j<dim;j++){
		psidata << j << "\t" << abs(psi(j))*abs(psi(j)) << endl;
	}

	psidata.close();

	double Tau=1.;
	V=0;
	H=Hamiltonian(V, dim, b, dE, epsmin);
	SH=Zeitentwicklungsop(V, dim, H, dTau);
	Zeitentwicklung(psi, SH, Tau, dTau, "V=0", "V=0_trans", dim, dE, epsmin);
	V=10;
	H=Hamiltonian(V, dim, b, dE, epsmin);
	SH=Zeitentwicklungsop(V, dim, H, dTau);
	Zeitentwicklung(psi, SH, Tau, dTau, "V=10", "V=10_trans",dim, dE, epsmin);
	V=30;
	H=Hamiltonian(V, dim, b, dE, epsmin);
	SH=Zeitentwicklungsop(V, dim, H, dTau);
	Zeitentwicklung(psi, SH, Tau, dTau, "V=30", "V=30_trans",dim, dE, epsmin);
	V=50;
	H=Hamiltonian(V, dim, b, dE, epsmin);
	SH=Zeitentwicklungsop(V, dim, H, dTau);
	Zeitentwicklung(psi, SH, Tau, dTau, "V=50", "V=50_trans",dim, dE, epsmin);


	cout << "Hamiltonian:  " << H(100,102);

	cout << "Ende der main-Funktion" << endl;
	return 0;
}


