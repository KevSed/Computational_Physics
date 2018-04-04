/*
 * Poisson.cpp
 *
 *  Created on: 07.06.2017
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

MatrixXd GaussSeidelStep(double rechter_rand, double diskretisierung, MatrixXd Phi, MatrixXd rho){
	int dim= rechter_rand/diskretisierung;
	MatrixXd M(dim,dim);
	M=Phi;
	for(int n=1;n<dim-1;n++){
			for(int m=1; m<dim-1; m++){
				M(n,m)=0.25*(M(n+1,m)+M(n-1,m)+M(n,m+1)+M(n,m-1)+diskretisierung*diskretisierung*rho(n,m));
			}
		}
	return M;
}

MatrixXd rho_disk(MatrixXd Ladungen, VectorXd qi, double rechter_rand, double diskretisierung){
	int dim= rechter_rand/diskretisierung;
	int AnzahlLadungen=Ladungen.cols();
	MatrixXd M(dim,dim);
	M.setZero();
	if(Ladungen.rows()!= 2){
		cout << "keine Ladungen vorhanden" << endl;
		return M;
	}
	for(int n=0;n<dim;n++){
		for(int m=0; m<dim; m++){
			for(int l=0; l<AnzahlLadungen; l++){
				if((n*diskretisierung)-Ladungen(0,l)==0 && (m*diskretisierung)-Ladungen(1,l)==0){M(n,m)=qi(l);}
			}
		}
	}
	return M;
}

MatrixXd GaussSeidel(double rechter_rand, double diskretisierung, MatrixXd Phi, MatrixXd rho, int N){
	int dim= rechter_rand/diskretisierung;
	MatrixXd M(dim,dim);
	M=Phi;
	for(int n=0;n<N;n++){
		M=GaussSeidelStep(rechter_rand, diskretisierung, M, rho);
	}

	return M;
}

MatrixXd elektrischesFeld(double rechter_rand, double diskretisierung, MatrixXd Phi){
	int dim= rechter_rand/diskretisierung;
	MatrixXd betragE(dim,dim);
	betragE.setZero();
	double Ex;
	double Ey;
	for(int n=1;n<dim-1;n++){
			for(int m=1; m<dim-1; m++){

				Ex=Phi(n+1,m)-Phi(n-1,m);
				Ex=Ex/(2*diskretisierung);
				Ey=Phi(n,m+1)-Phi(n,m-1);
				Ey=Ey/(2*diskretisierung);
				betragE(n,m)=sqrt(Ex*Ex+Ey*Ey);
			}
	}
	return betragE;

}




int main(){

	double rechter_rand=1.;
	double diskretisierung=0.05;
	int N=500000;

	int dim=rechter_rand/diskretisierung;
	MatrixXd Phi(dim,dim);

	for(int n=1;n<dim-1;n++){
			for(int m=1; m<dim-1; m++){
				Phi(n,m)=1;
			}
			Phi(0,n)=0;
			Phi(n,0)=0;
			Phi(n,dim-1)=0;
			Phi(dim-1,n)=0;
	}

	MatrixXd Ladungen(1,1);
	VectorXd qi(1,1);
	qi(0,0)=0;
	MatrixXd Rho=rho_disk(Ladungen,qi, rechter_rand, diskretisierung);
	MatrixXd Phi_fertig=GaussSeidel(rechter_rand, diskretisierung, Phi, Rho, N);

	cout << "Potential:" << endl;
	cout << Phi_fertig << endl << endl;

	MatrixXd betrE=elektrischesFeld(rechter_rand, diskretisierung, Phi_fertig);

	cout << "E-Feld:" << endl;
	cout << betrE << endl << endl;

	// ===============================
		ofstream phi_a_data;
		phi_a_data.open("phi_a.txt");
		phi_a_data.precision(10);//
	// ===============================
	// ===============================
		ofstream E_a_data;
		E_a_data.open("E_a.txt");
		E_a_data.precision(10);//
	// ===============================

	phi_a_data << Phi_fertig << endl;
	E_a_data << betrE << endl;

//Aufgabenteil b)---------------------------------------------------------------------------

	for(int n=1;n<dim-1;n++){
			for(int m=1; m<dim-1; m++){
				Phi(n,m)=1;
			}
			Phi(0,n)=0;
			Phi(n,0)=0;
			Phi(n,dim-1)=1;
			Phi(dim-1,n)=0;
	}

	MatrixXd Phi_b=GaussSeidel(rechter_rand, diskretisierung, Phi, Rho, N);

	cout << "Anfangs-Potential b:" << endl;
	cout << Phi << endl << endl;

	cout << "Potential b):" << endl;
	cout << Phi_b << endl << endl;

	MatrixXd betrEb=elektrischesFeld(rechter_rand, diskretisierung, Phi_b);

	cout << "E-Feld b):" << endl;
	cout << betrEb << endl << endl;

	// ===============================
		ofstream phi_b_data;
		phi_b_data.open("phi_b.txt");
		phi_b_data.precision(10);//
	// ===============================
	// ===============================
		ofstream E_b_data;
		E_b_data.open("E_b.txt");
		E_b_data.precision(10);//
	// ===============================

	phi_b_data << Phi_b << endl;
	E_b_data << betrEb << endl;

//Aufgabenteil c)---------------------------------------------------------------------------

	for(int n=1;n<dim-1;n++){
			for(int m=1; m<dim-1; m++){
				Phi(n,m)=1;
			}
			Phi(0,n)=0;
			Phi(n,0)=0;
			Phi(n,dim-1)=0;
			Phi(dim-1,n)=0;
	}
	MatrixXd Ladungen_c(2,1);
	Ladungen_c(0,0)=0.5;
	Ladungen_c(1,0)=0.5;
	VectorXd qi_c(1,1);
	qi_c(0,0)=1.;
	MatrixXd Rho_c=rho_disk(Ladungen_c,qi_c, rechter_rand, diskretisierung);

	cout << "Rho_c:" << endl;
	cout << Rho_c << endl << endl;

	MatrixXd Phi_c=GaussSeidel(rechter_rand, diskretisierung, Phi, Rho_c, N);


	cout << "Potential c):" << endl;
	cout << Phi_c << endl << endl;

	MatrixXd betrEc=elektrischesFeld(rechter_rand, diskretisierung, Phi_c);

	cout << "E-Feld c):" << endl;
	cout << betrEc << endl << endl;

	// ===============================
		ofstream phi_c_data;
		phi_c_data.open("phi_c.txt");
		phi_c_data.precision(10);//
	// ===============================
	// ===============================
		ofstream E_c_data;
		E_c_data.open("E_c.txt");
		E_c_data.precision(10);//
	// ===============================

	phi_c_data << Phi_c << endl;
	E_c_data << betrEc << endl;

	N=500001;
	MatrixXd Phi_c_probe=GaussSeidel(rechter_rand, diskretisierung, Phi, Rho_c, N);


	cout << "Probe:" << endl<< endl<< Phi_c_probe-Phi_c << endl;
	cout << "Analytisch....:" << endl<< endl;

//Analytische Loesung-----------------------------------------------------------------------------------------

	MatrixXd Phi_analytisch(dim,dim);
	Phi_analytisch.setZero();

	for(int j=1;j<dim-1;j++){
			for(int k=1; k<dim-1; k++){
				//cout << "(j,k):   " << j << "\t" << k << endl <<endl;
				int counter=0;
				int n=1;
				do{
					double term=2*(1-cos(n*M_PI))*sin(n*M_PI*j*diskretisierung)*sinh(n*M_PI*k*diskretisierung);
					//cout << "term 1:   " << term << endl;
					double nenner=n*M_PI*sinh(n*M_PI);
					//cout << "nenner:   " << nenner << endl;
					term=term/nenner;
					//cout << "term 2:   " << term << endl;
					Phi_analytisch(j,k)=Phi_analytisch(j,k)+term;
					n++;
					if(Phi_analytisch(j,k)!=0){counter++;}
				}while(counter<200);
			}
			Phi_analytisch(j,dim-1)=1;
	}

	cout << "Phi_analytisch:" << endl<< endl<< Phi_analytisch << endl;

	// ===============================
		ofstream phi_ana_data;
		phi_ana_data.open("phi_analytisch.txt");
		phi_ana_data.precision(10);//
	// ===============================

	phi_ana_data << Phi_analytisch << endl;


	cout << "Abweichung: " << Phi_analytisch-Phi_b << endl <<endl;

	cout << "Ende der main-Funktion" << endl;
	return 0;
}




