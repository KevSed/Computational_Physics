/*
 * LinKonGeneratoren.cpp
 *
 *  Created on: 04.07.2017
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
#include <random>

using namespace std;

//Aufgabe 1a. Funktion "Linear kongurente Generatoren"

void lkGenerator (int64_t r_0, int64_t a, int64_t c, int64_t m, int64_t N, string Name){

	// ===============================
	ofstream Data;
	Data.open(Name.c_str());
	Data.precision(10);//
	// ===============================

	int64_t r=r_0;
	for(int64_t n=0; n<N;n++){
		Data << (double)r/m << endl;
		r=(a*r+c)%m;
		if(r<0){cout << "Stop" << endl;}
	}
}

void marsenne (int64_t N, string Name){

	// ===============================
	ofstream Data;
	Data.open(Name.c_str());
	Data.precision(10);//
	// ===============================
	for(int64_t n=0; n<N;n++){
		mt19937 Generator (n);
		uniform_real_distribution<double> Gleichvert(0.0, 1.0);
		double randomNbr = Gleichvert(Generator);
		Data << randomNbr << endl;
	}

}

void marsenne_test (int64_t N, string Name, string Name2){

	// ===============================
	ofstream Data;
	Data.open(Name.c_str());
	Data.precision(10);//
	// ===============================
	for(int64_t n=0; n<N;n=n+2){
		mt19937 Generator1 (n);
		mt19937 Generator2 (n+1);
		uniform_real_distribution<double> Gleichvert(0.0, 1.0);
		double randomNbr1 = Gleichvert(Generator1);
		double randomNbr2 = Gleichvert(Generator2);
		Data << randomNbr1 << "\t" << randomNbr2 << endl;
	}

	// ===============================
	ofstream Data2;
	Data2.open(Name2.c_str());
	Data2.precision(10);//
	// ===============================
	for(int64_t n=0; n<N;n=n+3){
		mt19937 Generator1 (n);
		mt19937 Generator2 (n+1);
		mt19937 Generator3 (n+3);
		uniform_real_distribution<double> Gleichvert(0.0, 1.0);
		double randomNbr1 = Gleichvert(Generator1);
		double randomNbr2 = Gleichvert(Generator2);
		double randomNbr3 = Gleichvert(Generator3);
		Data2 << randomNbr1 << "\t" << randomNbr2 << "\t" << randomNbr3 << endl;
	}

}

//Funktioniert genauso wie lk Generator, nur dass jeweils zwei aufeinanderfolgende Werte gespeichert werden
void lkGenerator_Test (int64_t r_0, int64_t a, int64_t c, int64_t m, int64_t N, string Name1, string Name2){

	// ===============================
	ofstream Data;
	Data.open(Name1.c_str());
	Data.precision(10);//
	// ===============================

	int64_t r=r_0;
	int64_t ralt=r_0;
	for(int64_t n=0; n<N;n++){
		r=(a*ralt+c)%m;
		if((n%2)==0){
		Data << (double)ralt/m << "\t" << (double)r/m << endl;
		}
		ralt=r;
	}

	// ===============================
	ofstream Data2;
	Data2.open(Name2.c_str());
	Data2.precision(10);//
	// ===============================

	r=r_0;
	ralt=r_0;
	int64_t ralt2=r_0;
	for(int64_t n=0; n<N;n++){
		r=(a*ralt+c)%m;
		if((n%3)==1){
		Data2 << (double)ralt2/m << "\t" << (double)ralt/m << "\t" << (double)r/m<< endl;
		}
		ralt2=ralt;
		ralt=r;
	}
}

double Funk_p(double x){
	double n=8./3.;
	//double n=1;
	double p=sin(M_PI*x);
	p=p*p*p*p;
	return n*p;
}

void lkGenerator_Rueckweisung (int64_t r_0, int64_t a, int64_t c, int64_t m, int64_t N, string Name){

	// ===============================
	ofstream Data;
	Data.open(Name.c_str());
	Data.precision(10);//
	// ===============================

	int64_t r=r_0;
	for(int64_t n=0; n<N;n++){
		r=(a*r_0+c)%m;
		if((n%2)==0 && ((double)r/m)<Funk_p((double)r_0/m)){
		Data << ((double)r_0/m) << endl;
		}
		r_0=r;
	}
}

void mt_Rueckweisung (int64_t N, string Name){

	// ===============================
	ofstream Data;
	Data.open(Name.c_str());
	Data.precision(10);//
	// ===============================
	for(int64_t n=0; n<N;n=n+2){
		mt19937 Generator1 (n);
		mt19937 Generator2 (n+1);
		uniform_real_distribution<double> Gleichvert(0.0, 1.0);
		double randomNbr1 = Gleichvert(Generator1);
		double randomNbr2 = Gleichvert(Generator2);
		if(randomNbr2<Funk_p(randomNbr1)){
		Data << randomNbr1 << endl;
		}
	}
}


int main(){
	int64_t N=pow(10,5);
	cout << N << endl << endl;
	lkGenerator(1234,20,120,6075,N,"bi.txt");
	lkGenerator(123456789,65539,0,pow(2,31),N,"bii.txt");
	marsenne (N, "mt.txt");
	//lkGenerator_Test(1234,20,120,6075,N,"bi_test1.txt");
	lkGenerator_Test(1234,20,120,6075,4000,"bi_test1.txt","bi_test2.txt");
	lkGenerator_Test(123456789,65539,0,pow(2,31),6000,"bii_test1.txt", "bii_test2.txt");
	marsenne_test (10000, "mt_test.txt", "mt_test2.txt");

	lkGenerator_Rueckweisung(1234,20,120,6075,N,"bi_ruck.txt");
	lkGenerator_Rueckweisung(123456789,65539,0,pow(2,31),N,"bii_ruck.txt");
	mt_Rueckweisung (N, "mt_ruck.txt");

	// ===============================
	ofstream Data;
	Data.open("p.txt");
	Data.precision(10);//
	// ===============================
	for(double n=0; n<1;n=n+0.01){
		Data << n << "\t" << Funk_p(n) << endl;
	}
	Data.close();

	cout << "Ende der main-Funktion" << endl;
	return 0;
}
