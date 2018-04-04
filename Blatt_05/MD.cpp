/*
 * MD.cpp
 *
 *  Created on: 22.05.2017
 *      Author: mona
 */

#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <sstream>
#include <utility>
#include <math.h>
#include <fstream>
#include <functional>
#include "Eigen/Core"

using namespace std;
using namespace Eigen;

// Die Funktion ---Quadratgitter--- soll eine Matrix erstellen, deren Spaltenanzahl der Teilchenanzahl entspricht
// und deren Spalten die x und y Komponenten des Ortes des Teilchens enthalten
// Dabei entspricht L der Kantenlaenge der betrachteten Box, in der die Teilchen aequidistant positioniert werden.

MatrixXd Quadratgitter(int Teilchenzahl, double L){
	MatrixXd Teilchenvektor(2, Teilchenzahl);
	Vector2d Teilchenposition;
	int Teilchenzaehler =0;
	double Teilchenabstand = L/(sqrt(Teilchenzahl));

	// Hier werden die Teilchen gleichmaessig auf dem Gitter Positioniert.
	for(int i=0; i<sqrt(Teilchenzahl); i++){
		for(int j=0; j<sqrt(Teilchenzahl); j++){
				Teilchenposition << (Teilchenabstand/2.)*(1+2*i) , (Teilchenabstand/2.)*(1+2*j);
				for(int Komponente=0; Komponente<2; Komponente++){
					Teilchenvektor(Komponente, Teilchenzaehler)=Teilchenposition(Komponente);
				}
				Teilchenzaehler++;
			}
	}
	return Teilchenvektor;
}

//Die Funktion ---initGeschwindigkeit--- soll allen Teilchen zunaecht eine zufaellige Geschwindigkkeit geben,
//und dann die Schwerpunktsgeschwindigkeit des Systems von allen Teilchen abziehen.
//Die neue Schwerpunktsgeschwindigkeit ist also null
//Der Funktion wird außerdem die Temperatur T uebergeben, sodass die Geschwindigkeit einer Startenergie angepasst werden kann.

MatrixXd initGeschwindigkeit(int Teilchenzahl, double T){
	MatrixXd Teilchengeschw(2, Teilchenzahl);
	double xSchwerpunkt=0;
	double ySchwerpunkt=0;
	double xrand=0;
	double yrand=0;

	//Allen Teilchen eine zufaellige Geschwindigkeit geben
	for(int i=0; i<Teilchenzahl; i++){
		//nur Geschwindigkeiten kleiner als 1
		xrand=rand()* (1./RAND_MAX);
		yrand=rand()* (1./RAND_MAX);
		Teilchengeschw(0,i)=xrand;
		Teilchengeschw(1,i)=yrand;
		xSchwerpunkt=xSchwerpunkt+xrand;
		ySchwerpunkt=ySchwerpunkt+yrand;
	}

	//Schwerpunktsgeschwindigkeit berechnen
	xSchwerpunkt=(1./Teilchenzahl)*xSchwerpunkt;
	ySchwerpunkt=(1./Teilchenzahl)*ySchwerpunkt;

	//Schwerpunktsgeschwindigkeit von allen Teilchen abziehen
	for(int i=0; i<Teilchenzahl; i++){
		Teilchengeschw(0,i)=Teilchengeschw(0,i)-xSchwerpunkt;
		Teilchengeschw(1,i)=Teilchengeschw(1,i)-ySchwerpunkt;
	}

	//Isokinetischer Thermostat zur Reskalierung der Geschwindigkeiten
	//Berechnung der Momentanen Temperatur T0
	double Nf=3*Teilchenzahl-3; //Anzahl der Freiheitsgrade
	double T0=0;
	double alpha;
	double xs=0;
	double ys=0;

	for(int i=0; i<Teilchenzahl; i++){
			xs=Teilchengeschw(0,i);
			ys=Teilchengeschw(1,i);
			T0=T0+(xs*xs+ys*ys);
	}
	T0=(1./Nf)*T0;
	alpha=sqrt(T/T0);
	for(int i=0; i<Teilchenzahl; i++){
		Teilchengeschw(0,i)=Teilchengeschw(0,i)*alpha;
		Teilchengeschw(1,i)=Teilchengeschw(1,i)*alpha;
	}
	return Teilchengeschw;
}

//Die Funktion ---Randebdingung--- ueberprueft, ob ein Teichen sich innerhalb der Box befindet.
//Ist dies nicht der Fall, wird das Teilchen zurueck in die Box befoerdert

MatrixXd RB(MatrixXd Position, double L, int Teilchenzahl){
	for(int n=0; n< Teilchenzahl; n++){
		for(int i=0; i< 2; i++){
			if(Position(i,n)<0){Position(i,n)+=L;}
			if(Position(i,n)>L){Position(i,n)-=L;}
		}
	}
	return Position;
}

//Die Funktion ---LJForce--- giebt die Kraft aus, welche nach dem Lennard-Jones-Potential auf ein Teilchen wirkt,
//Welches relativ zu einem Anderen Teilchen um den Vektor r verschoben ist.
//Ausserdem soll die Kraft null sein, wenn die Teilchen um mehr als den kritischen radius von einander entfernt sind

Vector2d LJForce(Vector2d r, double rc){
	Vector2d LJF;
	double Abstand=r.norm();
	if (Abstand>rc ){
		LJF << 0,0;
		return LJF;
	}
	LJF=r*24;
	LJF=LJF*(2*pow(Abstand,-14)-pow(Abstand,-8));
	return LJF;
}

//Da ein Teilchen von jedem seiner Bildteilchen um L entfernt ist
//und die Bildteichen ebenfalls mindestens um L von einander entfernt sind
//erlaubt uns die eine Cutof-Distanz von rc<L/2 nur die Wechselwirkung mit dem naechten Teilchen zu berachten.
//Ist ein Teilchen/Bildteilchen von Teilchen 2 also naeher an Teilchen 1 als L/2
//So ist der Abstand zu allen anderen Teilchen zwangslaeufig groesser als L/2 und es gibt keine WW

Vector2d NaechstesTeilchen(Vector2d Teilchen1, Vector2d Teilchen2,double L, double rc){

	Vector2d v(L,0);
	Vector2d Bild1=Teilchen2+v;
	Vector2d Abstand=Teilchen1-Bild1;

	if(Abstand.norm()<rc){return Bild1;}

	Vector2d Bild2=Teilchen2-v;
	Abstand=Teilchen1-Bild2;
	if(Abstand.norm()<rc){return Bild2;}

	v << 0,L;
	Vector2d Bild3=Teilchen2+v;
	Abstand=Teilchen1-Bild3;
	if(Abstand.norm()<rc){return Bild3;}

	Vector2d Bild4=Teilchen2-v;
	Abstand=Teilchen1-Bild4;
	if(Abstand.norm()<rc){return Bild4;}

	v << L,L;
	Vector2d Bild5=Teilchen2+v;
	Abstand=Teilchen1-Bild5;
	if(Abstand.norm()<rc){return Bild5;}

	Vector2d Bild6=Teilchen2-v;
	Abstand=Teilchen1-Bild6;
	if(Abstand.norm()<rc){return Bild6;}

	v << L,-L;
	Vector2d Bild7=Teilchen2+v;
	Abstand=Teilchen1-Bild7;
	if(Abstand.norm()<rc){return Bild7;}

	Vector2d Bild8=Teilchen2-v;
	Abstand=Teilchen1-Bild8;
	if(Abstand.norm()<rc){return Bild8;}

	return Teilchen2;
}


//Die Funktion ---Kraft-- soll die Kraft berechnen, die auf jedes Teilchen wirkt
//Diese Kraft ist abhaengig von den Positionen aller anderen Teilchen, sowie aller Bildteilchen

MatrixXd Kraft(int Teilchenzahl, MatrixXd Position, double L, double rc){
	MatrixXd F(2, Teilchenzahl);
	F.setZero();
	Vector2d Abstandsvector;
	Vector2d Teilchen1;
	Vector2d Teilchen2;
	Vector2d Ftemp;

	// Dies funktioniert nur fuer rc<0.5L, weil sonst auch eine Wechselwirkung mit mehr als einer Kopie des Teilchens vorliegen kann
	for(int i=0; i<(Teilchenzahl-1); i++){
		for(int j=i+1; j<Teilchenzahl; j++){

			Teilchen1=Position.col(i);
			Teilchen2=Position.col(j);
			Teilchen2=NaechstesTeilchen(Teilchen1, Teilchen2, L, rc);
			//cout << "Teilchen1= " << Teilchen1 << endl;
			//cout << "Teilchen2= " << Teilchen2 << endl;
			Abstandsvector=Teilchen1-Teilchen2;
			if(Abstandsvector.norm()<rc){
			//cout << "Abstand zum naechsten Teilchen= " << Abstandsvector.norm() << endl;
			Ftemp=LJForce(Abstandsvector,L);
			//cout << "Ftemp für i= " << i << "und j= " << j << "ist " <<endl << Ftemp << endl;
			F(0,i)=F(0,i)+Ftemp(0);
			F(1,i)=F(1,i)+Ftemp(1);
			F(0,j)=F(0,j)-Ftemp(0);
			F(1,j)=F(1,j)-Ftemp(1);
			}
		}
	}
	return F;
}

//Funktion zur Berechnung der Schwerpunktsgeschwindigkeit

Vector2d Schwerpunktsgeschwindigkeit(MatrixXd Geschwindigkeit, int Teilchenzahl){
	double xSwpkt=0;
	double ySwpkt=0;
	Vector2d Swpkt;
	for(int i=0;i<Teilchenzahl;i++){
		xSwpkt=xSwpkt+Geschwindigkeit(0,i);
		ySwpkt=ySwpkt+Geschwindigkeit(1,i);
	}
	Swpkt << xSwpkt,ySwpkt;
	return Swpkt;
}

//Funktion zur Berechnung der kinetischen Energie

double kinetischeEnergie(MatrixXd Geschwindigkeit, int Teilchenzahl){
	double Ekin=0;
	for(int i=0;i<Teilchenzahl;i++){
		Ekin=Ekin+(Geschwindigkeit(0,i)*Geschwindigkeit(0,i)+Geschwindigkeit(1,i)*Geschwindigkeit(1,i));
	}
	Ekin=0.5*Ekin;
	return Ekin;
}

//Funktion zur Berechnung der potentiellen Energie

double potentielleEnergie(MatrixXd Position, int Teilchenzahl, double rc, double L){
	double Epot=0;
	double Abstand;
	Vector2d Teilchen1;
	Vector2d Teilchen2;
	Vector2d Abstandsvector;

	for(int i=0; i<Teilchenzahl; i++){
		for(int j=0; j<Teilchenzahl; j++){
			if(i!=j){
				Teilchen1=Position.col(i);
				Teilchen2=Position.col(j);
				Teilchen2=NaechstesTeilchen(Teilchen1, Teilchen2, L, rc);
				Abstandsvector=Teilchen1-Teilchen2;
				Abstand=Abstandsvector.norm();
				if(Abstand<rc){
					Epot=Epot+pow(Abstand,-12)-pow(Abstand,-6);
				}
			}
		}
	}
	Epot=2*Epot;
	return Epot;
}

VectorXd Paarkorrelation(int Nbin,double L, MatrixXd Position, int Teilchenzahl){

	double deltaR=L/(2*Nbin);
	VectorXd gvec(Nbin);
	gvec.setZero();
	Vector2d Teilchen1;
	Vector2d Teilchen2;
	Vector2d Abstandsvector;
	double Abstand;

	for(int i=0; i<Teilchenzahl; i++){
		for(int j=0; j<Teilchenzahl; j++){
			if(i!=j){
				Teilchen1=Position.col(i);
				Teilchen2=Position.col(j);
				Teilchen2=NaechstesTeilchen(Teilchen1, Teilchen2, L, L/2);
				Abstandsvector=Teilchen1-Teilchen2;
				Abstand=Abstandsvector.norm();
				if(Abstand<(L/2)){
					//cout << i << "\t" << Abstand << endl;
					gvec(floor(Abstand/deltaR))=gvec(floor(Abstand/deltaR))+1;
				}
			}
		}
	}
	gvec=gvec*0.5;
	return gvec;
}


//Die Funktion ---integrate-- soll einen Verlet-Algorithmus verwenden, um die Bewegungen der Teilchen zu integrieren

void integrate(MatrixXd Position, MatrixXd Geschwindigkeit, double N, double h, int Teilchenzahl, string Energie, string Paarkorr, double L, bool Energiebetrachtung, int Nbins, int teq){
	double rc=0.5*L;
	MatrixXd rtemp(2,Teilchenzahl);
	rtemp=Position;
	MatrixXd vtemp(2,Teilchenzahl);
	vtemp=Geschwindigkeit;
	MatrixXd Ftemp(2,Teilchenzahl);
	Ftemp=Kraft(Teilchenzahl, rtemp, L, rc);
	MatrixXd Ftempnext(2,Teilchenzahl);
	Vector2d Swpkt_v;
//#include <iostream>
	double Ekin;
	double Epot;
	double Tcurrent;
	double Nf=3*Teilchenzahl-3;
	double Tmean=0;
	VectorXd g(Nbins);
	g.setZero();
	int counter=0;
	double rho;
	double Delta_VL;
	double deltaR=L/(2*Nbins);

	// ===============================
	ofstream Data;
	Data.open(Energie.c_str());
	Data.precision(10);//
	// ===============================
	// ===============================
	ofstream gData;
	gData.open(Paarkorr.c_str());
	gData.precision(10);//
	// ===============================


	for(int i=0; i<N; i++){

		if (Energiebetrachtung == true){
		Swpkt_v=Schwerpunktsgeschwindigkeit(vtemp, Teilchenzahl);
		Ekin=kinetischeEnergie(vtemp, Teilchenzahl);
		Epot=potentielleEnergie(rtemp, Teilchenzahl, rc, L);
		Tcurrent=(2*Ekin)/Nf;
		Tmean=Tmean+Tcurrent;
		}

		for (int n=0; n<Teilchenzahl;n++){
			rtemp(0,n)=rtemp(0,n)+vtemp(0,n)*h+0.5*Ftemp(0,n)*h*h;
			rtemp(1,n)=rtemp(1,n)+vtemp(1,n)*h+0.5*Ftemp(1,n)*h*h;
		}

		rtemp=RB(rtemp,L,Teilchenzahl);
		Ftempnext=Kraft(Teilchenzahl, rtemp, L, rc);

		for (int n=0; n<Teilchenzahl;n++){
			vtemp(0,n)=vtemp(0,n)+0.5*h*(Ftempnext(0,n)+Ftemp(0,n));
			vtemp(1,n)=vtemp(1,n)+0.5*h*(Ftempnext(1,n)+Ftemp(1,n));
		}
		Ftemp=Ftempnext;

		if (Energiebetrachtung == true){
		Data << i*h << "\t" << Swpkt_v(0) << "\t"<< Swpkt_v(1)<< "\t" << Swpkt_v.norm();
		Data  << "\t"<< Ekin  << "\t" << Epot <<"\t" << Ekin+Epot;
		Data << "\t"<< Tcurrent << "\t"<< 1.284 << endl;
		}

		if (i*h>teq){
			counter=counter+1;
			g=g+Paarkorrelation(Nbins,L, rtemp,Teilchenzahl);
		}

	}
	g=(1./counter)*g;
	rho=(N)/(L*L);
	for(int l=1;l<=Nbins;l++){
		Delta_VL=M_PI*((l*deltaR)*(l*deltaR)-((l-1)*deltaR)*((l-1)*deltaR));
		g[l-1]=g[l-1]/(Delta_VL*rho*Teilchenzahl);
		gData << l*deltaR << "\t" << g[l-1] << endl;
	}

	cout << endl << endl << rtemp << endl;
	Data.close();
	gData.close();

}

int main(){

	int Teilchenzahl=16;
	double L=8;
	double h=0.01;
	MatrixXd Teilchenvektor=Quadratgitter(Teilchenzahl,L);
	cout << "Teilchenvektor:" << endl << Teilchenvektor << endl << endl;
	MatrixXd TeilchengeschwT1=initGeschwindigkeit(Teilchenzahl,1);
	MatrixXd TeilchengeschwT001=initGeschwindigkeit(Teilchenzahl,0.01);
	MatrixXd TeilchengeschwT100=initGeschwindigkeit(Teilchenzahl,100);
	cout << "Teilchengeschwindigkeit" << endl << TeilchengeschwT1 << endl << endl;


	integrate(Teilchenvektor, TeilchengeschwT1, 100000, h, Teilchenzahl, "EnergieT=0.txt", "gT=0.txt", L, true,3000,5);
	integrate(Teilchenvektor, TeilchengeschwT001, 100000, h, Teilchenzahl, "EnergieT=0.txt", "gT=001.txt", L, false,3000,5);
	integrate(Teilchenvektor, TeilchengeschwT100, 100000, h, Teilchenzahl, "EnergieT=0.txt", "gT=100.txt", L, false,3000,5);


	cout << "Ende der main-Funktion" << endl;
	return 0;
}
