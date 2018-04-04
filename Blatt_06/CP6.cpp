#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<stdlib.h>
#include<armadillo>
//Es wird das Package armadillo benutzt, weil man sich damit die Funktionen für die Operatorüberladung spart und
//einfach mit Matrizen rechnen kann

using namespace std;
using namespace arma;

//Festlegung von Delta, die Länge des Boxes und die Länge der kleinen Boxen
double delta = 0.05;
double L = 1.0;
double N = L / delta;

//a)
/*Es wird durch zwei for- Schleifen über alle Kästen in der Box durchgelaufen und dabei jeweils der dazugehörige
 * Wert von dem Potential mithilfe der Gaus-Seidel-Iteration berechnet */
mat gaussSeidel(mat phi,		// phi Startwert
				mat roh,		// roh Startwert
				string path,	// output pfad
				double N
				)
{
	ofstream A1(path);
	for (double j = 1.0; j < N-1; j++)
	{
		for (double l = 1.0; l < N-1; l++)
		{
			phi(j, l) = 0.25 * (phi(j + 1, l) + phi(j - 1, l) + phi(j, l + 1) + phi(j, l - 1)) + 0.25 * delta * delta * roh(j,l);
		}
	}
	for (double a = 0.0; a < N; a++){ //Aufrufen der Koordinaten von 0 bis 19 und der Matrix Phi und Speicherung in txt Datei
			for (double b = 0.0; b < N; b++){
				A1 << a << "\t" << b << "\t" << phi(a,b) << endl;
			}
		}
return phi;
}

//c)
/* Das Potential wird wieder mit Gaus-Seidel wie in der a) berechnet, nur jetzt solange bis alle Punkte innerhalb der Box
 * eine Genauigkeit von 10^-5 aufweisen. Dies geschieht durch eine while Schleife, bei der jeweils der aktuelle berechnete
 * Wert von dem Potential an einen Punkt j,l mit dem davorberechneten Wert verglichen. Ist dabei die Differenz kleiner
 * als 10^-5, wird der Counter c um eins erhöht. Erreicht der Counter die Anzahl der Plätze innerhalb der Box, wird
 * die while-Schleife beendet*/
mat gaussSeidelc(mat phi,		// phi Startwert
				 mat roh,		// roh Startwert
				 string path,	// output pfad
				 double N
				 )
{
	double c = 0.0; // Initialisierung von Counter c
	mat phi_n;
	phi_n.zeros(N,N);
	ofstream A2(path);
	while(c < (N-1)*(N-1)){ //Abruchbedingung für die while-Schleife
		for (double j = 1.0; j < N-1; j++)
		{
			for (double l = 1.0; l < N-1; l++)
			{
				phi_n(j, l) = 0.25 * (phi(j + 1, l) + phi(j - 1, l) + phi(j, l + 1) +
						phi(j, l - 1)) + 0.25 * delta * delta * roh(j,l);
				if(phi_n(j,l)-phi(j,l) <= pow(10,-5)){ //Vergleich von aktuellen phi(j,l) mit alteb phi(j,l)
					c++; //Erhöhung des Counters
				}
				phi(j,l) = phi_n(j,l);
			}
		}
	}
	for (double a = 0.0; a < N; a++){ //Speicherung von j,l und phi(j,l) in einer txt datei
		for (double b = 0.0; b < N; b++){
			A2 << a << "\t" << b << "\t" << phi(a,b) << endl;
		}
	}
	return phi;
}

/* Das in c berechnete Phi wird in diese Funktion eingegeben. In dieser Funktion wird wieder über alle Plätze in
 * dem Kasten gelaufen und dabei jeweils die x und die y Komponente von E-Feld mithilfe der Funktion des
 * Differenzenquotient von 2 symmetrischen Punkten (aus dem Skript) berechnet. h ist dabei delta, also der benachbarte
 * Punkt in der Box*/
void Efeld2(mat phi, double N, string path){
	double x;
	double y;
	ofstream A3(path);
	for(double j = 1.0; j < N-1; j++){
		for(double l = 1.0; l < N-1; l++){
			x = (-phi(j + 1, l) + phi(j - 1, l))/(2 * delta); //Berechnung der x-Komponente
			y = (-phi(j, l + 1) + phi(j, l - 1))/(2 * delta); //Berechnung der y-Komponente
			A3 << j << "\t" << l << "\t" << sqrt(x*x+y*y) << endl; //Speicherung von j, l und der Betrag des E-Feldes
																	//in eine txt Datei
		}
	}
}
//d)
/* Es wird über alle vier Ränder, gegen den Uhrzeigersinn, gegangen und dabei jeweils wieder E wie bei der c) berechnet.
 * Nur nun jeweils ein
 * Punkt vom Rand mit dem nächstliegenden Punkt in der Box. Das Ergebnis wird für jeden Punkt auf dem Rand aufsummiert.
 * Es muss auf die richtigen Vorzeichen durch den normalen Vektor geachtet werden.*/
double Linienintegral(mat phi, double N){
	double lin = 0.0;
		for(double i= 0.0; i < N - 1; i++){ //Durchgehen der Punkte von 0 bis 19
			lin += -(phi(0, i) - phi(1, i))/delta; //Rand von y = 0 und x = 0 bis 19
			lin += +(phi(i, N - 1) - phi(i, N -2))/delta; //Rand von x = 19 und y = 0 bis 19
			lin += -(phi(N - 1, i) - phi(N -2, i))/delta; //Rand von y = 19 und x = 0 bis 19
			lin += +(phi(i, 0) - phi(i, 1))/delta; //Rand von x = 0 und y = 0 bis 19
		}
	return lin;
}

//Berechnung des Winkels von E-Feld zu einem Rand. Dabei wird die y-Komponente (hier l) fest vorgegeben und über
//alle Punkte dieser y-Komponente innerhalb der Box gegangen und dabei wieder das E-Feld wie in der c berechnet
//Den Winkel bekommt man dann über den arctan(x/y)
void Winkel(mat phi, double l, double N, string path){
	ofstream A4(path);
	double x;
	double y;
	for(double j = 1.0; j < N-1 ; j++){ //Schleife über Punkte innerhalb der Box
		x = (-phi(j + 1, l) + phi(j - 1, l))/(2 * delta); //Berechnung x-Komponente des E-Feldes
		y = (-phi(j, l + 1) + phi(j, l - 1))/(2 * delta); //Berechnung der y-Komponente des E-Feldes
		A4 << j << "\t" << atan(x/y) << endl; //Speicherung von j und des dazugehörigen Winkel des E-Feldes
	}
}

int main()
{
	mat phi;
	phi.zeros(N,N); //Initialisierung von Matrix phi
	mat phi_new;
	phi_new.zeros(N,N); //Initialiserung von Matrix phi_new
	mat rho;
	rho.zeros(N,N); //Initialisierung von Matrix rho
	rho(5,5) = 1;  //mit Ladung 1 auf der Stelle 5,5
	rho(10,10) = -1; // und Ladung -1 auf der Stelle 10,10
	/*for(int i = 0; i <= 19; i++){ //Randbedingungen für b)
		phi(19,i) = 1.0;
	}*/
	/*for(int a = 0; a <= 10; a++){ //Potential für b wird 11 mal iteriert um einen etwas genaueren Wert zu bekommen,
	phi_new = gaussSeidel(phi, rho, "CP6_potential_b.txt", N); //in der Aufgabenstellung war kein Hinweis wie häufig
	phi = phi_new;
	}*/
	mat bla = gaussSeidelc(phi, rho,"CP6_Potential_e.txt", N ); //Berechnung von dem Potential
	//cout << Linienintegral(bla, N) << endl; //Berechnung des Linienintegral
	Efeld2(bla, N, "CP6_E_Feld_e.txt"); //Berechnung des E-Feldes
	//Winkel(bla, 5.0, N, "CP6_Winkel.txt"); //Berechnung des Winkels
	return 0;
}




