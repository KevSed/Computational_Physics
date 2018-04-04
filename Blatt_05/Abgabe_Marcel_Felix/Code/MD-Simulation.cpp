//============================================================================
// Name        : MD-Simulation.cpp
// Author      : Felix Brauers
//============================================================================

/*
 * Das Programm führt eine MD-Simulation durch
 * Die Teilchen befinden sich in einem Lennard-Jones-Potential
 * Es wird das Paket 'eigen' verwendet
 */

#define _USE_MATH_DEFINES

#include <iostream>
#include <stdio.h>
#include <time.h>
#include <random>
#include <complex>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <sstream>
#include <utility>
#include <math.h>
#include <string>
#include <fstream>
#include <Eigen/Dense>

using Eigen::MatrixXd;

using namespace std;

// Mersenne Twister Zufallsgenerator
// ============================================
random_device rd { };
mt19937 engine { rd() };
uniform_real_distribution<double> rdm_vel { -1.0, 1.0 };

double rdm_vel_generator() {
	return rdm_vel(engine);
}
// ============================================

// Konstanten festlegen
// ============================================
const double k_B = 1; // Boltzmann-Konstante

const int n_x = 4;
const int n_y = 4; // Erzeugt 4x4-Feld
const int N = n_x * n_y; // Gesamtzahl der Teilchen
const double L = 8; // Seitenlänge der betrachteten Fläche

const double h = 0.01; // Schrittweite für Verlet

const double N_H = 250; // Anzahl der Ringe, für Paarkorrelationsfunktion
const double Delta_r = 0.5 * L / N_H; // Abstand zwischen den Kreisen, für Paarkorrelationsfunktion

const double t_end = 250; // Abbruchzeit

double T0 = .01; // Anfangstemperatur
// ============================================

// *****************************************************
// BERECHNUNG DER ENERGIEN
// *****************************************************

// ============================================
// Berechne die potentielle Energie
double E_pot(MatrixXd r) {

	double E = 0; // Rückgabewert

	// Abstand zwischen zwei Teilchen
	// 1. Element: x-Abstand
	// 2. Element: y-Abstand
	vector<double> abs_ij(2);

	double abs_k = 0; // Betragsabstand zweier Teilchen
	double abs_quad = 0; // quadratischer Abstand in x- oder(!) y-Richtung

	for (int i = 0; i < N - 1; i++) {
		for (int j = i + 1; j < N; j++) {
			abs_quad = 0;

			for (int k = 0; k < 2; k++) {
				abs_ij[k] = r(i, k) - r(j, k);
				abs_k = abs(abs_ij[k]);

				// Beachte periodische Randbedingungen
				if (abs_k > 0.5 * L) {
					if (abs_ij[k] > 0) {
						abs_ij[k] -= L;
					} else {
						abs_ij[k] += L;
					}
				}
				abs_quad += abs_ij[k] * abs_ij[k];
			}
			if (sqrt(abs_quad) <= 0.5 * L) { // Beachte Cuf-Off
			// Das ist das Lennard-Jones-Potential, wenn Betragsabstand noch nicht gewurzelt ist
				E += 4 * (pow(abs_quad, -6) - pow(abs_quad, -3));
			}
		}
	}
	return E;
}
// ==========================================

// ==========================================
// Berechnung der kinetischen Energie
double E_kin(MatrixXd v) {
	double E = 0;
	for (int i = 0; i < N; i++) {
		E += 0.5 * (v(i, 0) * v(i, 0) + v(i, 1) * v(i, 1));
	}
	return E;
}
// =========================================

// =========================================
// Berechnung der Gesamtenergie
double E(MatrixXd r, MatrixXd v) {
	return E_kin(v) + E_pot(r);
}
// =========================================

// *****************************************************

// *****************************************************
// INITIALISIERUNG
// *****************************************************
MatrixXd Init() {
	// Wir erstellen eine Matrix mit Nx4 Einträgen
	// 1. + 2. Spalte: x- & y-Koordinate
	// 3. + 4. Spalte: x- & y-Geschwindigkeit
	MatrixXd conf(N, 4);

	// Anfangsorte
	// z.B. i = j = 0 --> linke obere Ecke --> x = 0.125 * L, y = 0.125 * L
	// i wird um 1 erhöht --> x rückt um 0.25*L nach rechts
	// x-Positionen sind 0.125*L , 0.375*L , 0.625*L, 0.875*L
	// genauso für j
	for (int i = 0; i < n_x; i++) {
		for (int j = 0; j < n_y; j++) {
			conf(i * n_x + j, 0) = 1. / (2 * n_x) * (1 + 2 * i) * L;
			conf(i * n_x + j, 1) = 1. / (2 * n_y) * (1 + 2 * j) * L;
		}
	}

	// Anfangsgeschwindigkeiten
	// Die Teilchen erhalten eine zufällige Geschwindigkeit zwischen -1 und 1
	for (int i = 0; i < N; i++) {
		for (int j = 2; j < 4; j++) {
			conf(i, j) = rdm_vel_generator();
		}
	}

	// Berechne Schwerpunktsgeschwindigkeit
	// Es werden die Geschwindigkeiten, getrennt in ihren Raumrichtungen aufaddiert
	double sum_x = 0, sum_y = 0;
	for (int i = 0; i < N; i++) {
		sum_x += conf(i, 2);
		sum_y += conf(i, 3);
	}

	// Setze Schwerpunktgeschwindigkeit auf 0
	// Schwerpunktgeschwindigkeit vorher (x-Richtung): v_x = (sum_i v_x,i)/N
	// Spg nachher:
	// v'_0 = (sum_i (v_i - v_0)) / N
	// = (sum_i v_i) / N - (sum_i (sum_i v_i) / N) / N
	// = v_0 - N * v_0 / N
	// = 0
	for (int i = 0; i < N; i++) {
		conf(i, 2) -= sum_x / N;
		conf(i, 3) -= sum_y / N;
	}

	return conf;
}
// ============================================

// ============================================
// Skaliere Geschwindigkeit auf bestimmte Temperatur
// Übergabeparameter sind
// Geschwindigkeitsmatrix v
// gewünschte Temperatur T
// Geschwindigkeiten werden um Faktor alpha skaliert: v_i --> alpha * v_i
// Dadurch ändert sich kinetische Energie:
// E_kin = 0.5*m*v*v --> E_kin' = alpha*alpha*E_kin
// Umstellen liefert:
// alpha = sqrt(E_kin' / E_kin) = sqrt(k_B * N_f * T / (2*E_kin))
// N_f = 2 * N - 2, Anzahl der Freiheitsgrade des Systems
// 2 Raumdimensionen pro Teilchen  2 Freiheitsgrade des erhaltenen Schwerpunktsimpulses
MatrixXd rescale_vel(MatrixXd v, double T) {
	double alpha = sqrt(k_B * (2 * N - 2) * T / (2 * E_kin(v)));
	return alpha * v;
}
// ===========================================
// *****************************************************

// *****************************************************
// MESSUNGEN
// *****************************************************

// ==========================================
// Berechnung der Schwerpunktsgeschwindigkeit
vector<double> v_SP(MatrixXd v) {

	// Rückgabevektor mit x- und y-Komponente
	vector<double> v_vec(2);

	for (int i = 0; i < N; i++) {
		for (int k = 0; k < 2; k++) {
			v_vec[k] += v(i, k);
		}
	}

	for (int k = 0; k < 2; k++) {
		v_vec[k] /= N;
	}

	return v_vec;
}
// ==========================================

// ==========================================
// Berechnung der momentanen Temperatur
// Formel: T(t) = 2 / (k_B * N_f) * sum_i E_kin,i
// N_f = 2*N - 2
double temp(MatrixXd v) {
	return 2 * E_kin(v) / (k_B * (2 * N - 2));
}
// ==========================================

// ==========================================
// Berechnung der Paarkorrelationsfunktion
// Idee der Paarkorrelationsfunktion:
// Bestimme den Abstand aller Teilchen untereinander
// Teile die möglichen Abstände in 'Bins' auf der Breite Delta_r
// Zu Messzeiten wird die Zahl der Abstände gezählt und den jeweiligen Ringen zugeordnet
// Die Gesamtzahl der Teilchen während der gesamten Simulation wird festgehalten
// Anhand der Paarkorrelationsfunktion kann die Phase des Systems bestimmt werden
vector<double> g_r(MatrixXd r, vector<double> g_r_val0) {

	vector<double> g_r_val = g_r_val0; // Vektor mit Werten zum Mitteln

	unsigned int cnt = 0; // Counter für Teilchen in einem Kreisring

	vector<double> distance(2); // Abstandsvektor mit Abstand in x- und y-Richtung zwischen zwei Teilchen

	double abs_distance = 0; // Abstandsvariable

	unsigned int pos = 0;

	for (double r_k = 0; r_k < 0.5 * L; r_k += Delta_r) {

		cnt = 0;

		for (int i = 0; i < N - 1; i++) {
			for (int j = i + 1; j < N; j++) {

				abs_distance = 0; // Setze Abstand zurück

				// Berechne nun den Abstand zwischen dem i-ten und dem j-ten Teilchen
				for (int k = 0; k < 2; k++) {
					distance[k] = r(i, k) - r(j, k);

					// Beachte die periodischen Randbedingungen
					if (abs(distance[k]) > 0.5 * L) {
						if (distance[k] > 0) {
							distance[k] -= L;
						} else {
							distance[k] += L;
						}
					}
					abs_distance += distance[k] * distance[k];
				}
				abs_distance = sqrt(abs_distance);

				// Überprüfe nun, ob der errechnete Abstand zwischen den Teilchen im betrachteten Kreisring liegt
				if (abs_distance >= r_k - 0.5 * Delta_r
						&& abs_distance <= r_k + 0.5 * Delta_r) {
					cnt += 1;
				}
			}
		}

		// Updaten, Zahl der Teilchen während der Simulation im pos-ten Ring hat sich um cnt erhöht
		g_r_val[pos] += cnt;

		// Betrachte dann den nächsten Ring
		pos += 1;

	}

	return g_r_val;

}
// *****************************************************

// *****************************************************
// VERLET-ALGORITHMUS
// *****************************************************
// Übergabeparameter sind
// r: Ortsmatrix
// v: Geschwindigkeitsmatrix
// ============================================
// Verlet-Algorithmus:
// r_n+1 = r_n + v_n*h + 0.5*a_n*h*h
// a_n+1 = a(r_n+1,t_n+1)
// v_n+1 = v_n + 0.5*h*(a_n+1 + a_n
void Verlet(MatrixXd r, MatrixXd v) {

	vector<double> distance(2); // Vektor, der Abstand in x- und y-Richtung zwischen zwei Teilchen bestimmt
	double abs_distance = 0; // Betragsabstand zwischen zwei Teilchen
	double F = 0; // Kraft-Variable
	vector<double> v_0(2); // x- und y-Variable für Schwerpunktgeschwindigkeit
	double E_k = 0; // Variable für kinetische Energie
	double E_p = 0; // Variable für potentielle Energie

	double T_avg = 0, E_kin_avg = 0, E_pot_avg = 0, E_ges_avg = 0;

	unsigned int cnt = 0; // Counter zum Abspeichern der Daten bei jedem 1000. Schleifendurchgang
	MatrixXd a(N, 2); // Erzeugung der Beschleunigungsmatrix

	// Variablen für Paarkorrelationsfunktion
	vector<double> paar_korr_vec(N_H); // Vektor, in jedem Eintrag wird die Summe der Teilchen im entsprechenden Kreisring abgespeichert

	//.txt, in der Werte für den Ort abgespeichert werden (eigentlich nur für gif)
	// ===============================
	ofstream data;
	data.open("MD.txt");
	data.precision(10);
	// ===============================

	//.txt, in der Werte für Energie abgespeichert werden
	// ===============================
	ofstream energy;
	energy.open("Energie.txt");
	energy.precision(10);
	energy << "Zeit\t E_pot\t E_kin\t E_ges" << "\n";
	// ===============================

	//.txt, in der Werte für die Schwerpunktsgeschwindigkeit abgespeichert werden
	// ===============================
	ofstream vSP;
	vSP.open("Schwerpunktsgeschwindigkeit.txt");
	vSP.precision(10);
	vSP << "Zeit\t v_x\t v_y" << "\n";
	// ===============================

	//.txt, in der Werte für Temperatur abgespeichert werden
	// ===============================
	ofstream temperatur;
	temperatur.open("Temperatur.txt");
	temperatur.precision(10);
	temperatur << "Zeit\t Temperatur" << "\n";
	temperatur << "0\t" << temp(v) << "\n";

	ofstream durchschnittsgroessen;
	durchschnittsgroessen.open("Durchschnittsgroessen.txt");
	durchschnittsgroessen.precision(10);
	// ===============================

	// Daten für Paarkorrelationsfunktion
	// ===============================
	ofstream paar_korr;
	paar_korr.open("Paarkorrelationsfunktion.txt");
	paar_korr.precision(10);
	paar_korr << "Radius\t Paarkorrelationsfunktion" << "\n";
	// ===============================

	for (double t = h; t < t_end; t += h) {

		if (cnt == 0) { // Berechne nur im 1. Schleifendurchlauf a_n, sonst mitnehmen aus der vorangegangenen Schleife
			// Berechnung mithilfe einer Doppelschleife
			// Beachtet wird, dass keine Doppeltzählung betrachtet wird
			// Berechne einfach Kraft und Gegenkraft gleichzeitig
			for (int i = 0; i < N - 1; i++) { // Zähle mit i und j alle Teilchen durch
				for (int j = i + 1; j < N; j++) { // Betrachtet wird dann die WeWe zwischen Teilchen i und Teilchen j

					abs_distance = 0;

					for (int k = 0; k < 2; k++) {
						distance[k] = r(i, k) - r(j, k); // Berechne Abstand in x- und y-Richtung

						if (abs(distance[k]) > 0.5 * L) {

							// Wertebereich distance in (-L/2 , L/2)
							// --> |distance| < 0.5*L
							// ist |distance| > 0.5*L, so ist der "Weg über den Rand kürzer" wegen periodischen Randbedingungen
							// --> korrigiere Distance um L, da sonst der "lange direkte Weg" gezählt werden würde

							if (distance[k] > 0)
								distance[k] -= L;
							else
								distance[k] += L;
						}

						// Berechne die Distanz zwischen den Teilchen (noch ohne die Wurzel zu ziehen)
						abs_distance += distance[k] * distance[k];
					}
					if (sqrt(abs_distance) <= 0.5 * L) { // Überprüfe Cutoff

					// Berechne die wirkende Kraft
					// Beachte die ungewurzelte Distanz
						F = 24
								* (2 * pow(abs_distance, -7)
										- pow(abs_distance, -4));

						for (int k = 0; k < 2; k++) {
							a(i, k) = a(i, k) + distance[k] * F; // Kraft von j auf i, beachte, dass Kraft mit Ortsrichtungskomponente ist, daher *distance[k]
							a(j, k) = a(j, k) - distance[k] * F; // (Gegen)kraft von i auf j
						}
					}
				}
			}
		}

		// Verlet-Algorithmus
		for (int i = 0; i < N; i++) {

			for (int k = 0; k < 2; k++) {

				// Berechne r_n+1 = r_n + v*h + 0.5*a*h*h
				r(i, k) = r(i, k) + h * v(i, k) + 0.5 * a(i, k) * h * h;

				// Periodische Randbedingungen prüfen, falls Teilchen aus Bild "rutscht"
				if (r(i, k) < 0)
					r(i, k) += L;
				else {
					if (r(i, k) > L)
						r(i, k) -= L;
				}

				// Geschwindigkeitsberechnung
				// Berechne v_n+1 = v_n + 0.5 * a_n+1 * h + 0.5 * a_h * h
				// Spalte auf:
				// ersten beiden Terme jetzt updaten, hinterer Term erst nach der Berechnung von a_n+1
				v(i, k) = v(i, k) + 0.5 * h * a(i, k);
			}

		}

		// Berechnung von a_n+1
		// Zunächst a wieder mit Nullen auffüllen
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < 2; j++) {
				a(i, j) = 0;
			}
		}

		// Berechnung von a_n+1 analog zu vorheriger Beschleunigungsrechnung
		// Daher unkommentiert
		for (int i = 0; i < N - 1; i++) {
			for (int j = i + 1; j < N; j++) {
				abs_distance = 0;
				for (int k = 0; k < 2; k++) {
					distance[k] = r(i, k) - r(j, k);
					if (abs(distance[k]) > 0.5 * L) {
						if (distance[k] > 0)
							distance[k] -= L;
						else
							distance[k] += L;
					}
					abs_distance = abs_distance + distance[k] * distance[k];
				}
				if (sqrt(abs_distance) <= 0.5 * L) {
					F =
							24
									* (2 * pow(abs_distance, -7)
											- pow(abs_distance, -4));
					for (int k = 0; k < 2; k++) {
						a(i, k) = a(i, k) + distance[k] * F;
						a(j, k) = a(j, k) - distance[k] * F;
					}
				}
			}
		}

		// Endgültige Geschwindigkeitsberechnung
		for (int i = 0; i < N; i++) {
			for (int k = 0; k < 2; k++) {
				v(i, k) += 0.5 * h * a(i, k);
			}
		}

		// data einschreiben, für .gif
		if (cnt % 50 == 0) {
			for (size_t i = 0; i < N; i++) {
				for (size_t k = 0; k < 2; k++) {
					data << r(i, k) << "\t";
				}
				data << r(i, 0) * r(i, 1) << endl;
			}
			data << endl;
		}

		// Energie, Schwerpunktgeschwindigkeit und Temperatur einschreiben
		if (cnt % 50 == 0) {
			E_k = E_kin(v);
			E_p = E_pot(r);
			v_0 = v_SP(v);

			energy << t << "\t" << E_p << "\t" << E_k << "\t" << E_p + E_k
					<< "\n";
			vSP << t << "\t" << v_0[0] << "\t" << v_0[1] << "\n";
		}

		v = rescale_vel(v,T0);

		if (cnt % 50 == 0) {
			temperatur << t << "\t" << temp(v) << "\n";
		}

		if (t > 50) { // Beginne mit Messung der Paarkorrelationsfunktion

			E_k = E_kin(v);
			E_p = E_pot(r);

			paar_korr_vec = g_r(r, paar_korr_vec);
			T_avg += temp(v);
			E_kin_avg += E_k;
			E_pot_avg += E_p;
			E_ges_avg += E_k + E_p;
		}

		cnt += 1;
	}

	data.close();
	temperatur.close();
	energy.close();
	vSP.close();

	// Berechne schließlich noch die Paarkorrelationsfunktion endgültig
	// Es muss noch gemittelt werden
	// g(r_l) = <p_l> / (N*rho*Delta V_l)
	// Delta V_l = 2*pi*((l*Delta_r)^2 - ((l-1)*Delta_r)^2) <-- beachte 2 Dimensionen, nicht 3
	// rho * N = 1 gesetzt, da es sich lediglich um einen konstanten (uninteressanten) Faktor handelt
	// <p_l> = p_l / Zeit --> Faktor 1 / (t_end - 50) / h
	int l = 0;

	for (double i = 0; i <= 0.5 * L; i += Delta_r) {

		if (i == 0) {
			paar_korr << 0 << "\t"
					<< paar_korr_vec[l] * 0.5
							/ (M_PI * (0.5 * Delta_r) * (0.5 * Delta_r))
							/ ((t_end - 50) / h) << "\n";
		} else {
			paar_korr << i << "\t"
					<< paar_korr_vec[l] * 0.5
							/ (M_PI
									* ((i + 0.5 * Delta_r) * (i + 0.5 * Delta_r)
											- (i - 0.5 * Delta_r)
													* (i - 0.5 * Delta_r)))
							/ ((t_end - 50) / h) << "\n";
		}
		l++;
	}

	durchschnittsgroessen << E_kin_avg / ((t_end - 50) / h) << "\t" << E_pot_avg / ((t_end - 50) / h)
			<< "\t" << E_ges_avg / ((t_end - 50) / h) << "\t" << T_avg / ((t_end - 50) / h)
			<< "\n";

	durchschnittsgroessen.close();

}
// *****************************************************

int main() {

	// Laufzeit starten
	clock_t check_clock;
	check_clock = clock();

	// Matrizen für Ort und Geschwindigkeit erstellen
	// r: Ortsmatrix
	// v: Geschwindigkeitsmatrix
	// beide Matrizen haben Dimension Nx2
	// Die Spalte nummeriert die Teilchen durch
	// Die beiden Zeilen stellen die x- und die y-Komponente dar
	// ============================================
	MatrixXd r(N, 2);
	MatrixXd v(N, 2);
	// ============================================

	// Initialisierung
	// Anfangspositionen + Anfangsgeschwindigkeit
	// ============================================
	MatrixXd matrix_init = Init();
	for (int i = 0; i < N; i++) {
		r(i, 0) = matrix_init(i, 0);
		r(i, 1) = matrix_init(i, 1);
		v(i, 0) = matrix_init(i, 2);
		v(i, 1) = matrix_init(i, 3);
	}
	// ============================================

	// Skalieren auf gewünschte Temperatur
	v = rescale_vel(v, T0);

	// Verlet-Aufruf
	Verlet(r, v);

	//Laufzeit berechnen
	check_clock = clock() - check_clock;
	cout << check_clock;

	return 0;
}
