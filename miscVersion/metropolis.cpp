
/*Montecarlo per generare le superfici di separazione*/
 
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <random>
#include <ctime>
using namespace std;

//Inserire il numero di iterazioni
#define N 100

int clusterize(int*** latt, int*** cluster, int boundary, int l1, int l2, int t, long double beta, int i, int j, int k, uniform_real_distribution<long double>*  distribution);
int surfaceCluster(int*** latt,  int*** cluster, int boundary, int l1, int l2, int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_real_distribution<long double>*  distribution);
void singleCluster(int*** latt,  int*** cluster, int boundary, int l1, int l2, int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_int_distribution<int>* distr_t, uniform_real_distribution<long double>*  distribution);
void hotStart(int*** latt, int l1, int l2 , int t);
void coldStart(int*** latt, int l1, int l2, int t, int spin);

random_device rd;
mt19937 mt(rd());

int main() {
	
	//Inserire i parametri delle simulazioni, verranno eseguite tutte le possibili combinazioni a meno di aggiunta di condizioni.
	const vector<long double> Beta = {0.2391};
	const vector<int> L1 = {4};
	const vector<int> L2 = {4};
	const vector<int> T = {4*3};
	
	int boundary = 1;

	for(int i = 0; i < Beta.size(); i++) {
		for(int j = 0; j < L1.size(); j++) {
			for (int h = 0; h < L2.size(); h++) {
				for(int k = 0; k < T.size(); k++) {

					//Questo è il posto ideale per inserire una eventuale parallelizzazione.
					
					long double beta = Beta[i];
					int l1 = L1[j];
					int l2 = L2[j];
					int t = T[k];
					
					//per ogni condizione data scrivo nuovo file il cui nome sarà nel formato beta_l1_l2_t
					ofstream file;
					ostringstream beta_str, l1_str, l2_str, t_str;
					
					beta_str << beta;
					l1_str << l1;
					l2_str << l2;
					t_str << t;
					
					file.open ("metropolis_" + beta_str.str() + "_" + l1_str.str() + "_" + l2_str.str() + "_" + t_str.str() + ".txt");
					
					//header file
					file << "#\t mag\t abs_mag\t b.c." << endl;

					//Inserire qui le eventuali condizioni per rimuovere le combinazioni non di interesse. Esempio
					//if(l1 != l2) continue;
					
  					uniform_int_distribution<int> distr_l1(0,l1-1);
  					uniform_int_distribution<int> distr_l2(0,l2-1);
					uniform_int_distribution<int> distr_t(0,t-1);
					uniform_real_distribution<long double> distribution(0,1);
					
					//Definisco la matrice tentativa
					uniform_int_distribution<int> distr_bcflip(0,1)
					uniform_int_distribution<int> distr_adjlink(0,5)
					
					int*** latt;
  					latt = (int***) malloc(l1 * sizeof(int**) );
					for(int a = 0 ; a < l1 ; a++) {
						latt[a] = (int**) malloc (l2 * sizeof(int*));
						for(int b = 0; b < l2; b++) {
						latt[a][b] = (int*) malloc (t * sizeof(int));
						}
					}

					
						
					//Inserire la condizione di partenza
					coldStart(latt, l1, l2, t, 1);
				
					int per = 0; 
					int aper = 0;
					int init = clock();
					for(long int rep = 0; rep < N; rep++) {
						
						int tot_mag = 0;
						int abs_tot_mag = 0;
						long double mag = 0;
						long double abs_mag = 0;
						
						cout << rep  << "\r";
						
						//Seleziono configurazione tentativa
						
						
						//ciclo stupido sul lattice per calcolare magnetizzazione (potrebbe essere inglobato altrove visto che tanto ci si passa sempre...)
						for(int a = 0 ; a < l1 ; a++){
							for(int b = 0 ; b < l2 ; b++){
								for(int c = 0 ; c < t ; c++){
									tot_mag += latt[a][b][c];
								}
							}
						}
						abs_tot_mag = abs(tot_mag);
						mag = (double)tot_mag/(l1*l2*t);
						abs_mag = (double)abs_tot_mag/(l1*l2*t);
						file << mag << "\t" << abs_mag << "\t" << boundary << endl;
					
					}
					int fine = clock() - init;

					double Aper = double(aper);
					double Per = double(per);
					double td = double(t);
					
					cout << "l1 = " << l1 << "\t" "l2 = " << l2 << "\t" << "t = " << t << "\t"  << "beta = " << beta << "\t" << "tempo = " << fine << endl;
					cout  << "Aper/Tot = " << Aper/(N) << "\t" << "Aper/Per = " << Aper/Per << "\t" << "F_s = " << -log(Aper) + log(Per) + log(td) << "\t" << "F_si = " << log(td) - log( 0.5 * log((1.0+ Aper/Per)/(1.0 - Aper/Per))) << endl;	
					
					for(int a = 0 ; a < l1 ; a++) {
						for(int b = 0; b < l2; b++) {
							free(latt[a][b]);
						}
						free(latt[a]);
					}
					free(latt);

					for(int a = 0 ; a < l1 ; a++) {
						for(int b = 0; b < l2; b++) {
							free(cluster[a][b]);
						}
						free(cluster[a]);
					}
					free(cluster);
					file.close();
	
				}
			}
		}
	}
				
	
	return 0;
}

void coldStart(int*** latt, int l1, int l2, int t, int spin) {

	for(int i = 0; i < l1; i++) 
		for(int j = 0; j < l2; j++) 
			for(int k = 0; k < t; k++) 
				latt[i][j][k] = spin;
}

void hotStart(int*** latt, int l1, int l2, int t) {

	random_device rd;
    	mt19937 mt(rd());
  	uniform_int_distribution<int> distribution(0,1);

	for(int i = 0; i < l1; i++) {
		for(int j = 0; j < l2; j++) {
			for(int k = 0; k < t; k++) {
				if(distribution(mt) == 0) latt[i][j][k] = 1;
				else latt[i][j][k] = -1;
			}
		}
	}

}
