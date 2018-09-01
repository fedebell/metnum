
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
	const vector<int> L1 = {30};
	const vector<int> L2 = {30};
	const vector<int> T = {30*3};
	
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
					
					file.open (beta_str.str() + "_" + l1_str.str() + "_" + l2_str.str() + "_" + t_str.str() + ".txt");
					
					//header file
					file << "#\t mag\t abs_mag\t b.c." << endl;

					//Inserire qui le eventuali condizioni per rimuovere le combinazioni non di interesse. Esempio
					//if(l1 != l2) continue;
					
  					uniform_int_distribution<int> distr_l1(0,l1-1);
  					uniform_int_distribution<int> distr_l2(0,l2-1);
					uniform_int_distribution<int> distr_t(0,t-1);
					uniform_real_distribution<long double> distribution(0,1);

					int*** latt;
  					latt = (int***) malloc(l1 * sizeof(int**) );
					for(int a = 0 ; a < l1 ; a++) {
						latt[a] = (int**) malloc (l2 * sizeof(int*));
						for(int b = 0; b < l2; b++) {
						latt[a][b] = (int*) malloc (t * sizeof(int));
						}
					}

			
					int*** cluster;
  					cluster = (int***) malloc(l1 * sizeof(int**) );
					for(int a = 0 ; a < l1 ; a++) {
						cluster[a] = (int**) malloc (l2 * sizeof(int*));
						for(int b = 0; b < l2; b++) {
							cluster[a][b] = (int*) malloc (t * sizeof(int));
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
						
						if(rep % 2 == 0){
							singleCluster(latt, cluster, boundary, l1, l2, t, beta, &distr_l1, &distr_l2, &distr_t, &distribution);
						}
						else
							boundary = surfaceCluster(latt, cluster, boundary, l1, l2, t, beta, &distr_l1, &distr_l2, &distribution);
						
						if(boundary == 1) per += 1;
						else aper += 1;
						
						
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

int clusterize(int*** latt, int*** cluster, int boundary, int l1, int l2, int t, long double beta, int precedente, int cross, int i, int j, int k, uniform_real_distribution<long double>*  distribution) {

	int size[3] = {l1, l2, t};
	int coord[3] = {i, j, k};
	
	int toClusterize[6][4];
	for(int i = 0; i < 6; i++) 
		for (int j = 0; j < 4; j++) 
			toClusterize[i][j] = -1;
	
	int flag = 0;
	int bound = +1;
	int newCross = -1;
	int b = 0;
	int count = 0;
	
	if(cluster[i][j][k] == (cross * precedente)) return 1;
	else if(cluster[i][j][k] == 0) {
		cluster[i][j][k] = precedente * (-cross);
		for(int d = 0; d < 3; d++) {
			for(int a = -1; a < 2; a = a + 2) {
			
				bound = +1;
				newCross = -1;
				coord[0] = i; coord[1] = j; coord[2] = k;
				
				coord[d] +=  a;
				if(coord[d] == size[d]) {
				 coord[d] = 0;
				 if(d == 2) { bound = boundary; newCross = 1; }
				}
				
				if(coord[d] == -1) {
				 coord[d] = size[d] - 1;
				 if(d == 2) { bound = boundary; newCross = 1; }
				}
				
				if(cluster[coord[0]][coord[1]][coord[2]] == 0) {
					if ((*distribution)(mt) < (1.0 - exp(-beta * (1.0 + bound*latt[i][j][k]*latt[coord[0]][coord[1]][coord[2]])))) {
						for(int i = 0; i<3; i++) toClusterize[count][i] = coord[i];
						toClusterize[count][3] = newCross;
						count++;
					}
				}
			}
		}
		for(int c = 0; c < count; c++) flag = clusterize(latt, cluster, boundary, l1, l2, t, beta, cluster[i][j][k], toClusterize[c][3], toClusterize[c][0], toClusterize[c][1], toClusterize[c][2], distribution) || flag;
	}
	
	return flag;

}

void singleCluster(int*** latt,  int*** cluster, int boundary, int l1, int l2, int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_int_distribution<int>* distr_t, uniform_real_distribution<long double>*  distribution) {

	for(int i = 0; i < l1; i++) {
		for(int j = 0; j < l2; j++) {
			for(int k = 0; k < t; k++) {
				cluster[i][j][k] = 0;
			}
		}
	}

	int i = (*distr_l1)(mt);
	int j = (*distr_l2)(mt);
	int k = (*distr_t)(mt);

	clusterize(latt, cluster, boundary, l1, l2, t, beta, -1, -1, i, j, k, distribution);	
	
	for(int i = 0; i < l1; i++) {
		for(int j = 0; j < l2; j++) {
			for(int k = 0; k < t; k++) {
				latt[i][j][k] = latt[i][j][k] * ((cluster[i][j][k] == 0) ? 1 : -1);
			}
		}
	}	
	
}

int surfaceCluster(int*** latt,  int*** cluster, int boundary, int l1, int l2, int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_real_distribution<long double>*  distribution) {

	int flag = 0;

	for(int i = 0; i < l1; i++) 
		for(int j = 0; j < l2; j++) 
			for(int k = 0; k < t; k++) 
				cluster[i][j][k] = 0;
	
	for(int i = 0; i < l1; i++) 
		for(int j = 0; j < l2; j++) 
			if(cluster[i][j][t-1] == 0) 
				if(clusterize(latt, cluster, boundary, l1, l2, t, beta, -1, -1, i, j, t-1, distribution) == 1) flag = + 1;
					
	if(flag == 0) {
		for(int i = 0; i < l1; i++) 
			for(int j = 0; j < l2; j++) 
				for(int k = 0; k < t; k++) 
					latt[i][j][k] = latt[i][j][k] * ((cluster[i][j][k] ==  0) ? +1 : cluster[i][j][k]);
		boundary = boundary * (-1);
	}

	return boundary;
}
