/*Montecarlo per generare le superfici di separazione*/
 
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <random>
using namespace std;

#define N 100000

int clusterize(int*** latt, int*** cluster, int boundary, unsigned int l1, unsigned int l2, unsigned int t, long double beta, unsigned int i, unsigned int j, unsigned int k, uniform_real_distribution<long double>*  distribution);
int surfaceCluster(int*** latt,  int*** cluster, int boundary, unsigned int l1, unsigned int l2, unsigned int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_real_distribution<long double>*  distribution);
void singleCluster(int*** latt,  int*** cluster, int boundary, unsigned int l1, unsigned int l2, unsigned int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_int_distribution<int>* distr_t, uniform_real_distribution<long double>*  distribution);
void hotStart(int*** latt, unsigned int l1, unsigned int l2 , unsigned int t);
void coldStart(int*** latt, unsigned int l1, unsigned int l2, unsigned int t, int spin);

random_device rd;
mt19937 mt(rd());

int main() {

	//E' necessario che questi vettori siano costanti perché si possano usare come dimensione di un ordinario array di c.
	//Questa cosa è in effetti non necessaria.	
	const vector<long double> Beta = {0.2391};
	const vector<unsigned int> L1 = {10};
	const vector<unsigned int> L2 = {10};
	const vector<unsigned int> T = {10};
	
	int boundary = 1;

	long double magn[2*N] = {0.0};

	for(int i = 0; i < Beta.size(); i++) {
		for(int j = 0; j < L1.size(); j++) {
			for (int h = 0; h < L2.size(); h++) {
				for(int k = 0; k < T.size(); k++) {


					//Questo è il posto ideale per inserire una eventuale parallelizzazione.
					long double beta = Beta[i];
					unsigned int l1 = L1[j];
					unsigned int l2 = L2[j];
					unsigned int t = T[k];

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
						
					coldStart(latt, l1, l2, t, 1);

					for(int i = 0; i<N; i++) magn[i] = 0.0;
				
					unsigned int per = 0; 
					unsigned int aper = 0;
					long double msquare = 0.0;
					for(long int rep = 0; rep < N; rep++) {
						msquare = 0.0;
						long double somma = 0.0;
						cout << rep << "\r";
						
						singleCluster(latt, cluster, boundary, l1, l2, t, beta, &distr_l1, &distr_l2, &distr_t, &distribution);
						
						if(boundary == 1) {
							per += 1;
							for(int i = 0; i < l1; i++) {
								for(int j = 0; j < l2; j++) {
									for(int k = 0; k < t; k++) {
										somma += latt[i][j][k];
									}
								}
							}
							magn[rep] = somma/(l1*l2*t);
						} else aper += 1;
						
						boundary = surfaceCluster(latt, cluster, boundary, l1, l2, t, beta, &distr_l1, &distr_l2, &distribution);
						//boundary = surfaceCluster(latt, cluster, boundary, l1, l2, t, beta, &distr_l1, &distr_l2, &distribution);
						
						if(boundary == 1) {
							per += 1;
							for(int i = 0; i < l1; i++) {
								for(int j = 0; j < l2; j++) {
									for(int k = 0; k < t; k++) {
										somma += latt[i][j][k];
									}
								}
							}
							magn[rep] = somma/(l1*l2*t);
						} else aper += 1;
					}
					double Aper = double(aper);
					double Per = double(per);
					double td = double(t);
					
					for(int i = 0; i<2*N; i++) msquare += magn[i] * magn[i];
					msquare = msquare/per;
					cout << "l1 = " << l1 << "\t" "l2 = " << l2 << "\t" << "t = " << t << "\t"  << "beta = " << beta << endl;
					//cout << msquare << string("---")  <<  string("---") << per <<  string("---") << aper << "---" << log(td) - log( 0.5 * log((1.0+ Aper/Per)/(1.0 - Aper/Per))) << endl;	
					//cout << msquare << string("---")  <<  string("---") << per <<  string("---") << aper << "---" << -(1/td) * log((1.0+ Aper/Per)/(1.0 - Aper/Per)) << endl;
					cout << msquare << string("---")  <<  string("---") << per <<  string("---") << aper << "---" << Aper/Per << "---" << -log(Aper) + log(Per) + log(td) << endl;	
					//cout << msquare << string("---")  <<  string("---") << per <<  string("---") << aper << "---" << Aper/Per << endl;		
					
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
	
				}
			}
		}
	}
				
	
	return 0;
}

void coldStart(int*** latt, unsigned int l1, unsigned int l2, unsigned int t, int spin) {

	for(int i = 0; i < l1; i++) {
		for(int j = 0; j < l2; j++) {
			for(int k = 0; k < t; k++) {
				latt[i][j][k] = spin;
			}
		}
	}

}

void hotStart(int*** latt, unsigned int l1, unsigned int l2, unsigned int t) {

	random_device rd;
    	mt19937 mt(rd());
  	uniform_int_distribution<int> distribution(0,1);

	for(int i = 0; i < l1; i++) {
		for(int j = 0; j < l2; j++) {
			for(int k = 0; k < t; k++) {
				if(distribution(mt) == 0) {
 					latt[i][j][k] = 1;
				} else latt[i][j][k] = -1;
			}
		}
	}

}

int clusterize(int*** latt, int*** cluster, int boundary, unsigned int l1, unsigned int l2, unsigned int t, long double beta, int precedente, int cross, unsigned int i, unsigned int j, unsigned int k, uniform_real_distribution<long double>*  distribution) {

	unsigned int size[3] = {l1, l2, t};
	unsigned int coord[3] = {i, j, k};
	int ret[6] = {0, 0, 0, 0, 0, 0};
	
	int toClusterize[6][4];
	
	//Inizializzazione della matrice
	for(int i = 0; i < 6; i++) 
		for (int j = 0; j < 4; j++) 
			toClusterize[i][j] = -1;
	
	long double p = 0.0;
	int flag = 0;
	int bound = +1;
	int newCross = -1;
	int b = 0;
	
	int count = 0;

	
	if(cluster[i][j][k] == (cross * precedente)) {
		/*cout << "Hello" << endl;
		
		
			for(int i = 0; i < l1; i++) {
			for(int j = 0; j < l2; j++) {
				for(int k = 0; k < t; k++) {
					cout << cluster[i][j][k] << " ";
				}
				cout << endl;
			}
			cout << endl;
			cout << endl;
			}*/
			
		return 1;
	
	}
	
	else if(cluster[i][j][k] == 0) {
	
		//Se non c'è cross cross vale -1, se c'è vale +1
		//cout << cluster[i][j][k] << precedente << cross << endl;
		cluster[i][j][k] = precedente * (-cross);
		//cout << cluster[i][j][k] << precedente << cross << endl;
		//cout << "Cluster:" << cluster[i][j][k] << endl;
		
		for(int d = 0; d < 3; d++) {
			for(int a = -1; a < 2; a = a + 2) {
				bound = +1;
				newCross = -1;
				coord[0] = i; coord[1] = j; coord[2] = k;
				b = coord[d];
				if((coord[d] + a == size[d]) || (coord[d] + a == -1)) { 
					b = size[d] - 1 - (coord[d] + a);
					if(d == 2) {
						bound = boundary;
						newCross = 1;
					}
				}

				coord[d] = a+b;
				
				//cout << "Cluster:" << cluster[i][j][k] << endl;
				
				if(cluster[coord[0]][coord[1]][coord[2]] == 0) {
				
					p = 1.0 - exp(-beta * (1.0 + bound*latt[i][j][k]*latt[coord[0]][coord[1]][coord[2]]));
			
					if ((*distribution)(mt) < p) {
					
						//Decide che questi sono da clusterizzare ma lo fa dopo sennò è un errore grave.
						//ret[a*(d+1)+3] = clusterize(latt, cluster, boundary, l1, l2, t, beta, cluster[i][j][k], newCross, coord[0], coord[1], coord[2], distribution);
						//cout << ret[a*(d+1)+3] << endl;
						
						toClusterize[count][0] = coord[0];
						toClusterize[count][1] = coord[1];
						toClusterize[count][2] = coord[2];
						toClusterize[count][3] = newCross;
						count++;
						
						//cout << "Cluster:" << cluster[i][j][k] << endl;
					
					}

				}

			}
		}
		
		for(int c = 0; c < count; c++) {
			//cout << "Cluster:" << cluster[i][j][k] << endl;
			
			
			/*if(cluster[i][j][k] == 0) {
				cout << "HELP: " << c << " - " <<  i << " - " << j << " - " << k << " - " <<  toClusterize[c][3] << " - " <<  toClusterize[c][0] << " - " <<  toClusterize[c][1] << " - " <<  toClusterize[c][2] << endl;
			}
			
			for(int i = 0; i < l1; i++) {
			for(int j = 0; j < l2; j++) {
				for(int k = 0; k < t; k++) {
					cout << cluster[i][j][k] << " ";
				}
				cout << endl;
			}
			cout << endl;
			cout << endl;
			}*/

			if((toClusterize[c][0]) != -1) ret[c] = clusterize(latt, cluster, boundary, l1, l2, t, beta, cluster[i][j][k], toClusterize[c][3], toClusterize[c][0], toClusterize[c][1], toClusterize[c][2], distribution);
		}
		
		for(int i = 0; i<6; i++) flag = (flag || ret[i]);
	}
	
	return flag;

}

void singleCluster(int*** latt,  int*** cluster, int boundary, unsigned int l1, unsigned int l2, unsigned int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_int_distribution<int>* distr_t, uniform_real_distribution<long double>*  distribution) {

	for(int i = 0; i < l1; i++) {
		for(int j = 0; j < l2; j++) {
			for(int k = 0; k < t; k++) {
				cluster[i][j][k] = 0;
			}
		}
	}

	unsigned int i = (*distr_l1)(mt);
	unsigned int j = (*distr_l2)(mt);
	unsigned int k = (*distr_t)(mt);

	//clusterize(latt, cluster, boundary, l1, l2, t, beta, -latt[i][j][k], -1, i, j, k, distribution);
	clusterize(latt, cluster, boundary, l1, l2, t, beta, -1, -1, i, j, k, distribution);

	
	
	for(int i = 0; i < l1; i++) {
		for(int j = 0; j < l2; j++) {
			for(int k = 0; k < t; k++) {
				//Penso che la prima condizione sia superflua. Ma vediamo con calma.
				if(cluster[i][j][k] == 1) cluster[i][j][k] = -1;
				if(cluster[i][j][k] == 0) cluster[i][j][k] = 1;
				latt[i][j][k] = latt[i][j][k] * cluster[i][j][k];
			}
		}
	}	
	
}

int surfaceCluster(int*** latt,  int*** cluster, int boundary, unsigned int l1, unsigned int l2, unsigned int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_real_distribution<long double>*  distribution) {

	int flag = 0;

	for(int i = 0; i < l1; i++) {
		for(int j = 0; j < l2; j++) {
			for(int k = 0; k < t; k++) {
				cluster[i][j][k] = 0;
			}
		}
	}
	
	
	for(int i = 0; i < l1; i++) {
		for(int j = 0; j < l2; j++) {
			if(cluster[i][j][t-1] == 0) {
				//if(clusterize(latt, cluster, boundary, l1, l2, t, beta, -latt[i][j][t-1], -1, i, j, t-1, distribution) == 1) flag = + 1;
				if(clusterize(latt, cluster, boundary, l1, l2, t, beta, -1, -1, i, j, t-1, distribution) == 1) flag = + 1;
			}
		}
	}
	
		/*int nonclusterized = 1;

	while(nonclusterized == 1) {

		unsigned int i = 0;
		unsigned int j = 0;

		do {
			i = (*distr_l1)(mt);
			j = (*distr_l2)(mt);
		} while(cluster[i][j][t-1] != 0);

		cluster[i][j][t-1] = -1;
		if(clusterize(latt, cluster, boundary, l1, l2, t, beta, i, j, t-1, distribution) == 1) flag = + 1;

		nonclusterized = 0;
		for(int i = 0; i < l1; i++) {
			for(int j = 0; j < l2; j++) {
				if(cluster[i][j][t-1] == 0) {
					nonclusterized = 1;
					break;
				}
			}
		}
	}*/
	
					
	if(flag == 0) {
	
		for(int i = 0; i < l1; i++) {
			for(int j = 0; j < l2; j++) {
				for(int k = 0; k < t; k++) {
					if(cluster[i][j][k] == 0) cluster[i][j][k] = +1;
					latt[i][j][k] = latt[i][j][k] * cluster[i][j][k];
				}
			}
		}

		boundary = boundary * (-1);
	}
	
	/*
	cout << endl;
	cout << endl;
	cout << "Boundary: " << ((boundary == +1) ? "periodic" : "Aperiodic") << endl;
	int count = 0;
	for(int i = 0; i < l1; i++) {
		for(int k = 0; k < t; k++) {
			count = 0;
			for(int j = 0; j < l2; j++) {
				count += ((latt[i][j][k] > 0) ? 1 : 0);
			}
			cout << ((count > (l2/2)) ? "o" : "-") << " ";
		}
		cout << endl;
	}*/
	

	return boundary;
}

	

	/*
	for(int i = 0; i < l; i++) {
		for(int j = 0; j < l; j++) {
			for(int k = 0; k < t; k++) {
				cout << latt[i][j][k] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}*/


