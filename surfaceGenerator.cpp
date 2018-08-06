/*MOntecarlo per generare le superfici di separazione*/

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <random>
using namespace std;

#define N 10000

int clusterize(int*** latt, int*** cluster, int boundary, unsigned int l , unsigned int t, long double beta, unsigned int i, unsigned int j, unsigned int k, uniform_real_distribution<long double>*  distribution);
int surfaceCluster(int*** latt,  int*** cluster, int boundary, unsigned int l , unsigned int t, long double beta, uniform_int_distribution<int>* distr_l, uniform_real_distribution<long double>*  distribution);
void singleCluster(int*** latt,  int*** cluster, int boundary, unsigned int l , unsigned int t, long double beta, uniform_int_distribution<int>* distr_l, uniform_int_distribution<int>* distr_t, uniform_real_distribution<long double>*  distribution);
void hotStart(int*** latt, unsigned int l , unsigned int t);
void coldStart(int*** latt, unsigned int l , unsigned int t, int spin);

random_device rd;
mt19937 mt(rd());


int main() {

	//E' necessario che questi vettori siano costanti perché si possano usare come dimensione di un ordinario array di c.
	//Questa cosa è in effetti non necessaria.	
	const vector<long double> Beta = {0.26};
	const vector<unsigned int> L = {10};
	const vector<unsigned int> T = {30};
	
	int boundary = 1;

	long double magn[N] = {0.0};

	for(int i = 0; i < Beta.size(); i++) {
		for(int j = 0; j < L.size(); j++) {
			for(int k = 0; k < T.size(); k++) {


				//Questo è il posto ideale per inserire una eventuale parallelizzazione.
				long double beta = Beta[i];
				unsigned int l = L[j];
				unsigned int t = T[k];

  				uniform_int_distribution<int> distr_l(0,l-1);
				uniform_int_distribution<int> distr_t(0,t-1);
				uniform_real_distribution<long double> distribution(0,1);

				//I lattici non sono implementai come vettori ordinari perché tanto sono a dimensione costante.
				int*** latt;
  				latt = (int***) malloc(l * sizeof(int**) );
				for(int a = 0 ; a < l ; a++) {
					latt[a] = (int**) malloc (l * sizeof(int*));
					for(int b = 0; b < l; b++) {
						latt[a][b] = (int*) malloc (t * sizeof(int));
					}
				}

			
				int*** cluster;
  				cluster = (int***) malloc(l * sizeof(int**) );
				for(int a = 0 ; a < l ; a++) {
					cluster[a] = (int**) malloc (l * sizeof(int*));
					for(int b = 0; b < l; b++) {
						cluster[a][b] = (int*) malloc (t * sizeof(int));
					}
				}
					
				coldStart(latt, l, t, 1);

				for(int i = 0; i<N; i++) magn[i] = 0.0;
				
				unsigned int per = 0; 
				unsigned int aper = 0;
				long double msquare = 0.0;
				for(long int rep = 0; rep < N; rep++) {
					msquare = 0.0;
					long double somma = 0.0;
					
					singleCluster(latt, cluster, boundary, l, t, beta, &distr_l, &distr_t, &distribution);
					if(rep%10 == 0) boundary = surfaceCluster(latt, cluster, boundary, l, t, beta, &distr_l, &distribution);
					if(boundary == 1) {
						per += 1;
						for(int i = 0; i < l; i++) {
							for(int j = 0; j < l; j++) {
								for(int k = 0; k < t; k++) {
									somma += latt[i][j][k];
								}
							}
						}
						magn[rep] = somma/(l*l*t);
					} else aper += 1;
					for(int i = 0; i<N; i++) msquare += magn[i] * magn[i];
					msquare = msquare/per;
					//cout << "l = " << l << "---" << "t = " << t << "---"  << "beta = " << beta << endl;
					cout << rep << "---" << msquare << string("---")  <<  string("---") << per <<  string("---") << aper << endl;	
				}
				
				
			}
		}
	}

	return 0;
}

void coldStart(int*** latt, unsigned int l , unsigned int t, int spin) {

	for(int i = 0; i < l; i++) {
		for(int j = 0; j < l; j++) {
			for(int k = 0; k < t; k++) {
				latt[i][j][k] = spin;
			}
		}
	}

}

void hotStart(int*** latt, unsigned int l , unsigned int t) {

	random_device rd;
    	mt19937 mt(rd());
  	uniform_int_distribution<int> distribution(0,1);

	for(int i = 0; i < l; i++) {
		for(int j = 0; j < l; j++) {
			for(int k = 0; k < t; k++) {
				if(distribution(mt) == 0) {
 					latt[i][j][k] = 1;
				} else latt[i][j][k] = -1;
			}
		}
	}

}

int clusterize(int*** latt, int*** cluster, int boundary, unsigned int l , unsigned int t, long double beta, unsigned int i, unsigned int j, unsigned int k, uniform_real_distribution<long double>*  distribution) {

	unsigned int size[3] = {l, l, t};
	unsigned int coord[3] = {i, j, k};
	long double p = 0.0;
	int flag = 0;
	int ret = 0;
	int bound = +1;
	int cross = -1;
	int b = 0;

	for(int d = 0; d < 3; d++) {
		for(int a = -1; a < 2; a = a + 2) {
			bound = +1;
			cross = -1;
			coord[0] = i; coord[1] = j; coord[2] = k;
			b = coord[d];
			if((coord[d] + a == size[d]) || (coord[d] + a == -1)) { 
				b = size[d] - 1 - (coord[d] + a);
				if(d == 2) {
					bound = boundary;
					cross = 1;
				}
			}

			coord[d] = a+b;
			p = 1.0 - exp(-beta * (1.0 + bound*latt[i][j][k]*latt[coord[0]][coord[1]][coord[2]]));
			
			if ((*distribution)(mt) < p) {
	
				if((cluster[coord[0]][coord[1]][coord[2]] * cluster[i][j][k]) == cross) flag = +1;
				if(cluster[coord[0]][coord[1]][coord[2]] == 0) {
					if(cross == 1) {
						cluster[coord[0]][coord[1]][coord[2]] = - cluster[i][j][k];
					} else cluster[coord[0]][coord[1]][coord[2]] = cluster[i][j][k];
					ret = clusterize(latt, cluster, boundary, l , t, beta, coord[0], coord[1], coord[2], distribution);
				}

			}

		}
	}
	
	return ret || flag;

}

void singleCluster(int*** latt,  int*** cluster, int boundary, unsigned int l , unsigned int t, long double beta, uniform_int_distribution<int>* distr_l, uniform_int_distribution<int>* distr_t, uniform_real_distribution<long double>*  distribution) {

	for(int i = 0; i < l; i++) {
		for(int j = 0; j < l; j++) {
			for(int k = 0; k < t; k++) {
				cluster[i][j][k] = 0;
			}
		}
	}

	unsigned int i = (*distr_l)(mt);
	unsigned int j = (*distr_l)(mt);
	unsigned int k = (*distr_t)(mt);

	cluster[i][j][k] = -1;
	clusterize(latt, cluster, boundary, l, t, beta, i, j, k, distribution);

	
	for(int i = 0; i < l; i++) {
		for(int j = 0; j < l; j++) {
			for(int k = 0; k < t; k++) {
				if(cluster[i][j][k] == 1) cluster[i][j][k] = -1;
				if(cluster[i][j][k] == 0) cluster[i][j][k] = 1;
				latt[i][j][k] = latt[i][j][k] * cluster[i][j][k];
			}
		}
	}	
	
}

int surfaceCluster(int*** latt,  int*** cluster, int boundary, unsigned int l , unsigned int t, long double beta, uniform_int_distribution<int>* distr_l, uniform_real_distribution<long double>*  distribution) {

	int flag = 0;

	for(int i = 0; i < l; i++) {
		for(int j = 0; j < l; j++) {
			for(int k = 0; k < t; k++) {
				cluster[i][j][k] = 0;
			}
		}
	}

	for(int i = 0; i < l; i++) {
		for(int j = 0; j < l; j++) {
			if(cluster[i][j][t-1] == 0) {
				cluster[i][j][t-1] = -1;
				if(clusterize(latt, cluster, boundary, l, t, beta, i, j, t-1, distribution) == 1) flag = + 1;
			}
		}
	}	

	if(flag == 0) {
	
		for(int i = 0; i < l; i++) {
			for(int j = 0; j < l; j++) {
				for(int k = 0; k < t; k++) {
					if(cluster[i][j][k] == 0) cluster[i][j][k] = +1;
					latt[i][j][k] = latt[i][j][k] * cluster[i][j][k];
				}
			}
		}

		boundary = boundary * (-1);
	}

	
	cout << endl;
	cout << endl;
	cout << "Boundary: " << ((boundary == +1) ? "periodic" : "Aperiodic") << endl;
	int count = 0;
	for(int i = 0; i < l; i++) {
		for(int k = 0; k < t; k++) {
			count = 0;
			for(int j = 0; j < l; j++) {
				count += ((latt[i][j][k] > 0) ? 1 : 0);
			}
			cout << ((count > (l/2)) ? "o" : "-") << " ";
		}
		cout << endl;
	}
	

	return boundary;
}
	
	/* int nonclusterized = 1;

	while(nonclusterized == 1) {

		unsigned int i = 0;
		unsigned int j = 0;

		do {
			i = (*distr_l)(mt);
			j = (*distr_l)(mt);
		} while(cluster[i][j][t-1] != 0);

		cluster[i][j][t-1] = -1;
		if(clusterize(latt, cluster, boundary, l, t, beta, i, j, t-1) == 1) flag = + 1;

		nonclusterized = 0;
		for(int i = 0; i < l; i++) {
			for(int j = 0; j < l; j++) {
				if(cluster[i][j][t-1] == 0) {
					nonclusterized = 1;
					break;
				}
			}
		}
	}*/

	

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


