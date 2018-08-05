#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <random>
using namespace std;

#define N 1000

void clusterize(int*** latt, int*** cluster, int boundary, unsigned int l , unsigned int t, long double beta, unsigned int i, unsigned int j, unsigned int k, int* flag, unsigned int* min);
int surfaceCluster(int*** latt, int boundary, unsigned int l , unsigned int t, long double beta);
void singleCluster(int*** latt, int boundary, unsigned int l , unsigned int t, long double beta);
void hotStart(int*** latt, unsigned int l , unsigned int t);
void coldStart(int*** latt, unsigned int l , unsigned int t, int spin);

int main() {

	//E' necessario che questi vettori siano costanti perché si possano usare come dimensione di un ordinario array di c.
	//Questa cosa è in effetti non necessaria.	
	const vector<long double> Beta = {0.2391};
	const vector<unsigned int> L = {4};
	const vector<unsigned int> T = {12};
	
	int boundary = 1;

	for(int i = 0; i < Beta.size(); i++) {
		for(int j = 0; j < L.size(); j++) {
			for(int k = 0; k < T.size(); k++) {

				//Questo è il posto ideale per inserire una eventuale parallelizzazione.
				long double beta = Beta[i];
				unsigned int l = L[j];
				unsigned int t = T[k];


				//I lattici non sono implementai come vettori ordinari perché tanto sono a dimensione costante.
				int*** latt;
  				latt = (int***) malloc(l * sizeof(int**) );
				for(int a = 0 ; a < l ; a++) {
					latt[a] = (int**) malloc (l * sizeof(int*));
					for(int b = 0; b < l; b++) {
						latt[a][b] = (int*) malloc (t * sizeof(int));
					}
				}
					
				hotStart(latt, l, t);
				long double magn[N] = {0.0};
				unsigned int per = 0; 
				unsigned int aper = 0;
				for(int rep = 0; rep < N; rep++) {
					long double somma = 0.0;
					
					//singleCluster(latt, boundary, l, t, beta);
					boundary = surfaceCluster(latt, boundary, l, t, beta);

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
					
					long double msquare = 0.0;
					for(int i = 0; i<N; i++) msquare += magn[i] * magn[i];
					msquare = msquare/per;	

					cout << msquare << string("---")  <<  per <<  string("---") <<  string("---") << aper << endl;
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

void clusterize(int*** latt, int*** cluster, int boundary, unsigned int l , unsigned int t, long double beta, unsigned int i, unsigned int j, unsigned int k, int* flag, unsigned int* min, int crossed) {

	unsigned int size[3] = {l, l, t};
	unsigned int coord[3] = {i, j, k};
	long double p = 0.0;
	int cross = 0;
	int bound = +1;
	int b = 0;

	//cout << coord[0] << " " << coord[1] << " " << coord[2] << endl;

	random_device rd;
    	mt19937 mt(rd());
  	uniform_real_distribution<long double> distribution(0,1);

	for(int d = 0; d < 3; d++) {
		for(int a = -1; a < 2; a = a + 2) {
			cross = 0;
			bound = +1;
			coord[0] = i; coord[1] = j; coord[2] = k;
			b = coord[d];
			//TODO: rivedere questo controllo, perché qualcosa non va.
			if((coord[d] + a == size[d]) || (coord[d] + a == -1)) { 
				b = size[d] - 1 - (coord[d] + a);
				if(d == 2) {
					bound = boundary;
					cross = 1;
				}
			}
			coord[d] = a+b;
			p = 1.0-exp(-beta * (1.0 + bound*latt[i][j][k]*latt[coord[0]][coord[1]][coord[2]]));
			if(distribution(mt) < p) {
				if(!crossed && (cluster[coord[0]][coord[1]][coord[2]] == -1)) *flag = 1; 
				if ((cluster[coord[0]][coord[1]][coord[2]] == 1)) {
					if((!crossed) and (coord[2] < *min)) *min = coord[2];
					cout << *min << endl;
					cluster[coord[0]][coord[1]][coord[2]] = -1;
					clusterize(latt, cluster, boundary, l , t, beta, coord[0], coord[1], coord[2], flag, min, crossed || cross);
				}
			}
		}
	}

}

void singleCluster(int*** latt, int boundary, unsigned int l , unsigned int t, long double beta) {

	random_device rd;
    	mt19937 mt(rd());
  	uniform_int_distribution<int> distr_l(0,l-1);
	uniform_int_distribution<int> distr_t(0,t-1);

	int flag = 0;
	unsigned int min = t-1;


	int*** cluster;
  	cluster = (int***) malloc(l * sizeof(int**) );
	for(int a = 0 ; a < l ; a++) {
		cluster[a] = (int**) malloc (l * sizeof(int*));
		for(int b = 0; b < l; b++) {
			cluster[a][b] = (int*) malloc (t * sizeof(int));
		}
	}

	for(int i = 0; i < l; i++) {
		for(int j = 0; j < l; j++) {
			for(int k = 0; k < t; k++) {
				cluster[i][j][k] = 1;
			}
		}
	}

	unsigned int i = distr_l(mt);
	unsigned int j = distr_l(mt);
	unsigned int k = distr_t(mt);

	cluster[i][j][k] = -1;
	clusterize(latt, cluster, boundary, l, t, beta, i, j, k, &flag, &min);


	for(int i = 0; i < l; i++) {
		for(int j = 0; j < l; j++) {
			for(int k = 0; k < t; k++) {
				latt[i][j][k] = latt[i][j][k] * cluster[i][j][k];
			}
		}
	}			
}

int surfaceCluster(int*** latt, int boundary, unsigned int l , unsigned int t, long double beta) {

	random_device rd;
    	mt19937 mt(rd());
  	uniform_int_distribution<int> distr_l(0,l-1);
	uniform_int_distribution<int> distr_t(0,t-1);

	int flag = 0;
	unsigned int min = t-1;

	int*** cluster;
  	cluster = (int***) malloc(l * sizeof(int**) );
	for(int a = 0 ; a < l ; a++) {
		cluster[a] = (int**) malloc (l * sizeof(int*));
		for(int b = 0; b < l; b++) {
			cluster[a][b] = (int*) malloc (t * sizeof(int));
		}
	}

	for(int i = 0; i < l; i++) {
		for(int j = 0; j < l; j++) {
			for(int k = 0; k < t; k++) {
				cluster[i][j][k] = 1;
			}
		}
	}


	//Scelgo sequenzialmente l'elemento iniziale dei cluster sul layer estremo	
	
	for(int i = 0; i < l; i++) {
		for(int j = 0; j < l; j++) {
			if(cluster[i][j][t-1] == 1) {
				cluster[i][j][t-1] = -1;
				clusterize(latt, cluster, boundary, l, t, beta, i, j, t-1, &flag, &min);
			}
		}
	}

	for(int i = 0; i < l; i++) {
		for(int j = 0; j < l; j++) {
			for(int k = 0; k < t; k++) {
				cout << cluster[i][j][k] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	
	cout << flag << " " << min << endl;
	/*
	int nonclusterized = 1;

	while(nonclusterized == 1) {

		unsigned int i = 0;
		unsigned int j = 0;

		do {
			i = distr_l(mt);
			j = distr_l(mt);
		} while(cluster[i][j][t-1] == -1);

		cluster[i][j][t-1] = -1;
		clusterize(latt, cluster, boundary, l, t, beta, i, j, t-1);

		nonclusterized = 0;
		for(int i = 0; i < l; i++) {
			for(int j = 0; j < l; j++) {
				if(cluster[i][j][t-1] == 1) {
					nonclusterized = 1;
					break;
				}
			}
		}
	}*/

	/*

	int sheet = 0;
	for(int k = 0; k < t; k++) {
		sheet = 1;
		for(int i = 0; i < l; i++) {
			for(int j = 0; j < l; j++) {
				if(cluster[i][j][k] == -1) {
					sheet = 0;
					break;
				}
			}
			if(sheet == 0) break;
		}	
		if(sheet ==  1) {
			for(int i = 0; i < l; i++) {
				for(int j = 0; j < l; j++) {
					
					for(int s = k+1; s < t; s++) {
						latt[i][j][s] = latt[i][j][s] * cluster[i][j][s];
						//latt[i][j][s] = latt[i][j][s] * (-1);
					}
				}
			}
			boundary = boundary * (-1);
			break;
		}
	}*/

	if(flag == 0) {
		for(int i = 0; i < l; i++) {
				for(int j = 0; j < l; j++) {
					for(int s = min; s < t; s++) {
						latt[i][j][s] = latt[i][j][s] * cluster[i][j][s];
						//latt[i][j][s] = latt[i][j][s] * (-1);
					}
				}
			}
		boundary = boundary * (-1);
	}

	return boundary;	
}

