/*Montecarlo per generare le superfici di separazione*/
 
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <random>
#include <iostream>
#include <sstream>
#include <fstream>
using namespace std;

#define N 1000

int clusterize(int*** latt, int*** cluster, int boundary, unsigned int l1, unsigned int l2, unsigned int t, long double beta, unsigned int i, unsigned int j, unsigned int k, uniform_real_distribution<long double>*  distribution);
int surfaceCluster(int*** latt,  int*** cluster, int boundary, unsigned int l1, unsigned int l2, unsigned int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_real_distribution<long double>*  distribution);
void singleCluster(int*** latt,  int*** cluster, int boundary, unsigned int l1, unsigned int l2, unsigned int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_int_distribution<int>* distr_t, uniform_real_distribution<long double>*  distribution);
void hotStart(int*** latt, unsigned int l1, unsigned int l2 , unsigned int t);
void coldStart(int*** latt, unsigned int l1, unsigned int l2, unsigned int t, int spin);
double energy(int*** latt, unsigned int l1, unsigned int l2, unsigned int t, int boundary);

random_device rd;
mt19937 mt(rd());

//NUmber of beta points
#define K 20

int main() {

	//E' necessario che questi vettori siano costanti perché si possano usare come dimensione di un ordinario array di c.
	//Questa cosa è in effetti non necessaria.	
	vector<long double> Beta(K);
	
	for(int i = 0; i < K; i++) Beta[i] = 0.1 + (1.5/K)*i; 
	
	const vector<unsigned int> L1 = {20};
	const vector<unsigned int> L2 = {20};
	const vector<unsigned int> T = {60};
	

	vector<double> surfaceEnergy(K);
	
	long double magn[N] = {0.0};

	
	for(int k = 0; k < T.size(); k++) {
		for(int j = 0; j < L1.size(); j++) {
			for (int h = 0; h < L2.size(); h++) {
				for(int i = 0; i < Beta.size(); i++) {
				
					//Questo è il posto ideale per inserire una eventuale parallelizzazione.
					long double beta = Beta[i];
					unsigned int l1 = L1[j];
					unsigned int l2 = L2[j];
					unsigned int t = T[k];

  					uniform_int_distribution<int> distr_l1(0,l1-1);
  					uniform_int_distribution<int> distr_l2(0,l2-1);
					uniform_int_distribution<int> distr_t(0,t-1);
					uniform_real_distribution<long double> distribution(0,1);


					for(int boundary = -1; boundary < 2; boundary += 2) {


						
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
				
					
						long double menergy = 0.0;
						for(long int rep = 0; rep < N; rep++) {
						
							cout << rep << "\r";
							singleCluster(latt, cluster, boundary, l1, l2, t, beta, &distr_l1, &distr_l2, &distr_t, &distribution);
							
						
						menergy += energy(latt, l1, l2, t, boundary);
						}

					
						menergy /= N;
					
					
						cout << ((boundary == -1) ? "Antiperiodic" : "Periodic") << "\t" << "l1 = " << l1 << "\t" "l2 = " << l2 << "\t" << "t = " << t << "\t"  << "beta = " << beta  << "\t" << "Energy = " << menergy << endl;
					
						if(boundary == -1)  surfaceEnergy[i] = -menergy;
						if(boundary == +1) surfaceEnergy[i] += menergy;
					
						//cout << msquare << string("---")  <<  string("---") << per <<  string("---") << aper << "---" << log(td) - log( 0.5 * log((1.0+ Aper/Per)/(1.0 - Aper/Per))) << endl;		
						//cout << msquare << string("---")  <<  string("---") << per <<  string("---") << aper << "---" << -(1/td) * log((1.0+ Aper/Per)/(1.0 - Aper/Per)) << endl;
						//cout << msquare << string("---")  <<  string("---") << per <<  string("---") << aper << "---" << log(Aper) - log(Per) + log(td) << endl;		
					
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
					
					/*for(int i = 0; i<K; i++) 
						cout << surfaceEnergy[i] << endl;*/
		
	
					vector<double> freeEnergy(K);
					freeEnergy[0] = 0.0;
	
					for(int i = 0; i < K; i++) {
						for(int j = 0; j <= i; j++)	{
							freeEnergy[i] += (Beta[1]-Beta[0])*surfaceEnergy[j];
						}	
					}
			
					ofstream myfile;
					
					std::stringstream ss;
					ss << "data" << l1 << "," << l2 << "," << t << ".txt";
					std::string fileName = ss.str();
					
  					myfile.open (fileName);
  					
  					for(int i = 0; i < K; i++) 
  						myfile << Beta[i] << "\t" << freeEnergy[i] << endl;
  						
  					myfile.close();
			
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

int clusterize(int*** latt, int*** latt, int*** cluster, int boundary, unsigned int l1, unsigned int l2, unsigned int t, long double beta, unsigned int i, unsigned int j, unsigned int k, uniform_real_distribution<long double>*  distribution) {

	unsigned int size[3] = {l1, l2, t};
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
					ret = clusterize(latt, cluster, boundary, l1, l2, t, beta, coord[0], coord[1], coord[2], distribution);
				}

			}

		}
	}
	
	return ret || flag;

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

	cluster[i][j][k] = -1;
	clusterize(latt, cluster, boundary, l1, l2, t, beta, i, j, k, distribution);
	
	/*cout << endl;
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

	
	for(int i = 0; i < l1; i++) {
		for(int j = 0; j < l2; j++) {
			for(int k = 0; k < t; k++) {
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
				cluster[i][j][t-1] = -1;
				if(clusterize(latt, cluster, boundary, l1, l2, t, beta, i, j, t-1, distribution) == 1) flag = + 1;
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


double energy(int*** latt, unsigned int l1, unsigned int l2, unsigned int t, int boundary) {

	double energy = 0.0;
	int bound = 1;
	unsigned int coord[3] = {0, 0, 0};
	unsigned int size[3] = {0, 0, 0};
	unsigned int b = 0;
	
	size[0] = l1; size[1] = l2; size[2] = t;
	
	for(int i = 0; i < l1; i++) {
		for(int j = 0; j < l2; j++) {
			for(int k = 0; k < t; k++) {
				for(int d = 0; d < 3; d++) {
					bound = 1;
					
					coord[0] = i; coord[1] = j; coord[2] = k;
					
					b = coord[d];
					if(coord[d] + 1 == size[d]) { 
						coord[d] = size[d] - (coord[d] + 1);
						if(d == 2) {
							bound = boundary;
						}
					} else coord[d] += 1;
					
					energy += bound * latt[coord[0]][coord[1]][coord[2]] * latt[i][j][k];
				}
			}
		}
	}
	
	return energy;
}

