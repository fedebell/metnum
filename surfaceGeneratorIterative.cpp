/*Montecarlo per generare le superfici di separazione*/
 
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <random>
#include <stack>
#include <ctime>
using namespace std;

//Inserire il numero di iterazioni
#define N 10000

int clusterize(int*** latt, int*** cluster, int boundary, int* size, long double beta, int i, int j, int k, uniform_real_distribution<long double>*  distribution);
int surfaceCluster(int*** latt,  int*** cluster, int boundary, int l1, int l2, int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_real_distribution<long double>*  distribution);
void singleCluster(int*** latt,  int*** cluster, int boundary, int l1, int l2, int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_int_distribution<int>* distr_t, uniform_real_distribution<long double>*  distribution);
void hotStart(int*** latt, int l1, int l2 , int t);
void coldStart(int*** latt, int l1, int l2, int t, int spin);

random_device rd;
mt19937 mt(rd());

typedef struct site {
	int c[3];
	int precedente;
	int cross;
} site;

int main() {
	
	//Inserire i parametri delle simulazioni, verranno eseguite tutte le possibili combinazioni a meno di aggiunta di condizioni.
	const vector<long double> Beta = {0.224};
	const vector<int> L1 = {10};
	const vector<int> L2 = {10};
	const vector<int> T = {30};
	
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
						cout << rep << "\r";
						
						//singleCluster(latt, cluster, boundary, l1, l2, t, beta, &distr_l1, &distr_l2, &distr_t, &distribution);
						boundary = surfaceCluster(latt, cluster, boundary, l1, l2, t, beta, &distr_l1, &distr_l2, &distribution);
						if(boundary == 1) per += 1;
						else aper += 1;
						
						boundary = surfaceCluster(latt, cluster, boundary, l1, l2, t, beta, &distr_l1, &distr_l2, &distribution);
						if(boundary == 1) per += 1;
						else aper += 1;
					}
					
					int fine = clock() - init;
					
					double Aper = double(aper);
					double Per = double(per);
					double td = double(t);
					
					cout << "l1 = " << l1 << "\t" "l2 = " << l2 << "\t" << "t = " << t << "\t"  << "beta = " << beta << "\t" << "tempo = " << fine << endl;
					cout  << "Aper/Tot = " << Aper/(2*N) << "\t" << "Aper/Per = " << Aper/Per << "\t" << "F_s = " << -log(Aper) + log(Per) + log(td) << "\t" << "F_si = " << log(td) - log( 0.5 * log((1.0+ Aper/Per)/(1.0 - Aper/Per))) << endl;	
					
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

int clusterize(int*** latt, int*** cluster, int boundary, int* size, long double beta, int i, int j, int k, uniform_real_distribution<long double>*  distribution) {

	stack<site> stack;
	
	site next = {{i, j, k}, -1, -1};
	stack.push(next);
	
	int flag = 0, value = 0, d = 0, a = 0, s = 0;

	while(!stack.empty()) {
		site current = stack.top();
		stack.pop();
		value = cluster[current.c[0]][current.c[1]][current.c[2]];
		if(value == (current.cross * current.precedente)) flag = 1;
		else if(value == 0) {
			value =  current.precedente * (-current.cross);
			cluster[current.c[0]][current.c[1]][current.c[2]] = value;
			next.precedente = value;
			for(d = 0; d < 3; d++) {
				for(a = -1; a < 2; a = a + 2) {
					next.cross = -1;
					for(s = 0; s<3; s++) next.c[s] = current.c[s]; 
					next.c[d] +=  a;
					if(next.c[d] == size[d]) { next.c[d] = 0; if(d == 2) next.cross = 1;}
					if(next.c[d] == -1) { next.c[d] = size[d] - 1; if(d == 2) next.cross = 1;}
					if((cluster[next.c[0]][next.c[1]][next.c[2]] == 0) && ((*distribution)(mt) < (1.0 - exp(-beta * (1.0 + ((next.cross == 1) ? boundary : 1) * latt[current.c[0]][current.c[1]][current.c[2]]*latt[next.c[0]][next.c[1]][next.c[2]]))))) stack.push(next);
				}
			}
		}
	}
	return flag;
}

void singleCluster(int*** latt,  int*** cluster, int boundary, int l1, int l2, int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_int_distribution<int>* distr_t, uniform_real_distribution<long double>*  distribution) {

	int size[3] = {l1, l2, t};

	for(int i = 0; i < l1; i++) 
		for(int j = 0; j < l2; j++) 
			for(int k = 0; k < t; k++) 
				cluster[i][j][k] = 0;
				
	int i = (*distr_l1)(mt);
	int j = (*distr_l2)(mt);
	int k = (*distr_t)(mt);

	clusterize(latt, cluster, boundary, size, beta, i, j, k, distribution);	
	
	for(int i = 0; i < l1; i++) 
		for(int j = 0; j < l2; j++) 
			for(int k = 0; k < t; k++) 
				latt[i][j][k] = latt[i][j][k] * ((cluster[i][j][k] == 0) ? 1 : -1);
}

int surfaceCluster(int*** latt,  int*** cluster, int boundary, int l1, int l2, int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_real_distribution<long double>*  distribution) {

	int size[3] = {l1, l2, t};
	int flag = 0;

	for(int i = 0; i < l1; i++) 
		for(int j = 0; j < l2; j++) 
			for(int k = 0; k < t; k++) 
				cluster[i][j][k] = 0;
	
	for(int i = 0; i < l1; i++) 
		for(int j = 0; j < l2; j++) 
			if(cluster[i][j][t-1] == 0) 
				if(clusterize(latt, cluster, boundary, size, beta, i, j, t-1, distribution) == 1) flag = + 1;
					
	if(flag == 0) {
		for(int i = 0; i < l1; i++) 
			for(int j = 0; j < l2; j++) 
				for(int k = 0; k < t; k++) 
					latt[i][j][k] = latt[i][j][k] * ((cluster[i][j][k] ==  0) ? +1 : cluster[i][j][k]);
		boundary = boundary * (-1);
	}

	return boundary;
}

