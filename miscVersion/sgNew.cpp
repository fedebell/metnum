/*Montecarlo per generare le superfici di separazione*/

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <stack>
#include <unistd.h>
#include <signal.h>
#include <sys/wait.h>
#include <cstring>
#include <cmath>
#include <sys/prctl.h>

using namespace std;

//Inserire il numero di iterazioni
#define N 1000000
//Percentuale iniziale rimossa dalla catena di Markov perché si considera termalizzazione.
#define frac 5
//Variabili che controllano il jackknife
#define start 10000
#define lenght 30
#define step 30

int clusterize(int*** latt, int*** cluster, int boundary, int l1, int l2, int t, long double beta, int i, int j, int k, uniform_real_distribution<long double>*  distribution);
int surfaceCluster(int*** latt,  int*** cluster, int boundary, int l1, int l2, int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_real_distribution<long double>*  distribution);
void singleCluster(int*** latt,  int*** cluster, int boundary, int l1, int l2, int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_int_distribution<int>* distr_t, uniform_real_distribution<long double>*  distribution);
void hotStart(int*** latt, int l1, int l2 , int t);
void coldStart(int*** latt, int l1, int l2, int t, int spin);
double jackknife(int* boundary, int size, int* blockDimentions, int len, double t);

typedef struct {
	int c[3];
	int precedente;
	int cross;
} site;

typedef struct {
	int c[3];
} siteS;


random_device rd;
mt19937 mt(rd());

int main(int argc, char *argv[]) {
	//Inseririmento di parametri da riga di comando
	
	vector<long double> Beta;
	vector<int> L;
	
	bool parall = false;
	double temp = 0.0;
	
	//Il primo comando da leggere è se fare le cose in parallelo o sequenzialmente.
	
	if(strcmp(argv[1], "-s") == 0)  parall = false;
	else if(strcmp(argv[1], "-p") == 0) parall = true;
	else {
		cerr << "Il primo argomento deve essere '-s' per l'esecuzione sequenziale o '-p' per l'esecuzione parallela." << endl;
		return -1;
	} 
	
	for (int i = 1; i < argc - 1 ; i++) {
		temp = atof(argv[i+1]);
		if(temp < 1.0) Beta.push_back(double(temp));
		else L.push_back(int(temp));
	}
	
	vector<pid_t> pids(Beta.size()*L.size());
	
	for (int i = 0; i < Beta.size()*L.size(); i++) pids[i] = 0; 
	
	int blockDimentions[lenght] = {0};
	for(int i = 0; i < lenght; i++) blockDimentions[i] = start+i*step;
	
	int boundary = 1;
	int measure[N-(N/100)*frac] = {0}; 
	
	if(parall) cout << "Per interrompere la simulazione digitare 1 (più invio). Al termine della simulazione per terminare il programma digitare 0 (più invio)." << endl;
			
	for(int i = 0; i < Beta.size(); i++) {
	
		long double beta = Beta[i];
	
		for(int j = 0; j < L.size(); j++) {
		
			if(parall) {
				pids[i*Beta.size() + j] = fork();
				
				if(pids[i*Beta.size() + j] == -1) {
					cerr << "fork() error" << endl;
					return -1;
				}
			
				if(pids[i*Beta.size() + j] != 0)  {
					//Aspetto 10 ms da un processo all'altro.
					usleep(10000);
					continue;
				}
			}
			
			int l1 = L[j];
			int l2 = L[j];
			int t = 3*L[j];
			
			char nameProcess[20];
			sprintf(nameProcess, "sim%Lf:%i", beta, l1);
			if(prctl(PR_SET_NAME, nameProcess , NULL, NULL, NULL) != 0) {
				cerr << "Error in renaming process: sim" << beta << ":" << l1 << endl;
				return -1;
			}	
			
			fstream file;
			file.open("data.txt");
			
			//Una di queste distribuzioni è inutile
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
			hotStart(latt, l1, l2, t);
				
			int per = 0; 
			int aper = 0;
			
			for(int rep = 0; rep < N; rep++) {
				
				//if(rep % 10000 == 0) cout << "beta = " << beta << " l = " << l1 << " rep = " << rep << endl; 
										
				if(rep%2 == 0) singleCluster(latt, cluster, boundary, l1, l2, t, beta, &distr_l1, &distr_l2, &distr_t, &distribution);
				else boundary = surfaceCluster(latt, cluster, boundary, l1, l2, t, beta, &distr_l1, &distr_l2, &distribution);

				if(rep >= (N/100)*frac) {
					
					if(boundary == 1) per += 1;
					else aper += 1;

					measure[rep-(N/100)*frac] = boundary;
				}
				
				
				if((rep > 1000) && (flag == false) && (boundary == -1)) flag = true;
				if(flag == true && count < 10) {
					for(int i = 0; i < l1; i++) 
						for(int j = 0; j < l2; j++) 
							for(int k = 0; k < t; k++) 
								file << latt[i][j][k] << endl;
					count++;
				}
			}
			
			double Aper = double(aper);
			double Per = double(per);
			double td = double(t);
					
			double F_mod = log(td) - log( 0.5 * log((1.0+ Aper/Per)/(1.0 - Aper/Per)));	
			
			//Jacknife
			double error = 0.0;
			error = jackknife(measure, N-(N/100)*frac, blockDimentions, lenght, td);
					
			//int fine = (clock() - init)/CLOCKS_PER_SEC;
					
			/*cout << "********************************************************************************************" << endl;
			cout << "l1 = " << l1 << "\t" "l2 = " << l2 << "\t" << "t = " << t << "\t"  << "beta = " << beta << "\t" << "tempo = " << fine << " s"  << endl;
			cout  << "Aper/Tot = " << Aper/(N) << "\t" << "Aper/Per = " << Aper/Per << "\t" << "F_s = " << -log(Aper) + log(Per) + log(td) << "\t" << "F_si = " << F_mod << "+/-" << error << endl;*/
			
+
			cout << beta << "\t" << l1 << "\t" << F_mod << "\t" << error << endl << flush;
				
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
			
			if(parall) return 0;
		}
	}
		
	//Ciclo per la rimozione dei sottoprocessi
	if(parall) {
		int killVar = -1;
		do {
			cin >> killVar;
			if(killVar == 1) {
				pid_t pid = 0;
					for(int i = 0; i < Beta.size()*L.size(); i++) {
					pid = pids[i];
					kill(pid, SIGTERM);
					bool died = false;
					for (int loop; !died && loop < 5; ++loop)
					{
    						int status;
    						pid_t id;
    						sleep(1);
    						if (waitpid(pid, &status, WNOHANG) == pid) died = true;
					}
		
					if (!died) kill(pid, SIGKILL);
				}
			
			}
		} while (killVar != 0 && killVar != 1);
	}
	
	return 0;
}


double jackknife(int* boundary, int size, int* blockDimentions, int len, double t) {

	int avgP = 0;
	int avgAP = 0;
	double davgP = 0.0;
	double davgAP = 0.0;
	
	double media = 0.0;
	double sum = 0.0;
	int Nb = 0;
	double F_mod = 0.0;

	double* block = 0;
	double* vars = (double*) malloc (len * sizeof(double));
	
	int i = 0;
	int j = 0;
	int k = 0;
	int count = 0;
	
	for(i = 0; i < len; i++) {
	
		Nb = size/blockDimentions[i];
		block = (double*) malloc (Nb * sizeof(double));
		
		for(j = 0; j < Nb; j++) {
		
			avgP = 0;
			avgAP = 0;
			
			for(k = 0; k < Nb; k++) {
			
				if(j == k) continue;
				
				for(count = 0; count < blockDimentions[i]; count++) {
					if( boundary[k*blockDimentions[i]+count] == 1 ) avgP += 1;
					else avgAP += 1;
				}
				
			}
			//cout << "Periodic average =" << avgP << endl;
			//cout << "AntiPeriodic average = " << avgAP << endl;
			davgP = double(avgP);
			davgAP = double(avgAP);
			F_mod = t - log( 0.5 * log( (1.0 + davgAP/davgP) / (1.0 - davgAP/davgP) ));
			
			//cout << "F_mod = " << F_mod << endl;
			
			block[j] = F_mod;
		}
		
		//Calcolo la varianza
		
		media = 0.0;
		for(j = 0; j < Nb; j++) media += block[j];	
		media /= Nb;
		
		//cout << "Media = " << media;
		
		sum = 0.0;
		for(j = 0; j < Nb; j++) sum += (block[j] - media)*(block[j] - media);
		sum /= (double(Nb) - 1.0)/double(Nb);
		
		vars[i] = sum;

		free(block);
	}
	
	//Eseguo la media di tutte le varianze in vars
	
	media = 0.0;
	for(int i = 0; i < len; i++) {
		media += vars[i];
		//cout << vars[i] << endl;
	}
	media /= double(len);
	
	free(vars);
	
	return sqrt(media);
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

int clusterize(int*** latt, int*** cluster, int boundary, int l1, int l2, int t, long double beta, int i, int j, int k, uniform_real_distribution<long double>*  distribution) {

	stack<site> stack;
	
	int size[3] = {l1, l2, t};
	
	site next = {{i, j, k}, -1, -1};
	stack.push(next);
	
	int flag = 0, value = 0, d = 0, a = 0, s = 0;
	double p = 0.0;

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
					else if(next.c[d] == -1) { next.c[d] = size[d] - 1; if(d == 2) next.cross = 1;}
					
					//Conviene o no togliere questo if... secondo me no.
					if(cluster[next.c[0]][next.c[1]][next.c[2]] == 0) {
						if( (((next.cross == 1) ? boundary : 1) * latt[current.c[0]][current.c[1]][current.c[2]]*latt[next.c[0]][next.c[1]][next.c[2]]) == 1) {
							p = (1.0 - exp(-2 * beta));
							if((*distribution)(mt) < p) stack.push(next);
						}
					}
					
				}
			}
		}
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

	clusterize(latt, cluster, boundary, l1, l2, t, beta, i, j, k, distribution);	
	
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
				if(clusterize(latt, cluster, boundary, l1, l2, t, beta, i, j, t-1, distribution) == 1) flag = + 1;
				
					
	if(flag == 0) {
		for(int i = 0; i < l1; i++) 
			for(int j = 0; j < l2; j++) 
				for(int k = 0; k < t; k++) 
					latt[i][j][k] = latt[i][j][k] * ((cluster[i][j][k] ==  0) ? +1 : cluster[i][j][k]);
		boundary = boundary * (-1);
	}

	return boundary;
}
