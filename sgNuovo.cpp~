/*Montecarlo per generare le superfici di separazione*/

#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <sstream>
#include <stack>
#include <unistd.h>
#include <signal.h>
#include <sys/wait.h>
#include <cstring>
#include <cmath>
#include <sys/prctl.h>
#include <iomanip>
#include <limits>

using namespace std;

//Inserire il numero di iterazioni
#define N 1000000
//Percentuale iniziale rimossa dalla catena di Markov perché si considera termalizzazione.
#define frac 10
//Variabili che controllano il jackknife
#define lenght 40

typedef struct {
	int c[3];
	int precedente;
	int cross;
} site;

typedef struct {
	double value;
	double error;
} F_mod;

int clusterize(int*** latt, int*** cluster, int boundary, int l1, int l2, int t, long double beta, int i, int j, int k, uniform_real_distribution<long double>*  distribution);
int surfaceCluster(int*** latt,  int*** cluster, int boundary, int l1, int l2, int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_real_distribution<long double>*  distribution);
void singleCluster(int*** latt,  int*** cluster, int boundary, int l1, int l2, int t, long double beta, uniform_int_distribution<int>* distr_l1, uniform_int_distribution<int>* distr_l2, uniform_int_distribution<int>* distr_t, uniform_real_distribution<long double>*  distribution);
void hotStart(int*** latt, int l1, int l2 , int t);
void coldStart(int*** latt, int l1, int l2, int t, int spin);
F_mod jackknife(int* boundary, int size, int* blockDimentions, int len, int t);


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
	
	//int blockDimentions[lenght] = {30, 32, 36, 40, 45, 48, 50, 60, 72, 75, 80, 90, 96, 100, 120, 125, 144, 150, 160, 180, 200, 225, 240, 250, 288, 300, 360, 375, 400, 450, 480, 500, 600, 625, 720, 750, 800, 900, 1000, 1125, 1200, 1250, 1440, 1500, 1800, 1875, 2000, 2250, 2400, 2500, 3000, 3125, 3600, 3750, 4000, 4500, 5000, 5625, 6000, 6250, 7200, 7500, 9000, 9375, 10000};
	
	int blockDimentions[lenght] = {300, 360, 375, 400, 450, 480, 500, 600, 625, 720, 750, 800, 900, 1000, 1125, 1200, 1250, 1440, 1500, 1800, 1875, 2000, 2250, 2400, 2500, 3000, 3125, 3600, 3750, 4000, 4500, 5000, 5625, 6000, 6250, 7200, 7500, 9000, 9375, 10000};
	
	int boundary = 1;
	int measure[N-(N/100)*frac] = {0}; 
	
	F_mod result = {0.0, 0.0};
	
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
			
			for(int rep = 0; rep < N; rep++) {
				
				//if(rep % 10000 == 0) cout << "beta = " << beta << " l = " << l1 << " rep = " << rep << endl; 
										
				if(rep % 2 == 0)singleCluster(latt, cluster, boundary, l1, l2, t, beta, &distr_l1, &distr_l2, &distr_t, &distribution);
				else boundary = surfaceCluster(latt, cluster, boundary, l1, l2, t, beta, &distr_l1, &distr_l2, &distribution);

				if(rep >= (N/100)*frac)
					measure[rep-(N/100)*frac] = boundary;
			}
			
					
			result = jackknife(measure, N-(N/100)*frac, blockDimentions, lenght, t);
					
			//int fine = (clock() - init)/CLOCKS_PER_SEC;
					
			/*cout << "********************************************************************************************" << endl;
			cout << "l1 = " << l1 << "\t" "l2 = " << l2 << "\t" << "t = " << t << "\t"  << "beta = " << beta << "\t" << " s"  << endl;
			cout  << "Aper/Tot = " << Aper/(N) << "\t" << "Aper/Per = " << Aper/Per << "\t" << "F_s = " << -log(Aper) + log(Per) + log(td) << "\t" << "F_si = " << F_mod << "+/-" << error << endl;*/
			
			cout << beta << "\t" << l1 << "\t" << result.value << setprecision(10) << "\t" << result.error << setprecision(10) << endl << flush;
			
			//per ogni condizione data scrivo nuovo file il cui nome sarà nel formato beta_l1_l2_t
			ofstream file;
			ostringstream beta_str, l1_str;
					
			beta_str << beta;
			l1_str << l1;

			file.open ("data" + beta_str.str() + ":" + l1_str.str() + ".txt");
			
			file << beta << "\t" << l1 << "\t" << result.value << setprecision(10) << "\t" << result.error << setprecision(10) << endl << flush;
				
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


F_mod jackknife(int* boundary, int size, int* blockDimentions, int len, int t) {

	int avgP = 0, avgAP = 0, Nb = 0;
	
	double media = 0.0, sum = 0.0, F_biased = 0.0;
	
	double* F_unbiased = (double*) malloc (len * sizeof(double));

	double* block = 0;
	double* vars = (double*) malloc (len * sizeof(double));
	
	F_mod finalValue = {0.0, 0.0};
	
	int i = 0, j = 0, k = 0, count = 0;
	
	for(i = 0; i < size; i++) {
		if(boundary[i] == 1) avgP += 1;
		else avgAP += 1;
	}
	
	F_biased = log((double) t) - log( 0.5 * log( (1.0 + ( (double) avgAP) / ( (double) avgP )) / (1.0 - ( (double) avgAP) / ( (double) avgP ) ) ) ) ;
	 
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

			block[j] = log((double) t) - log( 0.5 * log( ( 1.0 + ( (double) avgAP) / ( (double) avgP )) / (1.0 - ( (double) avgAP) / ( (double) avgP )) ) );
			//cout << ( (double) avgAP)/ ( (double) avgP )) )  << endl;
			
			//cout << "F_mod = " << F_mod << endl;
			
			
		}
		
		
		
		//Calcolo la varianza
		
		media = 0.0;
		for(j = 0; j < Nb; j++) {
			media += block[j];
			//cout << block[i] << endl;
		}
		
		//FIXME: Qui metterci lo stimatore unbiased
		//F_unbiased[i] = Nb * F_biased - ((( (double) Nb - 1.0) * (double) blockDimentions[i])/( (double) Nb * (double) blockDimentions[i])) * media;
			
		media /= (double) Nb;
		
		F_unbiased[i] = (Nb) * F_biased - ( Nb - 1) * media;  
		
		//cout << "Media = " << media;
		
		sum = 0.0;
		for(j = 0; j < Nb; j++) sum += (block[j] - media)*(block[j] - media);
		sum /= (double(Nb) - 1.0)/double(Nb);
		
		vars[i] = sqrt(sum);
		
		//if((i == 1) || (i == 20)) cout << Nb << "\t" << F_unbiased[i] << "\t" << F_biased << "\t" << media << "\t" << vars[i] << endl;

		free(block);
	}
	
	//Eseguo la media di tutte le varianze in vars
	
	media = 0.0;
	for(int i = 0; i < len; i++) {
		media += vars[i];
		//cout << vars[i] << endl;
	}
	media /= double(len);
	
	finalValue.error = media;
	
	//Eseguo la media di tutti i valori unbiased
	
	
	media = 0.0;
	for(int i = 0; i < len; i++) {
		media += F_unbiased[i];
		//cout << vars[i] << endl;
	}
	media /= double(len);
	
	//Eseguo la media di tutti i valori unbiased
	
	finalValue.value = media;
	
	free(F_unbiased);
	free(vars);
	
	return finalValue;
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
