int clusterize(int*** latt, int*** cluster, int boundary, unsigned int l1, unsigned int l2, unsigned int t, long double beta, int precedente, int cross, unsigned int i, unsigned int j, unsigned int k, uniform_real_distribution<long double>*  distribution) {

	stack<site> stack;
	unsigned int size[3] = {l1, l2, t};
	
	site init = {{i, j, k}, precedente, cross};
	site next = {0};
	
	stack.push(init);
	
	int flag = 0;
	
	int newCross = -1;
	long double p = 0.0;
	
	//cout << init.i << " " << init.j << " " << init.k << endl;
	
	int d = 0;
	int a = 0;
	int b = 0;
	
	while(!stack.empty()) {
	
		site current = stack.top();
		stack.pop();
		//cout << current.i << " " << current.j << " " << current.k << endl;
		
		//unsigned int coord[3] = {current.i, current.j, current.k};
	
		
		if(cluster[current.coord[0]][current.coord[1]][current.coord[2]] == (current.cross * current.precedente)) flag = 1;
		
		else if(precedente == 0) {
			cluster[current.coord[0]][current.coord[1]][current.coord[2]] = current.precedente * (-current.cross);
			next.precedente = cluster[current.coord[0]][current.coord[1]][current.coord[2]];
			for(d = 0; d < 3; d++) {
				for(a = -1; a < 2; a = a + 2) {
				
					next.cross = -1;
					//coord[0] = current.i; coord[1] = current.j; coord[2] = current.k;
					for(int i = 0; d < 3; d++) next.coord[i] = current.coord[i];
					
					/*coord[d] = coord[d] = a;
					if((coord[d] + a == size[d]) || (coord[d] + a == -1)) { 
						coord[d] = size[d] - 1 - coord[d];
						if(d == 2) { bound = boundary; newCross = 1; }
					}*/
					
					b = next.coord[d];
					if((next.coord[d] + a == size[d]) || (next.coord[d] + a == -1)) { 
						b = size[d] - 1 - (next.coord[d] + a);
						if(d == 2) next.cross = 1;
					}
					next.coord[d] = a+b;
					
					if(cluster[next.coord[0]][next.coord[1]][next.coord[2]] == 0) {
						p = 1.0 - exp(-beta * (1.0 + ((next.cross == 1) ? boundary : 1) * precedente * latt[next.coord[0]][next.coord[1]][next.coord[2]]));
						if ((*distribution)(mt) < p) {
							//next.i = coord[0]; next.j = coord[1]; next.k = coord[2]; 
							stack.push(next);
						}
					}
				}
			}
		}
	}
	return flag;
}

