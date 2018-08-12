int clusterize(int*** latt, int*** cluster, int boundary, unsigned int l1, unsigned int l2, unsigned int t, long double beta, int precendente, int cross, unsigned int i, unsigned int j, unsigned int k, uniform_real_distribution<long double>*  distribution) {

	unsigned int size[3] = {l1, l2, t};
	unsigned int coord[3] = {i, j, k};
	int ret[6] = {0, 0, 0, 0, 0, 0};
	long double p = 0.0;
	int flag = 0;
	int ret = 0;
	int bound = +1;
	int newCross = -1;
	int b = 0;

	
	if(cluster[i][j][k] == (cross * precendente)) return 1;
	
	else if(cluster[i][j][k] == 0) {
	
		//Se non c'è cross cross vale -1, se c'è vale +1
		cluster[i][j][k] = precedente * (-cross);
		
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
				
				if(cluster[coord[0]][coord[1]][coord[2]] == 0) {
				
					p = 1.0 - exp(-beta * (1.0 + bound*latt[i][j][k]*latt[coord[0]][coord[1]][coord[2]]));
			
					if ((*distribution)(mt) < p) {
						ret[a*(d+1)+3] = clusterize(latt, cluster, boundary, l1, l2, t, beta, cluster[i][j][k], newCross, coord[0], coord[1], coord[2], distribution);
					
					}

				}

			}
		}
		
		for(int i = 0; i<6; i++) flag = flag || ret[i];
	}
	
	return flag;

}
