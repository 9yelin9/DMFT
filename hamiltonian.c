// hamiltonian matrix 생성하기
// mu 추가하기
// i, j, h 만 저장하는 행렬 (0이 아닌 것만 저장)

#include "states.h"

#define t 1

int main() {
	int i, j, k, l, m, mu, **hamiltonian, **state, *h_state;
	double **H;

	// 저장 공간 할당
	hamiltonian = (int**)calloc(sizeof(int*), particle_num);
	for(i = 0; i < particle_num; i++) {
		hamiltonian[i] = (int*)calloc(sizeof(int), particle_num);
	}

	state = (int**)malloc(sizeof(int*) * state_num);
	for(i = 0; i < state_num; i++) {
		state[i] = (int*)malloc(sizeof(int) * particle_num);
	}

	H = (double**)calloc(sizeof(double*), state_num);
	for(i = 0; i < state_num; i++) {
		H[i] = (double*)calloc(sizeof(double), state_num);
	}

	// hamiltonian 정보 입력
	for(i = 0; i < particle_num; i++) {
		hamiltonian[i][i] = -1; // annihilation operator

		j = (i + 2) % particle_num;
		hamiltonian[i][j] = +1; // creation operator
	}

	generate_state(state); // N, S_Z를 기준으로 state list 생성하기

	// hamiltonian matrix 생성
	for(i = 0; i < state_num; i++) {
		for(j = 0; j < particle_num; j++) {

			for(l = 0; l < state_num; L
		
			h_state = (int*)malloc(sizeof(int) * particle_num);

			for(k = 0; k < particle_num; k++) {
				h_state[k] = state[i][k] + hamiltonian[j][k];
			}

			for(l = 0; l < state_num; l++) {
				for(k = 0; k < particle_num; k++) {
					if(h_state[k] != state[l][k]) break;
				}

				if(k == particle_num){
					for(m = 0; m < particle_num; m++) {
						if(hamiltonian[j][m] != 0) {
							if(state[i][m + 1] == 1) {
								H[i][l] = t; // 부호 바뀜
								break;
							}
							else {
								H[i][l] = -t; // 부호 안바뀜
								break;
							}
						}
					}
				}
			}
			free(h_state);
		}
	}	

	// hamiltonian matrix 출력

	printf("- hamiltonian matrix\n");

	for(i = 0; i < state_num + 1; i++) {
		printf("%-4d", i);
	}
	printf("\n");

	for(i = 0; i < state_num; i++) {
		
		printf("%-4d", i + 1);

		for(j = 0; j < state_num; j++) {
			printf("%-4.0f", H[i][j]);
		}
		printf("\n");
	}
	
	// 동적 할당 해제
	for(i = 0; i < particle_num; i++) {
		free(hamiltonian[i]);
	}
	free(hamiltonian);

	for(i = 0; i < state_num; i++) {
		free(state[i]);
		free(H[i]);
	}
	free(state);
	free(H);

	return 0;
}
