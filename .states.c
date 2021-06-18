// N, S_Z를 기준으로 state list 생성하기
// 2차원 배열 말고 1차원 배열로 함 해보기

#include "states.h"

void generate_state(int **state) {
	LINK head = NULL, cur = NULL, tail = NULL;
	int i, j, *spin, **bin, tmp, p_sum, s_sum;

	// 저장 공간 동적 할당
	spin = (int*)malloc(sizeof(int) * particle_num);

	bin = (int**)malloc(sizeof(int*) * state_num);
	for(i = 0; i < state_num; i++) {
		bin[i] = (int*)malloc(sizeof(int) * particle_num);
	}

	// spin 정보 입력(짝수 index : down, 홀수 index : up)
	for(i = 0; i < particle_num; i++) {
		if(i % 2 == 0) {
			spin[i] = -1;
		}
		else {
			spin[i] = 1;
		}
	}

	// binary list 생성
	for(i = 0; i < state_num; i++) {
		
		tmp = i;

		for(j = 0; j < particle_num; j++) {
			bin[i][j] = tmp % 2;
			tmp /= 2;
		}
	}

	// head node 설정
	head = (LINK)malloc(sizeof(NODE));
	head->index = 0;
	head->p_sum = 0;
	head->s_sum = 0;
	head->next = NULL;

	// N, S_Z를 기준으로 binary list index 정렬
	for(i = 1; i < state_num; i++) {

		p_sum = 0;
		s_sum = 0;

		for(j = 0; j < particle_num; j++) {
			p_sum += bin[i][j];
			s_sum += bin[i][j] * spin[j];
		}
		
		cur = (LINK)malloc(sizeof(NODE));
		cur->index = i;
		cur->p_sum = p_sum;
		cur->s_sum = s_sum;
		cur->next = NULL;

		tail = head;

		while(tail != NULL) {
			if(tail->next == NULL) {
				tail->next = cur;
				tail = cur;
				break;
			}
			else {
				if(cur->p_sum < tail->next->p_sum) {
					cur->next = tail->next;
					tail->next = cur;
					break;
				}
				else if(cur->p_sum == tail->next->p_sum && cur->s_sum < tail->next->s_sum) {
					cur->next = tail->next;
					tail->next = cur;
					break;
				}
				tail = tail->next;
			}
		}
	}

	// 정렬 결과 대입 후 state list 생성 및 출력
	i = 0;
	tail = head;

	printf("- state list\n");
	printf("%-10s%-10s%-10s%-10s\n", "index", "N", "S_Z", "state");

	while(tail != NULL) {

		printf("%-10d%-10d%-10d", i + 1, tail->p_sum, tail->s_sum);
		for(j = 0; j < particle_num; j++) {
			state[i][j] = bin[tail->index][j];
			printf("%d", state[i][j]);
		}
		printf("\n");
	
		i++;
		tail = tail->next;
	}
	printf("\n");

	// 동적 할당 해제
	free(spin);

	for(i = 0; i < state_num; i++) {
		free(bin[i]);
	}
	free(bin);
}
