#ifndef __STATES_H__
#define __STATES_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define particle_num (2 * 2) // spin up, down
#define state_num pow(2, particle_num)

struct linked_list {
        int index;
        int p_sum; // particle_sum
        int s_sum; // spin_sum
        struct linked_list *next;
};

typedef struct linked_list NODE;
typedef NODE *LINK;

void generate_state(int **state);

#endif
