#include <math.h>
#include <stdio.h>
#include "ip_lib.h"

int main () {
    ip_mat *mat1;
    mat1 = ip_mat_create(5, 6, 3, 1);
    /* test set*/
    set_val(mat1, 0,1,1,2);
    set_val(mat1, 0,1,2,4);
    set_val(mat1, 0,0,0, 3);
    set_val(mat1, 3,4,2,7);
    print_ipmat(mat1);
    /*test get*/
    printf("\nget val(3,4,2): %f\n", get_val(mat1, 3,4,2));
		
		ip_mat_init_random(mat1, 1, 1);

		print_ipmat(mat1);

		ip_mat_free(mat1);
    return 0;
}

