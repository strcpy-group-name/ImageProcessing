#include <math.h>
#include <stdio.h>
#include "ip_lib.h"

int main () {
    ip_mat *mat1;
    mat1 = ip_mat_create(5, 6, 3, 0);
    set_val(mat1, 0,1,0, 2);
    print_ipmat(mat1);
    printf("\n  get val: %f", get_val(mat1, 0,2,0));
    return 0;
}

