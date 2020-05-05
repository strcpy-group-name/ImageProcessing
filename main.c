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

    printf("\n");
    print_ipmat(mat1);
    /*test get*/
    printf("\nget val(3,4,2): %f\n", get_val(mat1, 3,4,2));

    /*test init_random*/
    ip_mat_init_random(mat1, 1, 1);

    printf("\n");
    print_ipmat(mat1);

    /*test ip_mat_free*/
    ip_mat_free(mat1);

    /* test subset*/
     printf("test subset\n");
    ip_mat *mat2, *subset2; 
    int rows, cols, channels;
    int i, j, k;
    int acc;

    rows = 5;
    cols = 6;
    channels = 3;
    acc = 0;
    mat2 = ip_mat_create(rows, cols, channels, 0);
    for (i = 0; i<rows; i++)
        for(j = 0; j< cols; j++)
            for(k=0; k<channels; k++, acc++)
                set_val(mat2,i,j,k,acc);
    printf("\n");
    print_ipmat(mat2);
    subset2 = ip_mat_subset(mat2, 2,4,1,4);
    printf("\n");
    print_ipmat(subset2);

    /*test ip_mat_mean*/
    ip_mat* a = ip_mat_create(2,2,3,6);
    ip_mat* b = ip_mat_create(2,2,3,3);

    ip_mat* c = ip_mat_mean(a, b);

    printf("\n");
    print_ipmat(c);

    ip_mat_free(a);
    ip_mat_free(b);
    ip_mat_free(c);

    /*test ip_mat_mul_scalar*/
    ip_mat* d = ip_mat_create(2,2,3,10);
    a = ip_mat_mul_scalar(d, 1.2f);

    printf("\n");
    print_ipmat(a);

    ip_mat_free(a);
    ip_mat_free(d);

    return 0;
}

