#include <math.h>
#include <stdio.h>
#include "ip_lib.h"

int main () {
    /* VARIABILI RICICLABILI NEI TEST*/
    /* indici e acc da poter usare per cicli vari nei test*/
    int ih, iw, ik;
    int acc;
    int val; /* valore*/
    /*3 dimensioni da poter usare nei test */
    int h, w, k;

    /* test ip_mat_create*/
    printf("\n TEST ip_mat_create \n");
    ip_mat *mat1;
    h=5;
    w=6;
    k=3;
    val = 2;
    printf("\n Matrice A %d x %d x %d : \n", h, w, k);
    mat1 = ip_mat_create(h, w, k, val);

    /* test set_val*/
    printf("\n TEST set_val \n");
    set_val(mat1, 0,1,1,2);
    set_val(mat1, 0,1,2,4);
    set_val(mat1, 0,0,0,3);
    set_val(mat1, 3,4,2,7);
    printf("\n Matrice A modificata con set: \n");
    ip_mat_print(mat1);

    /*test get_val*/
    printf("\n TEST get_val \n");
    printf("\n get val(3,4,2): %f\n", get_val(mat1, 3,4,2));

    /*test init_random*/
    printf("\n TEST ip_mat_init_random \n");
    ip_mat_init_random(mat1, 1, 1);

    printf("\n");
    ip_mat_print(mat1);

    /*test ip_mat_free*/
    printf("\n TEST ip_mat_free \n");
    ip_mat_free(mat1);

    /* test ip_mat_subset*/
    printf("\n TEST ip_mat_subset \n");
    ip_mat *mat2, *subset;    

    h = 5;
    w = 6;
    k = 3;
    acc = 0;
    mat2 = ip_mat_create(h, w, k, 0);
    for (ih = 0; ih<h; ih++)
        for(iw = 0; iw< w; iw++)
            for(ik=0; ik<k; ik++, acc++)
                set_val(mat2,ih,iw,ik,acc);

    printf("\n Matrice A: \n");
    ip_mat_print(mat2);

    subset = ip_mat_subset(mat2, 2,4,1,4);
    printf("\n Matrice subset di A: \n");
    ip_mat_print(subset);

    ip_mat_free(mat2);
    ip_mat_free(subset);

    /*test ip_mat_sum*/
    printf("\n TEST ip_mat_sum: \n");
    ip_mat *sum, *add1, *add2;
    
    h = 4;
    w = 2;
    k = 3;

    add1 = ip_mat_create(h,w,k,2);
    add2 = ip_mat_create(h,w,k,5);
    sum = ip_mat_sum(add1, add2);
    
    printf("\n A: \n");
    ip_mat_print(add1);

    printf("\n B: \n");
    ip_mat_print(add2);

    printf("\n A + B: \n");
    ip_mat_print(sum);
    printf("\n");
    ip_mat_free(sum);

    /*test ip_mat_sub*/
    printf("\n TEST ip_mat_sub: \n");
    ip_mat *sub;
    
    sub = ip_mat_sub(add1, add2);

    
    printf("\n A - B: \n");
    ip_mat_print(sub);
    printf("\n");
   
    ip_mat_free(sub);
    ip_mat_free(add1);
    ip_mat_free(add2);

    /*test ip_mat_mean*/
    printf("\n TEST ip_mat_mean: \n");
    ip_mat* a = ip_mat_create(2,2,3,6.0f);
    ip_mat* b = ip_mat_create(2,2,3,3.0f);

    ip_mat *c = ip_mat_mean(a, b);

    printf("\n Matrice A:\n");
    ip_mat_print(a);
    printf("\n Matrice B:\n");
    ip_mat_print(b);
    printf("\n Media di A e B\n");
    ip_mat_print(c);
    
    ip_mat_free(a);
    ip_mat_free(b);
    ip_mat_free(c);

    /*test ip_mat_mul_scalar*/
    printf("\n TEST ip_mat_mul_scalar: \n");
    ip_mat *d = ip_mat_create(2, 2, 3, 3.0f);
    ip_mat *e = ip_mat_mul_scalar(d, 3.0f);

    printf("\n Matrice A: \n");
    ip_mat_print(d);
    printf("\n A x 3:\n");
    ip_mat_print(e);

    ip_mat_free(d);
    ip_mat_free(e);

    return 0;
}

