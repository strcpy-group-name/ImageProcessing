#include <math.h>
#include <stdio.h>
#include "ip_lib.h"

int main () {
    int ih, iw, ik;
    int acc;
    int val; 
    int h, w, k;

    printf("\n TEST ip_mat_create \n");
    ip_mat *mat1;
    h=5;
    w=6;
    k=3;
    val = 2;
    printf("\n Matrice A %d x %d x %d : \n", h, w, k);
    mat1 = ip_mat_create(h, w, k, val);
    ip_mat_show(mat1);

    printf("\n TEST set_val \n");
    set_val(mat1, 0,1,1,2);
    set_val(mat1, 0,1,2,4);
    set_val(mat1, 0,0,0,3);
    set_val(mat1, 3,4,2,7);
    printf("\n Matrice A modificata con set: \n");
    ip_mat_show(mat1);

    printf("\n TEST get_val \n");
    printf("\n get val(3,4,2): %f\n", get_val(mat1, 3,4,2));

    printf("\n TEST ip_mat_init_random \n");
    ip_mat_init_random(mat1, 1, 1);

    printf("\n");
    ip_mat_show(mat1);
    
    printf("\n TEST ip_mat_free \n");
    ip_mat_free(mat1);

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
    ip_mat_show(mat2);

    subset = ip_mat_subset(mat2, 2,4,1,4);
    printf("\n Matrice subset di A: \n");
    ip_mat_show(subset);

    ip_mat_free(mat2);
    ip_mat_free(subset);

    printf("\n TEST ip_mat_sum: \n");
    ip_mat *sum, *add1, *add2;
    
    h = 4;
    w = 2;
    k = 3;

    add1 = ip_mat_create(h,w,k,2);
    add2 = ip_mat_create(h,w,k,5);
    sum = ip_mat_sum(add1, add2);
    
    printf("\n A: \n");
    ip_mat_show(add1);

    printf("\n B: \n");
    ip_mat_show(add2);

    printf("\n A + B: \n");
    ip_mat_show(sum);
    printf("\n");
    ip_mat_free(sum);

    printf("\n TEST ip_mat_sub: \n");
    ip_mat *sub;
    
    sub = ip_mat_sub(add1, add2);

    
    printf("\n A - B: \n");
    ip_mat_show(sub);
    printf("\n");
   
    ip_mat_free(sub);
    ip_mat_free(add1);
    ip_mat_free(add2);

    printf("\n TEST ip_mat_mean: \n");
    ip_mat* a = ip_mat_create(2,2,3,6.0f);
    ip_mat* b = ip_mat_create(2,2,3,3.0f);

    ip_mat *c = ip_mat_mean(a, b);

    printf("\n Matrice A:\n");
    ip_mat_show(a);
    printf("\n Matrice B:\n");
    ip_mat_show(b);
    printf("\n Media di A e B\n");
    ip_mat_show(c);
    
    ip_mat_free(a);
    ip_mat_free(b);
    ip_mat_free(c);

    printf("\n TEST ip_mat_mul_scalar: \n");
    ip_mat *d = ip_mat_create(2, 2, 3, 3.0f);
    ip_mat *e = ip_mat_mul_scalar(d, 3.0f);

    printf("\n Matrice A: \n");
    ip_mat_show(d);
    printf("\n A x 3:\n");
    ip_mat_show(e);

    ip_mat_free(d);
    ip_mat_free(e);

    printf("\n TEST ip_mat_padding: \n");
    ip_mat *no_pad = ip_mat_create(3,3,3,1.0f);
    ip_mat *pad = ip_mat_padding(no_pad, 2, 2);

    ip_mat_show(no_pad);
    ip_mat_show(pad);

    ip_mat_free(no_pad);
    ip_mat_free(pad);

    printf("\n TEST ip_mat_sum, ip_mat_show_stats, ip_mat_compute_stats: \n");
    ip_mat *aa = ip_mat_create(2, 2, 3, 3.0f);
    ip_mat *bb = ip_mat_create(2, 2, 3, 6.0f);
    set_val(aa, 1, 1, 1, 300.0f);
    ip_mat *cc = ip_mat_sum(aa, bb);
    printf("\n\nSOMMA:\n");
    ip_mat_show(cc);
    ip_mat_show_stats(cc);
    ip_mat_free(aa);
    ip_mat_free(bb);
    ip_mat_free(cc);

    printf("\n TEST ip_mat_copy: \n");
    ip_mat *a01 = ip_mat_create(3, 3, 3, 2.445f);
    set_val(a01, 1, 1, 1, 300.0f);
    ip_mat *a02 = ip_mat_copy(a01);
    printf("\nmatrice a1:");
    ip_mat_show(a01);
    printf("\nmatrice a2:");
    ip_mat_show(a02);
    ip_mat_free(a01);
    ip_mat_free(a02);

    printf("\n TEST ip_mat_concat: \n");
    ip_mat *a03 = ip_mat_create(3, 3, 3, 1.1f);
    ip_mat *a04 = ip_mat_create(3, 3, 3, 2.2f);
    set_val(a04, 1, 1, 1, 4.0f);
    ip_mat *a05 = ip_mat_concat(a03, a04, 0);
    printf("\nmatrice a05\n");
    ip_mat_show(a05);
    ip_mat *a06 = ip_mat_concat(a03, a04, 1);
    printf("\nmatrice a06\n");
    ip_mat_show(a06);
    ip_mat *a07 = ip_mat_concat(a03, a04, 2);
    printf("\nmatrice a07\n");
    ip_mat_show(a07);
    ip_mat_free(a03);
    ip_mat_free(a04);
    ip_mat_free(a05);
    ip_mat_free(a06);
    ip_mat_free(a07);

    printf("\n TEST ip_mat_brighten image: \n");
    Bitmap *image = bm_load("flower.bmp");
    ip_mat *img = bitmap_to_ip_mat(image);
    ip_mat *bright = ip_mat_brighten(img, 100.0f);
    Bitmap *bbmp = ip_mat_to_bitmap(bright);
    bm_save(bbmp, "flower_b.bmp");
    bm_free(image);
    bm_free(bbmp);
    ip_mat_free(bright);
    ip_mat_free(img);


    printf("\n TEST ip_mat_brighten darken image: \n");
    image = bm_load("flower.bmp");
    img = bitmap_to_ip_mat(image);
    bright = ip_mat_brighten(img, -100.0f);
    bbmp = ip_mat_to_bitmap(bright);
    bm_save(bbmp, "flower_d.bmp");
    bm_free(image);
    bm_free(bbmp);
    ip_mat_free(bright);
    ip_mat_free(img);

    printf("\n TEST ip_mat_blend image: \n");
    Bitmap *im1 = bm_load("mongolfiere.bmp");
    Bitmap *im2 = bm_load("flower2.bmp");
    ip_mat *b1 = bitmap_to_ip_mat(im1);
    ip_mat *b2 = bitmap_to_ip_mat(im2); 
    ip_mat *blend = ip_mat_blend(b1, b2, 0.5f);
    Bitmap *b3 = ip_mat_to_bitmap(blend);
    bm_save(b3, "blend.bmp");
    ip_mat_free(b1);
    ip_mat_free(b2);
    ip_mat_free(blend);
    bm_free(im1);
    bm_free(im2);
    bm_free(b3);
    

    printf("\n TEST ip_mat_padding image: \n");
    Bitmap *caf = bm_load("caf.bmp");
    ip_mat *caf_nopad = bitmap_to_ip_mat(caf);
    ip_mat *caf_pad = ip_mat_padding(caf_nopad, 100.0f, 100.0f);
    bbmp = ip_mat_to_bitmap(caf_pad);
    bm_save(bbmp, "caf_padding.bmp");
    bm_free(caf);
    bm_free(bbmp);
    ip_mat_free(caf_nopad);
    ip_mat_free(caf_pad);
    
    printf("\n TEST convoluzione average image: \n");
    caf = bm_load("caf.bmp");
    ip_mat *caf_noblur = bitmap_to_ip_mat(caf);
    ip_mat *ker_avg = create_average_filter(3,3,3);
    ip_mat *caf_blur = ip_mat_convolve(caf_noblur, ker_avg);
    bbmp = ip_mat_to_bitmap(caf_blur);
    bm_save(bbmp, "caf_average.bmp");
    bm_free(caf);
    bm_free(bbmp);
    ip_mat_free(caf_noblur);
    ip_mat_free(caf_blur);
    ip_mat_free(ker_avg);
    
    printf("\nTest ip_mat_to_gray image: \n");
    Bitmap *g_original_bmp = bm_load("flower.bmp");
    ip_mat *g_original = bitmap_to_ip_mat(g_original_bmp);
    ip_mat *gray = ip_mat_to_gray_scale(g_original);
    Bitmap *gray_bmp = ip_mat_to_bitmap(gray);
    bm_save(gray_bmp, "flower_gray.bmp");
    bm_free(g_original_bmp);
    bm_free(gray_bmp);
    ip_mat_free(gray);
    ip_mat_free(g_original);

    printf("\nTest ip_mat_to_gray_lum image: \n");
    g_original_bmp = bm_load("flower.bmp");
    g_original = bitmap_to_ip_mat(g_original_bmp);
    gray = ip_mat_to_gray_scale_lum_corr(g_original);
    gray_bmp = ip_mat_to_bitmap(gray);
    bm_save(gray_bmp, "flower_gray_lum.bmp");
    bm_free(g_original_bmp);
    bm_free(gray_bmp);
    ip_mat_free(gray);
    ip_mat_free(g_original);
    
    printf("\nTest ip_mat_to_gray_gamma image: \n");
    g_original_bmp = bm_load("flower.bmp");
    g_original = bitmap_to_ip_mat(g_original_bmp);
    gray = ip_mat_to_gray_scale_gamma_corr(g_original);
    gray_bmp = ip_mat_to_bitmap(gray);
    bm_save(gray_bmp, "flower_gray_gamma.bmp");
    bm_free(g_original_bmp);
    bm_free(gray_bmp);
    ip_mat_free(gray);
    ip_mat_free(g_original);

    return 0;
}

