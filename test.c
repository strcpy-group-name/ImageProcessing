#include <math.h>
#include <stdio.h>
#include "ip_lib.h"

int main()
{
    ip_mat* no_gauss, *ker_gauss, *ker_avg, *gauss, *no_blur, *blur, *b1, *b2
    ,*blend, *concatenazione, *img, *bright, *a08, *a09, *corruzione, *caf_nopad, *caf_pad
    ,*g_original, *gray, *mat1, *mat2, *subset, *sum, *add1, *add2, *sub, *a, *d, *no_pad, *a01,
    *a02, *a03, *a05, *a06, *a07, *a04, *aa, *bb, *cc, *pad, *e, *b, *c;
    Bitmap *im1, *im2, *b3, *bbmp, *image, *caf, *g_original_bmp, *gray_bmp;


    int ih, iw, ik;
    int acc;
    int val;
    int h, w, k;

    printf("\n TEST ip_mat_create \n");
    h=3;
    w=3;
    k=1;
    val = 2;
    printf("\n Matrice A %d x %d x %d : \n", h, w, k);
    mat1 = ip_mat_create(h, w, k, val);
    ip_mat_show(mat1);

    printf("\n TEST ip_mat_init_random \n");
    ip_mat_init_random(mat1, 1, 1);

    printf("\n");
    ip_mat_show(mat1);

    printf("\n TEST ip_mat_free \n");
    ip_mat_free(mat1);

    printf("\n TEST ip_mat_subset \n");

    h = 5;
    w = 6;
    k = 3;
    acc = 0;
    mat2 = ip_mat_create(h, w, k, 0);
    for (ih = 0; ih < h; ih++)
        for (iw = 0; iw < w; iw++)
            for (ik = 0; ik < k; ik++, acc++)
                set_val(mat2, ih, iw, ik, acc);

    printf("\n Matrice A: \n");
    ip_mat_show(mat2);

    subset = ip_mat_subset(mat2, 2, 4, 1, 4);
    printf("\n Matrice subset di A: \n");
    ip_mat_show(subset);

    ip_mat_free(mat2);
    ip_mat_free(subset);

    printf("\n TEST ip_mat_sum: \n");
    h = 4;
    w = 2;
    k = 3;

    add1 = ip_mat_create(h, w, k, 2);
    add2 = ip_mat_create(h, w, k, 5);
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

    sub = ip_mat_sub(add1, add2);

    printf("\n A - B: \n");
    ip_mat_show(sub);
    printf("\n");

    ip_mat_free(sub);
    ip_mat_free(add1);
    ip_mat_free(add2);

    printf("\n TEST ip_mat_mean: \n");
    a = ip_mat_create(2, 2, 3, 6.0f);
    b = ip_mat_create(2, 2, 3, 3.0f);

    c = ip_mat_mean(a, b);

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
    d = ip_mat_create(2, 2, 3, 3.0f);
    e = ip_mat_mul_scalar(d, 3.0f);

    printf("\n Matrice A: \n");
    ip_mat_show(d);
    printf("\n A x 3:\n");
    ip_mat_show(e);

    ip_mat_free(d);
    ip_mat_free(e);

    printf("\n TEST ip_mat_padding: \n");
    no_pad = ip_mat_create(3, 3, 3, 1.0f);
    pad = ip_mat_padding(no_pad, 2, 2);

    ip_mat_show(no_pad);
    ip_mat_show(pad);

    ip_mat_free(no_pad);
    ip_mat_free(pad);

    printf("\n TEST ip_mat_sum, ip_mat_show_stats, ip_mat_compute_stats: \n");
    aa = ip_mat_create(2, 2, 3, 3.0f);
    bb = ip_mat_create(2, 2, 3, 6.0f);
    set_val(aa, 1, 1, 1, 300.0f);
    cc = ip_mat_sum(aa, bb);
    printf("\n\nSOMMA:\n");
    ip_mat_show(cc);
    ip_mat_show_stats(cc);
    ip_mat_free(cc);
    ip_mat_free(aa);
    ip_mat_free(bb);

    printf("\n TEST ip_mat_copy: \n");
    a01 = ip_mat_create(3, 3, 3, 2.445f);
    set_val(a01, 1, 1, 1, 300.0f);
    a02 = ip_mat_copy(a01);
    printf("\nmatrice a1:");
    ip_mat_show(a01);
    printf("\nmatrice a2:");
    ip_mat_show(a02);
    ip_mat_free(a01);
    ip_mat_free(a02);

    printf("\n TEST ip_mat_concat: \n");
    a03 = ip_mat_create(3, 3, 3, 1.1f);
    a04 = ip_mat_create(3, 3, 3, 2.2f);
    set_val(a04, 1, 1, 1, 4.0f);
    a05 = ip_mat_concat(a03, a04, 0);
    printf("\nmatrice a05 altezza\n");
    ip_mat_show(a05);
    a06 = ip_mat_concat(a03, a04, 1);
    printf("\nmatrice a06 largezza\n");
    ip_mat_show(a06);
    a07 = ip_mat_concat(a03, a04, 2);
    printf("\nmatrice a07 canali\n");
    ip_mat_show(a07);
    ip_mat_free(a03);
    ip_mat_free(a04);
    ip_mat_free(a05);
    ip_mat_free(a06);
    ip_mat_free(a07);

    printf("\nTest ip_mat_corrupt\n");
    a08 = ip_mat_create(3, 3, 3, 50.0f);
    a09 = ip_mat_corrupt(a08, 10.0f);
    printf("Matrice a08:\n");
    ip_mat_show(a08);
    printf("Matrice a09 (a08 con corrupt):\n");
    ip_mat_show(a09);
    ip_mat_free(a08);
    ip_mat_free(a09);


    printf("\n TEST ip_mat_brighten image: \n");
    image = bm_load("flower.bmp");
    img = bitmap_to_ip_mat(image);
    bright = ip_mat_brighten(img, 100.0f);
    clamp(bright, 0.f, 255.f);
    bbmp = ip_mat_to_bitmap(bright);
    bm_save(bbmp, "flower_b.bmp");
    bm_free(image);
    bm_free(bbmp);
    ip_mat_free(bright);
    ip_mat_free(img);

    printf("\n TEST ip_mat_brighten darken image: \n");
    image = bm_load("flower.bmp");
    img = bitmap_to_ip_mat(image);
    bright = ip_mat_brighten(img, -100.0f);
    clamp(bright, 0.f, 255.f);
    bbmp = ip_mat_to_bitmap(bright);
    bm_save(bbmp, "flower_d.bmp");
    bm_free(image);
    bm_free(bbmp);
    ip_mat_free(bright);
    ip_mat_free(img);

    printf("\n TEST ip_mat_concat image: \n");
    image = bm_load("flower.bmp");
    img = bitmap_to_ip_mat(image);
    concatenazione = ip_mat_concat(img, img, 1);
    bbmp = ip_mat_to_bitmap(concatenazione);
    bm_save(bbmp, "flower_unito.bmp");
    bm_free(image);
    bm_free(bbmp);
    ip_mat_free(concatenazione);
    ip_mat_free(img);

    printf("\n TEST ip_mat_corrupt image: \n");
    image = bm_load("flower.bmp");
    img = bitmap_to_ip_mat(image);
    corruzione = ip_mat_corrupt(img, 255.0f);
    bbmp = ip_mat_to_bitmap(corruzione);
    bm_save(bbmp, "flower2_corrotto.bmp");
    bm_free(image);
    bm_free(bbmp);
    ip_mat_free(corruzione);
    ip_mat_free(img);

    printf("\n TEST ip_mat_blend image: \n");
    im1 = bm_load("mongolfiere.bmp");
    im2 = bm_load("flower2.bmp");
    b1 = bitmap_to_ip_mat(im1);
    b2 = bitmap_to_ip_mat(im2);
    blend = ip_mat_blend(b1, b2, 0.5f);
    b3 = ip_mat_to_bitmap(blend);
    bm_save(b3, "blend.bmp");
    ip_mat_free(b1);
    ip_mat_free(b2);
    ip_mat_free(blend);
    bm_free(im1);
    bm_free(im2);
    bm_free(b3);

    printf("\n TEST ip_mat_padding image: \n");
    caf = bm_load("flower.bmp");
    caf_nopad = bitmap_to_ip_mat(caf);
    caf_pad = ip_mat_padding(caf_nopad, 100.0f, 100.0f);
    bbmp = ip_mat_to_bitmap(caf_pad);
    bm_save(bbmp, "flower_padding.bmp");
    bm_free(caf);
    bm_free(bbmp);
    ip_mat_free(caf_nopad);
    ip_mat_free(caf_pad);

    printf("\n TEST convoluzione average image: \n");
    caf = bm_load("flower.bmp");
    no_blur = bitmap_to_ip_mat(caf);
    ker_avg = create_average_filter(5, 5, 3);
    blur = ip_mat_convolve(no_blur, ker_avg);
    bbmp = ip_mat_to_bitmap(blur);
    bm_save(bbmp, "flower_average.bmp");
    bm_free(caf);
    bm_free(bbmp);
    ip_mat_free(no_blur);
    ip_mat_free(blur);
    ip_mat_free(ker_avg);

    printf("\nTest ip_mat_to_gray image: \n");
    g_original_bmp = bm_load("flower.bmp");
    g_original = bitmap_to_ip_mat(g_original_bmp);
    gray = ip_mat_to_gray_scale(g_original);
    gray_bmp = ip_mat_to_bitmap(gray);
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

    printf("\n TEST convoluzione gaussian image: \n");
    caf = bm_load("flower.bmp");
    no_gauss = bitmap_to_ip_mat(caf);
    ker_gauss = create_gaussian_filter(3, 3, 3, 5.0f);
    gauss = ip_mat_convolve(no_gauss, ker_gauss);
    clamp(gauss, 0.0f, 255.0f);
    bbmp = ip_mat_to_bitmap(gauss);
    bm_save(bbmp, "flower_gauss.bmp");
    bm_free(caf);
    bm_free(bbmp);
    ip_mat_free(no_gauss);
    ip_mat_free(gauss);
    ip_mat_free(ker_gauss);

    printf("\n TEST convoluzione edge image: \n");
    caf = bm_load("flower.bmp");
    no_gauss = bitmap_to_ip_mat(caf);
    ker_gauss = create_edge_filter();
    gauss = ip_mat_convolve(no_gauss, ker_gauss);
    clamp(gauss, 0.0f, 255.0f);
    bbmp = ip_mat_to_bitmap(gauss);
    bm_save(bbmp, "flower_edge.bmp");
    bm_free(caf);
    bm_free(bbmp);
    ip_mat_free(no_gauss);
    ip_mat_free(gauss);
    ip_mat_free(ker_gauss);

    printf("\n TEST convoluzione emboss image: \n");
    caf = bm_load("fullmoon.bmp");
    no_gauss = bitmap_to_ip_mat(caf);
    ker_gauss = create_emboss_filter();
    gauss = ip_mat_convolve(no_gauss, ker_gauss);
    clamp(gauss, 0.0f, 255.0f);
    bbmp = ip_mat_to_bitmap(gauss);
    bm_save(bbmp, "fullmoon_emboss.bmp");
    bm_free(caf);
    bm_free(bbmp);
    ip_mat_free(no_gauss);
    ip_mat_free(gauss);
    ip_mat_free(ker_gauss);


    return 0;
}