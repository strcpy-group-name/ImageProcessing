/*
   Created by Sebastiano Vascon on 23/03/20.
   */

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"
#include <math.h>

void ip_mat_show(ip_mat *t)
{
    unsigned int i, l, j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n", t->h, t->w, t->k);
    for (l = 0; l < t->k; l++)
    {
        printf("Slice %d\n", l);
        for (i = 0; i < t->h; i++)
        {
            for (j = 0; j < t->w; j++)
            {
                printf("%f ", get_val(t, i, j, l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

/* DA IMPLEMENTARE O DÀ ERRORE */
/*
   void ip_mat_show_stats(ip_mat *t)
   {
   unsigned int k;

   compute_stats(t);

   for (k = 0; k < t->k; k++)
   {
   printf("Channel %d:\n", k);
   printf("\t Min: %f\n", t->stat[k].min);
   printf("\t Max: %f\n", t->stat[k].max);
   printf("\t Mean: %f\n", t->stat[k].mean);
   }
   }
   */

ip_mat *bitmap_to_ip_mat(Bitmap *img)
{
    unsigned int i = 0, j = 0;

    unsigned char R, G, B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat *out = ip_mat_create(h, w, 3, 0);

    for (i = 0; i < h; i++) /* rows */
    {
        for (j = 0; j < w; j++) /* columns */
        {
            bm_get_pixel(img, j, i, &R, &G, &B);
            set_val(out, i, j, 0, (float)R);
            set_val(out, i, j, 1, (float)G);
            set_val(out, i, j, 2, (float)B);
        }
    }

    return out;
}

Bitmap *ip_mat_to_bitmap(ip_mat *t)
{

    Bitmap *b = bm_create(t->w, t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++) /* rows */
    {
        for (j = 0; j < t->w; j++) /* columns */
        {
            bm_set_pixel(b, j, i, (unsigned char)get_val(t, i, j, 0),
                    (unsigned char)get_val(t, i, j, 1),
                    (unsigned char)get_val(t, i, j, 2));
        }
    }
    return b;
}

float get_val(ip_mat *a, unsigned int i, unsigned int j, unsigned int k)
{
    if (a && i < a->h && j < a->w && k < a->k)
    { /* j>=0 and k>=0 and i>=0 is non sense*/
        int n_channels;
        n_channels = a->k;
        return a->data[i][j*n_channels+k];    /* modificato per linearizzazione da data[i][j][k] a quello attuale*/
    }
    else
    {
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat *a, unsigned int i, unsigned int j, unsigned int k, float v)
{
    if (a && i < a->h && j < a->w && k < a->k)
    {
        int n_channels;
        n_channels = a->k;
        a->data[i][j*n_channels+k] = v;      /* modificato per linearizzazione da data[i][j][k] a quello attuale*/
    }
    else
    {
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random()
{
    float y1 = ((float)(rand()) + 1.) / ((float)(RAND_MAX) + 1.);
    float y2 = ((float)(rand()) + 1.) / ((float)(RAND_MAX) + 1.);
    return cos(2 * PI * y2) * sqrt(-2. * log(y1));
}
/* PARTE 1 */

ip_mat *ip_mat_create(unsigned int h, unsigned int w, unsigned int k, float v){
    ip_mat *mat;
    int ih, iw, ik;
    mat = (ip_mat *)malloc(sizeof(ip_mat));
    if(!mat) exit(1);
    mat -> w = w;
    mat -> h = h;
    mat -> k = k;
    mat->stat = (stats *) malloc(sizeof(stats));
    if(!mat->stat) exit(1);
    mat->data = (float **)malloc(sizeof(float *)*h);
    if(!mat->data) exit(1);
    mat->data[0] = (float *)malloc(sizeof(float)*w*k*h);
    if(!mat->data[0]) exit(1);
    
    for(ih=0; ih<h; ih++){
        if(ih!=0)
            mat->data[ih] = mat->data[ih-1]+w*k;
        for(iw=0;iw<w;iw++)
            for(ik=0; ik<k; ik++)
                set_val(mat, ih, iw, ik, v);
    }
    return mat;
}

void ip_mat_free(ip_mat *a)
{
    free(a->data[0]);
    free(a->data);
    free(a->stat);
    free(a);
}

void ip_mat_init_random(ip_mat *t, float mean, float var)
{
    if(t)
    {
        int ih, iw, ik;

        for (ih = 0; ih < t->h; ih++)
        {
            for (iw = 0; iw < t->w; iw++)
                for (ik = 0; ik < t->k; ik++)
                {
                    float random = get_normal_random();
                    float gaussian = (1 / sqrt(2 * PI * var * var)) * pow(E, -(((random - mean) * (random - mean)) / (2 * var * var)));
                    set_val(t, ih, iw, ik, gaussian);
                }
        }
    }
}

ip_mat *ip_mat_subset(ip_mat *t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end){
    
    if(t)
    {
        ip_mat *subset;
        int w, h, k;
        int ih, iw, ik;
        /* x comodità */
        h = row_end - row_start;
        w = col_end - col_start;
        k = t->k;
        subset = ip_mat_create(h, w, k, 0);
        for (ih = 0; ih < h; ih++)
        {
            for (iw = 0; iw < w; iw++)
                for (ik = 0; ik < k; ik++)
                {
                    float val = get_val(t, row_start + ih, col_start + iw, ik);
                    set_val(subset, ih, iw, ik, val);
                }
        }
        return subset;
    }
    else
        return NULL;
}

ip_mat *ip_mat_sum(ip_mat *a, ip_mat *b){
    if(a && b && a->h == b->h && a->w == b->w && a->k == b->k){
        ip_mat *c;
        int ih, iw, ik;

        c = ip_mat_create(a->h, a->w, a->k, 0);

        for(ih = 0; ih < c->h; ih++)
            for(iw = 0; iw < c->w; iw++)
                for(ik = 0; ik < c->k; ik++){
                    float sum = get_val(a, ih, iw, ik) + get_val(b, ih, iw, ik);
                    set_val(c, ih, iw, ik, sum);
                }
        return c;
    }
    else
        return NULL;
}

ip_mat *ip_mat_sub(ip_mat *a, ip_mat *b){
    if(a && b && a->h == b->h && a->w == b->w && a->k == b->k){
        ip_mat *c;
        int ih, iw, ik;

        c = ip_mat_create(a->h, a->w, a->k, 0);

        for(ih = 0; ih < c->h; ih++)
            for(iw = 0; iw < c->w; iw++)
                for(ik = 0; ik < c->k; ik++){
                    float sum = get_val(a, ih, iw, ik) - get_val(b, ih, iw, ik);
                    set_val(c, ih, iw, ik, sum);
                }
        return c;
    }
    else
        return NULL;
}

ip_mat *ip_mat_mean(ip_mat *a, ip_mat *b)
{
    /*TODO: in caso di matrici di dimensione diversa, estendere la matrice più piccola usando concat*/
    if(a && b && a->h == b->h && a->w == b->w && a->k == b->k)
    {
        int ih, iw, ik;
        ip_mat *c = ip_mat_create(a->h, a->w, a->k, 0);
        for(ih=0; ih<c->h; ih++){
            for(iw=0;iw<c->w;iw++)
                for(ik=0; ik < c->k; ik++)
                {
                    float v_a = get_val(a, ih, iw, ik);
                    float v_b = get_val(b, ih, iw, ik);
                    set_val(c, ih, iw, ik, (v_a+v_b)/2);
                }
        }
        return c;
    }
    else
        return NULL;
}

ip_mat *ip_mat_mul_scalar(ip_mat *a, float c)
{
    if(a)
    {
        int ih, iw, ik;
        ip_mat *mat = ip_mat_create(a->h, a->w, a->k, 0);

        for (ih = 0; ih < a->h; ih++)
        {
            for (iw = 0; iw < a->w; iw++)
                for (ik = 0; ik < a->k; ik++)
                {
                    float v_a = get_val(a, ih, iw, ik);
                    set_val(mat, ih, iw, ik, v_a * c);
                }
        }
        return mat;
    }
    else
        return NULL;
}

ip_mat *ip_mat_add_scalar(ip_mat *a, float c)
{
    if(a)
    {
        int ih, iw, ik;
        ip_mat *mat = ip_mat_create(a->h, a->w, a->k, 0);

        for (ih = 0; ih < a->h; ih++)
        {
            for (iw = 0; iw < a->w; iw++)
                for (ik = 0; ik < a->k; ik++)
                {
                    float v_a = get_val(a, ih, iw, ik);
                    set_val(mat, ih, iw, ik, v_a + c);
                }
        }
        return mat;
    }
    else
        return NULL;
}

/* PARTE 2 */

void ip_mat_adjust_to_rgb(ip_mat *a)
{
    int ih, iw, ik;
    for (ih = 0; ih < a->h; ih++)
    {
        for (iw = 0; iw < a->w; iw++)
            for (ik = 0; ik < a->k; ik++)
            {
                float v_a = get_val(a, ih, iw, ik);
                if(v_a > 255) v_a = 255;
                if(v_a < 0) v_a = 0;
                set_val(a, ih, iw, ik, v_a);
            }
    }
}

ip_mat *ip_mat_brighten(ip_mat *a, float bright)
{
    ip_mat *mat = ip_mat_add_scalar(a, bright);
    ip_mat_adjust_to_rgb(mat); 
    return mat;
}

void expand(ip_mat *a)
{
    if(a)
    {
        int ih, iw, ik;
        for (ih = 0; ih < a->h; ih++)
        {
            for (iw = 0; iw < a->w; iw++)
                for (ik = 0; ik < a->k; ik++)
                {
                    float v_a = get_val(a, ih, iw, ik) * 255.0f;
                    set_val(a, ih, iw, ik, v_a);
                }
        }
    }
}

ip_mat *ip_mat_blend(ip_mat *a, ip_mat *b, float alpha)
{
    if(a && b)
    {
        ip_mat *m_a = ip_mat_mul_scalar(a, alpha/255.0f);
        ip_mat *m_b = ip_mat_mul_scalar(b, (1 - alpha)/255.0f);
        ip_mat *mat = ip_mat_sum(m_a, m_b);
        ip_mat_free(m_a);
        ip_mat_free(m_b);
        expand(mat);
        return mat;
    }
    return NULL;
}


ip_mat *ip_mat_to_gray_scale(ip_mat *in){
    int ih, iw, ik;
    
    float val, acc;
    ip_mat *gray = ip_mat_create(in->h, in->w, in->k, 0.0f);
    for(ih=0; ih<(in->h); ih++)
        for(iw=0; iw<(in->w); iw++){
            
            for(ik=0, acc=0; ik<in->k; ik++){
                acc += get_val(in, ih, iw, ik);
            }
            val = acc/3;
            for(ik=0; ik<in->k; ik++)
                set_val(gray, ih, iw, ik, val);
        }
    return gray;
}

ip_mat *ip_mat_to_gray_scale_lum_corr(ip_mat *in){
    int ih, iw, ik;  
    float r, g, b, val;

    ip_mat *gray = ip_mat_create(in->h, in->w, in->k, 0.0f);
    for(ih=0; ih<(in->h); ih++)
        for(iw=0; iw<(in->w); iw++){
            r = get_val(in, ih, iw, 0);
            g = get_val(in, ih, iw, 1);
            b = get_val(in, ih, iw, 2);
            val = 0.3*r + 0.59*g + 0.11*b;
            for(ik=0; ik<in->k; ik++)
                set_val(gray, ih, iw, ik, val);
        }
    return gray;
}


float gamma_correction_exp_to_linear(float v){
    float c_srgb, c_linear;
    c_srgb = v/255.0f;
    c_linear = (c_srgb<=0.04045f)? c_srgb/12.92f: powf((c_srgb+0.055f)/1.055, 2.4);
    return c_linear;
}

ip_mat *ip_mat_to_gray_scale_gamma_corr(ip_mat *in){
    int ih, iw, ik;    
    float r_lin, g_lin, b_lin, y_lin, y_srgb;

    ip_mat *gray = ip_mat_create(in->h, in->w, in->k, 0.0f);
    for(ih=0; ih<(in->h); ih++)
        for(iw=0; iw<(in->w); iw++){
            r_lin = gamma_correction_exp_to_linear(get_val(in, ih, iw, 0));
            g_lin = gamma_correction_exp_to_linear(get_val(in, ih, iw, 1));
            b_lin = gamma_correction_exp_to_linear(get_val(in, ih, iw, 2));
            y_lin = 0.2126*r_lin + 0.7152*g_lin + 0.0722*b_lin;
            y_srgb = (y_lin<=0.0031308f)? 12.92f*y_lin : 1.055f*(powf(y_lin, 1.0f/2.4f)-0.055f);
            for(ik=0; ik<in->k; ik++)
                set_val(gray, ih, iw, ik, y_srgb*255.0f);
        }
    return gray;
}

/* --- Function implemented by our group --- */

/*
   void compute_stats(ip_mat *t);

   ip_mat *ip_mat_copy(ip_mat *in);

   ip_mat *ip_mat_concat(ip_mat *a, ip_mat *b, int dimensione);

   ip_mat *ip_mat_add_scalar(ip_mat *a, float c);

   ip_mat *ip_mat_to_gray_scale(ip_mat *in);

   ip_mat *ip_mat_blend(ip_mat *a, ip_mat *b, float alpha);

   ip_mat *ip_mat_brighten(ip_mat *a, float bright);

   ip_mat *ip_mat_corrupt(ip_mat *a, float amount);

   ip_mat *ip_mat_convolve(ip_mat *a, ip_mat *f);

   ip_mat *ip_mat_padding(ip_mat *a, int pad_h, int pad_w);

   ip_mat *create_sharpen_filter();

   ip_mat *create_edge_filter();

   ip_mat *create_emboss_filter();

   ip_mat *create_average_filter(int w, int h, int k);

   ip_mat *create_gaussian_filter(int w, int h, int k, float sigma);

   void rescale(ip_mat *t, float new_max);

   void clamp(ip_mat *t, float low, float high);
   */

  /* COSE DA CHIEDERE X IMPLEMENTAZIONI MIGLIORI
  * brighten in percentuale
  * grayscale con correzione gamma
  */
