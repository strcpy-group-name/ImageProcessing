/*
ID: 14
Gruppo: Cristian Nicolae Lupascu 880140, Marco Biondo 879994, Zanella Veronica 826585
*/

/*
   Created by Sebastiano Vascon on 23/03/20.
   */

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>

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

ip_mat *bitmap_to_ip_mat(Bitmap *img)
{
    unsigned int i = 0, j = 0;

    unsigned char R, G, B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat *out = ip_mat_create(h, w, 3, 0.f);

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
    {
        int index = k + j * a->k + i * a->k * a->w;
        return a->data[index];
    }
    else
    {
        printf("errore get_val: a era NULL oppure indici fuori intervallo\n");
        printf("%d %d %d\n", i, j, k);
        exit(1);
    }
}

void set_val(ip_mat *a, unsigned int i, unsigned int j, unsigned int k, float v)
{
    if (a && i < a->h && j < a->w && k < a->k)
    {

        int index = k + j * a->k + i * a->k * a->w;
        a->data[index] = v;
    }
    else
    {
        printf("Errore set_val: a era NULL oppure indici fuori intervallo\n");
        printf("%d %d %d\n", i, j, k);
        exit(1);
    }
}

/**** FUNZIONE AUSILIARIA ***/
/*
* Calcola gli indici di una matrice tridimensionale partendo dall'indice della matrice lineare
*/
void compute_indexes(unsigned int i, unsigned int *ih, unsigned int *iw, unsigned int *ik, unsigned int w, unsigned int k)
{
    if (w != 0 && k != 0)
    {
        *ik = i % k;
        *iw = (i / k) % w;
        *ih = i / (w * k);
    }
    else
    {
        printf("errore compute_indexes: divisione per zero\n");
        exit(1);
    }
}

float get_normal_random(float media, float std)
{

    float y1 = ((float)(rand()) + 1.) / ((float)(RAND_MAX) + 1.);
    float y2 = ((float)(rand()) + 1.) / ((float)(RAND_MAX) + 1.);
    float num = cos(2 * PI * y2) * sqrt(-2. * log(y1));

    return media + num * std;
}
/* PARTE 1 */

ip_mat *ip_mat_create(unsigned int h, unsigned int w, unsigned int k, float v)
{
    ip_mat *mat;
    unsigned int ih, iw, ik, i;
    mat = (ip_mat *)malloc(sizeof(ip_mat));
    if (!mat)
    {
        printf("Errore ip_mat_create: spazio non sufficiente a contenere la matrice\n");
        exit(1);
    }
    mat->w = w;
    mat->h = h;
    mat->k = k;
    mat->stat = (stats *)malloc(sizeof(stats) * k);
    if (!mat->stat)
    {
        printf("Errore ip_mat_create: spazio non sufficiente per contenere le statistiche\n");
        exit(1);
    }
    mat->data = (float *)malloc(sizeof(float) * h * w * k);
    if (!mat->data)
    {
        printf("Errore ip_mat_create: spazio non sufficiente a contenere i dati\n");
        exit(1);
    }
    for (i = 0; i < w * h * k; i++)
    {
        compute_indexes(i, &ih, &iw, &ik, w, k);
        set_val(mat, ih, iw, ik, v);
    }

    return mat;
}

void ip_mat_free(ip_mat *a)
{
    if (a)
    {
        free(a->stat);
        free(a->data);
        free(a);
    }
}

void ip_mat_init_random(ip_mat *t, float mean, float var)
{
    if (t)
    {
        unsigned int ih, iw, ik, i;
        float gaussian, random;

        for (i = 0; i < t->w * t->h * t->k; i++)
        {
            compute_indexes(i, &ih, &iw, &ik, t->w, t->k);

            random = get_normal_random(mean, var);
            gaussian = (1.f / sqrt(2.f * PI * var * var)) * exp(-(pow(random - mean, 2.f) / (2.f * var * var)));
            set_val(t, ih, iw, ik, gaussian);
        }
    }
    else
    {
        printf("Errore ip_mat_init_random: t era NULL\n");
        exit(1);
    }
}

ip_mat *ip_mat_subset(ip_mat *t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end)
{

    if (t)
    {
        ip_mat *subset;
        unsigned int w, h, k;
        unsigned int ih, iw, ik, i;
        float val;
        /* x comodità */
        h = row_end - row_start;
        w = col_end - col_start;
        k = t->k;
        subset = ip_mat_create(h, w, k, 0.f);
        for (i = 0; i < subset->w * subset->h * subset->k; i++)
        {
            compute_indexes(i, &ih, &iw, &ik, subset->w, subset->k);
            val = get_val(t, row_start + ih, col_start + iw, ik);
            set_val(subset, ih, iw, ik, val);
        }
        return subset;
    }
    else
    {
        printf("Errore ip_mat_subset: t era NULL\n");
        exit(1);
    }
}

ip_mat *ip_mat_sum(ip_mat *a, ip_mat *b)
{
    if (a && b && a->h == b->h && a->w == b->w && a->k == b->k)
    {
        ip_mat *c;
        unsigned int ih, iw, ik, i;
        float sum;

        c = ip_mat_create(a->h, a->w, a->k, 0.f);

        for (i = 0; i < a->w * a->h * a->k; i++)
        {
            compute_indexes(i, &ih, &iw, &ik, a->w, a->k);
            sum = get_val(a, ih, iw, ik) + get_val(b, ih, iw, ik);
            set_val(c, ih, iw, ik, sum);
        }
        return c;
    }
    else
    {
        printf("Errore ip_mat_sum: a o b erano NULL oppure di dimensione differente\\n");
        exit(1);
    }
}

ip_mat *ip_mat_sub(ip_mat *a, ip_mat *b)
{
    if (a && b && a->h == b->h && a->w == b->w && a->k == b->k)
    {
        ip_mat *c;
        unsigned int ih, iw, ik, i;
        float sum;

        c = ip_mat_create(a->h, a->w, a->k, 0.f);

        for (i = 0; i < a->w * a->h * a->k; i++)
        {
            compute_indexes(i, &ih, &iw, &ik, a->w, a->k);
            sum = get_val(a, ih, iw, ik) - get_val(b, ih, iw, ik);
            set_val(c, ih, iw, ik, sum);
        }
        return c;
    }
    else
    {
        printf("Errore ip_mat_sub: a o b erano NULL oppure di dimensione differente\\n");
        exit(1);
    }
}

ip_mat *ip_mat_mean(ip_mat *a, ip_mat *b)
{
    if (a && b && a->h == b->h && a->w == b->w && a->k == b->k)
    {
        unsigned int ih, iw, ik, i;
        float sum;
        ip_mat *c = ip_mat_create(a->h, a->w, a->k, 0.f);
        for (i = 0; i < a->w * a->h * a->k; i++)
        {
            compute_indexes(i, &ih, &iw, &ik, a->w, a->k);
            sum = get_val(a, ih, iw, ik) + get_val(b, ih, iw, ik);
            set_val(c, ih, iw, ik, sum / 2.f);
        }
        return c;
    }
    else
    {
        printf("Errore ip_mat_mean: a o b erano NULL oppure di dimensione differente\n");
        exit(1);
    }
}

ip_mat *ip_mat_mul_scalar(ip_mat *a, float c)
{
    if (a)
    {
        unsigned int ih, iw, ik, i;
        float sum;
        ip_mat *mat = ip_mat_create(a->h, a->w, a->k, 0.f);

        for (i = 0; i < a->w * a->h * a->k; i++)
        {
            compute_indexes(i, &ih, &iw, &ik, a->w, a->k);
            sum = get_val(a, ih, iw, ik) * c;
            set_val(mat, ih, iw, ik, sum);
        }
        return mat;
    }
    else
    {
        printf("Errore ip_mat_mul_scalar: a era NULL\n");
        exit(1);
    }
}

ip_mat *ip_mat_add_scalar(ip_mat *a, float c)
{
    if (a)
    {
        unsigned int ih, iw, ik, i;
        float sum;
        ip_mat *mat = ip_mat_create(a->h, a->w, a->k, 0.f);

        for (i = 0; i < a->w * a->h * a->k; i++)
        {
            compute_indexes(i, &ih, &iw, &ik, a->w, a->k);
            sum = get_val(a, ih, iw, ik) + c;
            set_val(mat, ih, iw, ik, sum);
        }
        return mat;
    }
    else
    {
        printf("Errore ip_mat_add_scalar: a era NULL\n");
        exit(1);
    }
}

void compute_stats(ip_mat *t)
{
    if (t)
    {
        unsigned int ih, iw, ik;
        for (ik = 0; ik < t->k; ik++) /*Per ogni canale calcolo e scrivo le statistiche*/
        {
            float min = get_val(t, 0, 0, ik);
            float max = get_val(t, 0, 0, ik);
            float tot = 0.f;
            for (ih = 0; ih < t->h; ih++)
            {
                for (iw = 0; iw < t->w; iw++)
                {
                    float current = get_val(t, ih, iw, ik);
                    if (current < min)
                        min = current;
                    if (current > max)
                        max = current;
                    tot += current;
                }
            }
            t->stat[ik].max = max;
            t->stat[ik].min = min;
            t->stat[ik].mean = tot / (t->h * t->w);
        }
    }
    else
    {
        printf("Errore compute_stats: t era NULL\n");
        exit(1);
    }
}

ip_mat *ip_mat_copy(ip_mat *in)
{
    if (in)
    {
        unsigned int h = in->h;
        unsigned int w = in->w;
        unsigned int k = in->k;
        unsigned int ih, iw, ik, i;
        float current;
        ip_mat *new = ip_mat_create(h, w, k, 0.f);
        for (i = 0; i < h * w * k; i++)
        {
            compute_indexes(i, &ih, &iw, &ik, in->w, in->k);
            current = get_val(in, ih, iw, ik);
            set_val(new, ih, iw, ik, current);
        }
        return new;
    }
    else
    {
        printf("Errore ip_mat_copy: in era NULL\n");
        exit(1);
    }
}

ip_mat *ip_mat_concat(ip_mat *a, ip_mat *b, int dimensione)
{
    if (a && b)
    {
        unsigned int ih, iw, ik, i;
        ip_mat *new = NULL;
        switch(dimensione)
        {
            case 0:
            {
                new = ip_mat_create(a->h + b->h, a->w, a->k, 0.f);
                for (i = 0; i < (a->h + b->h) * a->w * a->k; i++)
                {
                    float current = 0;
                    compute_indexes(i, &ih, &iw, &ik, a->w, a->k);
                    if (ih < a->h)
                        current = get_val(a, ih, iw, ik);
                    else
                        current = get_val(b, ih - a->h, iw, ik);
                    set_val(new, ih, iw, ik, current);
                }
                break;
            }
            case 1:
            {
                new = ip_mat_create(b->h, a->w + b->w, b->k, 0.f);
                for (i = 0; i < a->h * (a->w + b->w) * a->k; i++)
                {
                    float current;
                    compute_indexes(i, &ih, &iw, &ik, a->w + b->w, a->k);
                    if (iw < a->w)
                        current = get_val(a, ih, iw, ik);
                    else
                        current = get_val(b, ih, iw - a->w, ik);
                    set_val(new, ih, iw, ik, current);
                }
                break;
            }
            case 2:
            {
                new = ip_mat_create(b->h, b->w, a->k + b->k, 0.f);
                for (i = 0; i < a->h * a->w * (a->k + b->k); i++)
                {
                    float current;
                    compute_indexes(i, &ih, &iw, &ik, a->w, a->k + b->k);
                    if (ik < a->k)
                        current = get_val(a, ih, iw, ik);
                    else
                        current = get_val(b, ih, iw, ik - a->k);
                    set_val(new, ih, iw, ik, current);
                }
                break;
            }
            default: 
            {
                printf("Errore ip_mat_concat: valore di dimensione non esistente\n");
                exit(1);
            }
        }
        return new;
    }
    else
    {
        printf("Errore ip_mat_concat: a o b erano NULL\n");
        exit(1);
    }
}

/**** PARTE 2 ****/

ip_mat *ip_mat_brighten(ip_mat *a, float bright)
{
    if (a)
        return ip_mat_add_scalar(a, bright);
    else
    {
        printf("Errore ip_mat_brighten: a era NULL\n");
        exit(1);
    }
}

ip_mat *ip_mat_blend(ip_mat *a, ip_mat *b, float alpha)
{
    if (a && b)
    {
        ip_mat *m_a = ip_mat_mul_scalar(a, alpha);
        ip_mat *m_b = ip_mat_mul_scalar(b, (1 - alpha));
        ip_mat *mat = ip_mat_sum(m_a, m_b);
        ip_mat_free(m_a);
        ip_mat_free(m_b);
        return mat;
    }
    else
    {
        printf("Errore ip_mat_blend: a o b erano NULL\n");
        exit(1);
    }
}

ip_mat *ip_mat_to_gray_scale(ip_mat *in)
{
    if (in)
    {
        unsigned int ih, iw, ik, i;

        float val;
        ip_mat *gray = ip_mat_create(in->h, in->w, in->k, 0.0f);
        for (i = 0; i < in->h * in->w * in->k; i += 3)
        {
            compute_indexes(i, &ih, &iw, &ik, in->w, in->k);
            val = (get_val(in, ih, iw, ik) + get_val(in, ih, iw, ik + 1) + get_val(in, ih, iw, ik + 2)) / 3;

            set_val(gray, ih, iw, ik, val);
            set_val(gray, ih, iw, ik + 1, val);
            set_val(gray, ih, iw, ik + 2, val);
        }
        return gray;
    }
    else
    {
        printf("Errore ip_mat_to_gray_scale: in era NULL\n");
        exit(1);
    }
}

ip_mat *ip_mat_corrupt(ip_mat *a, float amount)
{
    if (a)
    {
        ip_mat *corrupt = ip_mat_create(a->h, a->w, a->k, 0.f);
        unsigned int i, ih, iw, ik;
        for (i = 0; i < a->h * a->w * a->k; i++)
        {
            compute_indexes(i, &ih, &iw, &ik, a->w, a->k);
            set_val(corrupt, ih, iw, ik, get_val(a, ih, iw, ik) + (get_normal_random(0, amount / 2.f)));
        }
        return corrupt;
    }
    else
    {
        printf("Errore ip_mat_corrupt: a era NULL\n");
        exit(1);
    }
}

/**** PARTE 3 ****/

ip_mat *ip_mat_padding(ip_mat *a, unsigned int pad_h, unsigned int pad_w)
{
    if (a)
    {
        ip_mat *nuova = ip_mat_create(a->h + 2 * pad_h, a->w + 2 * pad_w, a->k, 0.f);
        float v_a;
        unsigned int ih, iw, ik, i;
        for (i = 0; i < a->w * a->k * a->h; i++)
        {
            compute_indexes(i, &ih, &iw, &ik, a->w, a->k);
            v_a = get_val(a, ih, iw, ik);
            set_val(nuova, ih + pad_h, iw + pad_w, ik, v_a);
        }
        return nuova;
    }
    else
    {
        printf("Errore ip_mat_padding: a era NULL\n");
        exit(1);
    }
}

void clamp(ip_mat *t, float low, float high)
{
    if (t)
    {
        unsigned int ih, iw, ik, i;
        float v_a;
        for (i = 0; i < t->w * t->h * t->k; i++)
        {
            compute_indexes(i, &ih, &iw, &ik, t->w, t->k);
            v_a = get_val(t, ih, iw, ik);
            if (v_a > high)
                v_a = high;
            if (v_a < low)
                v_a = low;
            set_val(t, ih, iw, ik, v_a);
        }
    }
    else
    {
        printf("Errore clamp: t era NULL\n");
        exit(1);
    }
}

void rescale(ip_mat *t, float new_max)
{
    if (t)
    {
        unsigned int ih, iw, ik, i;
        float v_a;
        compute_stats(t);
        for (i = 0; i < t->h * t->w * t->k; i++)
        {
            compute_indexes(i, &ih, &iw, &ik, t->w, t->k);
            v_a = ((get_val(t, ih, iw, ik) - t->stat[ik].min) / (t->stat[ik].max - t->stat[ik].min)) * new_max;
            set_val(t, ih, iw, ik, v_a);
        }
    }
    else
    {
        printf("Errore rescale: t era NULL\n");
        exit(1);
    }
}

/**** FUNZIONE AUSILIARIA ***/
/*
* Calcola il valore della matrice in posizione (i,j,k) seguendo la regola della convoluzione
*/
float calculate_convolution(ip_mat *a, ip_mat *ker, int i, int j, int k)
{
    if (a && ker && a->k == ker->k)
    {
        unsigned int ih, iw, ik, idx;
        float acc = 0.0f;
        for (idx = k; idx < ker->h * ker->w * ker->k; idx += ker->k)
        {
            compute_indexes(idx, &ih, &iw, &ik, ker->w, ker->k);
            acc += get_val(a, ih + i, iw + j, k) * get_val(ker, ih, iw, k);
        }
        return acc;
    }
    else
    {
        printf("Errore calculate_convolution: a o ker erano NULL oppure il numero di canali non era lo stesso\n");
        exit(1);
    }
}

ip_mat *ip_mat_convolve(ip_mat *a, ip_mat *f)
{
    if (a && f && a->k == f->k)
    {
        unsigned int padh_amt = (f->h - 1) / 2;
        unsigned int padw_amt = (f->w - 1) / 2;
        unsigned int i = 0;
        ip_mat *pad_a = ip_mat_padding(a, padh_amt, padw_amt);
        ip_mat *mat = ip_mat_create(a->h, a->w, a->k, 0.0f);
        unsigned int ih, iw, ik;
        float val;
        
        for(i=0; i < a->w * a->h * a->k; i++){
            compute_indexes(i, &ih, &iw, &ik, a->w, a->k);
            val = calculate_convolution(pad_a, f, ih, iw, ik);
            set_val(mat, ih, iw, ik, val);
        }
        i=0;
        ip_mat_free(pad_a);
        return mat;
    }
    else
    {
        printf("Errore ip_mat_convolve: a o f erano NULL oppure il numero di canali non era lo stesso\n");
        exit(1);
    }
}

ip_mat *create_average_filter(unsigned int w, unsigned int h, unsigned int k)
{
    ip_mat *filter = ip_mat_create(h, w, k, 1 / (float)(w * h));
    return filter;
}

ip_mat *create_gaussian_filter(unsigned int w, unsigned int h, unsigned int k, float sigma)
{
    int cx, cy;
    unsigned int ih, iw, ik, i;
    int dx, dy;
    float val, acc;
    ip_mat *filter = ip_mat_create(h, w, k, 0.f);
    ip_mat *filter1;
    cx = (w - 1) / 2;
    cy = (h - 1) / 2;
    acc = 0;
    for (i = 0; i < w * h * k; i++)
    {
        compute_indexes(i, &ih, &iw, &ik, w, k);
        dx = ih - cx;
        dy = iw - cy;
        val = (1 / (2 * PI * sigma * sigma)) * exp(-(dx * dx + dy * dy) / (2 * sigma * sigma));
        set_val(filter, ih, iw, ik, val);
        acc += val;
    }
    acc /= 3;
    filter1 = ip_mat_mul_scalar(filter, 1 / acc);
    ip_mat_free(filter);
    return filter1;
}

ip_mat *create_edge_filter()
{
    ip_mat *filter = ip_mat_create(3, 3, 3, -1.0f);
    set_val(filter, 1, 1, 0, 8.0f);
    set_val(filter, 1, 1, 1, 8.0f);
    set_val(filter, 1, 1, 2, 8.0f);
    return filter;
}

ip_mat *create_emboss_filter()
{
    ip_mat *filter = ip_mat_create(3, 3, 3, 1.0f);
    int i;
    for (i = 0; i < 3; i++)
    {
        set_val(filter, 0, 0, i, -2.0f);
        set_val(filter, 1, 0, i, -1.0f);
        set_val(filter, 0, 1, i, -1.0f);
        set_val(filter, 0, 2, i, 0.0f);
        set_val(filter, 2, 0, i, 0.0f);
        set_val(filter, 2, 2, i, 2.0f);
    }
    return filter;
}

ip_mat *create_sharpen_filter()
{
    ip_mat *filter = ip_mat_create(3, 3, 3, 0.0f);
    int i;
    for (i = 0; i < 3; i++)
    {
        set_val(filter, 1, 0, i, -1.0f);
        set_val(filter, 0, 1, i, -1.0f);
        set_val(filter, 2, 1, i, -1.0f);
        set_val(filter, 1, 2, i, -1.0f);
        set_val(filter, 1, 1, i, 5.0f);
    }
    return filter;
}

/************************************************************/
/*                      FUNZIONI EXTRA                      */
/************************************************************/

/* scala di grigi con correzione approssimata per luminosità*/
ip_mat *ip_mat_to_gray_scale_lum_corr(ip_mat *in)
{
    if (in)
    {
        unsigned int ih, iw, ik;
        float r, g, b, val;

        ip_mat *gray = ip_mat_create(in->h, in->w, in->k, 0.0f);
        for (ih = 0; ih < (in->h); ih++)
            for (iw = 0; iw < (in->w); iw++)
            {
                r = get_val(in, ih, iw, 0);
                g = get_val(in, ih, iw, 1);
                b = get_val(in, ih, iw, 2);
                val = 0.299f * r + 0.587f * g + 0.114f * b;
                for (ik = 0; ik < (in->k); ik++)
                    set_val(gray, ih, iw, ik, val);
            }
        return gray;
    }
    else
    {
        printf("Errore ip_mat_to_gray_scale_lum_corr: in era NULL\n");
        exit(1);
    }
}

/**** FUNZIONE AUSILIARIA ***/
/*
* Linearizza un colore
*/
float gamma_correction_exp_to_linear(float v)
{
    float c_srgb, c_linear;
    c_srgb = v / 255.0f;
    c_linear = (c_srgb <= 0.04045f) ? c_srgb / 12.92f : powf((c_srgb + 0.055f) / 1.055f, 2.4f);
    return c_linear;
}

/* scala di grigi con gamma correction*/
ip_mat *ip_mat_to_gray_scale_gamma_corr(ip_mat *in)
{
    if (in)
    {
        unsigned int ih, iw, ik;
        float r_lin, g_lin, b_lin, y_lin, y_srgb;

        ip_mat *gray = ip_mat_create(in->h, in->w, in->k, 0.0f);
        for (ih = 0; ih < (in->h); ih++)
            for (iw = 0; iw < (in->w); iw++)
            {
                r_lin = gamma_correction_exp_to_linear(get_val(in, ih, iw, 0));
                g_lin = gamma_correction_exp_to_linear(get_val(in, ih, iw, 1));
                b_lin = gamma_correction_exp_to_linear(get_val(in, ih, iw, 2));
                y_lin = 0.2126f * r_lin + 0.7152f * g_lin + 0.0722f * b_lin;
                y_srgb = (y_lin <= 0.0031308f) ? 12.92f * y_lin : 1.055f * (powf(y_lin, 1.0f / 2.4f) - 0.055f);
                for (ik = 0; ik < in->k; ik++)
                    set_val(gray, ih, iw, ik, y_srgb * 255.0f);
            }
        return gray;
    }
    else
    {
        printf("Errore ip_mat_to_gray_scale_gamma_corr: in era NULL\n");
        exit(1);
    }
}

/*
* Copia w*h*k*sizeof(float) bytes dalla matrice in alla matrice di output, mantenendo le stesse dimensioni
*/
ip_mat *ip_mat_copy_mem(ip_mat *in)
{
    if(in)
    {
        ip_mat *mat = ip_mat_create(in->h, in->w, in->k, 0.f);
        memcpy(mat->data, in->data, sizeof(float)*in->w*in->h*in->k);
        return mat;
    }
    else
    {
        printf("Errore ip_mat_copy: in era NULL\n");
        exit(1);
    }
}