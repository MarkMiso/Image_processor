/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include "ip_lib.h"
#include "../bmp/bmp.h"

/**** GENERAZIONE ADT ****/

ip_mat * ip_mat_create(unsigned int h, unsigned int w, unsigned int k, float v) {
    
    ip_mat *ip;
    int i, j, n;

    ip = (ip_mat *) malloc(sizeof(ip_mat));

    if (ip) {
        ip -> w = w;
        ip -> h = h;
        ip -> k = k;

        ip -> stat = (stats *) malloc(k * sizeof(stats));
        for (i = 0; i < k; i++) {
            ip -> stat[i].min = v;
            ip -> stat[i].max = v;
            ip -> stat[i].mean = v;
        }

        ip -> data = (float ***) malloc(h * sizeof(float**));
        for (i = 0; i < h; i++) {
            ip -> data[i] = (float **) malloc(w * sizeof(float*));
            for (j = 0; j < w; j++) {
                ip -> data[i][j] = (float *) malloc(k * sizeof(float));
            }
        }

        for (i = 0; i < h; i++) {
            for (j = 0; j < w; j++) {
                for (n = 0; n < k; n++) {
                    ip -> data[i][j][n] = v;
                }
            }
        }
    } else {
        printf("ip_mat_create: malloc error");
        exit(1);
    }

    return ip;
}

void ip_mat_free(ip_mat *a) {
    int i, j;

    for (i = 0; i < a -> h; i++) {
        for (j = 0; j < a -> w; j++) {
            free(a -> data[i][j]);
        }
        free(a -> data[i]);
    }
    free(a -> data);

    free(a -> stat);

    free(a);
}

void compute_stats(ip_mat *t) {
    float accumulator, val;
    int  i, j, k;

    for (i = 0; i < t -> k; i++) {
        accumulator = 0;
        t -> stat[i].min = get_val(t, 0, 0, i);
        t -> stat[i].max = get_val(t, 0, 0, i);
        for (j = 0; j < t -> h; j++) {
            for (k = 0; k < t -> w; k++) {
                val = get_val(t, j, k, i);
                
                if (t -> stat[i].min > val)
                    t -> stat[i].min = val;
                
                if (t -> stat[i].max < val)
                    t -> stat[i].max = val;

                accumulator += val;
            }
        }
        t -> stat[i].mean = accumulator / ((t -> h) * (t -> w));
    }
}


void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0;i<t->h;i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", get_val(t,i,j,l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_show_stats(ip_mat * t){
    unsigned int k;

    compute_stats(t);

    for(k=0;k<t->k;k++){
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}

ip_mat * bitmap_to_ip_mat(Bitmap * img){
    unsigned int i=0,j=0;

    unsigned char R,G,B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }

    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t){

    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            bm_set_pixel(b, j,i, (unsigned char) get_val(t,i,j,0),
                    (unsigned char) get_val(t,i,j,1),
                    (unsigned char) get_val(t,i,j,2));
        }
    }
    return b;
}

float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(float mean, float std){
    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float num = cos(2*PI*y2)*sqrt(-2.*log(y1));

    return mean + num*std;
}

void ip_mat_init_random(ip_mat *t, float mean, float std) {
    int i, j, k;
    
    for (i = 0; i < t -> h; i++) {
        for (j = 0; j < t -> w; j++) {
            for (k = 0; k < t -> k; k++) {
                set_val(t, i, j, k, get_normal_random(mean, std));
            }
        }
    }
}

ip_mat * ip_mat_copy(ip_mat *in) {
    int i, j, k;
    ip_mat *out;

    out = ip_mat_create(in -> h, in -> w, in -> k, 0.0);

    for (i = 0; i < in -> h; i++) {
        for (j = 0; j < in -> w; j++) {
            for (k = 0; k < in -> k; k++) {
                set_val(out, i, j, k, get_val(in, i, j, k));
            }
        }
    }

    return out;
}

ip_mat * ip_mat_subset(ip_mat *t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end) {
    int i, io, j, jo, k;
    ip_mat *sub;

    if (row_start > row_end || col_start > col_end){
        printf("ip_mat_subset: Invalid input, starting row and column must be smaller then endig row and colum");
    }

    sub = ip_mat_create(row_end - row_start, col_end - col_start, t -> k, 0.0);

    io = 0;
    for (i = row_start; i < row_end; i++) {
        jo = 0;
        for (j = col_start; j < col_end; j++) {
            for (k = 0; k < t -> k; k++) {
                set_val(sub, io, jo, k, get_val(t, i, j, k));
            }
            jo++;
        }
        io++;
    }

    return sub;
}

ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione) {
    ip_mat *chain;
    unsigned int h, w, k, i, j, l;
 
    h = a -> h;
    w = a -> w;
    k = a -> k;
    
    if (dimensione == 0) {
        h += b -> h;
    } else if (dimensione == 1) {
        w += b -> w;
    } else if (dimensione == 2) {
        k += b -> k;
    } else {
        printf("ip_mat_concat: Invalid input");
        exit(1);
    }

    chain = ip_mat_create(h, w, k, 0.0);

    for (i = 0; i < a -> h; i++) {
        for (j = 0; j < a -> w; j++) {
            for (l = 0; l < a -> k; l++) {
                set_val(chain, i, j, l, get_val(a, i, j, l));
            }
        }
    }

    h -= b -> h;
    w -= b -> w;
    k -= b -> k;

    for (i = 0; i < b -> h; i++) {
        for (j = 0; j < b -> w; j++) {
            for (l = 0; l < b -> k; l++) {
                set_val(chain, i + h, j + w, l + k, get_val(b, i, j, l));
            }
        }
    }

    return chain;
}

int fit_to_size(ip_mat *in, int *h, int *w, int *k) {
    int out;

    if (*h != in -> h || *w != in -> w || *k != in -> k) {
        out = 1;
        if (*h > in -> h) *h = in -> h;
        if (*w > in -> w) *w = in -> w;
        if (*k > in -> k) *k = in -> k;
    } else {
        out = 0;
    }

    return out;
}

/**** OPERAZIONI MATEMATICHE ****/

ip_mat * ip_mat_sum(ip_mat *a, ip_mat *b) {
    ip_mat *sum;
    int i, j, k;
    if (a -> h == b -> h &&  a -> w == b -> w && a -> k == b -> k) {
        sum = ip_mat_create(a -> h, a -> w, a -> k, 0.0);

        for (i = 0; i < a -> h; i++) {
            for (j = 0; j < a -> w; j++) {
                for (k = 0; k < a -> k; k++) {
                    set_val(sum, i, j, k, get_val(a, i, j, k) + get_val(b, i, j, k));
                }
            }
        }
    } else {
        printf("ip_mat_sum: Invalid input");
        exit(1);
    }

    return sum;
}

ip_mat * ip_mat_sub(ip_mat *a, ip_mat *b) {
    ip_mat *sub;
    int i, j, k;
    if (a -> h == b -> h &&  a -> w == b -> w && a -> k == b -> k) {
        sub = ip_mat_create(a -> h, a -> w, a -> k, 0.0);

        for (i = 0; i < a -> h; i++) {
            for (j = 0; j < a -> w; j++) {
                for (k = 0; k < a -> k; k++) {
                    set_val(sub, i, j, k, get_val(a, i, j, k) - get_val(b, i, j, k));
                }
            }
        }
    } else {
        printf("ip_mat_sub: Invalid input");
        exit(1);
    }

    return sub;
}

ip_mat * ip_mat_mean(ip_mat *a, ip_mat *b) {
    ip_mat *mean;
    int i, j, l;
    int h = a -> h;
    int w = a -> w;
    int k = a -> k;

    if (fit_to_size(b, &h, &w, &k)) {
        printf("ip_mat_mean: Warning \nInput ip_mats have different dimentions the output ip_mat will be cropped to the smallest dimentions \n");
    }

    mean = ip_mat_create(h, w, k, 0.0);

    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            for (l = 0; l < k; l++) {
                set_val(mean, i, j, l, (get_val(a, i, j, l) + get_val(b, i, j, l)) / 2);
            }
        }
    }

    return mean;
}

ip_mat * ip_mat_mul_scalar(ip_mat *a, float c) {
    ip_mat *mult;
    int i, j, k;
    mult = ip_mat_create(a -> h, a -> w, a -> k, 0.0);

    for (i = 0; i < a -> h; i++) {
        for (j = 0; j < a -> w; j++) {
            for (k = 0; k < a -> k; k++) {
                set_val(mult, i, j, k, get_val(a, i, j, k) * c);
            }
        }
    }

    return mult;
}

ip_mat * ip_mat_add_scalar(ip_mat *a, float c) {
    ip_mat *sum;
    int i, j, k;
    sum = ip_mat_create(a -> h, a -> w, a -> k, 0.0);

    for (i = 0; i < a -> h; i++) {
        for (j = 0; j < a -> w; j++) {
            for (k = 0; k < a -> k; k++) {
                set_val(sum, i, j, k, get_val(a, i, j, k) + c);
            }
        }
    }

    return sum;
}

/**** OPERAZIONI SU IMMAGINI ****/

void ip_mat_clamp(ip_mat *in) {
    int i, j, k;
    float val;

    for (i = 0; i < in -> h; i++) {
        for (j = 0; j < in -> w; j++) {
            for (k = 0; k < in -> k; k++) {
                val = get_val(in, i, j, k);

                if (val < 0) {
                    set_val(in, i, j, k, 0);
                } else {
                    if (val > 255) {
                        set_val(in, i, j, k, 255);
                    }
                }
            }
        }
    }
}

ip_mat * ip_mat_to_gray_scale(ip_mat * in) {
    ip_mat *gray_scale;
    int i, j, k;
    float mean;

    gray_scale = ip_mat_create(in -> h, in -> w, in -> k, 0.0);

    for (i = 0; i < in -> h; i++) {
        for (j = 0; j < in -> w; j++) {
            mean = 0;

            for (k = 0; k < in -> k; k++) {
                mean += get_val(in, i, j, k);
            }

            mean = mean / (in -> k);

            for (k = 0; k < in -> k; k++) {
                set_val(gray_scale, i, j, k, mean);
            }
        }
    }

    return gray_scale;
}

ip_mat * ip_mat_blend(ip_mat *a, ip_mat *b, float alpha) {
    ip_mat *blend;
    int i, j, l, val;
    int h = a -> h;
    int w = a -> w;
    int k = a -> k;

    if (alpha < 0) {
        printf("ip_mat_blend: Warning \n");
        printf("Input alpha value out of range, capped at 0");
    }

    if (alpha > 1) {
        printf("ip_mat_blend: Warning \n");
        printf("Input alpha value out of range, capped at 1");
    }

    if (fit_to_size(b, &h, &w, &k)) {
        printf("ip_mat_blend: Warning \nInput ip_mats have different dimentions the output ip_mat will be cropped to the smallest dimentions \n");
    }

    blend = ip_mat_create(h, w, k, 0.0);

    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            for (l = 0; l < k; l ++) {
                val = (alpha * get_val(a, i, j, l)) + ((1 - alpha) * get_val(b, i, j, l));
                set_val(blend, i, j, l, val);
            }
        }
    }

    return blend;
}

ip_mat * ip_mat_brighten(ip_mat *a, float bright) {
    return ip_mat_add_scalar(a, bright);
}

ip_mat * ip_mat_corrupt(ip_mat *a, float amount) {
    ip_mat *noise;
    ip_mat *out;

    noise = ip_mat_create(a -> h, a -> w, a -> k, 0.0);
    ip_mat_init_random(noise, 0, amount / 2);

    out = ip_mat_sum(a, noise);
    ip_mat_free(noise);

    return out;
}
