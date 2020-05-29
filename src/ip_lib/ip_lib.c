/* 
 * Laboratorio di programmazione A.A 2019/2020
 *
 * Gruppo: 95
 * Partecipanti: Grassi Marco
 *
 * Tutor: Sebastiano Vascon
 *
 * */

#include <stdio.h>
#include "ip_lib.h"
#include "../bmp/bmp.h"

#define malloc_error(str) "\033[0;33m" str ": Error\n" "\033[0m" "Error during memory allocation\n"
#define null_error(str) "\033[0;33m" str ": Invalidi input\n" "\033[0m" "ip_lib not initialized\n"
#define dim_error(str) "\033[0;31m" str ": Invalid input\n" "\033[0m" "Input arguments must have same dimentions\n"
#define coordinate_error(str) "\033[0;31m" str ": Invalid input\n" "\033[0m" "Invalid input coordinates\n"
#define error(msg) printf(msg); exit(1)
#define check_null(t, str) if (t == NULL) {error(null_error(str));}
#define different_size(a, b) if (a -> h != b -> h || a -> w != b -> w || a -> k != b -> k)
#define check_different_size(a, b, str) different_size(a, b) {error(dim_error(str));}

/**** STRUTTURA DATI ****/

ip_mat * ip_mat_create(unsigned int h, unsigned int w, unsigned int k, float v) {
    ip_mat *ip;
    unsigned int i, j, n;

    ip = (ip_mat *) malloc(sizeof(ip_mat));

    if (!ip) {
        error(malloc_error("ip_mat_create"));
    }

    ip -> w = w;
    ip -> h = h;
    ip -> k = k;

    ip -> stat = (stats *) malloc(k * sizeof(stats));

    if (!ip -> stat) {
        error(malloc_error("ip_mat_create"));
    }

    for (i = 0; i < k; i++) {
        ip -> stat[i].min = v;
        ip -> stat[i].max = v;
        ip -> stat[i].mean = v;
    }

    ip -> data = (float ***) malloc(h * sizeof(float**));

    if (!ip -> data) {
        error(malloc_error("ip_mat_create"));
    }

    for (i = 0; i < h; i++) {
        ip -> data[i] = (float **) malloc(w * sizeof(float*));

        if (!ip -> data[i]) {
            error(malloc_error("ip_mat_create"));
        }

        for (j = 0; j < w; j++) {
            ip -> data[i][j] = (float *) malloc(k * sizeof(float));

            if (!ip -> data[i][j]) {
                error(malloc_error("ip_mat_create"));
            }
        }
    }

    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            for (n = 0; n < k; n++) {
                ip -> data[i][j][n] = v;
            }
        }
    }

    return ip;
}


void ip_mat_free(ip_mat *a) {
    unsigned int i, j;

    if (a) {
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
}

float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    check_null(a, "get_val");

    if(i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[i][j][k];
    }else{
        error(coordinate_error("get_val"));
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    check_null(a, "set_val");

    if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        error(coordinate_error("get_val"));
    }
}

void compute_stats(ip_mat *t) {
    float accumulator, val;
    unsigned int  i, j, k;

    check_null(t, "compute_stats");

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

float get_normal_random(float mean, float std){
    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float num = cos(2*PI*y2)*sqrt(-2.*log(y1));

    return mean + num*std;
}

void ip_mat_init_random(ip_mat *t, float mean, float std) {
    unsigned int i, j, k;

    check_null(t, "ip_mat_init_random");
    
    for (i = 0; i < t -> h; i++) {
        for (j = 0; j < t -> w; j++) {
            for (k = 0; k < t -> k; k++) {
                set_val(t, i, j, k, get_normal_random(mean, std));
            }
        }
    }
}

ip_mat * ip_mat_copy(ip_mat *in) {
    unsigned int i, j, k;
    ip_mat *out;

    check_null(in, "ip_mat_copy");

    out = ip_mat_create(in -> h, in -> w, in -> k, 0);

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
    unsigned int i, io, j, jo, k;
    unsigned int row, col;
    ip_mat *sub;

    check_null(t, "ip_mat_subset");

    if (row_start > row_end || col_start > col_end){
        error("\033[0;31m" "ip_mat_subset: Invalid input\n" "\033[0m" "Starting row and column must be smaller than endig row and colum\n");
    }

    row = row_end - row_start;
    col = col_end - col_start;

    if (t -> h < row || t -> w < col) {
        error("\033[0;31m" "ip_mat_subset: Invalid input\n" "\033[0m" "Subset dimentions must be smaller than input ip_mat dimentions\n");
    }

    sub = ip_mat_create(row, col, t -> k, 0);

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
    unsigned int h, w, k;
    unsigned int i, j, l;

    check_null(a, "ip_mat_concat");
    check_null(b, "ip_mat_concat");

    different_size(a, b) {
        printf("\033[0;33m" "ip_mat_concat: Warning\n" "\033[0m");
        printf("Input arguments have different dimentions, this might lead to errors\n");
    } 

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
        error("\033[0;31m" "ip_mat_concat: Invalid input\n" "\033[0m" "Third argument must have a value between 0 and 2 \n");
    }

    chain = ip_mat_create(h, w, k, 0);

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

void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    check_null(t, "ip_mat_show");
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0;i<t->h;i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f \t", get_val(t,i,j,l));
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
    unsigned int h;
    unsigned int w;
    ip_mat * out;

    check_null(img, "bitmap_to_ip_mat");

    h = img->h;
    w = img->w;

    out = ip_mat_create(h, w,3,0);

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
    unsigned int i, j;
    Bitmap *b ;

    check_null(t, "ip_mat_to_bitmap");

    b = bm_create(t->w,t->h);

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

/**** OPERAZIONI MATEMATICHE ****/

ip_mat * ip_mat_sum(ip_mat *a, ip_mat *b) {
    ip_mat *sum;
    unsigned int i, j, k;

    check_null(a, "ip_mat_sum");
    check_null(b, "ip_mat_sum");
    check_different_size(a, b, "ip_mat_sum");

    sum = ip_mat_create(a -> h, a -> w, a -> k, 0);

    for (i = 0; i < a -> h; i++) {
        for (j = 0; j < a -> w; j++) {
            for (k = 0; k < a -> k; k++) {
                set_val(sum, i, j, k, get_val(a, i, j, k) + get_val(b, i, j, k));
            }
        }
    }

    return sum;
}

ip_mat * ip_mat_sub(ip_mat *a, ip_mat *b) {
    ip_mat *sub;
    unsigned int i, j, k;

    check_null(a, "ip_mat_sub");
    check_null(b, "ip_mat_sub");
    check_different_size(a, b, "ip_mat_sub");

    sub = ip_mat_create(a -> h, a -> w, a -> k, 0);

    for (i = 0; i < a -> h; i++) {
        for (j = 0; j < a -> w; j++) {
            for (k = 0; k < a -> k; k++) {
                set_val(sub, i, j, k, get_val(a, i, j, k) - get_val(b, i, j, k));
            }
        }
    }

    return sub;
}

ip_mat * ip_mat_mul_scalar(ip_mat *a, float c) {
    ip_mat *mult;
    unsigned int i, j, k;

    check_null(a, "ip_mat_mul_scalar");

    mult = ip_mat_create(a -> h, a -> w, a -> k, 0);

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
    unsigned int i, j, k;

    check_null(a, "ip_mat_add_scalar");

    sum = ip_mat_create(a -> h, a -> w, a -> k, 0);

    for (i = 0; i < a -> h; i++) {
        for (j = 0; j < a -> w; j++) {
            for (k = 0; k < a -> k; k++) {
                set_val(sum, i, j, k, get_val(a, i, j, k) + c);
            }
        }
    }

    return sum;
}

ip_mat * ip_mat_mean(ip_mat *a, ip_mat *b) {
    ip_mat *mean;
    unsigned int i, j, l;
    unsigned int h = a -> h;
    unsigned int w = a -> w;
    unsigned int k = a -> k;

    check_null(a, "ip_mat_mean");
    check_null(b, "ip_mat_mean");
    check_different_size(a, b, "ip_mat_mean");

    mean = ip_mat_create(h, w, k, 0);

    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            for (l = 0; l < k; l++) {
                set_val(mean, i, j, l, (get_val(a, i, j, l) + get_val(b, i, j, l)) / 2);
            }
        }
    }

    return mean;
}

/**** OPERAZIONI SU IMMAGINI ****/


ip_mat * ip_mat_to_gray_scale(ip_mat * in) {
    ip_mat *gray_scale;
    unsigned int i, j, k;
    float mean;

    check_null(in, "ip_mat_to_gray_scale");

    gray_scale = ip_mat_create(in -> h, in -> w, in -> k, 0);

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
    unsigned int i, j, l;
    unsigned int h = a -> h;
    unsigned int w = a -> w;
    unsigned int k = a -> k;
    float val;

    check_null(a, "ip_mat_blend");
    check_null(b, "ip_mat_blend");

    if (alpha < 0) {
        printf("\033[0;33m" "ip_mat_blend: Warning\n" "\033[0m");
        printf("Input alpha value out of range, capped at 0\n");
        alpha = 0;
    }

    if (alpha > 1) {
        printf("\033[0;33m" "ip_mat_blend: Warning\n" "\033[0m");
        printf("Input alpha value out of range, capped at 1\n");
        alpha = 1;
    }

    check_different_size(a, b, "ip_mat_blend");

    blend = ip_mat_create(h, w, k, 0);

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
    check_null(a, "ip_mat_brighten");

    return ip_mat_add_scalar(a, bright);
}

ip_mat * ip_mat_corrupt(ip_mat *a, float amount) {
    ip_mat *noise;
    ip_mat *out;

    check_null(a, "ip_mat_corrupt");

    noise = ip_mat_create(a -> h, a -> w, a -> k, 0);
    ip_mat_init_random(noise, 0, amount / 2);

    out = ip_mat_sum(a, noise);
    ip_mat_free(noise);

    return out;
}

/**** CONVOLUZIONE E FILTRI ****/

ip_mat * ip_mat_convolve (ip_mat *a, ip_mat *f) {
    ip_mat *pad;
    ip_mat *out;
    int pad_h, pad_w;
    unsigned int i, j, k, l, p;
    float val;

    check_null(a, "ip_mat_padding");
    check_null(a, "ip_mat_padding");

    pad_h = (f -> h - 1) / 2;
    pad_w = (f -> w - 1) / 2;

    pad = ip_mat_padding(a, pad_h, pad_w);
    out = ip_mat_create(a -> h, a -> w, a -> k, 0);

    for (i = 0; i < out -> h; i++) {
        for (j = 0; j < out -> w; j++) {
            for (k = 0; k < out -> k; k++) {
                val = 0;

                for (l = 0; l < f -> h; l++) {
                    for (p = 0; p < f -> w; p++) {
                        val += get_val(pad, i + l, j + p, k) * get_val(f, l, p, k);
                    }
                }

                set_val(out, i, j, k, val);
            }
        }
    }
    ip_mat_free(pad);

    return out;
}

ip_mat * ip_mat_padding(ip_mat *a, unsigned int pad_h, unsigned int pad_w) {
    ip_mat *out;
    unsigned int i, j, k;

    check_null(a, "ip_mat_padding");

    out = ip_mat_create(a -> h + 2 * pad_h, a -> w + 2 * pad_w, a -> k, 0);

    for (i = 0; i < a -> h; i++) {
        for (j = 0; j < a -> w; j++) {
            for (k = 0; k < a -> k; k++) {
                set_val(out, i + pad_h, j + pad_w, k, get_val(a, i, j, k));
            }
        }
    }

    return out;
}

ip_mat * create_sharpen_filter() {
    ip_mat *sharpen;
    unsigned int k;
    
    sharpen = ip_mat_create(3, 3, 3, 0);

    for (k = 0; k < 3; k++) {
        set_val(sharpen, 0, 1, k, -1);
        set_val(sharpen, 1, 0, k, -1);
        set_val(sharpen, 1, 1, k, 5);
        set_val(sharpen, 1, 2, k, -1);
        set_val(sharpen, 2, 1, k, -1);
    }

    return sharpen;
}

ip_mat * create_edge_filter() {
    ip_mat *edge;
    unsigned int k;

    edge = ip_mat_create(3, 3, 3, -1);

    for (k = 0; k < 3; k++) {
        set_val(edge, 1, 1, k, 8);
    }

    return edge;
}

ip_mat * create_emboss_filter() {
    ip_mat *emboss;
    unsigned int k;

    emboss = ip_mat_create(3, 3, 3, 1);

    for (k = 0; k < 3; k++) {
        set_val(emboss, 0, 0, k, -2);
        set_val(emboss, 0, 1, k, -1);
        set_val(emboss, 0, 2, k, 0);
        set_val(emboss, 1, 0, k, -1);
        set_val(emboss, 2, 0, k, 0);
        set_val(emboss, 2, 2, k, 2);
    }

    return emboss;
}

ip_mat * create_average_filter(unsigned int h, unsigned int w, unsigned int k) {
    ip_mat *average;
    float c;

    if (h % 2 == 0) h--;
    if (w % 2 == 0) w--;

    c = 1.0 / (h * w);
    average = ip_mat_create(h, w, k, c);

    return average;
}

ip_mat * create_gaussian_filter(unsigned int h, unsigned int w, unsigned int k, float sigma) {
    ip_mat *gauss;
    ip_mat *norm_gauss;
    unsigned int i, j, l;
    int cx, cy, x, y;
    float val, exponent, counter;

    if (h % 2 == 0) h--;
    if (w % 2 == 0) w--;

    cx = (h - 1) / 2;
    cy = (w - 1) / 2;
    counter = 0;

    gauss = ip_mat_create(h, w, k, 0);

    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            x = i - cx;
            y = j - cy;
            exponent = ((x * x) + (y * y)) / (2 * sigma * sigma);
            val = 1 / (2 * PI * sigma * sigma) * exp(- exponent);
            counter += val;
            
            for (l = 0; l < k; l++) {
                set_val(gauss, i, j, l, val);
            }
        }
    }

    norm_gauss = ip_mat_mul_scalar(gauss, 1 / counter);
    ip_mat_free(gauss);

    return norm_gauss;
}

void rescale(ip_mat *t, float new_max) {
    unsigned int i, j, k;
    float val;

    check_null(t, "rescale");

    compute_stats(t);

    for (i = 0; i < t -> h; i++) {
        for (j = 0; j < t -> w; j++) {
            for (k = 0; k < t -> k; k++) {
                val = get_val(t, i, j, k);
                val = new_max * ((val - t -> stat[k].min) / (t -> stat[k].max - t -> stat[k].min));
                set_val(t, i, j, k, val);
            }
        }
    }
}

void clamp(ip_mat *t, float low, float high) {
    unsigned int i, j, k;
    float val;

    check_null(t, "clamp");

    for (i = 0; i < t -> h; i++) {
        for (j = 0; j < t -> w; j++) {
            for (k = 0; k < t -> k; k++) {
                val = get_val(t, i, j, k);

                if (val < low) {
                    set_val(t, i, j, k, low);
                } else {
                    if (val > high) {
                        set_val(t, i, j, k, high);
                    }
                }
            }
        }
    }
}
