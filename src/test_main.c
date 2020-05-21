#include <stdio.h>
#include <string.h>
#include "ip_lib/ip_lib.h"
#include "bmp/bmp.h"

int main () {
    
    /* copy test */

    /*
    ip_mat *t;
    t = ip_mat_create(5, 5, 1, 1.0);
    ip_mat *cp;
    cp = ip_mat_copy(t);

    ip_mat_init_random(t, 2, 3);
    ip_mat_show(t);
    ip_mat_show_stats(t);
    ip_mat_show(cp);
    ip_mat_show_stats(cp);

    ip_mat_free(t);
    ip_mat_free(cp);
    */

    /* subset test */

    /*
    ip_mat *t;
    t = ip_mat_create(4, 4, 2, 1.0);
    ip_mat *sub;

    ip_mat_init_random(t, 0, 1);
    ip_mat_show(t);

    sub = ip_mat_subset(t, 1, 2, 1, 2);
    ip_mat_show(sub);
    ip_mat_show_stats(sub);

    ip_mat_free(t);
    ip_mat_free(sub);
    */

    /* chain test */

    /*
    ip_mat *t;
    t = ip_mat_create(2, 2, 1, 1.0);
    ip_mat *cp;

    ip_mat_init_random(t, 0, 1);
    ip_mat_show(t);

    cp = ip_mat_copy(t);

    ip_mat *chain;
    chain = ip_mat_concat(t, cp, 2);
    ip_mat_show(chain);
    ip_mat_show_stats(chain);

    ip_mat_free(t);
    ip_mat_free(cp);
    ip_mat_free(chain);
    */

    /* sum/sub/mean test */
    /*

    ip_mat *t1;
    ip_mat *t2;
    ip_mat *sum;

    t1 = ip_mat_create(3, 3, 1, 1.1);
    t2 = ip_mat_create(3, 2, 2, 2.45);
    sum = ip_mat_mean(t1, t2);

    ip_mat_show(t1);
    ip_mat_show(t2);
    ip_mat_show(sum);
    ip_mat_show_stats(sum);

    ip_mat_free(t1);
    ip_mat_free(t2);
    ip_mat_free(sum);
    */

    /* scalar test */
    /*

    ip_mat *t1;
    ip_mat *sum;

    t1 = ip_mat_create(3, 3, 2, 1.1);
    sum = ip_mat_add_scalar(t1, 2.11);

    ip_mat_show(t1);
    ip_mat_show_stats(t1);
    ip_mat_show(sum);
    ip_mat_show_stats(sum);

    ip_mat_free(t1);
    ip_mat_free(sum);
    */

    /* crop test */
    /*
    
    ip_mat *t;
    
    t = ip_mat_create(3, 3, 2, 0.0);
    ip_mat_init_random(t, 0, 300);

    ip_mat_show(t);
    
    ip_mat_crop(t);

    ip_mat_show(t);
    */

    /* gray scale test */
    /*

    Bitmap *b = NULL;
    Bitmap *b_out = NULL;
    ip_mat *in_img = NULL;
    ip_mat *out_img = NULL;


    b = bm_load("caf.bmp");
    in_img = bitmap_to_ip_mat(b);
    

    bm_free(b);

    out_img = ip_mat_to_gray_scale(in_img);
    ip_mat_free(in_img);

    b_out = ip_mat_to_bitmap(out_img);
    ip_mat_free(out_img);

    bm_save(b_out, "gray_scale_output.bmp");
    bm_free(b_out);
    */

    /* blend test */
    /*

    Bitmap *b1 = NULL;
    Bitmap *b2 = NULL;
    Bitmap *b_out = NULL;
    ip_mat *img1 = NULL;
    ip_mat *img2 = NULL;
    ip_mat *out_img = NULL;


    b1 = bm_load("caf.bmp");
    b2 = bm_load("fullmoon.bmp");

    img1 = bitmap_to_ip_mat(b1);
    img2 = bitmap_to_ip_mat(b2);
    
    bm_free(b1);
    bm_free(b2);

    out_img = ip_mat_blend(img1, img2, 0.5);
    ip_mat_free(img1);
    ip_mat_free(img2);

    b_out = ip_mat_to_bitmap(out_img);
    ip_mat_free(out_img);

    bm_save(b_out, "Blend_output.bmp");
    bm_free(b_out);
    */

    /* brighten test */
    /*

    Bitmap *b = NULL;
    Bitmap *b_out = NULL;
    ip_mat *in_img = NULL;
    ip_mat *out_img = NULL;

    b = bm_load("caf.bmp");
    in_img = bitmap_to_ip_mat(b);

    bm_free(b);

    out_img = ip_mat_brighten(in_img, 100);
    ip_mat_free(in_img);
    ip_mat_clamp(out_img);

    b_out = ip_mat_to_bitmap(out_img);
    ip_mat_free(out_img);

    bm_save(b_out, "Brighten_output.bmp");
    bm_free(b_out);
    */

    /* gaussian noise test */

    Bitmap *b = NULL;
    Bitmap *b_out = NULL;
    ip_mat *in_img = NULL;
    ip_mat *out_img = NULL;

    b = bm_load("flower2.bmp");
    in_img = bitmap_to_ip_mat(b);

    bm_free(b);

    out_img = ip_mat_corrupt(in_img, 100);
    ip_mat_free(in_img);

    b_out = ip_mat_to_bitmap(out_img);
    ip_mat_free(out_img);

    bm_save(b_out, "Noise_output.bmp");
    bm_free(b_out);

    /* padding test */
    /*

    ip_mat *a;
    ip_mat *b;

    a = ip_mat_create(3, 3, 2, 1.0);
    b = ip_mat_padding(a, 1, 2);

    ip_mat_show(a);
    ip_mat_show(b);

    ip_mat_free(a);
    ip_mat_free(b);
    */

    /* convolve test */
    /*

    ip_mat *a;
    ip_mat *b;
    ip_mat *f;

    a = ip_mat_create(5, 5, 2, 1.0);
    f = ip_mat_create(3, 3, 2, 2.0);
    set_val(a, 2, 2, 1, 2.0);

    b = ip_mat_convolve(a, f);

    ip_mat_show(a);
    ip_mat_show(b);

    ip_mat_free(a);
    ip_mat_free(b);
    ip_mat_free(f);
    */

    /* sharpen/edge/emboss test */
    /*

    Bitmap *b = NULL;
    Bitmap *b_out = NULL;
    ip_mat *in_img = NULL;
    ip_mat *out_img = NULL;
    ip_mat *filter = NULL;

    b = bm_load("flower.bmp");
    in_img = bitmap_to_ip_mat(b);
    filter = create_emboss_filter();

    ip_mat_show(filter);

    bm_free(b);

    out_img = ip_mat_convolve(in_img, filter);
    ip_mat_free(in_img);
    ip_mat_clamp(out_img);

    b_out = ip_mat_to_bitmap(out_img);
    ip_mat_free(out_img);

    bm_save(b_out, "Emboss_output.bmp");
    bm_free(b_out);
    */

    /* average/gaussian test */
    /*

    Bitmap *b = NULL;
    Bitmap *b_out = NULL;
    ip_mat *in_img = NULL;
    ip_mat *out_img = NULL;
    ip_mat *filter = NULL;

    b = bm_load("flower.bmp");
    in_img = bitmap_to_ip_mat(b);
    filter = create_gaussian_filter(9, 9, in_img -> k, 5);

    ip_mat_show(filter);

    bm_free(b);

    out_img = ip_mat_convolve(in_img, filter);
    ip_mat_free(in_img);
    clamp(out_img, 0, 255);

    b_out = ip_mat_to_bitmap(out_img);
    ip_mat_free(out_img);

    bm_save(b_out, "Gaussian_output.bmp");
    bm_free(b_out);
    */

    /* clamp/rescale test */
    /*

    ip_mat *a;

    a = ip_mat_create(3, 3, 3, 0.0);
    ip_mat_init_random(a, 0, 10);
    ip_mat_show(a);
    ip_mat_show_stats(a);

    clamp(a, 0, 5);

    ip_mat_show(a);
    ip_mat_show_stats(a);

    ip_mat_free(a);
    */

    return 0;
}
