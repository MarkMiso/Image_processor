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

    t1 = ip_mat_create(3, 3, 2, 1.1);
    t2 = ip_mat_create(3, 3, 2, 2.45);
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
    
    return 0;
    
}
