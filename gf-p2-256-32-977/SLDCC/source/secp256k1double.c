#include "secp256k1.h"
#include "p25632977.h"

void secp256k1double(gej_secp256k1 *pp, const ge_secp256k1 *p) {
    gfe_p25632977 l1, l2, l3, x2, y2, y4, l1_2, xy2, ll2, nll2, nx, nl3, l2mx, yp;

    /* l1 = 3 * x2 */
    gfp25632977sqr(&x2, &p->x);
    gfp25632977mul(&l1, &(gfe_p25632977){3,0,0,0}, &x2);

    /* l2 = 4 * x1 * y2 */
    gfp25632977sqr(&y2, &p->y);
    gfp25632977mul(&xy2, &p->x, &y2);
    gfp25632977mul(&l2, &(gfe_p25632977){4,0,0,0}, &xy2);

    /* l3 = 8 * y4 */
    gfp25632977sqr(&y4, &y2);
    gfp25632977mul(&l3, &(gfe_p25632977){8,0,0,0}, &y4);

    /* z = 2 * y1 */
    gfp25632977mul(&pp->z, &(gfe_p25632977){2,0,0,0}, &p->y);

    /* x = l1_2 - 2 * l2 */
    gfp25632977mul(&ll2, &(gfe_p25632977){2,0,0,0}, &l2);
    gfp25632977negate(&nll2, &ll2);
    gfp25632977sqr(&l1_2, &l1);
    gfp25632977add(&pp->x, &l1_2, &nll2);

    /* y = l1 * (l2 - x3) - l3 */
    gfp25632977negate(&nl3, &l3);
    gfp25632977negate(&nx, &pp->x);
    gfp25632977add(&l2mx, &l2, &nx);
    gfp25632977mul(&yp, &l1, &l2mx);
    gfp25632977add(&pp->y, &yp, &nl3);

    if(pp->z.l[0] == 0 && pp->z.l[1] == 0 && pp->z.l[2] == 0 && pp->z.l[3] == 0) {
        pp->infinity = 1;
    } else {
        pp->infinity = 0;
    }
}

void secp256k1doublesub(gej_secp256k1 *pp, const ge_secp256k1 *p) {
    gfe_p25632977 l1, l2, l3, x2, y2, y4, l1_2, xy2, ll2, nll2, nx, nl3, l2mx, yp;

    /* l1 = 3 * x2 */
    gfp25632977sqr(&x2, &p->x);
    gfp25632977mul(&l1, &(gfe_p25632977){3,0,0,0}, &x2);

    /* l2 = 4 * x1 * y2 */
    gfp25632977sqr(&y2, &p->y);
    gfp25632977mul(&xy2, &p->x, &y2);
    gfp25632977mul(&l2, &(gfe_p25632977){4,0,0,0}, &xy2);

    /* l3 = 8 * y4 */
    gfp25632977sqr(&y4, &y2);
    gfp25632977mul(&l3, &(gfe_p25632977){8,0,0,0}, &y4);

    /* z = 2 * y1 */
    gfp25632977mul(&pp->z, &(gfe_p25632977){2,0,0,0}, &p->y);

    /* x = l1_2 - 2 * l2 */
    gfp25632977mul(&ll2, &(gfe_p25632977){2,0,0,0}, &l2);
    // gfp25632977negate(&nll2, &ll2);
    gfp25632977sqr(&l1_2, &l1);
    // gfp25632977add(&pp->x, &l1_2, &nll2);
    gfp25632977sub(&pp->x, &l1_2, &ll2);

    /* y = l1 * (l2 - x3) - l3 */
    // gfp25632977negate(&nl3, &l3);
    // gfp25632977negate(&nx, &pp->x);
    // gfp25632977add(&l2mx, &l2, &nx);
    gfp25632977sub(&l2mx, &l2, &pp->x);
    gfp25632977mul(&yp, &l1, &l2mx);
    // gfp25632977add(&pp->y, &yp, &nl3);
    gfp25632977add(&pp->y, &yp, &l3);

    if(pp->z.l[0] == 0 && pp->z.l[1] == 0 && pp->z.l[2] == 0 && pp->z.l[3] == 0) {
        pp->infinity = 1;
    } else {
        pp->infinity = 0;
    }
}