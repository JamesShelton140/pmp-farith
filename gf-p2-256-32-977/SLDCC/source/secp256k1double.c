#include "secp256k1.h"
#include "p25632977.h"

void secp256k1double(gej_secp256k1 *pp, const ge_secp256k1 *p) {
    gfe_p25632977 l1, l2, l3, x2, y2, y4, l1_2, xy2, ll2, l2mx, yp;

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
    gfp25632977sqr(&l1_2, &l1);
    gfp25632977sub(&pp->x, &l1_2, &ll2);

    /* y = l1 * (l2 - x3) - l3 */
    gfp25632977sub(&l2mx, &l2, &pp->x);
    gfp25632977mul(&yp, &l1, &l2mx);
    gfp25632977sub(&pp->y, &yp, &l3);

    if(pp->z.l[0] == 0 && pp->z.l[1] == 0 && pp->z.l[2] == 0 && pp->z.l[3] == 0) {
        pp->infinity = 1;
    } else {
        pp->infinity = 0;
    }
}

void secp256k1doublejacobian(gej_secp256k1 *pp, const gej_secp256k1 *p) {
    gfe_p25632977 l1, l2, l3, x2, y2, y4, l1_2, xy2, ll2, l2mx, yp, yz;

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

    /* z = 2 * y1 * z1 */
    gfp25632977mul(&yz, &p->y, &p->z);
    gfp25632977mul(&pp->z, &(gfe_p25632977){2,0,0,0}, &yz);

    /* x = l1_2 - 2 * l2 */
    gfp25632977mul(&ll2, &(gfe_p25632977){2,0,0,0}, &l2);
    gfp25632977sqr(&l1_2, &l1);
    gfp25632977sub(&pp->x, &l1_2, &ll2);

    /* y = l1 * (l2 - x3) - l3 */
    gfp25632977sub(&l2mx, &l2, &pp->x);
    gfp25632977mul(&yp, &l1, &l2mx);
    gfp25632977sub(&pp->y, &yp, &l3);

    if(pp->z.l[0] == 0 && pp->z.l[1] == 0 && pp->z.l[2] == 0 && pp->z.l[3] == 0) {
        pp->infinity = 1;
    } else {
        pp->infinity = 0;
    }
}


void secp256k1doublebernstein(gej_secp256k1 *pp, const ge_secp256k1 *p) {
    gfe_p25632977 xx, yy, yyyy, t0, t1, t2, t3, S, t4, M, t5, t6, T, t7, t8, t9;

    gfp25632977sqr(&xx, &p->x);
    gfp25632977sqr(&yy, &p->y);
    gfp25632977sqr(&yyyy, &yy);
    gfp25632977add(&t0, &p->x, &yy);
    gfp25632977sqr(&t1, &t0);
    gfp25632977sub(&t2, &t1, &xx);
    gfp25632977sub(&t3, &t2, &yyyy);
    gfp25632977mul(&S, &(gfe_p25632977){2,0,0,0}, &t3);
    gfp25632977mul(&t4, &(gfe_p25632977){3,0,0,0}, &xx);
    M = t4;
    gfp25632977sqr(&t5, &M);
    gfp25632977mul(&t6, &(gfe_p25632977){2,0,0,0}, &S);
    gfp25632977sub(&T, &t5, &t6);
    pp->x = T;
    gfp25632977sub(&t7, &S, &T);
    gfp25632977mul(&t8, &(gfe_p25632977){8,0,0,0}, &yyyy);
    gfp25632977mul(&t9, &M, &t7);
    gfp25632977sub(&pp->y, &t9, &t8);
    gfp25632977mul(&pp->z, &(gfe_p25632977){2,0,0,0}, &p->y);

    if(pp->z.l[0] == 0 && pp->z.l[1] == 0 && pp->z.l[2] == 0 && pp->z.l[3] == 0) {
        pp->infinity = 1;
    } else {
        pp->infinity = 0;
    }
}

// To be removed

// void secp256k1double(gej_secp256k1 *pp, const ge_secp256k1 *p) {
//     gfe_p25632977 l1, l2, l3, x2, y2, y4, l1_2, xy2, ll2, nll2, nx, nl3, l2mx, yp;

//     /* l1 = 3 * x2 */
//     gfp25632977sqr(&x2, &p->x);
//     gfp25632977mul(&l1, &(gfe_p25632977){3,0,0,0}, &x2);

//     /* l2 = 4 * x1 * y2 */
//     gfp25632977sqr(&y2, &p->y);
//     gfp25632977mul(&xy2, &p->x, &y2);
//     gfp25632977mul(&l2, &(gfe_p25632977){4,0,0,0}, &xy2);

//     /* l3 = 8 * y4 */
//     gfp25632977sqr(&y4, &y2);
//     gfp25632977mul(&l3, &(gfe_p25632977){8,0,0,0}, &y4);

//     /* z = 2 * y1 */
//     gfp25632977mul(&pp->z, &(gfe_p25632977){2,0,0,0}, &p->y);

//     /* x = l1_2 - 2 * l2 */
//     gfp25632977mul(&ll2, &(gfe_p25632977){2,0,0,0}, &l2);
//     gfp25632977negate(&nll2, &ll2);
//     gfp25632977sqr(&l1_2, &l1);
//     gfp25632977add(&pp->x, &l1_2, &nll2);

//     /* y = l1 * (l2 - x3) - l3 */
//     gfp25632977negate(&nl3, &l3);
//     gfp25632977negate(&nx, &pp->x);
//     gfp25632977add(&l2mx, &l2, &nx);
//     gfp25632977mul(&yp, &l1, &l2mx);
//     gfp25632977add(&pp->y, &yp, &nl3);

//     if(pp->z.l[0] == 0 && pp->z.l[1] == 0 && pp->z.l[2] == 0 && pp->z.l[3] == 0) {
//         pp->infinity = 1;
//     } else {
//         pp->infinity = 0;
//     }
// }