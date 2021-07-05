#include "secp256k1.h"
#include "p25632977.h"

void secp256k1add(gej_secp256k1 *pq, const ge_secp256k1 *p, const ge_secp256k1 *q) {
    gfe_p25632977 l3, l6, l7, l8, l9, ll6, ll3, l7ll3, twox3, lll3, l8lll3, l9l6, lastl;

    // l3 = xp - xq
    gfp25632977sub(&l3, &p->x, &q->x);

    // l6 = l4 - l5
    gfp25632977sub(&l6, &p->y, &q->y);

    // l7 = xp + xq
    gfp25632977add(&l7, &p->x, &q->x);

    // l8 = l4 + l5
    gfp25632977add(&l8, &p->y, &q->y);

    // z3 = l3
    pq->z = l3;

    // x3 = l6^2 - l7 * l3^2
    gfp25632977sqr(&ll6, &l6);
    gfp25632977sqr(&ll3, &l3);
    gfp25632977mul(&l7ll3, &l7, &ll3);
    gfp25632977sub(&pq->x, &ll6, &l7ll3);

    // l9 = l7 * l3^2 - 2 * x3
    gfp25632977add(&twox3, &pq->x, &pq->x);
    gfp25632977sub(&l9, &l7ll3, &twox3);

    // y3 = (l9 * l6 - l8 * l3^3) * 2^-1
    gfp25632977mul(&l9l6, &l9, &l6);
    gfp25632977mul(&lll3, &ll3, &l3);
    gfp25632977mul(&l8lll3, &l8, &lll3);
    gfp25632977sub(&lastl, &l9l6, &l8lll3);
    gfp25632977mul(&pq->y, &lastl, &twoinv);

    if(pq->z.l[0] == 0 && pq->z.l[1] == 0 && pq->z.l[2] == 0 && pq->z.l[3] == 0) {
        pq->infinity = 1;
    } else {
        pq->infinity = 0;
    }
}

void secp256k1addjacobian(gej_secp256k1 *pq, const gej_secp256k1 *p, const gej_secp256k1 *q) {
    gfe_p25632977 l1, l2, l3, l4, l5, l6, l7, l8, l9, zzq, zzp, zzzq, zzzp, zpzq, ll6, ll3, l7ll3, twox3, lll3, l8lll3, l9l6, lastl;

    if(p->infinity == 1) {
        pq->x = q->x;
        pq->y = q->y;
        pq->z = q->z;
        pq->infinity = q->infinity;
        return;
    }

    if(q->infinity == 1) {
        pq->x = p->x;
        pq->y = p->y;
        pq->z = p->z;
        pq->infinity = p->infinity;
        return;
    }

    // l1 = xp * zq^2
    gfp25632977sqr(&zzq, &q->z);
    gfp25632977mul(&l1, &p->x, &zzq);

    // l2 = xq * zp^2
    gfp25632977sqr(&zzp, &p->z);
    gfp25632977mul(&l2, &q->x, &zzp);

    // l3 = l1 - l2
    gfp25632977sub(&l3, &l1, &l2);

    // l4 = yp * zq^3
    gfp25632977mul(&zzzq, &zzq, &q->z);
    gfp25632977mul(&l4, &p->y, &zzzq);

    // l5 = yq * zp^3
    gfp25632977mul(&zzzp, &zzp, &p->z);
    gfp25632977mul(&l5, &q->y, &zzzp);

    // l6 = l4 - l5
    gfp25632977sub(&l6, &l4, &l5);

    // l7 = l1 + l2
    gfp25632977add(&l7, &l1, &l2);

    // l8 = l4 + l5
    gfp25632977add(&l8, &l4, &l5);

    // z3 = zp * zq * l3
    gfp25632977mul(&zpzq, &p->z, &q->z);
    gfp25632977mul(&pq->z, &zpzq, &l3);

    // x3 = l6^2 - l7 * l3^2
    gfp25632977sqr(&ll6, &l6);
    gfp25632977sqr(&ll3, &l3);
    gfp25632977mul(&l7ll3, &l7, &ll3);
    gfp25632977sub(&pq->x, &ll6, &l7ll3);

    // l9 = l7 * l3^2 - 2 * x3
    gfp25632977add(&twox3, &pq->x, &pq->x);
    gfp25632977sub(&l9, &l7ll3, &twox3);

    // y3 = (l9 * l6 - l8 * l3^3) * 2^-1
    gfp25632977mul(&l9l6, &l9, &l6);
    gfp25632977mul(&lll3, &ll3, &l3);
    gfp25632977mul(&l8lll3, &l8, &lll3);
    gfp25632977sub(&lastl, &l9l6, &l8lll3);
    gfp25632977mul(&pq->y, &lastl, &twoinv);

    if(pq->z.l[0] == 0 && pq->z.l[1] == 0 && pq->z.l[2] == 0 && pq->z.l[3] == 0) {
        pq->infinity = 1;
    } else {
        pq->infinity = 0;
    }
}