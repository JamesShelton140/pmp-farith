#include "secp256k1.h"
#include "p25632977.h"

void secp256k1_ge_from_gej(ge_secp256k1 *p, gej_secp256k1 *pj) {
    gfe_p25632977 z1, z2, z3;

    if (pj->infinity) {
        pj->x = (gfe_p25632977) {0,0,0,0};
        pj->y = (gfe_p25632977) {1,0,0,0};
        pj->z = (gfe_p25632977) {0,0,0,0};
        pj->infinity = 1;
        return;
    }

    gfp25632977inv(&z1, &pj->z); //z1 == z^-1
    gfp25632977sqr(&z2, &z1); //z2 == z1^2 == z^-2
    gfp25632977mul(&z3, &z2, &z1); //z3 == z1*z2 == z^-3

    gfp25632977mul(&p->x, &pj->x, &z2); //p.x == pj.x * z2 == x/z^2
    gfp25632977mul(&p->y, &pj->y, &z3); //p.y == pj.y * z3 == y/z^3
    p->infinity = 0;
}

void secp256k1_gej_from_ge(gej_secp256k1 *pj, ge_secp256k1 *p) {
    pj->x = p->x;
    pj->y = p->y;
    pj->z = (gfe_p25632977) {1,0,0,0};
    pj->infinity = p->infinity;
}