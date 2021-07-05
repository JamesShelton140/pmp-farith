#include "secp256k1.h"
#include "p25632977.h"
#include <stdio.h>

void secp256k1scalermult(gej_secp256k1 *nP, const gfe_p25632977 *n, const ge_secp256k1 *P) {
    gej_secp256k1 R0, R1;

    R0 = (gej_secp256k1){0,1,0,1};
    R1.x = P->x;
    R1.y = P->y;
    R1.infinity = P->infinity;
    if (R1.infinity) {
        R1.z = (gfe_p25632977){0,0,0,0};
    } else {
        R1.z = (gfe_p25632977){1,0,0,0};
    }
    
    printf("\n\n\nstarting loop\n\n\n");
    int i, bit, limb;
    uint64 mask, swap;
    for (i = 255; i >= 0; i--) {
        printf("i: %u\n",i);
        limb = i/64;printf("Limb: %u\n",limb);
        bit = i%64;printf("bit: %u\n",limb);
        mask = 1 << bit;printf("Mask: %16llX\n",mask);
        swap = mask & n->l[limb];printf("swap: %16llX\n\n",swap);
        
        // gfp25632977readbit(&bit, n, limb);

        if (swap = 0) {
            secp256k1addjacobian(&R1, &R0, &R1);
            secp256k1doublejacobian(&R0, &R0);
        } else {
            secp256k1addjacobian(&R0, &R0, &R1);
            secp256k1doublejacobian(&R1, &R1);
        }
    }
    
    nP->x = R0.x;
    nP->y = R0.y;
    nP->z = R0.z;
    nP->infinity = R0.infinity;

}