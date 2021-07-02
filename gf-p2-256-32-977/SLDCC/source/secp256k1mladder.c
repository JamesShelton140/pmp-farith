#include "secp256k1.h"
#include "p25632977.h"

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
    

    int bit, limb;
    for (int i = 255; i >= 0; i--) {
        limb = i/64;
        bit = i%64;
        gfp25632977readbit(&bit, &n, &limb);

        if (bit) {
            secp256k1addjacobian(&R1, &R0, &R1);
            secp256k1doublejacobian(&R0, &R0);
        } else {
            secp256k1addjacobian(&R0, &R0, &R1);
            secp256k1doublejacobian(&R1, &R1);
        }
    }
    
    nP = &R0;

}