#include "secp256k1.h"
#include "p25632977.h"
#include <stdio.h>

void print_felem(const gfe_p25632977 *);

// void secp256k1scalermult(gej_secp256k1 *nP, const gfe_p25632977 *n, const ge_secp256k1 *P) {
//     gej_secp256k1 R0, R1;

//     R0 = (gej_secp256k1){0,1,0,1};
//     R1.x = P->x;
//     R1.y = P->y;
//     R1.infinity = P->infinity;
//     if (R1.infinity) {
//         R1.z = (gfe_p25632977){0,0,0,0};
//     } else {
//         R1.z = (gfe_p25632977){1,0,0,0};
//     }
    
//     printf("\n\n\nstarting loop\n\n\n");
//     int i, bit, limb;
//     uint64 mask, swap;
//     for (i = 255; i >= 0; i--) {
//         limb = i/64;
//         bit = i%64;
//         mask = 1 << bit;
//         swap = mask & n->l[limb];
        
//         // gfp25632977readbit(&bit, n, limb);

//         if (swap = 0) {
//             secp256k1addjacobian(&R1, &R0, &R1);
//             secp256k1doublejacobian(&R0, &R0);
//         } else {
//             secp256k1addjacobian(&R0, &R0, &R1);
//             secp256k1doublejacobian(&R1, &R1);
//         }
//     }
    
//     nP->x = R0.x;
//     nP->y = R0.y;
//     nP->z = R0.z;
//     nP->infinity = R0.infinity;

// }

// Implement basic version of Joye ladder
void secp256k1scalermult(gej_secp256k1 *nP, const gfe_p25632977 *n, const ge_secp256k1 *P) {
    gej_secp256k1 R0, R1, R_temp;

    R0 = (gej_secp256k1){P->x,P->y,1,P->infinity};
    R1 = (gej_secp256k1){P->x,P->y,1,P->infinity};
    // R1.x = P->x;
    // R1.y = P->y;
    // R1.infinity = P->infinity;
    // if (R1.infinity) {
    //     R1.z = (gfe_p25632977){0,0,0,0};
    // } else {
    //     R1.z = (gfe_p25632977){1,0,0,0};
    // }
    
    printf("\n\n\nstarting loop\n\n\n");
    int i, bit, limb;
    uint64 mask, swap;
    for (i = 1; i <= 255; i++) {
        limb = i/64;
        bit = i%64;
        mask = 1 << bit;
        swap = mask & n->l[limb];printf("%16llX ",swap);
        
        // gfp25632977readbit(&bit, n, limb);

        if (swap = 0) {
            secp256k1doublejacobian(&R_temp, &R1);
            secp256k1addjacobian(&R1, &R_temp, &R0);
        } else {
            secp256k1doublejacobian(&R_temp, &R0);
            secp256k1addjacobian(&R0, &R_temp, &R1);
        }
    }
    
    nP->x = R0.x;
    nP->y = R0.y;
    nP->z = R0.z;
    nP->infinity = R0.infinity;

}

void print_felem(const gfe_p25632977 *e){

	uchar8  i;

	for (i=NLIMBS-1; i>0; --i) 
		printf("%16llX ",e->l[i]);
	printf("%16llX \n\n",e->l[0]);
}