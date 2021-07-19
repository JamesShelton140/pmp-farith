#include "secp256k1.h"
#include "p25632977.h"
#include <stdio.h>

void print_felem(const gfe_p25632977 *);

/*
* Elliptic curve scalar multiplication using basic Montgomery ladder
* (Spelling of name will be corrected)
*/

// void secp256k1scalermult(gej_secp256k1 *nP, const gfe_p25632977 *n, const ge_secp256k1 *P) {
//     gej_secp256k1 R0, R1;

		// Set R0 = 0/inf, R1 = P
//     R0 = (gej_secp256k1){0,1,0,1};
//     R1.x = P->x;
//     R1.y = P->y;
//     R1.infinity = P->infinity;
//     if (R1.infinity) {
//         R1.z = (gfe_p25632977){0,0,0,0};
//     } else {
//         R1.z = (gfe_p25632977){1,0,0,0};
//     }
    
//     int i, bit, limb;
//     uint64 mask, swap;
//     for (i = 255; i >= 0; i--) {
			// select bit at position i
//         limb = i/64;
//         bit = i%64;
//         mask = 1 << bit;
//         swap = mask & n->l[limb];
        
//         // gfp25632977readbit(&bit, n, limb);

//         if (swap == 0) {
				// R1 <- R0 + R1
				// R0 <- 2R0
//             secp256k1addjacobian(&R1, &R0, &R1);
//             secp256k1doublejacobian(&R0, &R0);
//         } else {
				// R0 = R0 + R1
				// R1 = 2R1
//             secp256k1addjacobian(&R0, &R0, &R1);
//             secp256k1doublejacobian(&R1, &R1);
//         }
//     }
    
        // Set nP = R0
//     nP->x = R0.x;
//     nP->y = R0.y;
//     nP->z = R0.z;
//     nP->infinity = R0.infinity;

// }

/*
* Elliptic curve scalar multiplication using basic Joye ladder
* (Spelling of name will be corrected)
*/

void secp256k1scalermult(gej_secp256k1 *nP, const gfe_p25632977 *n, const ge_secp256k1 *P) {
    gej_secp256k1 R0, R1, R_temp;

	// Set R0, R1 = P
    R0 = (gej_secp256k1){P->x,P->y,1,P->infinity};
    R1 = (gej_secp256k1){P->x,P->y,1,P->infinity};
    printf("R0.x: "); print_felem(&R0.x);
    printf("R1.x: "); print_felem(&R1.x);
    
    printf("\n\n\nstarting loop\n\n\n");
    int i, bit, limb;
    uint64 mask, swap;
    for (i = 1; i <= 255; i++) {
    	// select bit at position i
        limb = i/64;
        bit = i%64;
        mask = 1 << bit;
        swap = mask & n->l[limb];
        
        // gfp25632977readbit(&bit, n, limb);
        printf("\n\n\nSwap: %u\n\n",swap);
        if (swap == 0) {
        	// R1 <- 2R1 + R0
            secp256k1doublejacobian(&R_temp, &R1);
            secp256k1addjacobian(&R1, &R_temp, &R0);
        } else {
        	// R0 <- 2R0 + R1
            secp256k1doublejacobian(&R_temp, &R0);
            secp256k1addjacobian(&R0, &R_temp, &R1);
        }
        printf("R0.x: "); print_felem(&R0.x);
        printf("R1.x: "); print_felem(&R1.x);
    }
    
    // Set nP = R0
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
