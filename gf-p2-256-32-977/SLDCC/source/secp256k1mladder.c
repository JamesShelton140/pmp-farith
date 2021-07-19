#include "secp256k1.h"
#include "p25632977.h"
#include <stdio.h>

void print_felem(const gfe_p25632977 *);
void rewriten(gfe_p25632977 *, const gfe_p25632977 *);

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

// void secp256k1scalermult(gej_secp256k1 *nP, const gfe_p25632977 *n, const ge_secp256k1 *P) {
//     gej_secp256k1 R0, R1, R_temp;
//     gfe_p25632977 rn;

//     rewriten(&rn, n);

// 	// Set R0, R1 = P
//     R0 = (gej_secp256k1){P->x,P->y,1,P->infinity};
//     R1 = (gej_secp256k1){P->x,P->y,1,P->infinity};
    
//     int i, bit, limb;
//     uint64 mask, swap;
//     for (i = 1; i <= 255; i++) {
//     	// select bit at position i
//         limb = i/64;
//         bit = i%64;
//         mask = (uint64)1 << bit;
//         swap = mask & rn.l[limb];
        
//         // gfp25632977readbit(&bit, n, limb);
        
//         if (swap == 0) {
//         	// R1 <- 2R1 + R0
//             secp256k1doublejacobian(&R_temp, &R1);
//             secp256k1addjacobian(&R1, &R_temp, &R0);
//         } else {
//         	// R0 <- 2R0 + R1
//             secp256k1doublejacobian(&R_temp, &R0);
//             secp256k1addjacobian(&R0, &R_temp, &R1);
//         }

//         // if(swap != 0) {
//         //     printf("\n\n\ni: %u\n\n",i);
//         //     printf("limb: %u\n\n",limb);
//         //     printf("bit: %u\n\n",bit);
//         //     printf("l[limb]: %16llX\n\n",rn.l[limb]);
//         //     printf("mask: %16llX\n\n",mask);
//         //     printf("Swap: %16llX\n\n",swap);
//         //     printf("R0.x: "); print_felem(&R0.x);
//         //     printf("R1.x: "); print_felem(&R1.x);
//         // }
        
//     }
    
//     // Set nP = R0
//     nP->x = R0.x;
//     nP->y = R0.y;
//     nP->z = R0.z;
//     nP->infinity = R0.infinity;

// }

void rewriten(gfe_p25632977 *rn, const gfe_p25632977 *n) {
    gfe_p25632977 n1, n2, n3;

    // rn = 2 * ((n-1)/2 mod q) + 1
    gfp25632977sub(&n1, n, &(gfe_p25632977){1,0,0,0});
    printf("n1:\t\t");print_felem(&n1);

    gfp25632977mul(&n2, &n1, &twoinv);
    printf("n2:\t\t");print_felem(&n2);

    gfp25632977makeunique(&n2);
    printf("n2:\t\t");print_felem(&n2);

    gfp25632977add(&n3, &n2, &n2);
    printf("n3:\t\t");print_felem(&n3);

    gfp25632977add(rn, &n3, &(gfe_p25632977){1,0,0,0});
    printf("rn:\t\t");print_felem(rn);
}

/*
* Elliptic curve scalar multiplication using basic Joye ladder
* (Spelling of name will be corrected)
*/

void secp256k1scalermult(gej_secp256k1 *nP, const gfe_p25632977 *n, const ge_secp256k1 *P) {
    gej_secp256k1 R0, R1, R_temp;

	// Set R0 = 0, R1 = P
    R0 = (gej_secp256k1){(uint64)0,(uint64)1,(uint64)0,(uint64)1};
    R1 = (gej_secp256k1){P->x,P->y,1,P->infinity};
    
    int i, bit, limb;
    uint64 mask, swap;
    for (i = 0; i <= 255; i++) {
    	// select bit at position i
        limb = i/64;
        bit = i%64;
        mask = (uint64)1 << bit;
        swap = mask & n->l[limb];
        
        // gfp25632977readbit(&bit, n, limb);
        
        if (swap == 0) {
        	// R1 <- 2R1 + R0
            if(i == 0 || i == 1) {
                printf("\n\n\ni (swap 0): %u\n\n",i);
                printf("R0.x: "); print_felem(&R0.x);
                printf("R1.x: "); print_felem(&R1.x);
            }
            secp256k1doublejacobian(&R_temp, &R1);
            secp256k1addjacobian(&R1, &R_temp, &R0);
        } else {
        	// R0 <- 2R0 + R1
            printf("\n\n\ni: %u\n\n",i);
            printf("R0.x: "); print_felem(&R0.x);
            printf("R1.x: "); print_felem(&R1.x);
            secp256k1doublejacobian(&R_temp, &R0);
            secp256k1addjacobian(&R0, &R_temp, &R1);
        }

        if(i == 0 || i == 1) {
            printf("\n\n\ni: %u\n\n",i);
            printf("R0.x: "); print_felem(&R0.x);
            printf("R1.x: "); print_felem(&R1.x);
        }
        if(swap != 0) {
            printf("\n\n\ni: %u\n\n",i);
            printf("limb: %u\n\n",limb);
            printf("bit: %u\n\n",bit);
            printf("l[limb]: %16llX\n\n",n->l[limb]);
            printf("mask: %16llX\n\n",mask);
            printf("Swap: %16llX\n\n",swap);
            printf("R0.x: "); print_felem(&R0.x);
            printf("R1.x: "); print_felem(&R1.x);
        }
        
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
