#ifndef __secp256k1_H__
#define __secp256k1_H__

#include "datatype.h"

void secp256k1_ge_from_gej(ge_secp256k1 *, gej_secp256k1 *);
void secp256k1_gej_from_ge(gej_secp256k1 *, ge_secp256k1 *);

void secp256k1add(ge_secp256k1 *, const ge_secp256k1 *, const ge_secp256k1 *);
void secp256k1double(gej_secp256k1 *, const ge_secp256k1 *);

void secp256k1doublebernstein(gej_secp256k1 *, const ge_secp256k1 *);

void secp256k1scalermult(ge_secp256k1 *, const gfe_p25632977 *, const ge_secp256k1 *);

//                   sigr           , sigs           , privkey              , message              , nonce
void secp256k1sign(gfe_p25632977 *, gfe_p25632977 *, const gfe_p25632977 *, const gfe_p25632977 *, const gfe_p25632977 *);
//                     sigr           , sigs           , pubkey               , message
void secp256k1verify(gfe_p25632977 *, gfe_p25632977 *, const ge_secp256k1 *, const gfe_p25632977 *);

#endif