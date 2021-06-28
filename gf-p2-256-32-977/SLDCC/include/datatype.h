/*
+-----------------------------------------------------------------------------+
| This code corresponds to the Galois field F(2^256-2^32-977) from the paper  |
| "Efficient Arithmetic In (Pseudo-)Mersenne Prime Order Fields" authored by  |
| Kaushik Nath,  Indian Statistical Institute, Kolkata, India, and            |
| Palash Sarkar, Indian Statistical Institute, Kolkata, India.	              |
+-----------------------------------------------------------------------------+
| Copyright (c) 2018, Kaushik Nath and Palash Sarkar.                         |
|                                                                             |
| Permission to use this code is granted.                          	      |
|                                                                             |
| Redistribution and use in source and binary forms, with or without          |
| modification, are permitted provided that the following conditions are      |
| met:                                                                        |
|                                                                             |
| * Redistributions of source code must retain the above copyright notice,    |
|   this list of conditions and the following disclaimer.                     |
|                                                                             |
| * Redistributions in binary form must reproduce the above copyright         |
|   notice, this list of conditions and the following disclaimer in the       |
|   documentation and/or other materials provided with the distribution.      |
|                                                                             |
| * The names of the contributors may not be used to endorse or promote       |
|   products derived from this software without specific prior written        |
|   permission.                                                               |
+-----------------------------------------------------------------------------+
| THIS SOFTWARE IS PROVIDED BY THE AUTHORS ""AS IS"" AND ANY EXPRESS OR       |
| IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES   |
| OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.     |
| IN NO EVENT SHALL THE CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,      |
| INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT    |
| NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,   |
| DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY       |
| THEORY LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING |
| NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,| 
| EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                          |
+-----------------------------------------------------------------------------+
*/

#ifndef __DATATYPE_p25632977_H__
#define __DATATYPE_p25632977_H__

#define NLIMBS 4

typedef unsigned char uchar8;
typedef unsigned long long uint64;

typedef struct {
  uint64 l[NLIMBS]; 
}
gfe_p25632977;

//Group element of secp256k1 in affine coordinates (bitcoin-core/secp256k1)
typedef struct {
  gfe_p25632977 x;
  gfe_p25632977 y;
  int infinity; //indicates this represents the point at infinity
}
ge_secp256k1;

//Group element of secp256k1 in jacobian coordinates (bitcoin-core/secp256k1)
typedef struct {
    gfe_p25632977 x; /* actual X: x/z^2 */
    gfe_p25632977 y; /* actual Y: y/z^3 */
    gfe_p25632977 z;
    int infinity; //indicates this represents the point at infinity
}
gej_secp256k1;

#endif

