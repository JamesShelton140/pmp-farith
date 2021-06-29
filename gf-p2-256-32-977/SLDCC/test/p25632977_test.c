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

#include<stdio.h>
#include<math.h>
#include"datatype.h"
#include"p25632977.h"
#include"secp256k1.h"
#include"measure.h"

#define change_input(x,y,z)  {x.l[0] = y.l[0]^z.l[0];}
#define FILE stdout
void print_elem(const gfe_p25632977 *);

int main() {
	
	gfe_p25632977 e = {0xFFFFFFFEFFFFFC2E,-1,-1,-1};
	
	//gfe_p25632977 e = {0xFFFFFFFEFFFFFC2D,-1,-1,-1};
	//gfe_p25632977 e = {0xFFFFFFFEFFFFFC2C,-1,-1,-1};
	//gfe_p25632977 e = {0xFFFFFFFEFFFFFC2B,-1,-1,-1};
	//gfe_p25632977 e = {1,0,0,0};
	//gfe_p25632977 e = {0,0,1,0};
	//gfe_p25632977 e = {2,0,2,0};
	//gfe_p25632977 e = {3,2,5,5};

	gfe_p25632977 einv;
	
	fprintf(FILE,"\nThe field element is:\t\t"); print_elem(&e);  
	gfp25632977inv(&einv,&e);
	fprintf(FILE,"The found inverse is:\t\t"); print_elem(&einv);  

	gfe_p25632977 t;
	gfp25632977mul(&t,&e,&einv); gfp25632977makeunique(&t);
	fprintf(FILE,"The cross check value is:\t"); print_elem(&t);

	gfe_p25632977 ne;
	gfp25632977negate(&ne, &e);
	fprintf(FILE,"The found negation is:\t\t"); print_elem(&ne);

	gfe_p25632977 z;
	gfp25632977add(&z, &ne, &e); gfp25632977makeunique(&z);
	fprintf(FILE,"The cross check value is:\t"); print_elem(&z);

	gfp25632977makeunique(&ne);
	fprintf(FILE,"The reduced negation is:\t\t"); print_elem(&ne);
	gfp25632977add(&z, &ne, &e); gfp25632977makeunique(&z);
	fprintf(FILE,"The cross check value for reduced is:\t"); print_elem(&z);

	gfe_p25632977 x = {0x59f2815b16f81798,0x029bfcdb2dce28d9,0x55a06295ce870b07,0x79be667ef9dcbbac};
	gfe_p25632977 y = {0x9c47d08ffb10d4b8,0xfd17b448a6855419,0x5da4fbfc0e1108a8,0x483ada7726a3c465};
	ge_secp256k1 G = {x, y, 0};

	fprintf(FILE,"The point is:\n"); 
	fprintf(FILE,"x:\t\t"); print_elem(&G.x);
	fprintf(FILE,"y:\t\t"); print_elem(&G.y);

	gej_secp256k1 GGj;
	secp256k1double(&GGj, &G);
	fprintf(FILE,"The doubled point in projective coords is:\n"); 
	fprintf(FILE,"x:\t\t"); print_elem(&GGj.x);
	fprintf(FILE,"y:\t\t"); print_elem(&GGj.y);

	ge_secp256k1 GG;
	secp256k1_ge_from_gej(&GG, &GGj);
	fprintf(FILE,"The doubled point in affine coords is:\n"); 
	fprintf(FILE,"x:\t\t"); print_elem(&GG.x);
	fprintf(FILE,"y:\t\t"); print_elem(&GG.y);

	// fprintf(FILE,"Computing CPU-cycles. It will take some time. Please wait!\n\n");
	// MEASURE_TIME({gfp25632977mul(&t,&e,&e);change_input(e,t,e);});
	// fprintf(FILE,"CPU-cycles for a single field-multiplication is:%6.0lf\n\n", ceil(((get_median())/(double)(N))));
	// MEASURE_TIME({gfp25632977sqr(&t,&e);change_input(e,t,e);});
	// fprintf(FILE,"CPU-cycles for a single field-squaring is:%6.0lf\n\n", ceil(((get_median())/(double)(N))));
	// MEASURE_TIME({gfp25632977inv(&einv,&e);change_input(e,einv,e);});
	// fprintf(FILE,"CPU-cycles for a single field-inversion is:%6.0lf\n\n", ceil(((get_median())/(double)(N))));
	
	return 0;
}


void print_elem(const gfe_p25632977 *e){

	uchar8  i;

	for (i=NLIMBS-1; i>0; --i) 
		fprintf(FILE,"%16llX ",e->l[i]);
	fprintf(FILE,"%16llX \n\n",e->l[0]);
}

