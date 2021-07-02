#+-----------------------------------------------------------------------------+
#| This code corresponds to the Galois field F(2^256-2^32-977) from the paper  |
#| "Efficient Arithmetic In (Pseudo-)Mersenne Prime Order Fields" authored by  |
#| Kaushik Nath,  Indian Statistical Institute, Kolkata, India, and            |
#| Palash Sarkar, Indian Statistical Institute, Kolkata, India.	               |
#+-----------------------------------------------------------------------------+
#| Copyright (c) 2018, Kaushik Nath, Palash Sarkar.                            |
#|                                                                             |
#| Permission to use this code is granted.                          	       |
#|                                                                             |
#| Redistribution and use in source and binary forms, with or without          |
#| modification, are permitted provided that the following conditions are      |
#| met:                                                                        |
#|                                                                             |
#| * Redistributions of source code must retain the above copyright notice,    |
#|   this list of conditions and the following disclaimer.                     |
#|                                                                             |
#| * Redistributions in binary form must reproduce the above copyright         |
#|   notice, this list of conditions and the following disclaimer in the       |
#|   documentation and/or other materials provided with the distribution.      |
#|                                                                             |
#| * The names of the contributors may not be used to endorse or promote       |
#|   products derived from this software without specific prior written        |
#|   permission.                                                               |
#+-----------------------------------------------------------------------------+
#| THIS SOFTWARE IS PROVIDED BY THE AUTHORS ""AS IS"" AND ANY EXPRESS OR       |
#| IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES   |
#| OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.     |
#| IN NO EVENT SHALL THE CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,      |
#| INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT    |
#| NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,   |
#| DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY       |
#| THEORY LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING |
#| NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,| 
#| EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                          |
#+-----------------------------------------------------------------------------+
#*/

INCDRS     = -I../include/

SRCFLS = ../source/p25632977consts.s ../source/p25632977mul.s ../source/p25632977add.s ../source/p25632977negate.s ../source/p25632977sub.s ../source/p25632977nsqr.s ../source/p25632977makeunique.s ../source/p25632977inv.c ../source/secp256k1add.c ../source/secp256k1double.c ../source/secp256k1ecdsa.c ../source/secp256k1group.c ../source/secp256k1mladder.c ../source/p25632977readbit.s p25632977_test.c
         
OBJFLS = ../source/p25632977consts.o ../source/p25632977mul.o ../source/p25632977add.o ../source/p25632977negate.o ../source/p25632977sub.o ../source/p25632977nsqr.o ../source/p25632977makeunique.o ../source/p25632977inv.o ../source/secp256k1add.o ../source/secp256k1double.o ../source/secp256k1ecdsa.o ../source/secp256k1group.o ../source/secp256k1mladder.o ../source/p25632977readbit.o p25632977_test.o

EXE    = p25632977_test

CFLAGS = -march=skylake -mtune=skylake -m64 -O3 -funroll-loops -fomit-frame-pointer

CC     = gcc
LL     = gcc

$(EXE): $(OBJFLS)
	$(LL) -o $@ $(OBJFLS) -lm -no-pie 

.c.o:
	$(CC) $(INCDRS) $(CFLAGS) -o $@ -c $<

clean:
	-rm $(EXE)
	-rm $(OBJFLS)
