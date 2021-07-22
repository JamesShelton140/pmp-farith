// bit = n[bit]
// On call: %rdi = &bit, %rsi = &n, %rdx = limb

.p2align 5
.globl gfp25632977readbit
gfp25632977readbit:

push    %rbp
push    %rbx
push    %r12
push    %r13
push    %r14
push    %r15

movq    $8, %rax            ;Move 8 into rax
mulq    %rdx                ;Multiply rax and rdx
movq    (%rdi), %rcx        ;Move the value of bit into rcx
movq    $1, (%rdi)          ;Move 1 into bit
movq    (%rsi, %rax), %r8   ;Move the rax-limb of n(at rsi) into r8
shr     %cl, %r8            ;Shift r8 right by cl bits
andq    %r8, (%rdi)         ;Set bit (at rdi) to r8 & 1.

pop     %r15
pop     %r14
pop     %r13
pop     %r12
pop     %rbx
pop     %rbp

ret
