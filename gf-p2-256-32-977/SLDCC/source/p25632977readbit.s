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

movq    $8, %rax
mulq    %rdx
movq    (%rdi), %rcx
movq    $1, (%rdi)
movq    (%rsi, %rax), %r8
shr     %cl, %r8
andq    %r8, (%rdi)

pop     %r15
pop     %r14
pop     %r13
pop     %r12
pop     %rbx
pop     %rbp

ret
