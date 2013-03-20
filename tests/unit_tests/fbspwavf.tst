// Copyright (C) 2010 - H. Nahrstaedt
//
// shanwavf test



function y = sinc2(x)
y = ones(x);
k = find(x);
y(k) = sin(%pi*x(k))./(%pi*x(k));
endfunction

function [psi,x] = ref_fbspwavf(LB,UB,N,m,Fb,Fc);
x = linspace(LB,UB,N);  
psi = (Fb^0.5)*((sinc2(Fb*x/m).^m).*exp(2*%i*%pi*Fc*x));
endfunction


LB=-20;
UB=20;
N=1000;
Fb=1;
Fc=1.5;
M=2;

[psi,x] = ref_fbspwavf(LB,UB,N,M,Fb,Fc);
[PSI,X]=fbspwavf(LB,UB,N,M,Fc,Fb);

assert_checkalmostequal ( real(PSI) , real(psi) , %eps, 1e-9 );
assert_checkalmostequal ( imag(PSI) , imag(psi) , %eps, 1e-9 );
assert_checkalmostequal ( X , x , %eps, %eps );

LB=-20;
UB=20;
N=1000;
Fb=1.5;
Fc=1;
M=2;

[psi,x] = ref_fbspwavf(LB,UB,N,M,Fb,Fc);
[PSI,X]=fbspwavf(LB,UB,N,M,Fc,Fb);

assert_checkalmostequal ( real(PSI) , real(psi) , %eps, 1e-9 );
assert_checkalmostequal ( imag(PSI) , imag(psi) , %eps, 1e-9 );
assert_checkalmostequal ( X , x , %eps, %eps );

LB=-20;
UB=20;
N=1000;
Fb=1.5;
Fc=1.2;
M=3;

[psi,x] = ref_fbspwavf(LB,UB,N,M,Fb,Fc);
[PSI,X]=fbspwavf(LB,UB,N,M,Fc,Fb);

assert_checkalmostequal ( real(PSI) , real(psi) , %eps, 1e-9 );
assert_checkalmostequal ( imag(PSI) , imag(psi) , %eps, 1e-9 );
assert_checkalmostequal ( X , x , %eps, %eps );