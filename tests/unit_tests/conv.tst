// Copyright (C) 2010 - H. Nahrstaedt
//
// Utility Function Test
// Convolution

a=rand(1,100,'normal');
b=rand(1,51,'uniform');
c=conv(a,b);
c1=conv(a',b);
c2=conv(a,b');
c3=conv(a',b');

expected=convol(a,b);

assert_checkalmostequal ( c , expected , %eps*10 );
assert_checkalmostequal ( c1 , expected , %eps*10 );
assert_checkalmostequal ( c2 , expected , %eps*10 );
assert_checkalmostequal ( c3 , expected , %eps*10 );

