// Copyright (C) 2010 - H. Nahrstaedt
//
// filter Test
// dwtmode


a='asymh';
dwtmode(a);
ST=dwtmode('status','nodisp');
b=(ST==a);
x=rand(1,50,'normal');
[cA,cD]=dwt(x,'db2');
[caa,cdd]=dwt(x,'db2','mode',a);
x0=idwt(cA,cD,'db2',length(x));

assert_checkalmostequal ( cA , caa , %eps );
assert_checkalmostequal ( cD , cdd , %eps );
assert_checkalmostequal ( x , x0 , %eps*1000 );

assert_checktrue(b);


a='sp1';
dwtmode(a);
ST=dwtmode('status','nodisp');
b=(ST==a);
x=rand(1,50,'normal');
[cA,cD]=dwt(x,'db2');
[caa,cdd]=dwt(x,'db2','mode',a);
x0=idwt(cA,cD,'db2',length(x));

assert_checkalmostequal ( cA , caa , %eps );
assert_checkalmostequal ( cD , cdd , %eps );
assert_checkalmostequal ( x , x0 , %eps*1000 );

assert_checktrue(b);


dwtmode("symh");