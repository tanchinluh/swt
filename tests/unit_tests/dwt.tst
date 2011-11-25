// Copyright (C) 2010 - H. Nahrstaedt
//
// dwt1d  Test 

loadmatfile("-mat",get_swt_path()+"tests/unit_tests/Data.mat");

// dwt

// type 1 input
// haar
[cA,cD]=dwt(x1,'haar');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('haar');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );

[cA,cD]=dwt(x2,'haar');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'haar');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );

// db1
[cA,cD]=dwt(x1,'db1');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db1');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db1');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db1');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );

// db2
[cA,cD]=dwt(x1,'db2');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db2');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db2');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db2');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// db3
[cA,cD]=dwt(x1,'db3');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db3');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db3');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db3');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// db4
[cA,cD]=dwt(x1,'db4');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db4');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db4');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db4');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// db5
[cA,cD]=dwt(x1,'db5');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db5');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db5');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db5');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// db6
[cA,cD]=dwt(x1,'db6');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db6');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db6');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db6');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );

// db7
[cA,cD]=dwt(x1,'db7');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db7');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db7');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db7');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// db8
[cA,cD]=dwt(x1,'db8');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db8');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db8');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db8');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// db9
[cA,cD]=dwt(x1,'db9');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db9');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db9');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db9');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// db10
[cA,cD]=dwt(x1,'db10');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db10');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db10');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db10');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// coif1
[cA,cD]=dwt(x1,'coif1');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('coif1');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'coif1');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'coif1');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// coif2
[cA,cD]=dwt(x1,'coif2');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('coif2');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'coif2');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'coif2');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// coif3
[cA,cD]=dwt(x1,'coif3');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('coif3');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'coif3');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'coif3');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// coif4
[cA,cD]=dwt(x1,'coif4');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('coif4');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'coif4');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'coif4');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// coif5
[cA,cD]=dwt(x1,'coif5');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('coif5');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'coif5');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'coif5');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// sym4
[cA,cD]=dwt(x1,'sym4');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('sym4');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'sym4');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'sym4');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// sym5
[cA,cD]=dwt(x1,'sym5');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('sym5');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'sym5');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'sym5');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// sym6
[cA,cD]=dwt(x1,'sym6');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('sym6');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'sym6');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'sym6');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// sym7
[cA,cD]=dwt(x1,'sym7');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('sym7');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'sym7');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'sym7');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// sym8
[cA,cD]=dwt(x1,'sym8');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('sym8');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'sym8');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'sym8');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// sym9
[cA,cD]=dwt(x1,'sym9');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('sym9');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'sym9');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'sym9');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// sym10
[cA,cD]=dwt(x1,'sym10');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('sym10');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'sym10');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'sym10');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// bior1.1
[cA,cD]=dwt(x1,'bior1.1');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('bior1.1');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'bior1.1');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'bior1.1');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// bior1.3
[cA,cD]=dwt(x1,'bior1.3');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('bior1.3');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'bior1.3');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'bior1.3');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// bior1.5
[cA,cD]=dwt(x1,'bior1.5');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('bior1.5');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'bior1.5');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'bior1.5');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// bior2.2
[cA,cD]=dwt(x1,'bior2.2');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('bior2.2');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'bior2.2');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'bior2.2');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// bior2.4
[cA,cD]=dwt(x1,'bior2.4');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('bior2.4');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'bior2.4');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'bior2.4');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// bior2.6
[cA,cD]=dwt(x1,'bior2.6');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('bior2.6');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'bior2.6');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'bior2.6');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// bior2.8
[cA,cD]=dwt(x1,'bior2.8');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('bior2.8');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'bior2.8');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'bior2.8');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// bior3.1
[cA,cD]=dwt(x1,'bior3.1');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('bior3.1');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'bior3.1');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'bior3.1');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// bior3.3
[cA,cD]=dwt(x1,'bior3.3');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('bior3.3');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'bior3.3');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'bior3.3');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// bior3.5
[cA,cD]=dwt(x1,'bior3.5');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('bior3.5');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'bior3.5');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'bior3.5');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );



// bior3.7
[cA,cD]=dwt(x1,'bior3.7');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('bior3.7');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'bior3.7');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'bior3.7');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );



// bior3.9
[cA,cD]=dwt(x1,'bior3.9');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('bior3.9');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'bior3.9');
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'bior3.9');
caa=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// type 2 input
Lo_D=rand(1,20,'normal');
Hi_D=rand(1,20,'normal');
[cA,cD]=dwt(x1,Lo_D,Hi_D);
caa=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,Lo_D,Hi_D);
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,Lo_D,Hi_D);
caa=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// type 3 input
// symw
// haar
[cA,cD]=dwt(x1,'db10','mode','symw');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db10');
caa=dyaddown(wkeep(conv(wextend(1,'symw',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symw',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db10','mode','symw');
caa=dyaddown(wkeep(conv(wextend(1,'symw',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symw',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db10','mode','symw');
caa=dyaddown(wkeep(conv(wextend(1,'symw',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'symw',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// asymh
[cA,cD]=dwt(x1,'db10','mode','asymh');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db10');
caa=dyaddown(wkeep(conv(wextend(1,'asymh',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'asymh',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db10','mode','asymh');
caa=dyaddown(wkeep(conv(wextend(1,'asymh',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'asymh',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db10','mode','asymh');
caa=dyaddown(wkeep(conv(wextend(1,'asymh',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'asymh',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// asymw
[cA,cD]=dwt(x1,'db10','mode','asymw');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db10');
caa=dyaddown(wkeep(conv(wextend(1,'asymw',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'asymw',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db10','mode','asymw');
caa=dyaddown(wkeep(conv(wextend(1,'asymw',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'asymw',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db10','mode','asymw');
caa=dyaddown(wkeep(conv(wextend(1,'asymw',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'asymw',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );

// zpd
[cA,cD]=dwt(x1,'db10','mode','zpd');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db10');
caa=dyaddown(wkeep(conv(wextend(1,'zpd',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'zpd',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db10','mode','zpd');
caa=dyaddown(wkeep(conv(wextend(1,'zpd',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'zpd',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db10','mode','zpd');
caa=dyaddown(wkeep(conv(wextend(1,'zpd',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'zpd',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );

// sp0
[cA,cD]=dwt(x1,'db10','mode','sp0');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db10');
caa=dyaddown(wkeep(conv(wextend(1,'sp0',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'sp0',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db10','mode','sp0');
caa=dyaddown(wkeep(conv(wextend(1,'sp0',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'sp0',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db10','mode','sp0');
caa=dyaddown(wkeep(conv(wextend(1,'sp0',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'sp0',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// sp1
[cA,cD]=dwt(x1,'db10','mode','sp1');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db10');
caa=dyaddown(wkeep(conv(wextend(1,'sp1',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'sp1',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db10','mode','sp1');
caa=dyaddown(wkeep(conv(wextend(1,'sp1',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'sp1',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db10','mode','sp1');
caa=dyaddown(wkeep(conv(wextend(1,'sp1',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'sp1',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// ppd
[cA,cD]=dwt(x1,'db10','mode','ppd');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db10');
caa=dyaddown(wkeep(conv(wextend(1,'ppd',x1,length(Lo_D),'b'),Lo_D),(length(x1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'ppd',x1,length(Lo_D),'b'),Hi_D),(length(x1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db10','mode','ppd');
caa=dyaddown(wkeep(conv(wextend(1,'ppd',x2,length(Lo_D),'b'),Lo_D),(length(x2)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'ppd',x2,length(Lo_D),'b'),Hi_D),(length(x2)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db10','mode','ppd');
caa=dyaddown(wkeep(conv(wextend(1,'ppd',s1,length(Lo_D),'b'),Lo_D),(length(s1)+length(Lo_D)-1)));
cdd=dyaddown(wkeep(conv(wextend(1,'ppd',s1,length(Lo_D),'b'),Hi_D),(length(s1)+length(Lo_D)-1)));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );


// per
[cA,cD]=dwt(x1,'db10','mode','per');
[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db10');
caa=dyaddown(wkeep(conv(wextend(1,'per',x1,length(Lo_D),'b'),Lo_D),(length(x1))));
cdd=dyaddown(wkeep(conv(wextend(1,'per',x1,length(Lo_D),'b'),Hi_D),(length(x1))));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(x2,'db10','mode','per');
caa=dyaddown(wkeep(conv(wextend(1,'per',x2,length(Lo_D),'b'),Lo_D),(length(x2))));
cdd=dyaddown(wkeep(conv(wextend(1,'per',x2,length(Lo_D),'b'),Hi_D),(length(x2))));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
[cA,cD]=dwt(s1,'db10','mode','per');
caa=dyaddown(wkeep(conv(wextend(1,'per',s1,length(Lo_D),'b'),Lo_D),(length(s1))));
cdd=dyaddown(wkeep(conv(wextend(1,'per',s1,length(Lo_D),'b'),Hi_D),(length(s1))));
assert_checkalmostequal ( caa , cA , %eps );
assert_checkalmostequal ( cdd , cD , %eps );
