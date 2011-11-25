// Copyright (C) 2010 - H. Nahrstaedt
//
// wfilter  Test 

[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('haar');
w=dbwavf('db1');
[lo_d,hi_d,lo_r,hi_r]=orthfilt(w);
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );

[Lo_D,Hi_D]=wfilters('haar','d');
[Lo_R,Hi_R]=wfilters('haar','r');
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );

[Lo_D,Lo_R]=wfilters('haar','l');
[Hi_D,Hi_R]=wfilters('haar','h');
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );


[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db1');
w=dbwavf('db1');
[lo_d,hi_d,lo_r,hi_r]=orthfilt(w);
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );

[Lo_D,Hi_D]=wfilters('db1','d');
[Lo_R,Hi_R]=wfilters('db1','r');
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );

[Lo_D,Lo_R]=wfilters('db1','l');
[Hi_D,Hi_R]=wfilters('db1','h');
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );


[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db2');
w=dbwavf('db2');
[lo_d,hi_d,lo_r,hi_r]=orthfilt(w);
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );

[Lo_D,Hi_D]=wfilters('db2','d');
[Lo_R,Hi_R]=wfilters('db2','r');
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );

[Lo_D,Lo_R]=wfilters('db2','l');
[Hi_D,Hi_R]=wfilters('db2','h');
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );


[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('db3');
w=dbwavf('db3');
[lo_d,hi_d,lo_r,hi_r]=orthfilt(w);
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );

[Lo_D,Hi_D]=wfilters('db3','d');
[Lo_R,Hi_R]=wfilters('db3','r');
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );

[Lo_D,Lo_R]=wfilters('db3','l');
[Hi_D,Hi_R]=wfilters('db3','h');
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );


[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('bior1.1');
[rf,df]=biorwavf('bior1.1');
[lo_d,hi_d,lo_r,hi_r]=biorfilt(df,rf);
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );

[Lo_D,Hi_D]=wfilters('bior1.1','d');
[Lo_R,Hi_R]=wfilters('bior1.1','r');
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );

[Lo_D,Lo_R]=wfilters('bior1.1','l');
[Hi_D,Hi_R]=wfilters('bior1.1','h');
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );


[Lo_D,Hi_D,Lo_R,Hi_R]=wfilters('bior1.5');
[rf,df]=biorwavf('bior1.5');
[lo_d,hi_d,lo_r,hi_r]=biorfilt(df,rf);
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );

[Lo_D,Hi_D]=wfilters('bior1.5','d');
[Lo_R,Hi_R]=wfilters('bior1.5','r');
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );

[Lo_D,Lo_R]=wfilters('bior1.5','l');
[Hi_D,Hi_R]=wfilters('bior1.5','h');
assert_checkalmostequal ( lo_d , Lo_D , %eps );
assert_checkalmostequal ( lo_r , Lo_R , %eps );
assert_checkalmostequal ( hi_d , Hi_D , %eps );
assert_checkalmostequal ( hi_r , Hi_R , %eps );


clear Lo_D;
clear Lo_R;
clear Hi_R;
clear Hi_D;
clear lo_d;
clear lo_r;
clear hi_r;
clear hi_d;
