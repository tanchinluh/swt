// Copyright (C) 2010 - H. Nahrstaedt
//
// dwt1d  Test 
loadmatfile("-mat",get_swt_path()+"tests/unit_tests/Data.mat");

// appcoef
[c,l]=wavedec(s1,3,'sym10');
[cA1,cD1]=dwt(s1,'sym10');
[cA2,cD2]=dwt(cA1,'sym10');
[cA3,cD3]=dwt(cA2,'sym10');
[Lo_R,Hi_R]=wfilters('sym10','r');
cAA2=idwt(cA3,cD3,'sym10',length(cA2));
cAA1=idwt(cAA2,cD2,'sym10',length(cA1));
ca31=appcoef(c,l,'sym10');
ca32=appcoef(c,l,'sym10',3);
ca33=appcoef(c,l,Lo_R,Hi_R);
ca34=appcoef(c,l,Lo_R,Hi_R,3);
ca21=appcoef(c,l,'sym10',2);
ca22=appcoef(c,l,Lo_R,Hi_R,2);
ca11=appcoef(c,l,'sym10',1);
ca12=appcoef(c,l,Lo_R,Hi_R,1);

assert_checkalmostequal ( ca31 , cA3 , %eps );
assert_checkalmostequal ( ca32 , cA3 , %eps );
assert_checkalmostequal ( ca33 , cA3 , %eps );
assert_checkalmostequal ( ca21 , cAA2 , %eps );
assert_checkalmostequal ( ca22 , cAA2 , %eps );
assert_checkalmostequal ( ca11 , cAA1 , %eps );
assert_checkalmostequal ( ca12 , cAA1 , %eps );
