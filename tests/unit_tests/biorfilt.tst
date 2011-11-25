// Copyright (C) 2010 - H. Nahrstaedt
//
// filter Test
// biorfilt


df=rand(1,50,'normal');
rf=rand(1,50,'normal');
df_lo_r=df;
df_lo_d=wrev(df_lo_r);
df_hi_r=qmf(df_lo_r);
df_hi_d=wrev(df_hi_r);
rf_lo_r=rf;
rf_lo_d=wrev(rf_lo_r);
rf_hi_r=qmf(rf_lo_r);
rf_hi_d=wrev(rf_hi_r);
[Lo_D,Hi_D,Lo_R,Hi_R]=biorfilt(df,rf);

assert_checkalmostequal ( Hi_R , df_hi_r , %eps );
assert_checkalmostequal ( Lo_D , df_lo_d , %eps );
assert_checkalmostequal ( Hi_D , rf_hi_d , %eps );
assert_checkalmostequal ( Lo_R , rf_lo_r , %eps );
