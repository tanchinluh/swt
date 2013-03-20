// Copyright (C) 2010 - H. Nahrstaedt
//
// dwt1d  Test 

loadmatfile("-mat",get_swt_path()+"tests/unit_tests/Data.mat");

// wenergy
[c,l]=wavedec(s1,3,'db10');
[cA1,cD1]=dwt(s1,'db10');
[cA2,cD2]=dwt(cA1,'db10');
[cA3,cD3]=dwt(cA2,'db10');
ea=sum(cA3.*cA3);
ed1=sum(cD1.*cD1);
ed2=sum(cD2.*cD2);
ed3=sum(cD3.*cD3);
energy=sum(c.*c);
ea=ea*100/energy;
ed1=ed1*100/energy;
ed2=ed2*100/energy;
ed3=ed3*100/energy;
[Ea,Ed]=wenergy(c,l);
e=sum(abs(ea-Ea))+sum(abs(Ed-[ed3 ed2 ed1]));

assert_checkalmostequal ( ea , Ea , %eps , 1e-10);
assert_checkalmostequal ( Ed , [ed3 ed2 ed1] , %eps, 1e-10 );
