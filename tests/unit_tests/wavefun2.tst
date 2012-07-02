// Copyright (C) 2012 - H. Nahrstaedt
//
// wavefun2  Test 
ITER=4;
// dwt
// type 1 input
// haar

[S,W1,W2,W3,XYVAL]=wavefun2('haar',ITER);
assert_checkalmostequal ( sum(S)/(2^ITER)^2 , 1 , %eps, 1e-7 );
assert_checkalmostequal ( sum(W1.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W1) , 0 , %eps, 1e-10 );
assert_checkalmostequal ( sum(W2.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W2) , 0 , %eps, 1e-10 );
assert_checkalmostequal ( sum(W3.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W3) , 0 , %eps, 1e-10 );

//db family
db_fam={"db1","db2", "db3", "db4", "db5", "db6","db7", "db8", "db9", "db10", "db11",... 
"db12", "db13", "db14", "db15", "db16","db17", "db18", "db19", "db20"}; 
// db1
for i=1:max(size(db_fam))
  [S,W1,W2,W3,XYVAL]=wavefun2(db_fam(i),ITER);
    assert_checkalmostequal ( sum(S)/(2^ITER)^2 , 1 , %eps, 1e-7 );
    assert_checkalmostequal ( sum(W1.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
    assert_checkalmostequal ( sum(W1) , 0 , %eps, 1e-8 );
    assert_checkalmostequal ( sum(W2.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
    assert_checkalmostequal ( sum(W2) , 0 , %eps, 1e-8 );
    assert_checkalmostequal ( sum(W3.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
    assert_checkalmostequal ( sum(W3) , 0 , %eps, 1e-8 );
end;

//coif family
coif_fam={"coif1","coif2","coif3","coif4","coif5"};

for i=1:max(size(coif_fam))

  [S,W1,W2,W3,XYVAL]=wavefun2(coif_fam(i),ITER);
    assert_checkalmostequal ( sum(S)/(2^ITER)^2 , 1 , %eps, 1e-7 );
    assert_checkalmostequal ( sum(W1.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
    assert_checkalmostequal ( sum(W1) , 0 , %eps, 1e-7 );
    assert_checkalmostequal ( sum(W2.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
    assert_checkalmostequal ( sum(W2) , 0 , %eps, 1e-7 );
    assert_checkalmostequal ( sum(W3.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
    assert_checkalmostequal ( sum(W3) , 0 , %eps, 1e-7 );
end;

//symlets family
sym_fam={"sym2", "sym3", "sym4", "sym5", "sym6","sym7", "sym8", "sym9", "sym10", "sym11",... 
"sym12", "sym13", "sym14", "sym15", "sym16","sym17", "sym18", "sym19", "sym20"}; 

for i=1:max(size(sym_fam))
  [S,W1,W2,W3,XYVAL]=wavefun2(sym_fam(i),ITER);
    assert_checkalmostequal ( sum(S)/(2^ITER)^2 , 1 , %eps, 1e-7 );
    assert_checkalmostequal ( sum(W1.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
    assert_checkalmostequal ( sum(W1) , 0 , %eps, 1e-7 );
    assert_checkalmostequal ( sum(W2.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
    assert_checkalmostequal ( sum(W2) , 0 , %eps, 1e-7 );
    assert_checkalmostequal ( sum(W3.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
    assert_checkalmostequal ( sum(W3) , 0 , %eps, 1e-7 );
end;




//beylkin
  [S,W1,W2,W3,XYVAL]=wavefun2("beylkin",ITER);
assert_checkalmostequal ( sum(S)/(2^ITER)^2 , 1 , %eps, 1e-7 );
assert_checkalmostequal ( sum(W1.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W1) , 0 , %eps, 1e-5 );
assert_checkalmostequal ( sum(W2.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W2) , 0 , %eps, 1e-5 );
assert_checkalmostequal ( sum(W3.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W3) , 0 , %eps, 1e-5 );
  
  //vaidyanathan
    [S,W1,W2,W3,XYVAL]=wavefun2("vaidyanathan",ITER);
assert_checkalmostequal ( sum(S)/(2^ITER)^2 , 1 , %eps, 1e-7 );
assert_checkalmostequal ( sum(W1.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
// assert_checkalmostequal ( sum(W1) , 0 , %eps, 1e-5 );
assert_checkalmostequal ( sum(W2.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
// assert_checkalmostequal ( sum(W2) , 0 , %eps, 1e-5 );
assert_checkalmostequal ( sum(W3.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
// assert_checkalmostequal ( sum(W3) , 0 , %eps, 1e-5 );
  //dmey
    [S,W1,W2,W3,XYVAL]=wavefun2("dmey",ITER);
assert_checkalmostequal ( sum(S)/(2^ITER)^2 , 1 , %eps, 1e-1 );
 assert_checkalmostequal ( sum(W1.^2)/(2^ITER)^2 , 1 , %eps, 1e-1  );
// assert_checkalmostequal ( sum(W1) , 0 , %eps, 1e-1 );
 assert_checkalmostequal ( sum(W2.^2)/(2^ITER)^2 , 1 , %eps, 1e-1  );
// assert_checkalmostequal ( sum(W2) , 0 , %eps, 1e-1 );
assert_checkalmostequal ( sum(W3.^2)/(2^ITER)^2 , 1 , %eps, 1e-1  );
// assert_checkalmostequal ( sum(W3) , 0 , %eps, 1e-1 );
  //bath
  bath_fam={"bath4.0", "bath4.1", "bath4.2", "bath4.3", "bath4.4", "bath4.5",...
"bath4.6", "bath4.7", "bath4.8", "bath4.9", "bath4.10", ...
"bath4.11", "bath4.12", "bath4.13", "bath4.14", "bath4.15", ...
"bath6.0", "bath6.1", "bath6.2", "bath6.3", "bath6.4", ...
"bath6.5", "bath6.6", "bath6.7", "bath6.8", "bath6.9", ...
"bath6.10", "bath6.11", "bath6.12", "bath6.13", "bath6.14", ...
"bath6.15"};
for i=1:max(size(bath_fam))

  [S,W1,W2,W3,XYVAL]=wavefun2(bath_fam(i),ITER);
assert_checkalmostequal ( sum(S)/(2^ITER)^2 , 1 , %eps, 1e-2 );
assert_checkalmostequal ( sum(W1.^2)/(2^ITER)^2 , 1 , %eps, 1e-2  );
assert_checkalmostequal ( sum(W1) , 0 , %eps, 1e-1 );
assert_checkalmostequal ( sum(W2.^2)/(2^ITER)^2 , 1 , %eps, 1e-2  );
assert_checkalmostequal ( sum(W2) , 0 , %eps, 1e-1 );
assert_checkalmostequal ( sum(W3.^2)/(2^ITER)^2 , 1 , %eps, 1e-2  );
assert_checkalmostequal ( sum(W3) , 0 , %eps, 1e-1 );
end;
  //legd
legd_fam={"legd1", "legd2", "legd3", "legd4", "legd5", "legd6", "legd7", "legd8", "legd9"};

for i=1:max(size(legd_fam))

  [S,W1,W2,W3,XYVAL]=wavefun2(legd_fam(i),ITER);
assert_checkalmostequal ( sum(S)/(2^ITER)^2 , 1 , %eps, 1e-4 );
// assert_checkalmostequal ( sum(W1.^2)/(2^ITER)^2 , 1 , %eps, 1e-4  );
assert_checkalmostequal ( sum(W1) , 0 , %eps, 1e-5 );
// assert_checkalmostequal ( sum(W2.^2)/(2^ITER)^2 , 1 , %eps, 1e-4  );
assert_checkalmostequal ( sum(W2) , 0 , %eps, 1e-5 );
// assert_checkalmostequal ( sum(W3.^2)/(2^ITER)^2 , 1 , %eps, 1e-4  );
assert_checkalmostequal ( sum(W3) , 0 , %eps, 1e-5 );
end;

  //fa
  [S,W1,W2,W3,XYVAL]=wavefun2("fa1",ITER);
assert_checkalmostequal ( sum(S)/(2^ITER)^2 , 1 , %eps, 1e-7 );
assert_checkalmostequal ( sum(W1.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W1) , 0 , %eps, 1e-5 );
assert_checkalmostequal ( sum(W2.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W2) , 0 , %eps, 1e-5 );
assert_checkalmostequal ( sum(W3.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W3) , 0 , %eps, 1e-5 );
   [S,W1,W2,W3,XYVAL]=wavefun2("fa2",ITER);
assert_checkalmostequal ( sum(S)/(2^ITER)^2 , 1 , %eps, 1e-7 );
assert_checkalmostequal ( sum(W1.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W1) , 0 , %eps, 1e-5 );
assert_checkalmostequal ( sum(W2.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W2) , 0 , %eps, 1e-5 );
assert_checkalmostequal ( sum(W3.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W3) , 0 , %eps, 1e-5 );
  //  ksq
    [S,W1,W2,W3,XYVAL]=wavefun2("ksq1",ITER);
assert_checkalmostequal ( sum(S)/(2^ITER)^2 , 1 , %eps, 1e-7 );
assert_checkalmostequal ( sum(W1.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W1) , 0 , %eps, 1e-5 );
assert_checkalmostequal ( sum(W2.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W2) , 0 , %eps, 1e-5 );
assert_checkalmostequal ( sum(W3.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W3) , 0 , %eps, 1e-5 );

      [S,W1,W2,W3,XYVAL]=wavefun2("ksq2",ITER);
assert_checkalmostequal ( sum(S)/(2^ITER)^2 , 1 , %eps, 1e-7 );
assert_checkalmostequal ( sum(W1.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W1) , 0 , %eps, 1e-5 );
assert_checkalmostequal ( sum(W2.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W2) , 0 , %eps, 1e-5 );
assert_checkalmostequal ( sum(W3.^2)/(2^ITER)^2 , 1 , %eps, 1e-7  );
assert_checkalmostequal ( sum(W3) , 0 , %eps, 1e-5 );
  
  
  