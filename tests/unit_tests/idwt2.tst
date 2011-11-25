// Copyright (C) 2010 - H. Nahrstaedt
//
// dwt2d  Test 

loadmatfile("-mat",get_swt_path()+"tests/unit_tests/Data.mat");
clear row_low;
clear row_hi;	
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;


// idwt2
//loadmatfile("-mat","Data.mat");
d1=rand(64,65,'normal');
// haar
[cA,cH,cV,cD]=dwt2(d1,'haar');
[Lo_R,Hi_R]=wfilters('haar','r');
dd0=idwt2(cA,cH,cV,cD,'haar');
dd1=idwt2(cA,cH,cV,cD,'haar',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));

assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );

clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;

// db1
[cA,cH,cV,cD]=dwt2(d1,'db1');
[Lo_R,Hi_R]=wfilters('db1','r');
dd0=idwt2(cA,cH,cV,cD,'db1');
dd1=idwt2(cA,cH,cV,cD,'db1',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear row_low_r;
clear row_hi_r;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;


// db2
[cA,cH,cV,cD]=dwt2(d1,'db2');
[Lo_R,Hi_R]=wfilters('db2','r');
dd0=idwt2(cA,cH,cV,cD,'db2');
dd1=idwt2(cA,cH,cV,cD,'db2',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;

// db3
[cA,cH,cV,cD]=dwt2(d1,'db3');
[Lo_R,Hi_R]=wfilters('db3','r');
dd0=idwt2(cA,cH,cV,cD,'db3');
dd1=idwt2(cA,cH,cV,cD,'db3',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;

// db4
[cA,cH,cV,cD]=dwt2(d1,'db4');
[Lo_R,Hi_R]=wfilters('db4','r');
dd0=idwt2(cA,cH,cV,cD,'db4');
dd1=idwt2(cA,cH,cV,cD,'db4',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;

// db5
[cA,cH,cV,cD]=dwt2(d1,'db5');
[Lo_R,Hi_R]=wfilters('db5','r');
dd0=idwt2(cA,cH,cV,cD,'db5');
dd1=idwt2(cA,cH,cV,cD,'db5',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;

// db6
[cA,cH,cV,cD]=dwt2(d1,'db6');
[Lo_R,Hi_R]=wfilters('db6','r');
dd0=idwt2(cA,cH,cV,cD,'db6');
dd1=idwt2(cA,cH,cV,cD,'db6',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;

// db7
[cA,cH,cV,cD]=dwt2(d1,'db7');
[Lo_R,Hi_R]=wfilters('db7','r');
dd0=idwt2(cA,cH,cV,cD,'db7');
dd1=idwt2(cA,cH,cV,cD,'db7',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;

// db8
[cA,cH,cV,cD]=dwt2(d1,'db8');
[Lo_R,Hi_R]=wfilters('db8','r');
dd0=idwt2(cA,cH,cV,cD,'db8');
dd1=idwt2(cA,cH,cV,cD,'db8',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// db9
[cA,cH,cV,cD]=dwt2(d1,'db9');
[Lo_R,Hi_R]=wfilters('db9','r');
dd0=idwt2(cA,cH,cV,cD,'db9');
dd1=idwt2(cA,cH,cV,cD,'db9',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;



// db10
[cA,cH,cV,cD]=dwt2(d1,'db10');
[Lo_R,Hi_R]=wfilters('db10','r');
dd0=idwt2(cA,cH,cV,cD,'db10');
dd1=idwt2(cA,cH,cV,cD,'db10',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// coif1
[cA,cH,cV,cD]=dwt2(d1,'coif1');
[Lo_R,Hi_R]=wfilters('coif1','r');
dd0=idwt2(cA,cH,cV,cD,'coif1');
dd1=idwt2(cA,cH,cV,cD,'coif1',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// coif2
[cA,cH,cV,cD]=dwt2(d1,'coif2');
[Lo_R,Hi_R]=wfilters('coif2','r');
dd0=idwt2(cA,cH,cV,cD,'coif2');
dd1=idwt2(cA,cH,cV,cD,'coif2',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// coif3
[cA,cH,cV,cD]=dwt2(d1,'coif3');
[Lo_R,Hi_R]=wfilters('coif3','r');
dd0=idwt2(cA,cH,cV,cD,'coif3');
dd1=idwt2(cA,cH,cV,cD,'coif3',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// coif4
[cA,cH,cV,cD]=dwt2(d1,'coif4');
[Lo_R,Hi_R]=wfilters('coif4','r');
dd0=idwt2(cA,cH,cV,cD,'coif4');
if exists('dd0')
dd1=idwt2(cA,cH,cV,cD,'coif4',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;
end

// coif5
[cA,cH,cV,cD]=dwt2(d1,'coif5');
[Lo_R,Hi_R]=wfilters('coif5','r');
dd0=idwt2(cA,cH,cV,cD,'coif5');
if exists('dd0')
dd1=idwt2(cA,cH,cV,cD,'coif5',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;
end

// sym4
[cA,cH,cV,cD]=dwt2(d1,'sym4');
[Lo_R,Hi_R]=wfilters('sym4','r');
dd0=idwt2(cA,cH,cV,cD,'sym4');
if exists('dd0')
dd1=idwt2(cA,cH,cV,cD,'sym4',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;
end

// sym5
[cA,cH,cV,cD]=dwt2(d1,'sym5');
[Lo_R,Hi_R]=wfilters('sym5','r');
dd0=idwt2(cA,cH,cV,cD,'sym5');
if exists('dd0')
dd1=idwt2(cA,cH,cV,cD,'sym5',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;
end

// sym6
[cA,cH,cV,cD]=dwt2(d1,'sym6');
[Lo_R,Hi_R]=wfilters('sym6','r');
dd0=idwt2(cA,cH,cV,cD,'sym6');
if exists('dd0')
dd1=idwt2(cA,cH,cV,cD,'sym6',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;
end

// sym7
[cA,cH,cV,cD]=dwt2(d1,'sym7');
[Lo_R,Hi_R]=wfilters('sym7','r');
dd0=idwt2(cA,cH,cV,cD,'sym7');
if exists('dd0')
dd1=idwt2(cA,cH,cV,cD,'sym7',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;
end

// sym8
[cA,cH,cV,cD]=dwt2(d1,'sym8');
[Lo_R,Hi_R]=wfilters('sym8','r');
dd0=idwt2(cA,cH,cV,cD,'sym8');
dd1=idwt2(cA,cH,cV,cD,'sym8',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// sym9
[cA,cH,cV,cD]=dwt2(d1,'sym9');
[Lo_R,Hi_R]=wfilters('sym9','r');
dd0=idwt2(cA,cH,cV,cD,'sym9');
dd1=idwt2(cA,cH,cV,cD,'sym9',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// sym10
[cA,cH,cV,cD]=dwt2(d1,'sym10');
[Lo_R,Hi_R]=wfilters('sym10','r');
dd0=idwt2(cA,cH,cV,cD,'sym10');
dd1=idwt2(cA,cH,cV,cD,'sym10',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// bior1.1
[cA,cH,cV,cD]=dwt2(d1,'bior1.1');
[Lo_R,Hi_R]=wfilters('bior1.1','r');
dd0=idwt2(cA,cH,cV,cD,'bior1.1');
dd1=idwt2(cA,cH,cV,cD,'bior1.1',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;

// bior1.3
[cA,cH,cV,cD]=dwt2(d1,'bior1.3');
[Lo_R,Hi_R]=wfilters('bior1.3','r');
dd0=idwt2(cA,cH,cV,cD,'bior1.3');
dd1=idwt2(cA,cH,cV,cD,'bior1.3',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// bior1.5
[cA,cH,cV,cD]=dwt2(d1,'bior1.5');
[Lo_R,Hi_R]=wfilters('bior1.5','r');
dd0=idwt2(cA,cH,cV,cD,'bior1.5');
dd1=idwt2(cA,cH,cV,cD,'bior1.5',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// bior2.2
[cA,cH,cV,cD]=dwt2(d1,'bior2.2');
[Lo_R,Hi_R]=wfilters('bior2.2','r');
dd0=idwt2(cA,cH,cV,cD,'bior2.2');
dd1=idwt2(cA,cH,cV,cD,'bior2.2',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// bior2.4
[cA,cH,cV,cD]=dwt2(d1,'bior2.4');
[Lo_R,Hi_R]=wfilters('bior2.4','r');
dd0=idwt2(cA,cH,cV,cD,'bior2.4');
dd1=idwt2(cA,cH,cV,cD,'bior2.4',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// bior2.6
[cA,cH,cV,cD]=dwt2(d1,'bior2.6');
[Lo_R,Hi_R]=wfilters('bior2.6','r');
dd0=idwt2(cA,cH,cV,cD,'bior2.6');
dd1=idwt2(cA,cH,cV,cD,'bior2.6',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// bior2.8
[cA,cH,cV,cD]=dwt2(d1,'bior2.8');
[Lo_R,Hi_R]=wfilters('bior2.8','r');
dd0=idwt2(cA,cH,cV,cD,'bior2.8');
dd1=idwt2(cA,cH,cV,cD,'bior2.8',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// bior3.1
[cA,cH,cV,cD]=dwt2(d1,'bior3.1');
[Lo_R,Hi_R]=wfilters('bior3.1','r');
dd0=idwt2(cA,cH,cV,cD,'bior3.1');
dd1=idwt2(cA,cH,cV,cD,'bior3.1',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// bior3.3
[cA,cH,cV,cD]=dwt2(d1,'bior3.3');
[Lo_R,Hi_R]=wfilters('bior3.3','r');
dd0=idwt2(cA,cH,cV,cD,'bior3.3');
dd1=idwt2(cA,cH,cV,cD,'bior3.3',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// bior3.5
[cA,cH,cV,cD]=dwt2(d1,'bior3.5');
[Lo_R,Hi_R]=wfilters('bior3.5','r');
dd0=idwt2(cA,cH,cV,cD,'bior3.5');
dd1=idwt2(cA,cH,cV,cD,'bior3.5',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// bior3.7
[cA,cH,cV,cD]=dwt2(d1,'bior3.7');
[Lo_R,Hi_R]=wfilters('bior3.7','r');
dd0=idwt2(cA,cH,cV,cD,'bior3.7');
dd1=idwt2(cA,cH,cV,cD,'bior3.7',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;


// bior3.9
[cA,cH,cV,cD]=dwt2(d1,'bior3.9');
[Lo_R,Hi_R]=wfilters('bior3.9','r');
dd0=idwt2(cA,cH,cV,cD,'bior3.9');
dd1=idwt2(cA,cH,cV,cD,'bior3.9',size(d1));
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R);
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1));
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,ss);
x1=wkeep(row_low_r+row_hi_r,size(d1));
assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );
clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;

// type 2
[cA,cH,cV,cD]=dwt2(d1,'bior3.9');
[Lo_R,Hi_R]=wfilters('bior3.9','r');
dd0=idwt2(cA,cH,cV,cD,'bior3.9','mode','symw');
if exists('dd0')
dd1=idwt2(cA,cH,cV,cD,'bior3.9',size(d1),'mode','zpd');
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,'mode','sp0');
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1),'mode','ppd');

assert_checkalmostequal ( dd0 , dd2 , %eps );
assert_checkalmostequal ( dd1 , dd3 , %eps );
clear cA;
clear cH;
clear cV;
clear cD;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
end
// type 3

[cA,cH,cV,cD]=dwt2(d1,'bior3.9','mode','per');
[Lo_R,Hi_R]=wfilters('bior3.9','r');
dd0=idwt2(cA,cH,cV,cD,'bior3.9','mode','per');
if exists('dd0')
dd1=idwt2(cA,cH,cV,cD,'bior3.9',size(d1),'mode','per');
dd2=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,'mode','per');
dd3=idwt2(cA,cH,cV,cD,Lo_R,Hi_R,size(d1),'mode','per');
[r,c]=size(cA);
for i=1:c,
   col_low_low(:,i)=(conv(dyadup(cA(:,i)),Lo_R))';
   col_low_hi(:,i)=(conv(dyadup(cH(:,i)),Hi_R))';
   col_hi_low(:,i)=(conv(dyadup(cV(:,i)),Lo_R))';
   col_hi_hi(:,i)=(conv(dyadup(cD(:,i)),Hi_R))';
end
row_low=col_low_low+col_low_hi;
row_hi=col_hi_low+col_hi_hi;
[rc,cc]=size(row_low);
for i=1:rc,
   row_low_r(i,:)=conv(dyadup(row_low(i,:)),Lo_R);
   row_hi_r(i,:)=conv(dyadup(row_hi(i,:)),Hi_R);
end
ss=2*size(cA)-length(Lo_R)+2;
x0=wkeep(row_low_r+row_hi_r,2*size(cA));
x1=wkeep(row_low_r+row_hi_r,size(d1));

assert_checkalmostequal ( x0 , dd0 , %eps );
assert_checkalmostequal ( x1 , dd1 , %eps );
assert_checkalmostequal ( x0 , dd2 , %eps );
assert_checkalmostequal ( x1 , dd3 , %eps );

clear x0;
clear x1;
clear i;
clear rc;
clear cc;
clear row_low;
clear row_hi;
clear col_low_low;
clear col_low_hi;
clear col_hi_low;
clear col_hi_hi;
clear r;
clear c;
clear Lo_R;
clear Hi_R;
clear dd0;
clear dd1;
clear dd2;
clear dd3;
clear row_low_r;
clear row_hi_r;
end
