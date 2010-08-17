mode(7);clc; figure(); clf;
//        generate sine wave
x=[1:512];
y=sin(2*%pi*x/64);
//        increase memory size
stacksize(10000000);
//        continous wavelet transform
coef=cwt(y,1:128,'DOG');
mesh(coef);
//        pseduo color display, SIVP needed
//        colormap obtaining
//cmap=jetcolormap(512);
//        matrix normalization
//n=wcodemat(abs(coef),512);
//        encode matrix
figure();clf;cwtplot(coef,1:128);
//w=ind2rgb(n,cmap);
//        display the matrix
//imshow(w); 


