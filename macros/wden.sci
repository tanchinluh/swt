function [XD,CXD,LXD] = wden(C,L,TPTR,SORH,SCAL,N,wname)
//Automatic 1-D de-noising
//Calling Sequence
//[XD,CXD,LXD] = wden(X,TPTR,SORH,SCAL,N,wname)
//[XD,CXD,LXD] = wden(C,L,TPTR,SORH,SCAL,N,wname)
//Parameters
//X: input vector
//C: coefficent array
//L: length array
//wname: wavelet name
//N: decompostion level
//SCAL: threshold rescaling
//: 'one' for no rescaling
//: 'sln' for rescaling using a single estimation of level noise based on first-level coefficients
//: 'mln' for rescaling done using level-dependent estimation of level noise
//SORH: ('s' or 'h') soft or hard thresholding
//TPTR: threshold selection rule
//:'rigrsure' uses the principle of Stein's Unbiased Risk.
//: 'heursure' is an heuristic variant of the first option.
//: 'sqtwolog' for universal threshold
//: 'minimaxi' for minimax thresholding
//CXD: de-noised coefficent array
//LXD: de-noised length array
//XD: de-noised signal
//Description
//wden performs an automatic de-noising process of a one-dimensional signal using wavelets.
//Examples
//snr = 3; init = 2055615866; 
//[xref,x] = wnoise(3,11,snr,init);
//lev = 5;
//xd = wden(x,'heursure','s','one',lev,'sym8');
// 
// See also
//thselect
//wthresh
//wavedec 
//Authors
//Holger Nahrstaedt


[nargout,nargin]=argn(0);
  if (nargin < 6 | nargin>7),
	error ( '6 or 7 parameters are required' ); 
  end
 
  if (nargin==6)
    X=C;
    wname=N;
    N=SCAL;
    SCAL=SORH;
    SORH=TPTR;
    TPTR=L;
    clear C,L;
    if (min(size(X))>1)
      error("This is a 1-D denoising funtion. X must be a vector!");
    end
    [C,L] = wavedec(X,N,wname);
  end
D=list();
for n=1:N
 D(n)=detcoef(C,L,n);
end;

  THR = thselect(X,TPTR);

select convstr(SCAL),
   case 'one'
      sigma=ones(1,N);
   case 'sln'
       sigma=ones(1,N)*wnoisest(C,L,1);
      // THR = thselect(D(1),TPTR);
   case 'mln'
       sigma=wnoisest(C,L,1:N);
end
 // THR = thselect(X,TPTR);

  CXD=C;
  LXD=L;
 i=0;
for n=N:-1:1
  i=i+1;

  CXD(sum(L(1:i))+(1:L(i+1)))=wthresh(D(n),SORH,thselect(D(n)/sigma(n),TPTR)*sigma(n));
end



 XD=waverec(CXD,LXD,wname);

 endfunction
 
