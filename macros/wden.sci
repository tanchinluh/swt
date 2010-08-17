function [XD,CXD,LXD] = wden(C,L,TPTR,SORH,SCAL,N,wname)
// Noisy wavelet test data
//  Usage
//X = wnoise(FUN,N)
//[X,XN] = wnoise(FUN,N,SQRT_SNR)
//[X,XN] = wnoise(FUN,N,SQRT_SNR,INIT)

//  Inputs
//   FUN = 1     or 	'blocks'
//    C,L   output of the wavedec function   
//  Outputs
//    STDC  estimation of sigma
//
//  Description
//    The estimator used is Median Absolute Deviation / 0.6745, well suited for zero mean Gaussian white noise in the de-noising one-dimensional model 
//

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
 
