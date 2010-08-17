function THR = thselect(X,TPTR)

//thselect is a one-dimensional de-noising oriented function.

//THR = thselect(X,TPTR) returns threshold X-adapted value using selection rule defined by string TPTR.

//Available selection rules are

//    * TPTR = 'rigrsure', adaptive threshold selection using principle of Stein's Unbiased Risk Estimate.
//    * TPTR = 'heursure', heuristic variant of the first option.
//    * TPTR = 'sqtwolog', threshold is sqrt(2*log(length(X))).
//    * TPTR = 'minimaxi', minimax thresholding.

  [nargout,nargin]=argn(0);
  if (nargin < 2),
	error ( '2 parameters are required' ); 
  end

 if ~(convstr(TPTR)=='rigrsure' | convstr(TPTR)=='heursure' | convstr(TPTR)=='sqtwolog' | convstr(TPTR)=='minimaxi')
	    error("TPTR must be rigrsure, heursure, sqtwolog or minimaxi!");
 end;

  select convstr(TPTR),
    case 'rigrsure'
      THR  = ValSUREThresh(X);
    case 'heursure'
	[n,j] = dyadlength(X);
	magic = sqrt(2*log(n));
	eta = (norm(X).^2 - n)/n;
	crit = j^(1.5)/sqrt(n);
	if eta < crit
		THR=magic;
	else
		THR  = min(ValSUREThresh(X), magic);
	end
    case 'sqtwolog'
	[n,j] = dyadlength(X);
	THR=sqrt(2*log(n));
    case 'minimaxi'
	lamlist = [0  0 0 0 0 1.27 1.474 1.669 1.860 2.048 2.232 2.414 2.594 2.773 2.952 3.131 3.310 3.49 3.67 3.85 4.03 4.21];
	[n,j] = dyadlength(X);
	THR = lamlist(j);
   else
      error("wrong parameter for threshhold selection!");
  end

endfunction


function THR = ValSUREThresh(X)
// ValSUREThresh -- Adaptive Threshold Selection Using Principle of SURE
//  Usage 
//    thresh = ValSUREThresh(y)
//  Inputs 
//    y        Noisy Data with Std. Deviation = 1
//  Outputs 
//    thresh   Value of Threshold
//
//  Description
//    SURE referes to Stein's Unbiased Risk Estimate.
//
	X=X(:);
	//a = mtlb_sort(abs(X)).^2 ;
	a = gsort(abs(X),'g','i').^2;
	b = cumsum(a,'m');
	n = length(X);
	c = linspace(n-1,0,n);
	s = b+c(:).*a;
	risk = (n - ( 2 .* (1:n ))  + s')'/n;
	[guess,ibest] = min(risk);
	THR = sqrt(a(ibest));

//
// Copyright (c) 1993-5.  Jonathan Buckheit, David Donoho and Iain Johnstone
//
    
endfunction

function [n,J] = dyadlength(x)
// dyadlength -- Find length and dyadic length of array
//  Usage
//    [n,J] = dyadlength(x)
//  Inputs
//    x    array of length n = 2^J (hopefully)
//  Outputs
//    n    length(x)
//    J    least power of two greater than n
//
//

  n = length(x) ;
  J = ceil(log(n)/log(2));
//   if 2^J ~= n ,
//       disp('Warning in dyadlength: n != 2^J')
//   end
    
    endfunction
 