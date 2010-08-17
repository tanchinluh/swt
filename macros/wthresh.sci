function Y = wthresh(X,SORH,T)
// Soft or hard thresholding
//returns the soft (if SORH = 's') or hard (if SORH = 'h') thresholding o
//  Usage 
//    Y = wthresh(X,SORH,T)
//  Inputs 
//    X     Noisy Data 
//    SORH  'h' or 's'
//    T     Threshold
//  Outputs 
//    Y   


      [nargout,nargin]=argn(0);
      if (nargin < 3),
      error ( '3 parameters are required' ); 
      end
    if (convstr(SORH) ~= 'h' &  convstr(SORH) ~= 's') then
      error(' SORH must be either h or s');
    end;

	if (convstr(SORH)=='h') then
	  Y   = X .* (abs(X) > T);
        else
	  res = (abs(X) - T);
	  res = (res + abs(res))/2;
	  Y   = sign(X).*res;
	end;
endfunction
