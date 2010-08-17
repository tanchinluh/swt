function [THR,SORH,KEEPAPP,CRIT] = ddencmp(IN1,IN2,X)
// Default values for de-noising or compression
//  Usage
//   [THR,SORH,KEEPAPP,CRIT] = ddencmp(IN1,IN2,X)
//   [THR,SORH,KEEPAPP] = ddencmp(IN1,'wv',X)
//   [THR,SORH,KEEPAPP,CRIT] = ddencmp(IN1,'wp',X)
//  Inputs
//    IN1 is 'den' for de-noising or 'cmp' for compression.
//    IN2 is 'wv' for wavelet or 'wp' for wavelet packet.
//  Outputs
//    
//
//  Description
//     
//

  [nargout,nargin]=argn(0);
  if (nargin < 3),
	error ( '3 parameters are required' ); 
  end

  if ~(convstr(IN1)=='den' | convstr(IN1)=='cmp')
	    error("IN1 must be den or cmp!");
 end;
  if ~(convstr(IN2)=='wv' | convstr(IN2)=='wp')
	    error("IN2 must be wv or wp!");
 end;
 CRIT='threshold';
n=prod(size(X));

  if (convstr(IN1)=='den' & convstr(IN2)=='wv') then
	SORH='s';
	KEEPAPP=1;
       [C,L]=wavedec(X,1,'db3');
       s=wnoisest(C,L,1);
        s=median(abs(X))/0.6745;
	THR = sqrt(2*log(n)) * s;
  elseif (convstr(IN1)=='cmp' & convstr(IN2)=='wv') then
        SORH='h';
	KEEPAPP=1;
	[C,L]=wavedec(X,1,'db3');
	THR= median(abs(C(L(1)+1:$)));
	if (THR==0),
	  THR= 0.05 * max(abs(C(L(1)+1:$)));
	end;
  elseif (convstr(IN1)=='den' & convstr(IN2)=='wp') then
    SORH='h';
    KEEPAPP=1;
     CRIT='sure';
    THR = sqrt(2*log(n*log(n)/log(2)));
 elseif (convstr(IN1)=='cmp' & convstr(IN2)=='wp') then
    SORH='h';
    KEEPAPP=1;
    THR = median(abs(X));
  end;
    
 endfunction
 
