function [F]=scal2frq(A,wname,DELTA)

  [nargout,nargin]=argn(0);
  if (nargin < 2),
	error ( '2 parameters are required !' ); 
  end

  if (nargin<3)
    DELTA=1;
  end;

Fc=centfrq(wname);
F=Fc ./ (A .* DELTA);

//RECFREQ=sin(2*%pi*FREQ*XVAL);

endfunction