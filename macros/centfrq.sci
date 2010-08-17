//function [FREQ,XVAL,RECFREQ]=centfrq(wname,ITER)
function [FREQ]=centfrq(wname,ITER)

  [nargout,nargin]=argn(0);
  if (nargin < 1),
	error ( '1 parameter is required !' ); 
  end

  if (nargin<2)
    ITER=10;
  end;

clear XVAL;
clear psi;
clear phi;
try
[phi,psi,XVAL]=wavefun(wname,ITER);
catch 
end


if (~exists('XVAL'))
[psi,XVAL]=wavefun(wname,ITER);
end

if (~exists('XVAL'))
[phi1,psi,phi2,psi2,XVAL]=wavefun(wname,ITER);

end



sample_rate=1/(XVAL(2)-XVAL(1));
N=size(XVAL,'*');
f=sample_rate*(0:(N/2))/N;
n=size(f,'*');

y=fft(real(psi));
[value,place]=max(abs(y(1:n)));

//second max?
moremax=find(ceil(abs(y(1:n)))==ceil(real(value)));
if(size(moremax,'*')>1)
  FREQ=(f(moremax(1))+f(moremax($)))/2;
else
FREQ=f(place);
end;

//RECFREQ=sin(2*%pi*FREQ*XVAL);

endfunction