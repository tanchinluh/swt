function cwtplot(coef,scales,cbar,f)



[nargout,nargin]=argn(0);
  if (nargin < 2),
    error ( '2 parameters are required !' ); 
  end

if (nargin<3)
  cbar=%t;
end

if (nargin< 4)
  f=figure();
  scf(f);clf(f);
end;


x=1:size(coef,2);
f.color_map =pinkcolormap(64);
grayplot(x,scales,abs(coef)');
if cbar
colorbar(0,max(abs(coef)));
end;
ylabel('Scale');
xlabel('time (or space)');
//xtitle('Continuous Transform, absolute coefficients');
a=gca();a.tight_limits="on";
endfunction