function E = wentropy(X,T,P)



      [nargout,nargin]=argn(0);
      if (nargin < 2),
      error ( '2 parameters are required' ); 
      end
      if (nargin < 3),
	  P=0;
      end
 
      select(convstr(T))
	case "shannon" 
	 E= - sum(sum(X.^2 .* log(%eps+X.^2)));
	case "log energy"
	 E=  sum(sum( log(%eps+X.^2)));
        case "threshold"
	 if (P<0)
	  error("P must be >=0!");
	 end;
	 E=sum(abs(X)>P);
	case "sure"
	 if (P<0)
	  error("P must be >=0!");
	 end;
	 E=sum(abs(X)>P)-sum(abs(X)>=P)+sum(min(X.^2,P^2));
        case "norm"
	 if (P<1)
	  error("P must be >=1!");
	 end;
	 E=sum( abs(X).^P );
        case "risk"
	  if (P<0)
	  error("P must be >=0!");
	 end;
	 E=sum(min(X.^2,P^2));
        else
         error(' T must be shannon, log energy, treshold, sure, norm or risk ');
      end;
       
endfunction
