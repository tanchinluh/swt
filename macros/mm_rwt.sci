function [maxmap] = mm_rwt(rwt_in,par)
// MM_RWT -- Modulus Maxima of a Real Wavelet Transform
//  Usage
//    maxmap = MM_RWT(rwt,par)
//  Inputs
//    rwt    Output of RWT
//    par    optional. If present, keep thresholds only
//           above a certain value. default = 1000
//  Outputs
//    maxmap binary array indicating presence of max or not
//
//  Description
//    Used to calculate fractal exponents etc. 
//
	if argn(2)<2 then
		par = 1000;
	end 
	sz = size(rwt_in);
	nscale = sz(2);
	n = sz(1);
	
	maxmap = mtlb_zeros(sz);

	t      = 1:n;
	tplus  = rshift(t);
	tminus = lshift(t);
	rwt_in    = abs(rwt_in);
		
	for k=1:nscale,
	     localmax =  rwt_in(t,k) > rwt_in(tplus,k) & rwt_in(t,k) > rwt_in(tminus,k) ;
		y =  localmax(t) .* rwt_in(:,k);
		maxy = mtlb_max(y);
		maxmap(t,k) = (y >= maxy/par);
	end
	
// Written by Maureen Clerc and Jerome Kalifa, 1997
// clerc@cmapx.polytechnique.fr, kalifa@cmapx.polytechnique.fr

   endfunction
    
 
 
//
//  Part of Wavelab Version 850
//  Built Tue Jan  3 13:20:39 EST 2006
//  This is Copyrighted Material
//  For Copying permissions see COPYING.m
//  Comments? e-mail wavelab@stat.stanford.edu 
