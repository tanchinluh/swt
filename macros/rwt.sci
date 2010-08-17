function [rwt_out] = rwt(x,nvoice,wavelet,oct,rscale)
// RWT -- Real Wavelet Transform
//  Usage
//    rwt_out = rwt(x,nvoice,wavelet)
//  Inputs
//    x        signal, dyadic length n=2^J, real-valued
//    nvoice   number of voices/octave
//    wavelet  string 'Gauss', 'DerGauss','Sombrero', 'Morlet'
//  Outputs
//    rwt_out      matrix n by n where
//             nscale = nvoice .* noctave
//
//  Description
//     see sections 4.3.1 and 4.3.3 of Mallat's book  
//
	if argn(2)<4 then
		oct = 2;
		rscale = 4;
	end	
// preparation
	x = ShapeAsRow(x);
	n = length(x);
	xhat = mtlb_fft(x);
	xi   = [ (0: (n/2)) (((-n/2)+1):-1) ] .* (2*%pi/n);
	
// root
	omega0 = 5;
	
//	noctave = floor(log2(n))-2;
//	noctave = floor(log2(n))-1;
	noctave = floor(log2(n))-oct;
	nscale  = nvoice .* noctave;
	
	rwt_out = zeros(n,nscale);

	kscale  = 1;
//	scale   = 4;
//	scale = 16;

	for jo = 1:noctave,
	    for jv = 1:nvoice,
		   qscale = rscale .* (2^(jv/nvoice));
		   omega =  n .* xi ./ qscale ;
		   if mtlb_strcmp(wavelet,'Gauss') then
				rwindow = exp(-omega.^2 ./2);
		elseif mtlb_strcmp(wavelet,'DerGauss') then
                                rwindow = %i.*omega.*exp(-omega.^2 ./2);
		   elseif mtlb_strcmp(wavelet,'Sombrero') then
				rwindow = (omega.^2) .* exp(-omega.^2 ./2);
		   elseif mtlb_strcmp(wavelet,'Morlet') then
				rwindow = exp(-(omega - omega0).^2 ./2) - exp(-(omega.^2 + omega0.^2)/2);
		   end
//(  Renormalisation comme dans le bouquin )
		   rwindow = rwindow ./ sqrt(qscale);
		   rwhat = rwindow .* xhat;
		   w    = mtlb_ifft(rwhat);
		   rwt_out(1:n,kscale) = real(w)';
		   kscale = kscale+1;
		end
		rscale  = rscale .*2;
    end 
// kscale = 1 -> low frequencies
// the matrix rwt is ordered from low to high frequencies 

// Written by Maureen Clerc and Jerome Kalifa, 1997
// clerc@cmapx.polytechnique.fr, kalifa@cmapx.polytechnique.fr

   
    
endfunction 
 
//
//  Part of Wavelab Version 850
//  Built Tue Jan  3 13:20:39 EST 2006
//  This is Copyrighted Material
//  For Copying permissions see COPYING.m
//  Comments? e-mail wavelab@stat.stanford.edu 
