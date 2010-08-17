function [cwt_out] = cwt_wavelab(x,nvoice,wavelet,oct,rscale)
// CWT -- Continuous Wavelet Transform
//  Usage
//    cwt_out = CWT_Wavelab(x,nvoice,wavelet,oct,rscale)
//  Inputs
//    x        signal, dyadic length n=2^J, real-valued
//    nvoice   number of voices/octave
//    wavelet  string 'Gauss', 'DerGauss','Sombrero', 'Morlet'
//    octave   Default=2
//    rscale    Default=4
//  Outputs
//    cwt_out      matrix n by nscale where
//             nscale = nvoice .* noctave
//
//  Description
//    
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
	
	cwt_out = zeros(n,nscale);

	kscale  = 1;
//	scale   = 4;
//	scale = 16;

	for jo = 1:noctave,
	    for jv = 1:nvoice,
		   qscale = rscale .* (2^(jv/nvoice));
		   omega =  n .* xi ./ qscale ;
		   if strcmp(wavelet,'Gauss') then
				rwindow = exp(-omega.^2 ./2);
		elseif strcmp(wavelet,'DerGauss') then
                                rwindow = %i.*omega.*exp(-omega.^2 ./2);
		   elseif strcmp(wavelet,'Sombrero') then
				rwindow = (omega.^2) .* exp(-omega.^2 ./2);
		   elseif strcmp(wavelet,'Morlet') then
				rwindow = exp(-(omega - omega0).^2 ./2) - exp(-(omega.^2 + omega0.^2)/2);
		   end
		   // Renormalization
		   rwindow = rwindow ./ sqrt(qscale);
		   wxhat = rwindow .* xhat;
		   w    = mtlb_ifft(wxhat);
		   cwt_out(1:n,kscale) = real(w)';
		   kscale = kscale+1;
		end
		rscale  = rscale .*2;
    end 
// kscale = 1 -> low frequencies
// the matrix cwt_out is ordered from low to high frequencies 
//
// Originally created for WaveLab.701.
//
// Modified by Maureen Clerc and Jerome Kalifa, 1997
// clerc@cmapx.polytechnique.fr, kalifa@cmapx.polytechnique.fr
    
endfunction    
 
 
//
//  Part of Wavelab Version 850
//  Built Tue Jan  3 13:20:39 EST 2006
//  This is Copyrighted Material
//  For Copying permissions see COPYING.m
//  Comments? e-mail wavelab@stat.stanford.edu 
