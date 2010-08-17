function [STDC] = wnoisest(C,L,S)
// Estimate noise of 1-D wavelet coefficients
//  Usage
//   STDC = wnoisest(C,L,S) 
//   STDC = wnoisest(C,L,S) returns estimates of the detail coefficients' standard deviation for levels contained in the input vector S. [
//  Inputs
//    x     1-d signal
//    C,L   output of the wavedec function   
//  Outputs
//    STDC  estimation of sigma
//
//  Description
//    The estimator used is Median Absolute Deviation / 0.6745, well suited for zero mean Gaussian white noise in the de-noising one-dimensional model 
//
	maxLevel=length(L)-2;
	if(maxLevel<length(S))
	    error("C,L does not contain so much levels. reduce S!");
	end;
	
	STDC=zeros(1,length(S));
	for level=S

	 //STDC(level) = median(abs(C(sum(l(1:(maxLevel-level+1)))+1:sum(l(1:(maxLevel-level+2))))))/.6745;
         STDC(level)= median(abs(detcoef(C,L,level)))/.6745;
	//STDC(level) = median(abs(C(sum(L(1:level))+(1:L(level+1)))))/.6745;
	end;
    
 endfunction
 
