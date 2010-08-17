/*
 * -------------------------------------------------------------------------
 * cowt.h -- Declarations for complex wavelet functions.
 * SWT - Scilab wavelet toolbox
 * Copyright (C) 2005-2006  Roger Liu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * -------------------------------------------------------------------------
 */

#include "swt_common.h"

/*********************************************
 * Function Prototype
 ********************************************/
extern 
void
cowavedec (double *sigIn, int sigInLength, double *sigOutR,
	   double *sigOutI, int sigOutLength, 
	   double *lowDTree1S1, double *hiDTree1S1,
	   double *lowDTree2S1, double *hiDTree2S1,
	   double *lowDTree1S2, double *hiDTree1S2, 
	   double *lowDTree2S2, double *hiDTree2S2,
	   int filterLen, int *waveDecLengthArray, 
	   int lengthArrayLengh, int stride, extend_method extMethod);

extern 
void
cowaverec (double *sigInR, double *sigInI, int sigInLength, 
	   double *sigOut, int sigOutLength, 
	   double *lowRTree1S1, double *hiRTree1S1,
	   double *lowRTree2S1, double *hiRTree2S1,
	   double *lowRTree1S2, double *hiRTree1S2, 
	   double *lowRTree2S2, double *hiRTree2S2,
	   int filterLen, int *waveDecLengthArray,
	   int lengthArraylength, int stride, 
	   extend_method extMethod);

extern void
cowavedec2 (double *matrixIn, int matrixInRow, int matrixInCol,
	    double *lowDTree1S1, double *hiDTree1S1,
	    double *lowDTree1S2, double *hiDTree1S2,
	    int filterLen, int *pLen, double *coef, 
	    int sigOutLength, int stride, extend_method extMethod);

extern void
cowavedec2a (double *matrixIn, int matrixInRow, int matrixInCol,
	     double *lowDTree1S1R, double *hiDTree1S1R,
	     double *lowDTree1S1C, double *hiDTree1S1C,
	     double *lowDTree1S2R, double *hiDTree1S2R,
	     double *lowDTree1S2C, double *hiDTree1S2C,
	     int filterLen, int *pLen, double *coef, 
	     int sigOutLength, int stride, extend_method extMethod);

extern void
cowaverec2 (double *coef, int sigInLength, 
	    double *lowRTree1S1, double *hiRTree1S1,
	    double *lowRTree1S2, double *hiRTree1S2,
	    int filterLen, double *matrixOut, int matrixOutRow, 
	    int matrixOutCol, int *pLen, int stride, 
	    extend_method extMethod);

extern void
cowaverec2a (double *coef, int sigInLength, 
	     double *lowRTree1S1R, double *hiRTree1S1R,
	     double *lowRTree1S1C, double *hiRTree1S1C,
	     double *lowRTree1S2R, double *hiRTree1S2R,
	     double *lowRTree1S2C, double *hiRTree1S2C,
	     int filterLen, double *matrixOut, int matrixOutRow, 
	     int matrixOutCol, int *pLen, int stride, 
	     extend_method extMethod);

extern void
copmd (double *matrixInR, double *matrixInI, int sigInLength,
      int InRow, int InCol, double *matrixOutR, double *matrixOutI);

extern void
copmr (double *matrixInR, double *matrixInI, int sigInLength,
      int InRow, int InCol, double *matrixOutR, double *matrixOutI);
