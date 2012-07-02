/*
 * -------------------------------------------------------------------------
 * cwt.c -- continuous wavelet transform
 * SWT - Scilab wavelet toolbox
 * Copyright (C) 2005-2006  Roger Liu
 * Copyright (C) 20010-2012  Holger Nahrstaedt
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
#include "cwt.h"


cwt_identity ci[] = {
    {"sinus", REAL, SINUS, PSI_ONLY, -0.5, 0.5, 1, sinus},
	{"poisson", REAL, POISSON, PSI_ONLY,-10, 10, 1, poisson},
	{"mexh", REAL, MEXICAN_HAT,PSI_ONLY, -5, 5, 1.0, mexihat},
	{"morl",REAL,  MORLET, PSI_ONLY, -4, 4, 1.0, morlet},
	{"DOG", REAL, DOGAUSS, PSI_ONLY, -5, 5, 0.6455109, DOGauss},
	{"gaus1", REAL, GAUSS, PSI_ONLY, -5, 5, 1, Gaus1},
	{"gaus2", REAL, GAUSS, PSI_ONLY, -5, 5, 1, Gaus2},
	{"gaus3", REAL, GAUSS, PSI_ONLY, -5, 5, 1, Gaus3},
	{"gaus4", REAL, GAUSS, PSI_ONLY, -5, 5, 1, Gaus4},
	{"gaus5", REAL, GAUSS, PSI_ONLY, -5, 5, 1, Gaus5},
	{"gaus6", REAL, GAUSS, PSI_ONLY, -5, 5, 1, Gaus6},
	{"gaus7", REAL, GAUSS, PSI_ONLY, -5, 5, 1, Gaus7},
	{"gaus8", REAL, GAUSS, PSI_ONLY, -5, 5, 1, Gaus8},
	{"cmor",COMPLEX, CMORLET, PSI_ONLY, -8, 8, 1, cmorlet_packet},
	{"shan",  COMPLEX, SHANNON, PSI_ONLY, -20, 20, 1, shanwavf_packet},
	{"fbsp", COMPLEX, FBSP, PSI_ONLY, -20, 20, 1, fbspwavf_packet},
	{"cauchy", COMPLEX, CAUCHY, PSI_ONLY, -5, 5, 1, cauchy_packet},
	{"cgau1", COMPLEX, CGAUSS, PSI_ONLY, -5, 5, 1, cgau1_packet},
	{"cgau2", COMPLEX, CGAUSS, PSI_ONLY, -5, 5, 1, cgau2_packet},
	{"cgau3", COMPLEX, CGAUSS, PSI_ONLY, -5, 5, 1, cgau3_packet},
	{"cgau4", COMPLEX, CGAUSS, PSI_ONLY, -5, 5, 1, cgau4_packet},
	{"cgau5", COMPLEX, CGAUSS, PSI_ONLY, -5, 5, 1, cgau5_packet},
	{"cgau6", COMPLEX, CGAUSS, PSI_ONLY, -5, 5, 1, cgau6_packet},
	{"cgau7", COMPLEX, CGAUSS, PSI_ONLY, -5, 5, 1, cgau7_packet},
	{"cgau8", COMPLEX, CGAUSS, PSI_ONLY, -5, 5, 1, cgau8_packet}
};
int cwtIdentityNum = sizeof(ci)/sizeof(cwt_identity);

cwt_family cif[] = {
    {"sinus", "REAL", "SINUS"},
    {"poisson", "REAL", "POISSON"},
    {"mexh", "REAL", "MEXICAN_HAT"},
    {"morl","REAL",  "MORLET"},
    {"DOG", "REAL", "DOGAUSS"},
    {"cmor","COMPLEX", "CMORLET"},
    {"shan",  "COMPLEX", "SHANNON"},
    {"fbsp", "COMPLEX", "FBSP"},
    {"cauchy", "COMPLEX", "CAUCHY"},
    {"gaus", "REAL", "GAUSS"},    
    {"cgau", "COMPLEX", "CGAUSS"}
};
int cwtFamilyNum = sizeof(cif)/sizeof(cwt_family);

/*void haar(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;

	for(count=0;count<sigInLength;count++)
	{
		if (x[count]<0)
			psi[count] = -1/ys;
		else
			psi[count] = 1/ys;
	}
	return;
}*/

/*-------------------------------------------
 * Sinus Scale Filter Generation
 *-----------------------------------------*/

void sinus(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;

	for(count=0;count<sigInLength;count++)
		psi[count] = -1*sin(2*PI*x[count])/ys;
	return;
}

/*-------------------------------------------
 * Poisson Scale Filter Generation
 *-----------------------------------------*/

void poisson(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2;

    for(count=0;count<sigInLength;count++)
	{
		x2 = x[count] * x[count] * 4.0 / 9.0;
		psi[count] = (1-x2)/(PI*(1+x2)*ys);
	}
	return;
}

/*-------------------------------------------
 * Mexican Hat Scale Filter Generation
 *-----------------------------------------*/

void mexihat(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2, con;

    con = 2*sqrt(sqrt(PI))/sqrt(3);
	for(count=0;count<sigInLength;count++)
	{
		x2 = x[count] * x[count];
		psi[count] = (1-x2)*exp(-x2/2)*con/ys;
	}
	return;
}

/*-------------------------------------------
 * Morlet Scale Filter Generation
 *-----------------------------------------*/

void morlet(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2;

    for(count=0;count<sigInLength;count++)
	{
		x2 = x[count] * x[count];
		psi[count] = cos(5*x[count])*exp(-x2/2)/ys;
	}
	return;
}

/*-------------------------------------------
 * DoG Scale Filter Generation
 *-----------------------------------------*/

void DOGauss(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2;
	
    for(count=0;count<sigInLength;count++)
	{
		x2 = x[count] * x[count];
		psi[count] = exp(-x2/2)/ys-exp(-x2/8)/(2*ys);
	}
	return;
}

/*-------------------------------------------
 * Gauss Scale Filter Generation
 *-----------------------------------------*/

void Gauss(double *x, int sigInLength, double *psi, int sigOutLength, int n, double ys)
{
	switch (n) {
		case 1:
			Gaus1(x, sigInLength, psi, sigOutLength, ys);
			break;
		case 2:
			Gaus2(x, sigInLength, psi, sigOutLength, ys);
			break;
		case 3:
			Gaus3(x, sigInLength, psi, sigOutLength, ys);
			break;
		case 4:
			Gaus4(x, sigInLength, psi, sigOutLength, ys);
			break;
		case 5:
			Gaus5(x, sigInLength, psi, sigOutLength, ys);
			break;
		case 6:
			Gaus6(x, sigInLength, psi, sigOutLength, ys);
			break;
		case 7:
			Gaus7(x, sigInLength, psi, sigOutLength, ys);
			break;
		case 8:
			Gaus8(x, sigInLength, psi, sigOutLength, ys);
			break;
		default:
			break;
	}

	
	return;
}

void Gaus1(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2;
	
	for(count=0;count<sigInLength;count++)
	{
		x2 = x[count] * x[count];
		psi[count] = -2*x[count]*exp(-x2)/sqrt(sqrt(PI/2));
	}
	return;
}


void Gaus2(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2;
	
	for(count=0;count<sigInLength;count++)
	{
		x2 = x[count] * x[count];
		psi[count] = 2*(2*x2-1)*exp(-x2)/sqrt(3*sqrt(PI/2));
	}
	return;
}


void Gaus3(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2,x3;
	
	for(count=0;count<sigInLength;count++)
	{
		x2 = x[count] * x[count];
		x3 = x2 * x[count];
		psi[count] = -4*(2*x3-3*x[count])*exp(-x2)/sqrt(15*sqrt(PI/2));
	}
	return;
}


void Gaus4(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2,x4;
	
	for(count=0;count<sigInLength;count++)
	{
	 	x2 = x[count] * x[count];
		x4 = x2 * x2;
		psi[count] = 4*(-12*x2+4*x4+3)*exp(-x2)/sqrt(105*sqrt(PI/2));
	}
	return;
}


void Gaus5(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2,x3,x5;
	
	for(count=0;count<sigInLength;count++)
	{
         x2 = x[count] * x[count];
	 	 x3 = x2* x[count];
		 x5 = x3 * x2;
		psi[count] = 8*(-4*x5+20*x3-15*x[count])*exp(-x2)/sqrt(105*9*sqrt(PI/2));
	}
	return;
}


void Gaus6(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2,x4,x6;
	
	for(count=0;count<sigInLength;count++)
	{
        x2 = x[count] * x[count];
		x4 = x2 * x2;
		x6 = x4 * x2;
		psi[count] = 8*(8*x6-60*x4+90*x2-15)*exp(-x2)/sqrt(105*9*11*sqrt(PI/2));
	}
	return;
}


void Gaus7(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2,x3,x5,x7;
	
	for(count=0;count<sigInLength;count++)
	{
        x2 = x[count] * x[count];
	    x3 = x2 * x[count];
		x5 = x3 * x2;
		x7 = x5 * x2;
		psi[count] = 16*(-8*x7+84*x5-210*x3+105*x[count])*exp(-x2)/sqrt(105*9*11*13*sqrt(PI/2));
	}
	return;
}


void Gaus8(double *x, int sigInLength, double *psi, int sigOutLength, double ys)
{
	int count;
	double x2,x4,x6,x8;
	
	for(count=0;count<sigInLength;count++)
	{
        x2 = x[count] * x[count];
    	x4 = x2 * x2;
	    x6 = x4 * x2;
		x8 = x6 * x2;
		psi[count] = 16*(16*x8-224*x6+840*x4-840*x2+105)*exp(-x2)/sqrt(105*9*11*13*15*sqrt(PI/2));
	}
	return;
}

/*-------------------------------------------
 * Complex Gauss Scale Filter Generation
 *-----------------------------------------*/
void cgauss(double *x, int sigInLength, int p, double *psir, double *psii, int sigOutLength, double ys)
{

    switch (p) {
		case 1:
			cgau1(x, sigInLength, psir, psii, sigOutLength, ys);
			break;
		case 2:
			cgau2(x, sigInLength, psir, psii, sigOutLength, ys);
			break;
		case 3:
			cgau3(x, sigInLength, psir, psii, sigOutLength, ys);
			break;
		case 4:
			cgau4(x, sigInLength, psir, psii, sigOutLength, ys);
			break;
		case 5:
			cgau5(x, sigInLength, psir, psii, sigOutLength, ys);
			break;
		case 6:
			cgau6(x, sigInLength, psir, psii, sigOutLength, ys);
			break;
		case 7:
			cgau7(x, sigInLength, psir, psii, sigOutLength, ys);
			break;
		case 8:
			cgau8(x, sigInLength, psir, psii, sigOutLength, ys);
			break;
		default:
			break;
	}
	return;
}

void cgau1(double *x, int sigInLength, 
			double *psir, double *psii, 
			int sigOutLength, double ys)
{
    int count;
	double x2,cosx,sinx;

	for(count=0;count<sigInLength;count++)
	{
       x2 = x[count] * x[count];
	   cosx = cos(x[count]);
	   sinx = sin(x[count]);
	   psir[count]=(-2*x[count]*cosx-sinx)*exp(-x2)/sqrt(2*sqrt(PI/2));
	   psii[count]=(2*x[count]*sinx-cosx)*exp(-x2)/sqrt(2*sqrt(PI/2));
	}
	return;
}

void cgau1_packet(double *x, int sigInLength, 
			double *f, int sigOutLength, double ys)
{
	cgau1(x, sigInLength,f, f+sigInLength, sigOutLength,1);
    return;
}

void cgau2(double *x, int sigInLength, 
			double *psir, double *psii, 
			int sigOutLength, double ys)
{
    int count;
	double x2,cosx,sinx;

	for(count=0;count<sigInLength;count++)
	{
       x2 = x[count] * x[count];
	   cosx = cos(x[count]);
	   sinx = sin(x[count]);
	   psir[count]=(4*x2*cosx+4*x[count]*sinx-3*cosx)*exp(-x2)/sqrt(10*sqrt(PI/2));
	   psii[count]=(-4*x2*sinx+4*x[count]*cosx+3*sinx)*exp(-x2)/sqrt(10*sqrt(PI/2));
	}
	return;
}

void cgau2_packet(double *x, int sigInLength, 
			double *f, int sigOutLength, double ys)
{
	cgau2(x, sigInLength,f, f+sigInLength, sigOutLength,1);
    return;
}

void cgau3(double *x, int sigInLength, 
			double *psir, double *psii, 
			int sigOutLength, double ys)
{
    int count;
	double x2,cosx,sinx,x3;

	for(count=0;count<sigInLength;count++)
	{
       x2 = x[count] * x[count];
	   x3 = x[count] * x2;
	   cosx = cos(x[count]);
	   sinx = sin(x[count]);
	   psir[count]=(-8*x3*cosx-12*x2*sinx+18*x[count]*cosx+7*sinx)*exp(-x2)/sqrt(76*sqrt(PI/2));
	   psii[count]=(8*x3*sinx-12*x2*cosx-18*x[count]*sinx+7*cosx)*exp(-x2)/sqrt(76*sqrt(PI/2));
	}
	return;
}

void cgau3_packet(double *x, int sigInLength, 
			double *f, int sigOutLength, double ys)
{
	cgau3(x, sigInLength,f, f+sigInLength, sigOutLength,1);
    return;
}

void cgau4(double *x, int sigInLength, 
			double *psir, double *psii, 
			int sigOutLength, double ys)
{
    int count;
	double x2,cosx,sinx,x3,x4;

	for(count=0;count<sigInLength;count++)
	{
       x2 = x[count] * x[count];
	   x3 = x[count] * x2;
	   x4 = x2 * x2;
	   cosx = cos(x[count]);
	   sinx = sin(x[count]);
	   psir[count]=(16*x4*cosx+32*x3*sinx-72*x2*cosx-56*x[count]*sinx+25*cosx)*exp(-x2)/sqrt(774*sqrt(PI/2));
	   psii[count]=(-16*x4*sinx+32*x3*cosx+72*x2*sinx-56*x[count]*cosx-25*sinx)*exp(-x2)/sqrt(774*sqrt(PI/2));
	}
	return;
}

void cgau4_packet(double *x, int sigInLength, 
			double *f, int sigOutLength, double ys)
{
	cgau4(x, sigInLength,f, f+sigInLength, sigOutLength,1);
    return;
}

void cgau5(double *x, int sigInLength, 
			double *psir, double *psii, 
			int sigOutLength, double ys)
{
    int count;
	double x2,cosx,sinx,x3,x4,x5;

	for(count=0;count<sigInLength;count++)
	{
       x2 = x[count] * x[count];
	   x3 = x[count] * x2;
	   x4 = x2 * x2;
	   x5 = x2 * x3;
	   cosx = cos(x[count]);
	   sinx = sin(x[count]);
	   psir[count]=(-32*x5*cosx-80*x4*sinx+240*x3*cosx+280*x2*sinx-250*x[count]*cosx-81*sinx)*exp(-x2)/sqrt(29746*sqrt(PI/2));
	   psii[count]=(32*x5*sinx-80*x4*cosx-240*x3*sinx+280*x2*cosx+250*x[count]*sinx-81*cosx)*exp(-x2)/sqrt(29746*sqrt(PI/2));
	}
	return;
}

void cgau5_packet(double *x, int sigInLength, 
			double *f, int sigOutLength, double ys)
{
	cgau5(x, sigInLength,f, f+sigInLength, sigOutLength,1);
    return;
}

void cgau6(double *x, int sigInLength, 
			double *psir, double *psii, 
			int sigOutLength, double ys)
{
    int count;
	double x2,cosx,sinx,x3,x4,x5, x6;

	for(count=0;count<sigInLength;count++)
	{
       x2 = x[count] * x[count];
	   x3 = x[count] * x2;
	   x4 = x2 * x2;
	   x5 = x2 * x3;
	   x6 = x3 * x3;
	   cosx = cos(x[count]);
	   sinx = sin(x[count]);
	   psir[count]=(64*x6*cosx+192*x5*sinx-720*x4*cosx-1120*x3*sinx+1500*x2*cosx+972*x[count]*sinx-331*cosx)*exp(-x2)/sqrt(147322*sqrt(PI/2));
	   psii[count]=(-64*x6*sinx+192*x5*cosx+720*x4*sinx-1120*x3*cosx-1500*x2*sinx+972*x[count]*cosx+331*sinx)*exp(-x2)/sqrt(147322*sqrt(PI/2));
	}
	return;
}

void cgau6_packet(double *x, int sigInLength, 
			double *f, int sigOutLength, double ys)
{
	cgau6(x, sigInLength,f, f+sigInLength, sigOutLength,1);
    return;
}

void cgau7(double *x, int sigInLength, 
			double *psir, double *psii, 
			int sigOutLength, double ys)
{
    int count;
	double x2,cosx,sinx,x3,x4,x5, x6, x7;

	for(count=0;count<sigInLength;count++)
	{
       x2 = x[count] * x[count];
	   x3 = x[count] * x2;
	   x4 = x2 * x2;
	   x5 = x2 * x3;
	   x6 = x3 * x3;
	   x7 = x4 * x3;
	   cosx = cos(x[count]);
	   sinx = sin(x[count]);
	   psir[count]=(-128*x7*cosx-448*x6*sinx+2016*x5*cosx+3920*x4*sinx-7000*x3*cosx-6804*x2*sinx+4364*x[count]*cosx+1303*sinx)*exp(-x2)/sqrt(2099965*sqrt(PI/2));
	   psii[count]=(128*x7*sinx-448*x6*cosx-2016*x5*sinx+3920*x4*cosx+7000*x3*sinx-6804*x2*cosx-4364*x[count]*sinx+1303*cosx)*exp(-x2)/sqrt(2099965*sqrt(PI/2));
	}
	return;
}

void cgau7_packet(double *x, int sigInLength, 
			double *f, int sigOutLength, double ys)
{
	cgau7(x, sigInLength,f, f+sigInLength, sigOutLength,1);
    return;
}

void cgau8(double *x, int sigInLength, 
			double *psir, double *psii, 
			int sigOutLength, double ys)
{
    int count;
	double x2,cosx,sinx,x3,x4,x5, x6, x7, x8;

	for(count=0;count<sigInLength;count++)
	{
       x2 = x[count] * x[count];
	   x3 = x[count] * x2;
	   x4 = x2 * x2;
	   x5 = x2 * x3;
	   x6 = x3 * x3;
	   x7 = x4 * x3;
	   x8 = x4 * x4;
	   cosx = cos(x[count]);
	   sinx = sin(x[count]);
	   psir[count]=(256*x8*cosx+1024*x7*sinx-5376*x6*cosx-12544*x5*sinx+28000*x4*cosx+36288*x3*sinx-36532*x2*cosx-20578*x[count]*sinx+5667*cosx)*exp(-x2)/sqrt(42102009.0*sqrt(PI/2));
	   psii[count]=(-256*x8*sinx+1024*x7*cosx+5736*x6*sinx-12544*x5*cosx-28000*x4*sinx+36288*x3*cosx+36532*x2*sinx-20578*x[count]*cosx-5667*sinx)*exp(-x2)/sqrt(42102009.0*sqrt(PI/2));
	}
	return;
}

void cgau8_packet(double *x, int sigInLength, 
			double *f, int sigOutLength, double ys)
{
	cgau8(x, sigInLength,f, f+sigInLength, sigOutLength,1);
    return;
}


/*-------------------------------------------
 * Complex Morlet Scale Filter Generation
 *-----------------------------------------*/

void cmorlet(double *x, int sigInLength, 
			double fb, double fc, double *psir, double *psii, 
			int sigOutLength, double ys)
{
	int count;
	double con, econ, x2;

	con = 1/sqrt(PI*fb);
	for(count=0;count<sigInLength;count++)
	{
		x2=x[count]*x[count];
		econ = exp(-x2/fb); 
		psir[count] = cos(2*PI*fc*x[count])*econ*con/ys;
		psii[count] = sin(2*PI*fc*x[count])*econ*con/ys;
	}
	return;
}

void cmorlet_packet(double *x, int sigInLength, 
			double *f, int sigOutLength, double ys)
{
	int count;
	double con, econ, x2;

	con = 1/sqrt(PI);
	for(count=0;count<sigInLength;count++)
	{
		x2=x[count]*x[count];
		econ = exp(-x2); 
		f[count] = cos(2*PI*x[count])*econ*con/ys;
		f[count+sigInLength] = sin(2*PI*x[count])*econ*con/ys;
	}
	return;
}

/*-------------------------------------------
 * Shannon Scale Filter Generation
 *-----------------------------------------*/

void shanwavf(double *x, int sigInLength, 
			double fb, double fc, double *psir, double *psii, 
			int sigOutLength, double ys)
{
	int count;
	double con, econ;

	con = sqrt(fb);
	for(count=0;count<sigInLength;count++)
	{
		if (x[count] != 0)
			econ = sin(x[count]*fb*PI)/(x[count]*fb*PI); 
		else
			econ = 1;
		psir[count] = cos(2*PI*fc*x[count])*econ*con/ys;
		psii[count] = sin(2*PI*fc*x[count])*econ*con/ys;
	}
	return;
}

void shanwavf_packet(double *x, int sigInLength, 
			double *f, int sigOutLength, double ys)
{
	int count;
	double con, econ;

	con = 1;
	for(count=0;count<sigInLength;count++)
	{
		if (x[count] != 0)
			econ = sin(x[count]*PI)/(x[count]*PI); 
		else
			econ = 1;
		f[count] = cos(2*PI*x[count])*econ*con/ys;
		f[count+sigInLength] = sin(2*PI*x[count])*econ*con/ys;
	}
	return;
}

/*-------------------------------------------
 * Frequency B-Spline Scale Filter Generation
 *-----------------------------------------*/
void fbspwavf(double *x, int sigInLength,int m, 
			double fb, double fc, double *psir, double *psii, 
			int sigOutLength, double ys)
{
	int count, i;
	double con, econ;

	con = sqrt(fb);
	for(count=0;count<sigInLength;count++)
	{
		if (x[count] != 0)
			econ = sin(x[count]*fb*PI/m)/(x[count]*fb*PI/m); 
		else
			econ = 1;
        for(i=0;i<m;i++)
			econ*=econ;

		psir[count] = cos(2*PI*fc*x[count])*econ*con/ys;
		psii[count] = sin(2*PI*fc*x[count])*econ*con/ys;
	}
	return;
}

void fbspwavf_packet(double *x, int sigInLength, 
			double *f, int sigOutLength, double ys)
{
	int count;
	double con, econ;

	con = 1;
	for(count=0;count<sigInLength;count++)
	{
		if (x[count] != 0)
			econ = sin(x[count]*PI)/(x[count]*PI); 
		else
			econ = 1;
     	f[count] = cos(2*PI*x[count])*econ*con/ys;
		f[count+sigInLength] = sin(2*PI*x[count])*econ*con/ys;
	}
	return;
}

/*-------------------------------------------
 * Cauchy Scale Filter Generation
 *-----------------------------------------*/

void cauchy(double *x, int sigInLength, 
			double fb, double fc, double *psir, double *psii, 
			int sigOutLength, double ys)
{
	int count;
	double con;

	
	for(count=0;count<sigInLength;count++)
	{
		con=2*PI*x[count]*fc/fb;
		//econ = exp(-x2/fb); 
		psir[count] = (1.0-con*con)/((1.0+con*con)*(1.0+con*con)*ys);
		psii[count] = -2.0*con/((1.0+con*con)*(1.0+con*con)*ys);
	}
	return;
}

void cauchy_neo(double *x, int sigInLength, double *psir, double *psii, int sigOutLength, double ys)
{
    int count;
	double x2, con, econ;

	for (count=0;count<sigInLength;count++)
	{
		x2 = x[count]*x[count];
		con = 1-x2;
		econ = (1-x2)*(1-x2)+4*x2;
		psir[count]=con/(2*PI*econ*ys);
		psii[count]=2*x2/(2*PI*econ*ys);
	}
	return;
}

void cauchy_packet(double *x, int sigInLength, 
			double *f, int sigOutLength, double ys)
{
	//int count;
	//double con;

	//con = 1/sqrt(PI);
	/*for(count=0;count<sigInLength;count++)
	{
		con=PI*x[count];
		
		f[count] = (1.0-con*con)/((1.0+con*con)*(1.0+con*con)*ys);;
		f[count+sigInLength] = -2.0*con/((1.0+con*con)*(1.0+con*con)*ys);
	}*/
    cauchy_neo(x, sigInLength, f, f+sigInLength, sigOutLength, ys);
	return;
}

/*-------------------------------------------
 * Meyer Filter Generation
 *-----------------------------------------*/
void meyeraux(double x, double *y)
{
    double x4, x5, x6, x7;


	x4 = x*x*x*x;
	x5 = x4*x;
	x6 = x5*x;
	x7 = x6*x;

	*y = 35*x4-84*x5+70*x6-20*x7;
	return;
}


/*-------------------------------------------
 * CWT Utility
 *-----------------------------------------*/
void cwt_fun_parser(char *wname, int *ind)
{
    int count;

	*ind = -1;
	for(count=0;count<cwtIdentityNum;count++)
	{
		if (strcmp(wname,ci[count].wname) == 0)
		{
			*ind = count;
			break;
		}
	}
	return;
}


void full_range_scalef (char *wname, double *f, int sigOutLength)
{
    int level, ind, family, member, count, s1, s2, l;
	double one, *lowfltr, *hifltr, *trange;
	char d[2]="d";
    Func syn_fun, ana_fun;
    swt_wavelet pWaveStruct;

    level = 10;
	one = 1.0;
	wavelet_fun_parser (wname, &ind);
    if ((ind!=-1) && (wi[ind].rOrB==ORTH))
    {
      wavelet_parser(wname,&family,&member);
	  syn_fun = wi[ind].synthesis;
      (*syn_fun)(member, &pWaveStruct);
	  cwt_upcoef_len_cal (1, pWaveStruct.length, level, 
	       &s1, &s2);
	  l=1;
	  //l=(int)(floor((sigOutLength-s1)/2));
	  for(count=0;count<sigOutLength;count++)
		  f[count] = 0;
	  cwt_upcoef (&one, 1, pWaveStruct.pLowPass,
		  pWaveStruct.pHiPass, pWaveStruct.length, &(f[l]), 
	      s1, s1, d, level);
	  if ((family==COIFLETS) || (family==SYMLETS) || (family==DMEY))
	  {
		  for(count=0;count<sigOutLength;count++)
		      f[count] = -1*f[count];
	  }
	  for(count=0;count<sigOutLength;count++)
		      f[count] = f[count]*pow(sqrt(2),level);
	  filter_clear();
	  return;
	 }

	if ((ind!=-1) && (wi[ind].rOrB==BIORTH))
    {
      wavelet_parser(wname,&family,&member);
	  ana_fun = wi[ind].analysis;
      (*ana_fun)(member, &pWaveStruct);
	  cwt_upcoef_len_cal (1, pWaveStruct.length, level, 
	       &s1, &s2);
	  //l=(int)(floor((sigOutLength-s1)/2));
	  l=1;
   	  for(count=0;count<sigOutLength;count++)
	       f[count]=0;
	  lowfltr = malloc(pWaveStruct.length*sizeof(double));
	  hifltr = malloc(pWaveStruct.length*sizeof(double));
      wrev(pWaveStruct.pLowPass, pWaveStruct.length, lowfltr, pWaveStruct.length);
	  qmf_wrev(lowfltr,pWaveStruct.length,hifltr,pWaveStruct.length);
      cwt_upcoef (&one, 1, lowfltr, hifltr, pWaveStruct.length, &(f[l]), 
	      s1, s1, d, level);
	  free(lowfltr);
	  free(hifltr);
	  filter_clear();
      for(count=0;count<sigOutLength;count++)
	  	   f[count] = -1*f[count]*pow(sqrt(2),level);
	  return;
   }
   
	cwt_fun_parser(wname, &ind);
  if ((ind!=-1) && (ci[ind].realOrComplex==REAL))
  {
	  trange = malloc(sigOutLength*sizeof(double));
	  linspace(ci[ind].lb, ci[ind].ub, sigOutLength, trange, sigOutLength);
	  (*(ci[ind].scalef))(trange,sigOutLength,f,sigOutLength,ci[ind].cpsi);
	  free(trange);
	  return;
  }

  if ((ind!=-1) && (ci[ind].realOrComplex==COMPLEX))
  {
	  trange = malloc(sigOutLength*sizeof(double)/2);
	  linspace(ci[ind].lb, ci[ind].ub, sigOutLength/2, trange, sigOutLength/2);
	  (*(ci[ind].scalef))(trange,sigOutLength/2,f,sigOutLength,ci[ind].cpsi);
	  free(trange);
	  return;
  }

	return;
}

void cwt_len_cal (int sigInLength, int scale, int *sigOutLength, double *delta)
{
    *sigOutLength = scale;
	*delta = (double)(sigInLength-1)/(double)(scale-1);
	return;
}

void scale_real (double *f, int sigInLength, double delta, double *fout, int sigOutLength)
{
    int count;
	for(count=0;count<sigOutLength;count++)
	{
		fout[count] = f[(int)(floor(count*delta))];
		//sciprint("%d\n",count*delta+1);
	}
	return;
}

/*void scale_complex (double *f, int sigInLength, int delta, double *fout, int sigOutLength)
{

	return;
}*/

void cwt_conv_real (double *sigIn, int sigInLength, double *f, int filterLen, double *sigOut, int sigOutLength)
{
	int len;
    double *fTemp, *buf;

	len = sigInLength+filterLen-1;
	buf = malloc(len*sizeof(double));
    fTemp = malloc(filterLen*sizeof(double));

	wrev(f, filterLen, fTemp, filterLen);
	conv(sigIn,sigInLength,buf,len,fTemp,filterLen);
	free(fTemp);
    //for (i=1;i<len;i++)
	//	buf[i]=buf[i]-buf[i-1];
	wkeep_1D_center(buf,len,sigOut,sigOutLength);
	free(buf);
	return;
}

void cwt_iconv_real (double *sigIn, int sigInLength, double *f, int filterLen, double *sigOut, int sigOutLength)
{
	int len;
    double *fTemp, *buf;

	len = sigInLength+filterLen-1;
	buf = malloc(len*sizeof(double));
    fTemp = malloc(filterLen*sizeof(double));

	wrev(f, filterLen, fTemp, filterLen);
	//iconv(sigIn,sigInLength,buf,len,fTemp,filterLen);
	i_conv(sigIn,sigInLength,buf,len,fTemp,filterLen);
	free(fTemp);
	wkeep_1D_center(buf,len,sigOut,sigOutLength);
	free(buf);
	return;
}

void cwt_conv_complex (double *sigIn, int sigInLength, double *fr, double *fi, int filterLen, 
					   double *sigOutR, double *sigOutI, int sigOutLength)
{
    cwt_conv_real(sigIn,sigInLength,fr,filterLen,sigOutR,sigOutLength);
	cwt_conv_real(sigIn,sigInLength,fi,filterLen,sigOutI,sigOutLength);
	return;
}

void cwt_conv_complex_complex (double *a, double *b, int sigInLength,double *c, double *d,
							   int filterLen, double *sigOutR, double *sigOutI, int sigOutLength)
{
	int count;
    double *ac, *bd, *bc, *ad;
	ac = malloc(sigOutLength*sizeof(double));
	bd = malloc(sigOutLength*sizeof(double));
	bc = malloc(sigOutLength*sizeof(double));
	ad = malloc(sigOutLength*sizeof(double));

	cwt_conv_real(a,sigInLength,c,filterLen,ac,sigOutLength);
	cwt_conv_real(b,sigInLength,d,filterLen,bd,sigOutLength);
	cwt_conv_real(b,sigInLength,c,filterLen,bc,sigOutLength);
	cwt_conv_real(a,sigInLength,d,filterLen,ad,sigOutLength);

    for(count=0;count<sigOutLength;count++)
	{
		sigOutR[count]=ac[count]-bd[count];
		sigOutI[count]=bc[count]+ad[count];
	}
	free(ac);
	free(bd);
	free(bc);
	free(ad);

	return;
}


void
cwt_upcoef_len_cal (int sigInLength, int filterLen, int stride, 
		int *sigOutLength, int *sigOutLengthDefault)
{
  int count;
  *sigOutLength = sigInLength;
  *sigOutLengthDefault = sigInLength;
      for(count=0;count<stride;count++)
      {
	// original version
	*sigOutLengthDefault = 2*(*sigOutLengthDefault) + filterLen - 1;
	*sigOutLength = 2*(*sigOutLength) + filterLen - 2;
	
      }
  return;
}

void
cwt_upcoef (double *sigIn, int sigInLength, double *lowRe,double *hiRe, 
	int filterLen, double *sigOut, int sigOutLength, 
	int defaultLength, char *coefType, int step)
{
  int count, sigInLengthTemp, leng;
  double *sigInTemp, *sigOutTemp;

  // works with wavefun, cwt
  sigInLengthTemp = 2 * sigInLength + filterLen - 2;
   

  
  //sigInLengthTemp = 2 * sigInLength + filterLen - 1;
  sigInTemp = (double *) malloc(defaultLength*sizeof(double));
  
  if (strcmp(coefType,"a")==0)
  {
	  //sciprint("recognized\n");
// 	  printf("sigInLength %d, filterLen%d, sigInLengthTemp %d\n",sigInLength,filterLen,sigInLengthTemp);
	  idwt_approx_neo (sigIn, sigInLength, lowRe, filterLen, 
		 sigInTemp, sigInLengthTemp);
// 	  sciprint("recognized\n");
  }
  else
    idwt_detail_neo (sigIn, sigInLength, hiRe, filterLen, 
		 sigInTemp, sigInLengthTemp); 

  if (step > 1)
    {
      sigOutTemp = (double *) malloc(defaultLength*sizeof(double));
      for(count=0;count<defaultLength;count++)
	sigOutTemp[count] = 0;
      leng = sigInLengthTemp;
      for(count=0;count<(step-1);count++) //for cwt
	{
	  //printf("leng %d, filterLen%d, leng*2-filterLen+2 %d\n",leng,filterLen,leng*2-filterLen+2);
	  // original version
	  idwt_approx_neo (sigInTemp, leng, lowRe, filterLen,
	               sigOutTemp, leng*2+filterLen-2);
	  leng = leng*2+filterLen-2;
	  
	  verbatim_copy (sigOutTemp, leng, sigInTemp, leng);
	}
      sigInLengthTemp = leng;
      free(sigOutTemp);
    }

 
  wkeep_1D_center (sigInTemp, sigInLengthTemp, sigOut, sigOutLength);
  free(sigInTemp);
  return;
}

/*void real_scale (double lb, double ub, double scale, int length, double *f, double ys, RWScaleFunc w)
{
    int ns, count;
	double *x, *ft;

    ns = (int)(floor(ub*scale)-ceil(lb*scale));
	x = malloc(ns*sizeof(double));
	ft = malloc(ns*sizeof(double));
    for(count=0;count<ns;count++)
		x[count]=ceil(lb*scale)+count;
	(*w)(x, ns, ft, ns, ys);
	free(x);
    wextend_1D_center (ft, ns, f, length, ZPD);
	free(ft);
	return;
}*/


//void complex_scale(int lb, int ub, double scale, int length, double *fr, double *fi, double ys, WScaleFunc w)
//{
//	int ns, count;
//	double *x, *ftr, *fti;
//
  //  ns = (int)(floor(ub*scale)-ceil(lb*scale));
//	x = malloc(ns*sizeof(double));
//	ftr = malloc(ns*sizeof(double));
//	fti = malloc(ns*sizeof(double));
//    for(count=0;count<ns;count++)
//		x[count]=ceil(lb*scale)+count;
//	(*w)(x, ns, ft, ns, ys);
//	free(x);
//    wextend_1D_center (ft, ns, fr, length, ZPD);
//	free(ft);
//	return;
//}


//void cwt (double *sigIn, int sigInLength, 
//		  double *matrixOut, int matrixOutRow, 
//		  int matrixCol, char *wname, int *scale)
//{
//	int count;
//    double *scaleFunc;

//	cwt_fun_parser(wname,&ind);

//	scaleFunc = malloc(sigInLength*sizeof(double));
//	if (ci[count].realOrComplex == 0)
//	{
//	    for(count=0:count<matrixOutRow;count++)
//	    {
//		    ci[count].scalef(ci[ind].lb, ci[ind],ub, scale[count], sigInLength, scaleFunc,sqrt(scale[count]));
//		    iconv(sigIn,sigInLength,matrixOut+count*sigInLength,sigInLength,scaleFunc,sigInLength);
//      }
//	}
//	else
//	{
//
//	}
//    free(scaleFunc);
//
//	return;
//}