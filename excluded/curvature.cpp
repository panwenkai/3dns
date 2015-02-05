/*//  ____________________________________________________________________________
// |   3d_Curvature.h    Curvature header
// |   (C) 1998  Columbia University, MSME
// |____________________________________________________________________________
//______________________________________________________________________________
// Public Functions
//______________________________________________________________________________
void Curvature (const int i, const int j, const int k, double &sCurv, double &dTGT);
*///  ____________________________________________________________________________
// |   3d_Curvature.cpp    Curvature of interface routines
// |     v1.0  For 3DPC v3   A.B.Limanov     10-29-99
// |   (C) 1999  Columbia University
// |____________________________________________________________________________
#include <conio.h>
#include "3DNS.h"
#include "phase.h"
#include "math.h"
#include "curvature.h"
#define GT 1E-8 // Capillary constant GT = Hlc/Slc [1/cm] for c-Si
#define ADIABAT 1 //      ----|----               angle ninty degrees
#define FLOW 2    //Interface boundary condition; same angle 
#define PERIODIC 3  //      ----|----               first depends on last
#define I_EDGE ADIABAT
#define J_EDGE ADIABAT //ADIBATIC NEUMANN BOUNDARY CONDITION
#define K_EDGE ADIABAT //PERIODIC
#define dia(a,b) sqrt(a * a + b * b)
#define pi PI
//______________________________________________________________________________
// Private functions
//______________________________________________________________________________
void DirectTo_kij100(int indi, int indj,           double &line_s, double &curv, int &kneib);
void DirectTo_ijk100(int indi, int indj,           double &line_s, double &curv, int &kneib);
void DirectTo_jki100(int indi, int indj,           double &line_s, double &curv, int &kneib);
void DirectTo_ijk110(int indi, int indd, int indk, double &line_s, double &curv, int &kneib);
void DirectTo_jki110(int indi, int indd, int indk, double &line_s, double &curv, int &kneib);
void DirectTo_kij110(int indi, int indd, int indk, double &line_s, double &curv, int &kneib);
void FindBoundary100(int i, int j, int k, int indi, int indj, int ijk,
                     double dX, double xFrac, int i_first, int i_last, int i_mode,
                     double dY, double yFrac, int j_first, int j_last, int j_mode, int &kneib);
void FindBoundary110(int i, int j, int k, int indi, int indd, int indj, int ijk,
                     double dX, double xFrac, int i_first, int i_last, int i_mode,
                                              int j_first, int j_last, int j_mode,
                     double dZ, double KPos, int k_first, int k_last, int k_mode, int &kneib);
void FindNeighbor100 (int il, int jl, int kl,           int ijk, int i_first, int i_last,
                                                                 int j_first, int j_last, int &kneib);
void FindNeighbor110 (int il, int jl, int kl, int indd, int ijk, int i_first, int i_last,
                                                                 int j_first, int j_last,
                                                                 int k_first, int k_last, int &kneib);
void ChekNeibNumber(int &kneibSw, int &kneib, int nodeD, int in, int jn, int kn, int uncode);
double Curv2d (int ci, int cj, double DX, double DY, double &line_s, double &alp12);
int Find2NeighborsFrom3();
int Si_bottom (void);
double Mediana(double ma, double mb, double mc);
/*
Note neib stands for 'neighbor'
*/
//______________________________________________________________________________
// Private variables
//______________________________________________________________________________
unsigned short j_LAST = 0; //last+1 Si column; local
int uncode, uncode1, uncode2;
unsigned short i, j, k; //modified to unsigned short
int cri, crj, crk, indd, nodeA, nodeB, nodeC;
int in1, jn1, kn1,   in2, jn2, kn2,   in3, jn3, kn3,   in4, jn4, kn4; 
double x1, yy1,   x2, y2;
double alp1, alp2,   alpX, alpY, alpZ;  //alpX is an average one for X-norm.plane
int source;

/*****************
    Si_bottom: 
        Iterates through all J-nodes until it finds a J-node where phase
        transitions are forbidden. Returns that node.
    
    Possible complications: 
        What if the nodes in the J-region isn't continuous?
        What about for back-side irradiation? It seems that j counts from 0.
********************/
int Si_bottom (void)
{
	int j = J_FIRST;
	for (; j < J_LAST; j++)
		if ( j > 0 && (Geometry.canChangeJ[j] != Geometry.canChangeJ[j-1]))
			return(j); 
    return(J_LAST);
		ErrorMsg("can't find last Si film row");
	return(j);
}

/*******************************
    Curvature:
        Find neighboring slush cells.

		ii [int]: i-index of the slush node's direction
		jj [int]: j-index of the slush node's direction
		kk [int]: k-index of the slush node's direction
		sCurv [double] (pass-by-reference): ??? curvature of the slush cell ???
		dTGT [double] (pass-by-reference): ??? no idea wtf this is...???
*********************************/
void Curvature (const int ii, const int jj, const int kk, double &sCurv, double &dTGT)
{
	int kneib1, kneib2, kneib3;
	double curvS, curv1, curv2, curv3;
	double line1_s, line2_s, line3_s;
	alp1 = 0.0; alp2 = 0.0; 
	alpX = -1.0; alpY = -1.0; alpZ = -1.0;   //initial (no/a) values
	in2 = -1;                          //in order to find circle-boundary point
	kneib1 = 0; kneib2 = 0; kneib3 = 0;
	i = ii; j = jj; k = kk;
	cri = crystalI[SlushDirection[i][j][k]]; 
	crj = crystalJ[SlushDirection[i][j][k]]; 
	crk = crystalK[SlushDirection[i][j][k]];
	if (j_LAST == 0) 
    {
		j_LAST = Si_bottom(); 
		printf("I_F=%d j_F=%d K_F=%d I_L=%d j_L=%d K_L=%d\n", I_FIRST, J_FIRST, K_FIRST, I_LAST, j_LAST, K_LAST);
	}
	nodeA = IJKToIndex(i, j, k);

	if (abs(cri) == 1 && crj == 0 && crk == 0)//------- i=1 -----------------------(100)-set
	{
		DirectTo_kij100(crk, cri, line1_s, curv1, kneib1);
		if (kneib1 != 2) 
			line1_s = A(DelZ);
		DirectTo_ijk100(cri, crj, line2_s, curv2, kneib2);
		if (kneib2 != 2) 
			line2_s = A(DelY);
		curvS = curv1 + curv2;
		sCurv = line1_s * line2_s;
 	}
	else if ((cri == 0) && (abs(crj) == 1) && (crk == 0))// j=1
	{
		DirectTo_ijk100(cri, crj, line1_s, curv1, kneib1);
		if (kneib1 != 2) line1_s = A(DelX);
		DirectTo_jki100(crj, crk, line2_s, curv2, kneib2);
		if (kneib2 != 2) line2_s = A(DelZ);
		curvS = curv1 + curv2;
		sCurv = line1_s * line2_s;
	}
	else if ((cri == 0) && (crj == 0) && (abs(crk) == 1))// k=1
	{
		DirectTo_jki100(crj, crk, line1_s, curv1, kneib1);
		if (kneib1 != 2) line1_s = A(DelY);
		DirectTo_kij100(crk, cri, line2_s, curv2, kneib2);
		if (kneib2 != 2) line2_s = A(DelX);
		curvS = curv1 + curv2;
		sCurv = line1_s * line2_s;
	}
	else if ((crk == 0) && (abs(cri) == 1) && (abs(crj) == 1))//(110)---- k=0 ---------
	{
		DirectTo_ijk100(cri,                            crj, line1_s, curv1, kneib1);
		if (kneib1 != 2) line1_s = dia(A(DelX),A(DelY));
//		printf("i%2d j%2d k%2d  %2d %2d %2d  k1%1d x1=%5.2f y1=%5.2f x2=%5.2f y2=%5.2f DT1=%5.3f L1=%5.2f dT=%5.3f \n", 		i, j, k, cri, crj, crk, kneib1, x1*1e7,yy1*1e7,x2*1e7,y2*1e7, (1.1e-8 * 1685 * curv1), line1_s * 1e7, 1685 - A(T)); getche();
		DirectTo_ijk110(cri, crj, crk, line2_s, curv2, kneib2);
		curvS = curv1 + curv2;
		sCurv = line1_s * line2_s;
//		printf("i%2d j%2d k%2d  %2d %2d %2d  k1%1d x1=%5.2f y1=%5.2f x2=%5.2f y2=%5.2f DT1=%5.3f L1=%5.2f ! DT1=%5.3f L1=%5.2f \n", 		i, j, k, cri, crj, crk, kneib2, x1*1e7,yy1*1e7,x2*1e7,y2*1e7, (1.1e-8 * 1685 * curv2), line2_s * 1e7, (1.1e-8 * 1685 * curvS), sCurv * 1e14); getche();
	}
	else if ((crj == 0) && (abs(cri) == 1) && (abs(crk) == 1))//(101)---- j=0
	{
		DirectTo_kij100(crk,                            cri, line1_s, curv1, kneib1);
		if (kneib1 != 2) line1_s = dia(A(DelZ),A(DelX));
		DirectTo_kij110(crk, cri, crj, line2_s, curv2, kneib2);
		curvS = curv1 + curv2;
		sCurv = line1_s * line2_s;
	}
	else if ((cri == 0) && (abs(crj) == 1) && (abs(crk) == 1))//(011)---- i=0
	{
		DirectTo_jki100(crj,                            crk, line1_s, curv1, kneib1);
		if (kneib1 != 2) line1_s = dia(A(DelY),A(DelZ));
		DirectTo_jki110(crj, crk, cri, line2_s, curv2, kneib2);
		curvS = curv1 + curv2;
		sCurv = line1_s * line2_s;
	}
	else if ((abs(cri) == 1) && (abs(crj) == 1) && (abs(crk) == 1))//(111)-set
	{
		DirectTo_ijk110(cri, crj, crk, line1_s, curv1, kneib1);
		DirectTo_jki110(crj, crk, cri, line2_s, curv2, kneib2);
		DirectTo_kij110(crk, cri, crj, line3_s, curv3, kneib3);
		curvS = 0.6666667 * (curv1 + curv2 + curv3); 
		if ((alpZ <= alpY) && (alpZ <= alpX)) 	
		{	
			printf("XY %2d %2d %2d  alpX=%7.4f alpY=%7.4f alpZ=%7.4f\n", i, j, k, alpX*180/pi, alpY*180/pi, alpZ*180/pi);
			kneib2 = 0;
			DirectTo_ijk100(cri, crj, line2_s, curv2, kneib2);
			if (kneib2 != 2) line2_s = A(DelXY);
		}
		else if ((alpX <= alpY) && (alpX <= alpZ)) 	
		{	
			printf("YZ %2d %2d %2d  alpX=%7.4f alpY=%7.4f alpZ=%7.4f\n", i, j, k, alpX*180/pi, alpY*180/pi, alpZ*180/pi);
			line1_s = line2_s; kneib1 = kneib2; kneib2 = 0;
			DirectTo_jki100(crj, crk, line2_s, curv2, kneib2);
			if (kneib2 != 2) line2_s = A(DelYZ);
		}
		else if ((alpY <= alpX) && (alpY <= alpZ)) 	
		{	
			printf("ZX %2d %2d %2d  alpX=%7.4f alpY=%7.4f alpZ=%7.4f\n", i, j, k, alpX*180/pi, alpY*180/pi, alpZ*180/pi);
			line1_s = line3_s; kneib1 = kneib3; kneib2 = 0;
			DirectTo_kij100(crk, cri, line2_s, curv2, kneib2);
			if (kneib2 != 2) line2_s = A(DelXY);
		}
		else	ErrorMsg("Cant find SURFACE for 111 plane");
		sCurv = line1_s * line2_s;
	}
	else
	{
		printf("Can't direct around cell of i%d j%d k%d cri%d crj%d crk%d, dir%d\n", 
			i, j, k, cri, crj, crk, SlushDirection[i][j][k]);
		curvS =  0.0;
		sCurv = pow((A(DelX) * A(DelY) * A(DelZ)),0.66666);
	}
	dTGT = (1.1e-8 * 1685 * curvS);
	return;
};//end of function

//______________________________________________________________________________
//      Curvature2d
//    Find curvature using three points of curcle in vector-ortogonal plane
//______________________________________________________________________________
double Curv2d (int ci, int cj, double DX, double DY, double &line_s, double &alp12)
{
	double X0Cros, Y0Cros, diaIncl, curv2d;
	int signC = 1; source = -1;
	line_s = 0;
	alp1 = 0.0; alp2 = 0.0;
	alp1 = (yy1 == 0.0) ? 0.49999 * pi : -atan(x1/yy1);
	alp2 = ( y2 == 0.0) ? 0.49999 * pi : -atan(x2/y2);
	alp12 = fabs(0.5 * (alp1 + alp2));
	X0Cros = (yy1 != y2) ? (x1*y2 - x2*yy1)/(y2 - yy1) : 0.0;   //Crossing of X-axis by line {{x1,yy1},{x2,y2}}
	Y0Cros = ( x1 != x2) ? (x1*y2 - x2*yy1)/(x1 -  x2) : 0.0;   //Crossing of Y-axis by line {{x1,yy1},{x2,y2}}
	diaIncl = atan(DX/DY);
	if (alp1 == alp2)
	{
		source = 0; curv2d = 0.0;
		if (fabs(alp1) >= diaIncl) line_s = fabs(DX/sin(alp1)); else line_s = DY/cos(alp1);
	}
	else if ((((alp1 >= 0) && (alp2 <= 0)) || ((alp1 <= 0) && (alp2 >= 0))) && (((yy1>=0) && (y2>=0)) || ((yy1<=0) && (y2<=0))))
	{
		source = 1;
		curv2d = (cos(alp1) + cos(alp2))/DX;
		line_s = (pi - fabs(alp1) - fabs(alp2))/curv2d;
	}
	else if ((alp12 >= diaIncl) && (((alp1 >= 0) && (alp2 >= 0)) || ((alp1 <= 0) && (alp2 <= 0))))
	{
		source = 2;
		curv2d = fabs(cos(alp1) - cos(alp2))/DX;
		line_s = fabs(fabs(alp1) - fabs(alp2))/curv2d;
	}
	else if ((((alp1 >= 0) && (alp2 <= 0)) || ((alp1 <= 0) && (alp2 >= 0))) && (((x1>=0) && (x2>=0)) || ((x1<=0) && (x2<=0))))
	{
		source = 3;
		curv2d = (fabs(sin(alp1)) + fabs(sin(alp2)))/DY;
		line_s = (fabs(alp1) + fabs(alp2))/curv2d;
	}
	else if ((alp12 <= diaIncl) && (((alp1 >= 0) && (alp2 >= 0)) || ((alp1 <= 0) && (alp2 <= 0))))
	{
		source = 4;
		curv2d = fabs((fabs(sin(alp1)) - fabs(sin(alp2))))/DY;
		line_s = fabs(fabs(alp1) - fabs(alp2))/curv2d;
	}
	else
	{	printf("Cant find curv-type %d %d %d ci=%2d cj=%2d alp1=%6.3f alp2=%6.3f diaIncl=%5.2f\n",
			i, j, k, ci, cj, alp1*180/pi, alp2*180/pi, diaIncl*180/pi);
		printf("x1=%8.5f y1=%8.5f x2=%8.5f y2=%8.5f ______\n", x1*1e7,yy1*1e7,x2*1e7,y2*1e7);
		getche();
	};
	if (((Y0Cros > 0.0) && (cj == 1)) || ((Y0Cros < 0.0) && (cj == -1)) ||
		((X0Cros > 0.0) && (ci == 1)) || ((X0Cros < 0.0) && (ci == -1)))
			signC = -1;
	return (signC * curv2d);
};//endfunction

//______________________________________________________________________________
//      FindBoundary100
//    Find boundary points within LOCAL X-Y-system
 //______________________________________________________________________________
void FindBoundary100(int i, int j, int k, int indi, int indj, int ijk,
                     double dX, double xFrac, int i_first, int i_last, int i_mode,
                     double dY, double yFrac, int j_first, int j_last, int j_mode, int &kneib)
{
	int in, jn;
	if ((j == j_first) || (j == j_last-1)) // Y-axis  j priority over i results in k-curcle priority over j-flow --
    {
		switch(j_mode)
        {
        	case PERIODIC:
            {
  		    	jn = (j == j_first) ? j_last-1 : j_first;
				for (in = i-1; in <= i+1; in++)
					if ((in >= i_first) && (in < i_last))
                    {
           				switch (ijk)
						{
							case 0: {in2 = in; jn2 = jn; kn2 = k;}; break;
							case 1: {in2 = k; jn2 = in; kn2 = jn;}; break;
							case 2: {in2 = jn; jn2 = k; kn2 = in;}; break;
						};//end switch ijk
                        nodeC = IJKToIndex(in2, jn2, kn2);
						if (C(State) == SLUSH)
                        {
							//if (Sim.curTime > 40E-9) printf("B0Y  %d %d %d ", in2,jn2,kn2);
							kneib++;	return;
                        };
					};
			}; //endcaseCURCLE
			break;
			case FLOW:
			{
         		if (((j == j_first) && (yy1/x1 > 0) &&
					(((indi == -1) && (indj == 0)) || ((indi == 0) && (indj ==  1)) || ((indi == -1) && (indj ==  1)))) ||
					((j == j_first) && (yy1/x1 < 0) &&
            		(((indi ==  1) && (indj == 0)) || ((indi == 0) && (indj ==  1)) || ((indi ==  1) && (indj ==  1)))))
						{x2 = 0; y2 = -yFrac * dY; kneib++; //if (Sim.curTime > 40E-9) printf("B0F");
				return;}
				else if (((j == j_last-1) && (yy1/x1 > 0) &&
					(((indi ==  1) && (indj == 0)) || ((indi == 0) && (indj == -1)) || ((indi ==  1) && (indj == -1)))) ||
					((j == j_last-1) && (yy1/x1 < 0) &&
					(((indi == -1) && (indj == 0)) || ((indi == 0) && (indj == -1)) || ((indi == -1) && (indj == -1)))))
               			{x2 = 0; y2 = (1 - yFrac) * dY; kneib++; //if (Sim.curTime > 40E-9) printf("B0F");
				return;}
				else
            		{x2 = -x1; y2 = -yy1; kneib++; //if (Sim.curTime > 40E-9) printf("B0F");
				return;};
			};//endcaseFLOW
			break;
			case ADIABAT:
			{
         		y2 = (j == j_first) ? -dY * yFrac : dY * (1 - yFrac);
				x2 = 0; kneib++; //if (Sim.curTime > 40E-9) printf("B0C");
				return;
			};
			break;
		};//endswitch_b_mode
	};//end if j
	if ((i == i_first) || (i == i_last-1)) // X-axis --------------------------------------------
    {
		switch(i_mode)
        {
        	case PERIODIC:
            {
  	      		in = (i == i_first) ? i_last-1 : i_first;
				for (jn = j-1; jn <= j+1; jn++)
					if ((jn >= j_first) && (jn < j_last))
                    {
           				switch (ijk)
						{
							case 0: {in2 = in; jn2 = jn; kn2 = k;}; break;
							case 1: {in2 = k; jn2 = in; kn2 = jn;}; break;
							case 2: {in2 = jn; jn2 = k; kn2 = in;}; break;
						};//end switch ijk
                        nodeC = IJKToIndex(in2, jn2, kn2);
						if (C(State) == SLUSH)
                        {
							kneib++;
							//if (Sim.curTime > 40E-9) printf("B0X  %d %d %d ", in2,jn2,kn2);
							return;
                        };
					};
			}; //endcaseCURCLE
			break;
			case FLOW:
			{
				if (((i == i_first) && 	(yy1/x1 > 0) &&
            		(((indi ==  1) && (indj == 0)) || ((indi == 0) && (indj == -1)) || ((indi ==  1) && (indj == -1)))) ||
					((i == i_first) && (yy1/x1 < 0) &&
					(((indi ==  1) && (indj == 0)) || ((indi == 0) && (indj ==  1)) || ((indi ==  1) && (indj ==  1)))))
						{y2 = 0; x2 =     -xFrac  * dX; kneib++; //if (Sim.curTime > 40E-9) printf("B0F");
				return;}
				else if (((i == i_last-1) && (yy1/x1 > 0) &&
					(((indi == -1) && (indj == 0)) || ((indi == 0) && (indj ==  1)) || ((indi == -1) && (indj ==  1)))) ||
					((i == i_last-1) && (yy1/x1 < 0) &&
            		(((indi == -1) && (indj == 0)) || ((indi == 0) && (indj == -1)) || ((indi == -1) && (indj == -1)))))
						{y2 = 0; x2 = (1 - xFrac) * dX; kneib++; //if (Sim.curTime > 40E-9) printf("B0F");
				return;}
				else
            			{x2 = -x1; y2 = -yy1; kneib++; //if (Sim.curTime > 40E-9) printf("B0F");
				return;};
			};//endcaseFLOW
			break;
			case ADIABAT:
			{
				x2 = (i == i_first) ? -dX * xFrac : dX * (1 - xFrac);
				y2 = 0; kneib++; //if (Sim.curTime > 40E-9) printf("B0C");
				return;
			};
			break;
		};//endswitch_b_mode
	};//end if i
};//end function
//______________________________________________________________________________
//      FindBoundary110
//    Find boundary points within LOCAL system X-Z, where X is daigonal of base square(i,j) and
//    Z is cube edge k of the base     // ind1 = 1 if daigonal is right and -1 is left one
//______________________________________________________________________________
void FindBoundary110(int i, int j, int k, int indi, int indd, int indk, int ijk,
                     double dX, double xFrac, int i_first, int i_last, int i_mode,
                     						  int j_first, int j_last, int j_mode,
                     double dZ, double KPos, int k_first, int k_last, int k_mode, int &kneib)
{
	int in, jn, kn;
	if ((k == k_first) || (k == k_last-1)) //Z-axis (k) in LOCAL system -----------------------------------
	{
		switch(k_mode)
		{
			case PERIODIC:
			{
				kn = (k == k_first) ? k_last-1 : k_first;
				for (in = i-1, jn = j-indd; in <= i+1, jn <= j+indd; in++, jn+=indd)
					if ((in >= i_first) && (in < i_last) && (jn >= j_first) && (jn < j_last))
                    {
           				switch (ijk)
						{
							case 0: {in2 = in; jn2 = jn; kn2 = kn;}; break;
							case 1: {in2 = kn; jn2 = in; kn2 = jn;}; break;
							case 2: {in2 = jn; jn2 = kn; kn2 = in;}; break;
						};//end switch ijk
                        nodeC = IJKToIndex(in2, jn2, kn2);
						if (C(State) == SLUSH)
                        {
							//if (Sim.curTime > 40E-9) printf("B1X  %d %d %d ", in2,jn2,kn2);
				kneib++;	return;
                        };
					};
			}; //endcaseCURCLE
			break;
			case FLOW:
			{
         		if (((k == k_first) && (yy1/x1 > 0) &&
				(((indi == -1) && (indk == 0)) || ((indi == 0) && (indk ==  1)) || ((indi == -1) && (indk ==  1)))) ||
			    	((k == k_first) && (yy1/x1 < 0) &&
            	(((indi ==  1) && (indk == 0)) || ((indi == 0) && (indk ==  1)) || ((indi ==  1) && (indk ==  1)))))
					{x2 = 0; y2 =     -KPos  * dZ; kneib++; //if (Sim.curTime > 40E-9) printf("B1F");
				return;}
				else if (((k == k_last-1) && (yy1/x1 > 0) &&
				(((indi ==  1) && (indk == 0)) || ((indi == 0) && (indk == -1)) || ((indi ==  1) && (indk == -1)))) ||
					((k == k_last-1 ) && (yy1/x1 < 0) &&
				(((indi == -1) && (indk == 0)) || ((indi == 0) && (indk == -1)) || ((indi == -1) && (indk == -1)))))
					{x2 = 0; y2 = (1 - KPos) * dZ; kneib++; //if (Sim.curTime > 40E-9) printf("B1F");
				return;}
				else
				{
            		x2 = -x1; y2 = -yy1; kneib++; //if (Sim.curTime > 40E-9) printf("B1F");
					return;
				};
			};//endcaseFLOW
			break;
			case ADIABAT:
			{
				y2 = (i == i_first) ? -dZ * KPos : dZ * (1 - KPos);
				x2 = 0; kneib++; //if (Sim.curTime > 40E-9) printf("B1C"); return;
			};
			break;
		};//endswitch_b_mode
	};//endif k
	if ((j == j_first) || (j == j_last-1)) // second base component j of LOCAL X-axis-(diagonal in base) --------------
    {
		switch(j_mode)
		{
        	case PERIODIC:
			{
  	        	in, jn = (j == j_first) ? i-indd, j_last-1 : i+indd, j_first;
				for (kn = k-1; kn <= k+1; kn++)
					if ((kn >= k_first) && (kn < k_last))
                    {
           				switch (ijk)
						{
							case 0: {in2 = in; jn2 = jn; kn2 = kn;}; break;
							case 1: {in2 = kn; jn2 = in; kn2 = jn;}; break;
							case 2: {in2 = jn; jn2 = kn; kn2 = in;}; break;
						};//end switch ijk
                        nodeC = IJKToIndex(in2, jn2, kn2);
						if (C(State) == SLUSH)
                        {
							//if (Sim.curTime > 40E-9) printf("B1Z  %d %d %d ", in2,jn2,kn2);
							kneib++;	return;
                        };
					};
			}; //endcaseCURCLE
			break;
			case FLOW:
			{
				if (((j == j_first) && (yy1/x1 > 0) &&
                (((indi ==  indd) && (indk == 0)) || ((indi == 0) && (indk == -indd)) || ((indi ==  indd) && (indk == -indd)))) ||
                   ((j == j_first) && (yy1/x1 < 0) &&
                (((indi == -indd) && (indk == 0)) || ((indi == 0) && (indk == -indd)) || ((indi == -indd) && (indk == -indd)))))
						{y2 = 0; x2 = (indd == 1) ?   -xFrac  * dX : (1-xFrac) * dX; kneib++; //if (Sim.curTime > 40E-9) printf("B1F");
				return;}
				else if (((j == j_last-1) && (yy1/x1 > 0) &&
				(((indi == -indd) && (indk == 0)) || ((indi == 0) && (indk ==  indd)) || ((indi == -indd) && (indk ==  indd)))) ||
					((j == j_last-1) && (yy1/x1 < 0) &&
            	(((indi ==  indd) && (indk == 0)) || ((indi == 0) && (indk ==  indd)) || ((indi ==  indd) && (indk ==  indd)))))
						{y2 = 0; x2 = (indd == 1) ? (1-xFrac) * dX :     xFrac * dX; kneib++; //if (Sim.curTime > 40E-9) printf("B1F");
				return;}
				else
				{
            		x2 = -x1; y2 = -yy1; kneib++; //if (Sim.curTime > 40E-9) printf("B1F");
					return;
				};
			};//endcaseFLOW
			break;
			case ADIABAT:
			{
				x2 = (i == i_first) ? -dX * xFrac : dX * (1 - xFrac);
				y2 = 0; kneib++; //if (Sim.curTime > 40E-9) printf("B1C"); return;
			};
			break;
		};//endswitch_b_mode
	};//end if j
	if ((i == i_first) || (i == i_last-1)) // first base component i of LOCAL X-axis-(diagonal in base) ---------------
    {
    	switch(i_mode)
		{
        	case PERIODIC:
			{
        		jn, in = (i == i_first) ? j-indd, i_last-1 : j+indd, i_first;
				for (kn = k-1; k <= k+1; kn++)
					if ((kn >= k_first) && (kn < k_last) && (jn >= j_first) && (jn < j_last))
                    {
           				switch (ijk)
						{
							case 0: {in2 = in; jn2 = jn; kn2 = kn;}; break;
							case 1: {in2 = kn; jn2 = in; kn2 = jn;}; break;
							case 2: {in2 = jn; jn2 = kn; kn2 = in;}; break;
						};//end switch ijk
                        nodeC = IJKToIndex(in2, jn2, kn2);
						if (C(State) == SLUSH)
                        {
							//if (Sim.curTime > 40E-9) printf("B1Y  %d %d %d ", in2,jn2,kn2);
							kneib++;	return;
                        };
					};
			}; //endcaseCURCLE
            break;
			case FLOW:
			{
         		if (((i == i_first) && (yy1/x1 > 0) &&
            	(((indi ==  1) && (indk == 0)) || ((indi == 0) && (indk == -1)) || ((indi ==  1) && (indk == -1)))) ||
					((i == i_first) && (yy1/x1 < 0) &&
                (((indi ==  1) && (indk == 0)) || ((indi == 0) && (indk ==  1)) || ((indi ==  1) && (indk ==  1)))))
                		{y2 = 0; x2 = -xFrac * dX; kneib++;// if (Sim.curTime > 40E-9) printf("B1F"); 
				return;}
				else if (((i == i_last-1) && (yy1/x1 > 0) &&
			    (((indi == -1) && (indk == 0)) || ((indi == 0) && (indk ==  1)) || ((indi == -1) && (indk ==  1)))) ||
					((i == i_last-1) && (yy1/x1 < 0) &&
            	(((indi == -1) && (indk == 0)) || ((indi == 0) && (indk == -1)) || ((indi == -1) && (indk == -1)))))
               		{y2 = 0; x2 = (1-xFrac) * dX; kneib++; //if (Sim.curTime > 40E-9) printf("B1F");
				return;}
				else
				{
            		x2 = -x1; y2 = -yy1; kneib++; //if (Sim.curTime > 40E-9) printf("B1F"); return;
				};
			};//endcaseFLOW
			break;
			case ADIABAT:
			{
				x2 = (i == i_first) ? -dX * xFrac : dX * (1 - xFrac);
				y2 = 0; kneib++; //if (Sim.curTime > 40E-9) printf("B1C"); return;
			};
			break;
		  };//endswitch_b_mode
	};//endif i
};//end function
//______________________________________________________________________________
//      FindNeighbor100              Find neighbor slush cells
//______________________________________________________________________________
void FindNeighbor100(int il, int jl, int kl, int ijk, int i_first, int i_last,
                                          int j_first, int j_last, int &kneib)
{
	int kneibSw = 0;
	int in, jn, nodeD;kneib = 0; in2 = -1;
	int StateTest;
	int ind, jnd, knd;
	for (jn = jl-1; jn <= jl+1; jn +=2)
	{
		StateTest = SLUSH;
		for (in = il-1; in <= il+1; in++)
		{
			if ((in >= i_first) && (in < i_last) && (jn >= j_first) && (jn < j_last))
			{
				switch (ijk)
				{
					case 0: {ind = in; jnd = jn; knd = kl;}; break;
					case 1: {ind = kl; jnd = in; knd = jn;}; break;
					case 2: {ind = jn; jnd = kl; knd = in;}; break;
				};//end switch ijk
				nodeD = IJKToIndex(ind, jnd, knd);
				if (D(State) == SLUSH)
				{
					uncode = 0; ChekNeibNumber(kneibSw, kneib, nodeD, ind, jnd, knd, uncode);
				}
				else
				{
					if ((in > il-1) && (in > i_first) && (StateTest != D(State)) && 
						((D(State) == LIQUID) || (D(State) == SOLID)) && ((StateTest == LIQUID) || (StateTest == SOLID)))
					{
						uncode = 1; ChekNeibNumber(kneibSw, kneib, nodeD,ind, jnd, knd, uncode); //if (Sim.curTime > 40E-9) printf("L/S");
					}
				};
				StateTest = D(State);
			};
		};
	};
	for (in = il-1; in <= il+1; in +=2)
	{
		StateTest = SLUSH;
		for (jn = jl-1; jn <= jl+1; jn++)
		{
			if ((in >= i_first) && (in < i_last) && (jn >= j_first) && (jn < j_last))
			{
			switch (ijk)
				{
					case 0: {ind = in; jnd = jn; knd = kl;}; break;
					case 1: {ind = kl; jnd = in; knd = jn;}; break;
					case 2: {ind = jn; jnd = kl; knd = in;}; break;
				};//end switch ijk
				nodeD = IJKToIndex(ind, jnd, knd);
				if (((jn != jl-1) && (jn != jl+1)) && (D(State) == SLUSH))
				{
					uncode = 0; ChekNeibNumber(kneibSw, kneib, nodeD, ind, jnd, knd, uncode);
				}
				else
				{
					if ((jn > jl-1) && (jn > j_first) && (StateTest != D(State)) &&
					   ((D(State) == LIQUID) || (D(State) == SOLID)) && ((StateTest == LIQUID) || (StateTest == SOLID)))
					{
						uncode = 2; ChekNeibNumber(kneibSw, kneib, nodeD, ind, jnd, knd, uncode); //if (Sim.curTime > 40E-9) printf("L/S");
					};
				};
				StateTest = D(State);
			};
		};
	};
	if (kneibSw == 3)
	{			
		kneib = Find2NeighborsFrom3();
		nodeB = IJKToIndex(in1, jn1, kn1);
		nodeC = IJKToIndex(in2, jn2, kn2);
		if (kneib == 2)
			InfoMsg("100 2 neighbor from 3 ones have     found around i%d j%d k%d  %d %d %d \n", i, j, k, cri, crj, crk);
		else
			InfoMsg("100 2 neighbor from 3 ones have not found around i%d j%d k%d  %d %d %d\n", i, j, k, cri, crj, crk);
	}
	if (kneibSw > 3)
	{
		InfoMsg("100 %d neighbors in wrong wall i%d j%d k%d  %d %d %d   %d %d %d  %d %d %d  %d %d %d  %d %d %d\n",kneibSw, i, j, k, cri, crj, crk, 
				crystalI[SlushDirection[in1][jn1][kn1]], crystalJ[SlushDirection[in1][jn1][kn1]], crystalK[SlushDirection[in1][jn1][kn1]],
				crystalI[SlushDirection[in2][jn2][kn2]], crystalJ[SlushDirection[in2][jn2][kn2]], crystalK[SlushDirection[in2][jn2][kn2]],
				crystalI[SlushDirection[in3][jn3][kn3]], crystalJ[SlushDirection[in3][jn3][kn3]], crystalK[SlushDirection[in3][jn3][kn3]],
				crystalI[SlushDirection[in4][jn4][kn4]], crystalJ[SlushDirection[in4][jn4][kn4]], crystalK[SlushDirection[in4][jn4][kn4]]); //getche();
		if ((SlushDirection[in1][jn1][kn1] == SlushDirection[in2][jn2][kn2]) &&
			(SlushDirection[in3][jn3][kn3] == SlushDirection[in2][jn2][kn2]))
		{
//			  if (Sim.curTime > 40E-9) 
			{printf("Repair wall 100:i%d j%d k%d  %d %d %d    new:%d %d %d\n",i, j, k, cri, crj, crk, 
				crystalI[SlushDirection[in1][jn1][kn1]], crystalJ[SlushDirection[in1][jn1][kn1]], crystalK[SlushDirection[in1][jn1][kn1]]); //getche();
			};
			A(SlushDirection) = SlushDirection[in2][jn2][kn2];
			PositionInterface (A(SlushDirection), A(FractionSolid), A(IPos), A(JPos), A(KPos));
		};	
	};
}; //endfunction

//______________________________________________________________________________
//      FindNeighbor110              Find naighbor slush cells
//______________________________________________________________________________
void FindNeighbor110 (int il, int jl, int kl, int indd, int ijk, int i_first, int i_last,
                                                                 int j_first, int j_last,
                                                                 int k_first, int k_last, int &kneib)
{
	int kneibSw = 0;
	int in, jn, kn, nodeD;
	kneib = 0; in2 = -1;
	int StateTest;
	int ind, jnd, knd;
	for (kn = kl-1; kn <= kl+1; kn +=2)
	{
		StateTest = SLUSH;
		for (in = il-1, jn = jl-indd; in <= il+1; in++, jn+=indd)
		{
			if ((in >= i_first) && (in < i_last) && (jn >= j_first) && (jn < j_last) && (kn >= k_first) && (kn < k_last))
			{
				switch (ijk)
				{
					case 0: {ind = in; jnd = jn; knd = kn;}; break;
					case 1: {ind = kn; jnd = in; knd = jn;}; break;
					case 2: {ind = jn; jnd = kn; knd = in;}; break;
				};//end switch ijk
				nodeD = IJKToIndex(ind, jnd, knd);
				if (D(State) == SLUSH)
				{
					uncode = 0; ChekNeibNumber(kneibSw, kneib, nodeD, ind, jnd, knd, uncode);
				}
				else
				{
					if ((in > il-1) && (in > i_first) && (jn > j_first) && (jn < j_last-1) && (StateTest != D(State)) &&
						((D(State) == LIQUID) || (D(State) == SOLID)) && ((StateTest == LIQUID) || (StateTest == SOLID)))
					{
						uncode = 1; ChekNeibNumber(kneibSw, kneib, nodeD, ind, jnd, knd, uncode); 
					}
				};
				StateTest = D(State);
			};
		};
	};
	for (in = il-1, jn = jl-indd; in <= il+1; in +=2, jn+=2*indd)
	{
		StateTest = SLUSH;
		for (kn = kl-1; kn <= kl+1; kn++)
		{
			if ((in >= i_first) && (in < i_last) && (jn >= j_first) && (jn < j_last) && (kn >= k_first) && (kn < k_last))
			{
			switch (ijk)
				{
					case 0: {ind = in; jnd = jn; knd = kn;}; break;
					case 1: {ind = kn; jnd = in; knd = jn;}; break;
					case 2: {ind = jn; jnd = kn; knd = in;}; break;
				};//end switch ijk
				nodeD = IJKToIndex(ind, jnd, knd);
				if (((kn != kl-1) && (kn != kl+1)) && (D(State) == SLUSH))
				{
					uncode = 0; ChekNeibNumber(kneibSw, kneib, nodeD, ind, jnd, knd, uncode);
				}
				else
				{
					if ((kn > kl-1) && (kn > k_first) && (StateTest != D(State)) &&      //???????????? @
						((D(State) == LIQUID) || (D(State) == SOLID)) && ((StateTest == LIQUID) || (StateTest == SOLID)))
					{
						uncode = 2; ChekNeibNumber(kneibSw, kneib, nodeD, ind, jnd, knd, uncode); //if (Sim.curTime > 40E-9) printf("L/S");
					};
				};
				StateTest = D(State);
			};
		};//end for k
	};//end for i,j
	if (kneibSw == 3)
	{			
		kneib = Find2NeighborsFrom3();
		nodeB = IJKToIndex(in1, jn1, kn1);
		nodeC = IJKToIndex(in2, jn2, kn2);
		if (kneib == 2)
			InfoMsg("110 2 neighbor from 3 ones have     found around i%d j%d k%d  %d %d %d \n", i, j, k, cri, crj, crk);
		else
		InfoMsg("110 2 neighbor from 3 ones have not found around i%d j%d k%d  %d %d %d \n", i, j, k, cri, crj, crk);
	}
	if (kneibSw > 3)
	{
		InfoMsg("110 %d neighbors in wrong wall i%d j%d k%d  %d %d %d  1:%d %d %d 2:%d %d %d 3:%d %d %d\n",kneibSw, i, j, k, cri, crj, crk, 
				crystalI[SlushDirection[in1][jn1][kn1]], crystalJ[SlushDirection[in1][jn1][kn1]], crystalK[SlushDirection[in1][jn1][kn1]],
				crystalI[SlushDirection[in2][jn2][kn2]], crystalJ[SlushDirection[in2][jn2][kn2]], crystalK[SlushDirection[in2][jn2][kn2]],
				crystalI[SlushDirection[in3][jn3][kn3]], crystalJ[SlushDirection[in3][jn3][kn3]], crystalK[SlushDirection[in3][jn3][kn3]],
				crystalI[SlushDirection[in4][jn4][kn4]], crystalJ[SlushDirection[in4][jn4][kn4]], crystalK[SlushDirection[in4][jn4][kn4]]); //getche();
		if ((SlushDirection[in1][jn1][kn1] == SlushDirection[in2][jn2][kn2]) &&
			(SlushDirection[in3][jn3][kn3] == SlushDirection[in2][jn2][kn2]))
		{
			{
				printf("Repair wall 110:i%d j%d k%d  %d %d %d    new:%d %d %d\n",i, j, k, cri, crj, crk, 
				crystalI[SlushDirection[in1][jn1][kn1]], crystalJ[SlushDirection[in1][jn1][kn1]], crystalK[SlushDirection[in1][jn1][kn1]]); //getche();
			};
	
			A(SlushDirection) = SlushDirection[in2][jn2][kn2];
			PositionInterface (A(SlushDirection), A(FractionSolid), A(IPos), A(JPos), A(KPos));
		};	
	};
}; //endfunction

//______________________________________________________________________________
//      ChekNeibNumber           Hold neighbor number and set "B" and "C" ones
//______________________________________________________________________________
void ChekNeibNumber(int &kneibSw, int &kneib, int nodeD, int in, int jn, int kn, int uncode)
{
	switch (++kneibSw)
	{
		case 1:
			{in1 = in; jn1 = jn; kn1 = kn; nodeB = nodeD; kneib = kneibSw; uncode1 = uncode;};
   		break;
		case 2:
			{in2 = in; jn2 = jn; kn2 = kn; nodeC = nodeD; kneib = kneibSw; uncode2 = uncode;};
 		break;
		case 3:
			{in3 = in; jn3 = jn; kn3 = kn;};
		break;
		default:
			{in4 = in; jn4 = jn; kn4 = kn; kneib = kneibSw;}; 
		break;
	}; //endswitch
};//end of function

//______________________________________________________________________________
//      Find2NeighborsFrom3              Find 2 right naighbors from 3 ones
//______________________________________________________________________________
int Find2NeighborsFrom3()
{
	if ((i == in1) || (j == jn1))
	{
		if ((in1 == in2) || (jn1 == jn2))
		{
			in2 = in3; jn2 = jn3; kn2 = kn3;
		};
		return(2);
	}
	else if ((i == in2) || (j == jn2))
	{
		if ((in2 == in1) || (jn2 == jn1))
		{
			in1 = in3; jn1 = jn3; kn1 = kn3;
		};
		return(2);
	}
	else if ((i == in3) || (j == jn3))
	{
		if ((in3 == in2) || (jn3 == jn2))
		{
			in2 = in3; jn2 = jn3; kn2 = kn3;
		}
		else
		{
			in1 = in3; jn1 = jn3; kn1 = kn3;
		};
		return(2);
	}
	else
		return(3);
}; //end of function

//______________________________________________________________________________
//    DirectTo_ijk001      Direct finding of two neighbors within
//                         local ijk (001) plain with respect to initial ijk_basis
//______________________________________________________________________________
void DirectTo_ijk100(int indi, int indj, double &line_s, double &curv, int &kneib)
{
	double iBFr, jBFr, iCFr, jCFr;
	FindNeighbor100(i, j, k, 0, I_FIRST, I_LAST, J_FIRST, j_LAST, kneib);
	if (kneib != 0)
	{
        switch (uncode1)
		{
			case 0: {iBFr  = B(IPos); jBFr = B(JPos);}; break;
			case 1: {iBFr  =     0.0; jBFr =     0.5;}; break;
			case 2: {iBFr  =     0.5; jBFr =     0.0;}; break;
		};
		x1 =  B(DelX) * iBFr - A(DelX) * A(IPos) + Geometry.xMap[in1] - Geometry.xMap[i];
		yy1 = B(DelY) * jBFr - A(DelY) * A(JPos) + Geometry.yMap[jn1] - Geometry.yMap[j];
	};
	if (kneib == 2)
	{
        switch (uncode2)
		{
			case 0: {iCFr  = C(IPos); jCFr = C(JPos);}; break;
			case 1: {iCFr  =     0.0; jCFr =     0.5;}; break;
			case 2: {iCFr  =     0.5; jCFr =     0.0;}; break;
		};
		x2 = C(DelX) * iCFr - A(DelX) * A(IPos) + Geometry.xMap[in2] - Geometry.xMap[i];
		y2 = C(DelY) * jCFr - A(DelY) * A(JPos) + Geometry.yMap[jn2] - Geometry.yMap[j];
	};
	if (kneib == 1)
	{
	  	FindBoundary100(i, j, k, indi, indj, 0, A(DelX), A(IPos), I_FIRST, I_LAST, I_EDGE,
                                                A(DelY), A(JPos), J_FIRST, j_LAST, J_EDGE, kneib);
	  	if ((kneib == 2) && (in2 != -1))
     	{
			if ((in2 == I_FIRST) || (in2 == I_LAST -1))
			{
				x2 = (in2 == I_FIRST) ? A(DelX) : -A(DelX);
				x2 = C(DelX) * C(IPos) - A(DelX) * A(IPos) + x2;
				y2 = C(DelY) * C(JPos) - A(DelY) * A(JPos) + Geometry.yMap[jn2] - Geometry.yMap[j];
			}
			else if ((jn2 == J_FIRST) ||(jn2 == j_LAST -1))
			{
				y2 = (jn2 == J_FIRST) ? A(DelY) : -A(DelY);
				x2 = C(DelX) * C(IPos) - A(DelX) * A(IPos) + Geometry.xMap[in2] - Geometry.xMap[i];
     			y2 = C(DelY) * C(JPos) - A(DelY) * A(JPos) + y2;
	   		};
		};
	};
	if (kneib != 2)  curv = 0;
	else  curv = Curv2d(indi, indj, A(DelX), A(DelY), line_s, alpZ);
};//end of function
//______________________________________________________________________________
//    DirectTo_jki100      Direct finding of two neighbors within
//                         local jki (100) plain with respect to initial ijk_basis
//______________________________________________________________________________
void DirectTo_jki100(int indi, int indj, double &line_s, double &curv, int &kneib)
{
	double jBFr, kBFr, jCFr, kCFr;
	FindNeighbor100(j, k, i, 1, J_FIRST, j_LAST, K_FIRST, K_LAST, kneib);
	if (kneib != 0)
	{
        switch (uncode1)
		{
			case 0: {jBFr  = B(JPos); kBFr = B(KPos);}; break;
			case 1: {jBFr  =     0.0; kBFr =     0.5;}; break;
			case 2: {jBFr  =     0.5; kBFr =     0.0;}; break;
		};
		x1 =  B(DelY) * jBFr - A(DelY) * A(JPos) + Geometry.yMap[jn1] - Geometry.yMap[j];
		yy1 = B(DelZ) * kBFr - A(DelZ) * A(KPos) + Geometry.zMap[kn1] - Geometry.zMap[k];
	};
	if (kneib == 2)
	{
        switch (uncode2)
		{
			case 0: {jCFr  = C(JPos); kCFr = C(KPos);}; break;
			case 1: {jCFr  =     0.0; kCFr =     0.5;}; break;
			case 2: {jCFr  =     0.5; kCFr =     0.0;}; break;
		};
		x2 = C(DelY) * jCFr - A(DelY) * A(JPos) + Geometry.yMap[jn2] - Geometry.yMap[j];
		y2 = C(DelZ) * kCFr - A(DelZ) * A(KPos) + Geometry.zMap[kn2] - Geometry.zMap[k];
	};
	if (kneib == 1)
	{
	  	FindBoundary100(j, k, i, indi, indj, 1, A(DelY), A(JPos), J_FIRST, j_LAST, J_EDGE,
                                                A(DelZ), A(KPos), K_FIRST, K_LAST, K_EDGE, kneib);
	  	if ((kneib == 2) && (in2 != -1))
     	{
			if ((kn2 == K_FIRST) ||(kn2 == K_LAST-1))
			{
				y2 = (kn2 == K_FIRST) ? A(DelZ) : -A(DelZ);
				x2 = C(DelY) * C(JPos) - A(DelY) * A(JPos) + Geometry.yMap[jn2] - Geometry.yMap[j];
				y2 = C(DelZ) * C(KPos) - A(DelZ) * A(KPos) + y2;
			}
			else if ((jn2 == J_FIRST) || (jn2 == j_LAST-1))
			{
		   		x2 = (jn2 == J_FIRST) ? A(DelY) : -A(DelY);
				x2 = C(DelY) * C(JPos) - A(DelY) * A(JPos) + x2;
				y2 = C(DelZ) * C(KPos) - A(DelZ) * A(KPos) + Geometry.zMap[kn2] - Geometry.zMap[k];
			};
		};
	};
	if (kneib != 2)  curv = 0;
	else  curv = Curv2d(indi, indj, A(DelY), A(DelY), line_s, alpY);
};//end of function
//______________________________________________________________________________
//    DirectTo_kij010      Direct finding of two neighbors
//                         into kij (010) plain with respect to initial ijk_basis
//______________________________________________________________________________
void DirectTo_kij100(int indi, int indj, double &line_s, double &curv, int &kneib)
{
	double kBFr, iBFr, kCFr, iCFr;
	FindNeighbor100(k, i, j, 2, K_FIRST, K_LAST, I_FIRST, I_LAST, kneib);
	if (kneib != 0)
	{
        switch (uncode1)
		{
			case 0: {kBFr  = B(KPos); iBFr = B(IPos);}; break;
			case 1: {kBFr  =     0.0; iBFr =     0.5;}; break;
			case 2: {kBFr  =     0.5; iBFr =     0.0;}; break;
		};
		x1 =  B(DelZ) * kBFr - A(DelZ) * A(KPos) + Geometry.zMap[kn1] - Geometry.zMap[k];
		yy1 = B(DelX) * iBFr - A(DelX) * A(IPos) + Geometry.xMap[in1] - Geometry.xMap[i];
	};
	if (kneib == 2)
	{
        switch (uncode2)
		{
			case 0: {kCFr  = C(KPos); iCFr = C(IPos);}; break;
			case 1: {kCFr  =     0.0; iCFr =     0.5;}; break;
			case 2: {kCFr  =     0.5; iCFr =     0.0;}; break;
		};
	  	x2 = C(DelZ) * kCFr - A(DelZ) * A(KPos) + Geometry.zMap[kn2] - Geometry.zMap[k];
		y2 = C(DelX) * iCFr - A(DelX) * A(IPos) + Geometry.xMap[in2] - Geometry.xMap[i];
	};
	if (kneib == 1)
	{
	  	FindBoundary100(k, i, j, indi, indj, 2, A(DelZ), A(KPos), K_FIRST, K_LAST, K_EDGE,
                                                A(DelX), A(IPos), I_FIRST, I_LAST, I_EDGE, kneib);
	  	if ((kneib == 2) && (in2 != -1))
     	{
			if ((kn2 == K_FIRST) || (kn2 == K_LAST-1))
			{
				x2 = (kn2 == K_FIRST) ? A(DelZ) : -A(DelZ);
				x2 = C(DelZ) * C(KPos) - A(DelZ) * A(KPos) + x2;
   				y2 = C(DelX) * C(IPos) - A(DelX) * A(IPos) + Geometry.xMap[in2] - Geometry.xMap[i];
			}
			else if ((in2 == I_FIRST) ||(in2 == I_LAST-1))
			{
				y2 = (in2 == I_FIRST) ? A(DelX) : -A(DelX);
				x2 = C(DelZ) * C(KPos) - A(DelZ) * A(KPos) + Geometry.zMap[kn2] - Geometry.zMap[k];
     			y2 = C(DelX) * C(IPos) - A(DelX) * A(IPos) + y2;
			};//end if
		};//end if
	};//endif
	if (kneib != 2)  curv = 0;
	else  curv = Curv2d(indi, indj, A(DelZ), A(DelX), line_s, alpX);
};//end of function
//______________________________________________________________________________
//    DirectTo_ijk110      Direct finding of two neighbors within
//                         local kij (101) plain with respect to initial ijk_basis
//______________________________________________________________________________
void DirectTo_ijk110(int indi, int indj, int indk, double &line_s, double &curv, int &kneib)
{
	indd = ((indi == indj) ? 1 : -1);
	double iBFr, jBFr, kBFr, iCFr, jCFr, kCFr;
	FindNeighbor110 ( i, j, k, indd, 0, I_FIRST, I_LAST, J_FIRST, j_LAST, K_FIRST, K_LAST, kneib);
	if (kneib != 0)
    {
		switch (uncode1)
		{
			case 0: {iBFr = B(IPos); jBFr = B(JPos); kBFr = B(KPos);}; break;
			case 1: {iBFr =     0.0; jBFr = (indd > 0) ? 0.0 : 1.0; kBFr =     0.5;}; break;
			case 2: {iBFr =     0.5; jBFr =     0.5; kBFr =     0.0;}; break;
		};
		x1 = B(DelXY) * iBFr - A(DelXY) * A(IPos);
		if (in1 > i) x1 = x1 + A(DelXY);  
		else if (in1 < i) x1 = x1 - B(DelXY);  
		yy1 = B(DelZ) * kBFr - A(DelZ) * A(KPos) +  Geometry.zMap[kn1] - Geometry.zMap[k];

//		printf("Ax=%7.5f Ay=%7.5f Ak=%7.5f  %2d %2d %2d  Bx=%6.5f By=%6.5f Bk=%6.5f B(St)=%1d UC=%d\n", A(IPos), A(JPos), A(KPos),  in1, jn1, kn1, iBFr, jBFr, kBFr, B(State), uncode1);	
	};
	if (kneib == 2)
	{
		switch (uncode2)
		{
			case 0: {iCFr = C(IPos); jCFr = C(JPos); kCFr = C(KPos);}; break;
			case 1: {iCFr =     0.0; jCFr = (indd > 0) ? 0.0 : 1.0; kCFr =     0.5;}; break;
			case 2: {iCFr =     0.5; jCFr =     0.5; kCFr =     0.0;}; break;
		};
		x2 = C(DelXY) * iCFr - A(DelXY) * A(IPos);
		if (in2 > i) x2 = x2 + A(DelXY);  
		else if (in2 < i) x2 = x2 - C(DelXY);  
		y2 = C(DelZ) * kCFr - A(DelZ) * A(KPos) + Geometry.zMap[kn2] - Geometry.zMap[k];
//		printf("Ax=%7.5f Ay=%7.5f Ak=%7.5f  %2d %2d %2d  Bx=%6.5f By=%6.5f Bk=%6.5f B(St)=%d UC=%d \n", A(IPos), A(JPos), A(KPos),  in2, jn2, kn2, iCFr, jCFr, kCFr, C(State), uncode2);	
	};
	if (kneib == 1)
	{
	  	FindBoundary110(i, j, k, indi, indd, indk, 0, A(DelXY),
                                                      dia(A(IPos), ((indd > 0) ? A(JPos) : (1 - A(JPos)))),
                                            I_FIRST, I_LAST, I_EDGE,
                                            J_FIRST, j_LAST, J_EDGE,
                         A(DelZ), A(KPos), K_FIRST, K_LAST, K_EDGE, kneib);
	  	if ((kneib == 2) && (in2 != -1))
     	{
			if ((kn2 == K_FIRST) || (kn2 == K_LAST-1))
			{
           		y2 = (kn2 == K_FIRST) ? A(DelZ) : - A(DelZ);
				x2 = C(DelXY) * dia(C(IPos),((indd > 0) ? C(JPos) : (1 - C(JPos)))) -
   					 A(DelXY) * dia(A(IPos),((indd > 0) ? C(JPos) : (1 - C(JPos)))) +
	   				 dia((Geometry.xMap[in2] - Geometry.xMap[i]),(Geometry.yMap[jn2] - Geometry.yMap[j])) * indd;
				y2 = C(DelZ) * C(KPos) - A(DelZ) * A(KPos) +    y2;
			}
			else if ((in2 == I_FIRST) ||(in2 == I_LAST-1))
			{
				x2 = (in2 == I_FIRST) ? A(DelX) : -A(DelX);
				x2 = C(DelXY) * dia(C(IPos),((indd > 0) ? C(JPos) : (1 - C(JPos)))) -
   				     A(DelXY) * dia(A(IPos),((indd > 0) ? C(JPos) : (1 - C(JPos)))) +
	   	             dia(x2,                                     (Geometry.yMap[jn2] - Geometry.yMap[j])) * indd;
 				y2 = C(DelZ) * C(KPos) - A(DelZ) * A(KPos) +    Geometry.zMap[kn2] - Geometry.zMap[k];
			}//end if
			else if ((jn2 == J_FIRST) ||(jn2 == j_LAST-1))
			{
           		x2 = (jn2 == J_FIRST) ? A(DelY) : -A(DelY);
				x2 = C(DelXY) * dia(C(IPos),((indd > 0) ? C(JPos) : (1 - C(JPos)))) -
   					 A(DelXY) * dia(A(IPos),((indd > 0) ? C(JPos) : (1 - C(JPos)))) +
	   				 dia((Geometry.xMap[in2] - Geometry.xMap[i]), x2) * indd;
				y2 = C(DelZ) * C(KPos) - A(DelZ) * A(KPos) +    Geometry.zMap[kn2] - Geometry.zMap[k];
	   		};//end if
		};//end if
	};//endif
	if (kneib != 2)
		{curv = 0; line_s = A(DelZ); alpZ = pi/4;}
	else 	curv = Curv2d(indi, indk, A(DelXY), A(DelZ), line_s, alpZ);
};//end of function
//______________________________________________________________________________
//    DirectTo_jki110      Direct finding of two neighbors within
//                         local jki (011) plain with respect to initial ijk_basis
//______________________________________________________________________________
void DirectTo_jki110(int indi, int indj, int indk, double &line_s, double &curv, int &kneib)
{
	indd = ((indi == indj) ? 1 : -1);
	double iBFr, jBFr, kBFr, iCFr, jCFr, kCFr;
	FindNeighbor110 ( j, k, i, indd, 1, J_FIRST, j_LAST, K_FIRST, K_LAST, I_FIRST, I_LAST, kneib);
	if (kneib != 0)
	{
		switch (uncode1)
		{
			case 0: {jBFr = B(JPos); kBFr = B(KPos); iBFr = B(IPos);}; break;
			case 1: {jBFr =     0.0; kBFr = (indd > 0) ? 0.0 : 1.0; iBFr =     0.5;}; break;
			case 2: {jBFr =     0.5; kBFr =     0.5; iBFr =     0.0;}; break;
		};
		x1 = B(DelYZ) * jBFr - A(DelYZ) * A(JPos);
		if (jn1 > j) x1 = x1 + A(DelYZ);  
		else if (jn1 < j) x1 = x1 - B(DelYZ);  
		yy1 = B(DelX) * iBFr - A(DelX) * A(IPos) + Geometry.xMap[in1] - Geometry.xMap[i];
//		printf("Ax=%7.5f Ay=%7.5f Ak=%7.5f  %2d %2d %2d  Bx=%4.3f By=%4.3f Bk=%4.3f B(St)=%d UC=%d\n", A(JPos), A(KPos), A(IPos),  in1, jn1, kn1, jBFr, kBFr, iBFr, B(State), uncode1);	
	};
	if (kneib == 2)
	{
		switch (uncode2)
		{
			case 0: {jCFr = C(JPos); kCFr = C(KPos); iCFr = C(IPos);}; break;
			case 1: {jCFr =     0.0; kCFr = (indd > 0) ? 0.0 : 1.0; iCFr =     0.5;}; break;
			case 2: {jCFr =     0.5; kCFr =     0.5; iCFr =     0.0;}; break;
		};
		x2 = C(DelYZ) * jCFr - A(DelYZ) * A(JPos);
		if (jn2 > j) x2 = x2 + A(DelYZ);  
		else if (jn2 < j) x2 = x2 - C(DelYZ);  
		y2 = C(DelX) * iCFr - A(DelX) * A(IPos) + Geometry.xMap[in2] - Geometry.xMap[i];
//		printf("Ax=%7.5f Ay=%7.5f Ak=%7.5f  %2d %2d %2d  Cx=%4.3f Cy=%4.3f Ck=%4.3f C(St)=%d UC=%d\n", A(JPos), A(KPos), A(IPos),  in2, jn2, kn2, jCFr, kCFr, iCFr, C(State), uncode2);	
	};
	if (kneib == 1)
	{
		FindBoundary110( j, k, i, indi, indd, indk, 1, A(DelYZ),
                                                     dia(A(JPos), ((indd > 0) ? A(KPos) : (1 - A(KPos)))),
                                            J_FIRST, j_LAST, J_EDGE,
                                            K_FIRST, K_LAST, K_EDGE,
                         A(DelX), A(IPos), I_FIRST, I_LAST, I_EDGE, kneib);
	  	if ((kneib == 2) && (in2 != -1))
     	{
			if ((kn2 == K_FIRST) || (kn2 == K_LAST-1))
			{
           		x2 = (kn2 == K_FIRST) ? A(DelZ) : - A(DelZ);
				x2 = C(DelYZ) * dia(C(JPos),((indd > 0) ? C(KPos) : (1 - C(KPos)))) -
   					 A(DelYZ) * dia(A(JPos),((indd > 0) ? A(KPos) : (1 - A(KPos)))) +
					 dia((Geometry.yMap[jn2] - Geometry.yMap[j]),x2) * indd;
 				y2 = C(DelX) * C(IPos) - A(DelX) * A(IPos) +    Geometry.xMap[in2] - Geometry.xMap[i]; 
			}
			else if ((in2 == I_FIRST) ||(in2 == I_LAST-1))
			{
				y2 = (in2 == I_FIRST) ? A(DelY) : -A(DelY);
				x2 = C(DelYZ) * dia(C(JPos),((indd > 0) ? C(KPos) : (1 - C(KPos)))) -
   					 A(DelYZ) * dia(A(JPos),((indd > 0) ? A(KPos) : (1 - A(KPos)))) +
	   				 dia((Geometry.yMap[jn2] - Geometry.yMap[j]),(Geometry.zMap[kn2] - Geometry.zMap[k])) * indd;
 				y2 = C(DelX) * C(IPos) - A(DelX) * A(IPos) +   y2;
			}
			else if ((jn2 == J_FIRST) ||(jn2 == j_LAST-1))
			{
				x2 = (jn2 == J_FIRST) ? A(DelY) : -A(DelY);
				x2 = C(DelYZ) * dia(C(JPos),((indd > 0) ? C(KPos) : (1 - C(KPos)))) -
   					 A(DelYZ) * dia(A(JPos),((indd > 0) ? A(KPos) : (1 - A(KPos)))) +
	   				 dia(x2,                                     (Geometry.zMap[kn2] - Geometry.zMap[k])) * indd;
 				y2 = C(DelX) * C(IPos) - A(DelX) * A(IPos) +    Geometry.xMap[in2] - Geometry.xMap[i];
			};//end if
		};//end if
	};//endif
	if (kneib != 2)
		{curv = 0; line_s = A(DelX); alpX = pi/4;}
	else 	curv = Curv2d(indi, indk, A(DelYZ), A(DelX), line_s, alpX);
};//end of function
//______________________________________________________________________________
//    DirectTo_kij110      Direct finding of two neighbors within
//                         local kij (101) plain with respect to initial ijk_basis
//______________________________________________________________________________
void DirectTo_kij110(int indi, int indj, int indk, double &line_s, double &curv, int &kneib)
{
	indd = ((indi == indj) ? 1 : -1);
	double iBFr, jBFr, kBFr, iCFr, jCFr, kCFr;
	FindNeighbor110 ( k, i, j, indd, 2, K_FIRST, K_LAST, I_FIRST, I_LAST, J_FIRST, j_LAST, kneib);
	if (kneib != 0)
	{
		switch (uncode1)
		{
			case 0: {kBFr = B(KPos); iBFr = B(IPos); jBFr = B(JPos);}; break;
			case 1: {kBFr =     0.0; iBFr = (indd > 0) ? 0.0 : 1.0; jBFr =     0.5;}; break;
			case 2: {kBFr =     0.5; iBFr =     0.5; jBFr =     0.0;}; break;
		};		//		switch (uncode)
		x1 = B(DelXZ) * kBFr - A(DelXZ) * A(KPos);
		if (kn1 > k) x1 = x1 + A(DelXZ);  
		else if (kn1 < k) x1 = x1 - B(DelXZ);  
		yy1 = B(DelY) * jBFr - A(DelY) * A(JPos) +   Geometry.yMap[jn1] - Geometry.yMap[j];
	};
	if (kneib == 2)
	{
		switch (uncode2)
		{
			case 0: {kCFr = C(KPos); iCFr = C(IPos); jCFr = C(JPos);}; break;
			case 1: {kCFr =     0.0; iCFr = (indd > 0) ? 0.0 : 1.0; jCFr =     0.5;}; break;
			case 2: {kCFr =     0.5; iCFr =     0.5; jCFr =     0.0;}; break;
		};
		x2 = C(DelXZ) * kCFr - A(DelXZ) * A(KPos);
		if (kn2 > k) x2 = x2 + A(DelXZ);  
		else if (kn2 < k) x2 = x2 - C(DelXZ);  
		y2 = C(DelY) * jCFr - A(DelY) * A(JPos) +    Geometry.yMap[jn2] - Geometry.yMap[j];
//		printf("Ax=%7.5f Ay=%7.5f Ak=%7.5f  %2d %2d %2d  Bx=%4.3f By=%4.3f Bk=%4.3f C(st)=%d\n", A(KPos), A(IPos), A(JPos), in2, jn2, kn2, kCFr, iCFr, jCFr, C(State));	
	};
	if (kneib == 1)
	{
	  	FindBoundary110( k, i, j, indi, indd, indk, 2, A(DelXZ),
                                                   dia(A(KPos), ((indd > 0) ? A(IPos) : (1 - A(IPos)))),
				   										K_FIRST, K_LAST, K_EDGE,
					      								I_FIRST, I_LAST, I_EDGE,
                           A(DelY), A(JPos),J_FIRST, j_LAST, J_EDGE, kneib);
	  	if ((kneib == 2) && (in2 != -1))
     	{
		if ((kn2 == K_FIRST) || (kn2 == K_LAST-1))
			{
           		x2 = (kn2 == K_FIRST) ? A(DelZ) : -A(DelZ);
				x2 = C(DelXZ) * dia(C(KPos),((indd > 0) ? C(IPos) : (1 - C(IPos)))) -
   					 A(DelXZ) * dia(A(KPos),((indd > 0) ? A(IPos) : (1 - A(IPos)))) +
	   	             dia(                                     x2, (Geometry.xMap[in2] - Geometry.xMap[i])) * indd;
 				y2 = C(DelY) * C(JPos) - A(DelY) * A(JPos) +     Geometry.yMap[jn2] - Geometry.yMap[j];
			}
			else if ((in2 == I_FIRST) ||(in2 == I_LAST-1))
			{
           		x2 = (in2 == I_FIRST) ? A(DelX) : -A(DelX);
				x2 = C(DelXZ) * dia(C(KPos),((indd > 0) ? C(IPos) : (1 - C(IPos)))) -
   					 A(DelXZ) * dia(A(KPos),((indd > 0) ? A(IPos) : (1 - A(IPos)))) +
	   				 dia((Geometry.zMap[kn2] - Geometry.zMap[k]), x2) * indd;
 				y2 = C(DelY) * C(JPos) - A(DelY) * A(JPos) +     Geometry.yMap[jn2] - Geometry.yMap[j];
			}//end if
			else if ((jn2 == J_FIRST) || (jn2 == j_LAST-1))
			{
           		y2 = (jn2 == J_FIRST) ? A(DelY) : -A(DelY);
				x2 = C(DelXZ) * dia(C(KPos),((indd > 0) ? C(IPos) : (1 - C(IPos)))) -
   					 A(DelXZ) * dia(A(KPos),((indd > 0) ? A(IPos) : (1 - A(IPos)))) +
	   				 dia((Geometry.zMap[kn2] - Geometry.zMap[k]), (Geometry.xMap[in2] - Geometry.xMap[i])) * indd;
				y2 = C(DelY) * C(JPos) - A(DelY) * A(JPos) +    y2;
			};//end if
		};//end if
	};//endif
	if (kneib != 2)
		{curv = 0; line_s = A(DelY); alpY = pi/4;}
	else
		curv = Curv2d(indi, indk, A(DelXZ), A(DelY), line_s, alpY);
};//end of function