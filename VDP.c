/*------------------------------------------------------------------------------*
 * File Name: VDP.c			 													*
 * Creation: 																	*
 * Purpose: OriginC Source C file												*
 * Copyright (c) ABCD Corp.	2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010		*
 * All Rights Reserved															*
 * 																				*
 * This file contains functions required to solve the VDP Resistance formula	*
 * 																				*
 * Modification Log:															*
 *------------------------------------------------------------------------------*/
 
////////////////////////////////////////////////////////////////////////////////////
// Including the system header file Origin.h should be sufficient for most Origin
// applications and is recommended. Origin.h includes many of the most common system
// header files and is automatically pre-compiled when Origin runs the first time.
// Programs including Origin.h subsequently compile much more quickly as long as
// the size and number of other included header files is minimized. All NAG header
// files are now included in Origin.h and no longer need be separately included.
//
// Right-click on the line below and select 'Open "Origin.h"' to open the Origin.h
// system header file.
#include <Origin.h>
////////////////////////////////////////////////////////////////////////////////////

//#pragma labtalk(0) // to disable OC functions for LT calling.
//#pragma labtalk(1,UserDefined)

////////////////////////////////////////////////////////////////////////////////////
// Include your own header files here.


////////////////////////////////////////////////////////////////////////////////////
// Start your functions here.

double VanDerPauw(double R1, double R2, double RS)
{
	return (exp(-pi*R1/RS) + exp(-pi*R2/RS));
}

double VDPderiv(double R1, double R2, double RS)
{
	return (pi*(R1*exp(-pi*R1/RS)+R2*exp(-pi*R2/RS))/pow(RS,2));
}

double VPDSheetRes(double VDPA, double VDPB, int numloops)
{
	double sheetres = pi*(VDPA+VDPB)/(2*ln(2));
	for (int i = 0; i < numloops; i++)
	{
		sheetres = sheetres + (1.0 - VanDerPauw(VDPA, VDPB, sheetres))/(VDPderiv(VDPA, VDPB, sheetres));
	}
	return sheetres;
}