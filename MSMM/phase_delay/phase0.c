#include "alis.h"


int main(int argc, char **argv)
{
	
	float wBox = 500; 
	float hBox = 450;
	float tSlot = 40;
	float a = atof(argv[1]); // lattice constant
	float wLength = atof(argv[2]);

	//set calculation space
	res calc_Res = {5, {10}, {1240}};
	dom TF  = {{INF}, {INF}, {-INF, wLength + hBox}};
	dom calc_Dom = {{-a/2, a/2}, {0, 0}, {-wLength * 2 - hBox, hBox + wLength *2.5}};
	sur sur_BC = {{PBC}, {PBC}, {PML}, {24}};
	
	world calc_Wld = createWorld(calc_Dom, calc_Res, sur_BC, "%s_a%03.0f_l%.0f", argv[0], a, wLength);//output folder	

	object agBox = {Box, {{-wBox/2, wBox/2}, {-INF, INF}, {-hBox/2, hBox/2}}};
	object tio2Slot = {Box, {{-tSlot/2, tSlot/2}, {-INF, INF}, {-hBox/2, hBox/2}}};
	
	object AgSlot = {Difference, {2}, objects{agBox, tio2Slot}};

	putObjects (calc_Wld, Air);
	planewave(calc_Wld, TF, Ppol, 0, 0, Sine, wLength, 10);
	// insert the object in calculation space.
	// insert the light source as plane wave, Ex polarization.
	
	//declare variables in for loop
	phaser Phs = createPhaser(calc_Wld, wLength);
	planewave(calc_Wld, TF, Ppol, 0, 0, Sine, wLength, 10);

	float mesR = 0, mesT = 0, mesA = 0, mesTotal = 0;
	float out=0;
	
	//slice test
	slice XY = createSliceXY(calc_Wld, 0);
	slice XZ = createSliceXZ(calc_Wld, 0);
    slice YZ = createSliceYZ(calc_Wld, 0);

	sliceSnap(calc_Wld, LogRI, XY, png(jet, 2), "/%%");
	sliceSnap(calc_Wld, LogRI, XZ, png(jet, 2), "/%%");
    sliceSnap(calc_Wld, LogRI, YZ, png(jet, 2), "/%%");

	writeTxt(calc_Wld, "/RTA_ratio", "Wavelength\tReflection\tTransmission\tAbsorption\tTotal\r\n");
	
	//for loop for getting data
	for (int n=0, N=100*wLength/calc_Wld->dt; timer(n, N); n++) {
		updateH(calc_Wld);
		mesR += -poyntingZ(calc_Wld, hBox/2+150);
		mesT += -poyntingZ(calc_Wld, -250);
		mesA += 0; //objectAbsorption(calc_Wld, AgSlot);
		mesTotal = mesR + mesT + mesA;
		
		updateE(calc_Wld);
		mesR += -poyntingZ(calc_Wld, hBox/2+150);
		mesT += -poyntingZ(calc_Wld, -250);
		mesA += 0; //objectAbsorption(calc_Wld, AgSlot);
		mesTotal = mesR + mesT + mesA;

		if (N-n < calc_Wld -> T) updatePhaser(calc_Wld, Phs);
		
		if ( !(n%( calc_Wld ->T))) {
			writeRow(calc_Wld, "/RTA_ratio", calc_Wld->dt*n/100, mesR/mesTotal, mesT/mesTotal, mesA/mesTotal, mesTotal);
			mesR = 0;
			mesA = 0;
			mesT = 0;
		}
		
		//output images
		if (N-n <= 2*calc_Wld->T) {
			sliceSnap(calc_Wld, Ex, XZ, 20, png(dkbr, -1), "/XZ-Ex/");
			sliceSnap(calc_Wld, Hy, XZ, 20, png(dkbr, -1), "/XZ-Hy/");
			
			sliceSnap(calc_Wld, Ez, XZ, 20, png(dkbr, -1), "/XZ-Ez/");
			sliceSnap(calc_Wld, Hy, XZ, 20, png(dkbr, -1), "/XZ-Hy/");
			
			sliceSnap(calc_Wld, Ez, XY, 20, png(dkbr, -1), "/XY-Ez/");
			sliceSnap(calc_Wld, Hy, XY, 20, png(dkbr, -1), "/XY-Hy/");

			sliceSnap(calc_Wld, EE, XZ, 20, png(jet, -1), "/XZ-EE/");
			sliceSnap(calc_Wld, HH, XZ, 20, png(jet, -1), "/XZ-HH/");
						
			}
		
		//for getting all electric field data for upper side of entire structure
		if (n > N-1)
		{
		for(float i = - wLength*2 - hBox; i <=wLength * 2.5 + hBox ; i=i+calc_Wld->dz)
			{
				writeRow(calc_Wld, "/Total_EE", i, get(calc_Wld, EE, 0, 0, i));
				writeRow(calc_Wld, "/Total_Ex", i, get(calc_Wld, Ex, 0, 0, i));
				writeRow(calc_Wld, "/Total_Ez", i, get(calc_Wld, Ez, 0, 0, i));
				if(i > wLength * 0.25 && i < wLength * 2.25) writeRow(calc_Wld, "/Above_EE", i, get(calc_Wld, EE, 0, 0, i));
				if(i > wLength * 0.25 && i < wLength * 2.25) writeRow(calc_Wld, "/Above_Ex", i, get(calc_Wld, Ex, 0, 0, i));
				if(i > wLength * 0.25 && i < wLength * 2.25) writeRow(calc_Wld, "/Above_Ez", i, get(calc_Wld, Ez, 0, 0, i));
			}		
		}
	}
	farFieldProfile(calc_Wld, Phs, AzimuthalMap, 360, 180, png(jet, 0), "/Log-");
	farFieldProfile(calc_Wld, Phs, AzimuthalMap, 360, 180, txt, "/Log-");
	farFieldTheta(calc_Wld, Phs, 1, "/far_theta");
}
