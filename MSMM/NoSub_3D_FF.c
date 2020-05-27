#include <alis.h>

int main(int argc, char **argv)
{
	float wBox = 500;
	float hBox = 450;
	float tSlot = 40;
	float a = atof(argv[1]); // domain scale 
	float wLength = atof(argv[2]); // wavelength

	res calc_Res = {2.5, {5}, {1240}}; // resolution of calculation
	dom TF = {{INF}, {INF}, {-INF, hBox/2 + 200}}; // set TF domain area
	dom calc_Dom = {{a/2}, {a/2}, {1500}}; // set calc. domain
	sur sur_BC = {{SYM, PML}, {SYM, PML}, {PML}, {24}}; //set B.C.
	world calc_Wld = createWorld(calc_Dom, calc_Res, sur_BC, "%s_a%.0f_l%.0f",argv[0], a, wLength);
	// define variables & define calcualtion space

	object agBox = {Box, {{-wBox/2, wBox/2}, {-wBox/2, wBox/2}, {-hBox/2, hBox/2}}};
	object slotAra = {Combination, {2}, objects{
		{Box, {{-tSlot/2, tSlot/2}, {-wBox/2, wBox/2}, {-hBox/2, hBox/2}}},
		{Box, {{-wBox/2, wBox/2}, {-tSlot/2, tSlot/2}, {-hBox/2, hBox/2}}}
		}
	};

	object AgSlot = {Difference, {2}, objects{agBox, slotAra}};

	putObjects (calc_Wld, Ag, AgSlot, n(2.0567), slotAra, Air);
	planewave(calc_Wld, TF, Ppol, 0, 0, Sine, wLength, 10);
	// insert the object in calculation space.
	// insert the light source as plane wave, Ex polarization.


	phaser Phs = createPhaser(calc_Wld, wLength);
	planewave(calc_Wld, TF, Ppol, 0, 0, Sine, wLength, 10);

	float mesR = 0, mesT = 0, mesA = 0, mesTotal = 0;

	slice XY = createSliceXY(calc_Wld, 0);
	slice XZ = createSliceXZ(calc_Wld, 0);
    slice YZ = createSliceYZ(calc_Wld, 0);
	// define domain slice

	sliceSnap(calc_Wld, LogRI, XY, png(jet, 2), "/%%");
	sliceSnap(calc_Wld, LogRI, XZ, png(jet, 2), "/%%");
    sliceSnap(calc_Wld, LogRI, YZ, png(jet, 2), "/%%");
	// export sturcture images as grayscale

	writeTxt(calc_Wld, "/RTA_ratio", "Wavelength\tReflection\tTransmission\tAbsorption\tTotal\r\n");

	for (int n=1, N=100*wLength/calc_Wld -> dt; timer(n, N); n++) {
		updateH(calc_Wld);
		mesR += poyntingZ(calc_Wld, hBox/2+150);
		mesT += -poyntingZ(calc_Wld, -250);
		mesA += objectAbsorption(calc_Wld, AgSlot);
		mesTotal = mesR + mesT + mesA;

		updateE(calc_Wld);
		mesR += poyntingZ(calc_Wld, hBox/2+150);
		mesT += -poyntingZ(calc_Wld, -250);
		mesA += objectAbsorption(calc_Wld, AgSlot);
		mesTotal = mesR + mesT + mesA;

		if (N-n < calc_Wld -> T) updatePhaser(calc_Wld, Phs);

		if ( !(n%( calc_Wld ->T))) {
			writeRow(calc_Wld, "/RTA_ratio", calc_Wld->dt*n/100, mesR, mesT, mesA, mesTotal);
			mesR = 0;
			mesA = 0;
			mesT = 0;
		}
		if (N-n <= 2*calc_Wld -> T) {
			sliceSnap(calc_Wld, Ex, XZ, 20, png(dkbr, -1), "/XZ-Ex/");
			sliceSnap(calc_Wld, Ex, YZ, 20, png(dkbr, -1), "/YZ-Ex/");
			sliceSnap(calc_Wld, Ex, XY, 20, png(dkbr, -1), "/XY-Ex/");
			sliceSnap(calc_Wld, EE, XZ, 20, png(jet, -1), "/XZ-EE/");
			sliceSnap(calc_Wld, EE, YZ, 20, png(jet, -1), "/YZ-EE/");
			sliceSnap(calc_Wld, EE, XY, 20, png(jet, -1), "/XY-EE/");
			//getting field data as png image or h5 or txt file.
		}
	}
	farFieldProfile(calc_Wld, Phs, AzimuthalMap, 360, 180, png(hot, 0), "/Log-");
	farFieldProfile(calc_Wld, Phs, AzimuthalMap, 360, 180, txt, "/Log-");
	farFieldTheta(calc_Wld, Phs, 1, "/far_theta");
}


