#include <alis.h>

int main(int argc, char **argv)
{
	float wBox = atof(argv[1]);
	float hBox = atof(argv[2]);
	float tSlot = atof(argv[3]);
	float a = 1500; // domain scale 
	float wLength = atof(argv[4]); // wavelength

	res calc_Res = {2.5, {5}, {1240}}; // resolution of calculation
	dom TF = {{wBox/2 + 200}, {0}, {hBox/2 + 200}}; // set TF domain area
	dom calc_Dom = {{a/2}, {0}, {a/2}}; // set calc. domain
	sur sur_BC = {{PML}, {PBC}, {PML}, {24}}; //set B.C.
	world calc_Wld = createWorld(calc_Dom, calc_Res, sur_BC, "%s_w%.0f_h%.0f_t%.0f_l%.0f",argv[0], wBox, hBox, tSlot,  wLength);
	// define variables & define calcualtion space

	object agBox = {Box, {{-wBox/2, wBox/2}, {-INF, INF}, {-hBox/2, hBox/2}}};
	object tio2Slot = {Box, {{-tSlot/2, tSlot/2}, {-INF, INF}, {-hBox/2, hBox/2}}};
	object subPDMS = {Combination, {2}, objects{
		{Box, {{-INF, INF}, {-INF, INF}, {-INF, -hBox/2}}},
		{Box, {{-INF, INF}, {-INF, INF}, {hBox/2, INF}}}}};

	object AgSlot = {Difference, {2}, objects{agBox, tio2Slot}};

	putObjects (calc_Wld, Ag, AgSlot, n(2.0567), tio2Slot, n(1.4039), subPDMS, Air);
	planewave(calc_Wld, TF, Ppol, 0, 0, Sine, wLength, 10);
	// insert the object in calculation space.
	// insert the light source as plane wave, Ex polarization.

	float sca = 0;

	slice XY = createSliceXY(calc_Wld, 0);
	slice XZ = createSliceXZ(calc_Wld, 0);
    slice YZ = createSliceYZ(calc_Wld, 0);
	// define domain slice

	sliceSnap(calc_Wld, LogRI, XY, png(jet, 2), "/%%");
	sliceSnap(calc_Wld, LogRI, XZ, png(jet, 2), "/%%");
    sliceSnap(calc_Wld, LogRI, YZ, png(jet, 2), "/%%");
	// export sturcture images as grayscale

	writeTxt(calc_Wld, "/scattering", "Wavelength\tScattering_out\r\n");

	for (int n=1, N=100*wLength/calc_Wld -> dt; timer(n, N); n++) {
		updateH(calc_Wld);
		sca += poyntingOut (calc_Wld, -wBox/2-250, wBox/2+250, -100, 100, -hBox/2-250, hBox/2+250);

		updateE(calc_Wld);
		sca += poyntingOut (calc_Wld, -wBox/2-250, wBox/2+250, -100, 100, -hBox/2-250, hBox/2+250);

		if ( n%calc_Wld -> T == 0 ) {
			writeRow(calc_Wld, "/scattering", "%f\t%f\r\n", calc_Wld -> dt*n/100 , sca/calc_Wld->T );
			sca = 0;
			//exporting scattering cross-section data as txt file to roof of the executable file folder.
		}
		if (N-n <= 2*calc_Wld -> T) {

			sliceSnap(calc_Wld, Ex, XZ, 20, png(dkbr, -1), "/XZ-Ex/");
			sliceSnap(calc_Wld, ScEx, XZ, 20, png(dkbr, -1), "/XZ-ScEx/");
			
			sliceSnap(calc_Wld, EE, XZ, 20, png(jet, -1), "/XZ-EE/");
			sliceSnap(calc_Wld, ScEE, XZ, 20, png(jet, -1), "/XZ-ScEE/");
			//getting field data as png image or h5 or txt file.
		}
	}
}


