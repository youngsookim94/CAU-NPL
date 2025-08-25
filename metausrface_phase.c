#include "alis.h"


int main(int argc, char **argv)
{
	
	float h = 260; // Si height
	float w = 520; // Si box x-width
	float a = atof(argv[1]); // lattice constant
	float h_SiO2 = 500; // SiO2 thickness (substrate)
	float Wavelength = atof(argv[2]);

	//set parameters
	float idx_Si = 3.4777;
	float idx_SiO2 = 1.440;
	float dAg = 300;
	
	//set calculation space
	res Res = {5, {10}, {1240}};
	dom TF  = {{INF}, {INF}, {-INF, Wavelength * 2.5 - 200}};
	dom Dom = {{-a/2, a/2}, {-a/2, a/2}, {-h-h_SiO2-dAg, Wavelength *2.5}};
	sur Sur = {{PBC}, {PBC}, {PML}, {24}};
	
	world W = createWorld(Dom, Res, Sur, "phase0_a%03.0f_lambda%04.f", a, Wavelength);//output folder	

	//set objects
//	object Si = {Box, {{-w/2, w/2}, {-w/2, w/2}, {-h, 0}}};
//	object SiO2 = {Box, {{-INF, INF}, {-INF, INF}, {-h - h_SiO2, -h}}};
//	object Ag_sub = {Box, {{-INF, INF}, {-INF, INF}, {-INF, -h - h_SiO2}}};
	
	//input objects and source
	putObjects(W, Air);
	planewave(W, TF, Ppol, 0, 0, Sine, Wavelength, 10);
	
	//declare variables in for loop
	float abs_Ag = 0;
	float out=0;
	float ref1=0;
	float tsm2=0;
	
	//slice test
	slice XZ = createSliceXZ(W, 0);
	slice XY = createSliceXY(W, 0);
	sliceSnap(W, LogRI, XZ, png(jet,2), "/XZ");
	sliceSnap(W, LogRI, XY, png(jet,2), "/XY");

	writeTxt(W, "/data", "Wavelength\tAbsorption\tReflection\tTransmission\tTotal\r\n");//output file
	
	//for loop for getting data
	for (int n=0, N=300*Wavelength/W->dt; timer(W->ID, n, N); n++) {
		updateH(W);
//		abs_Ag += objectAbsorption(W, Ag_sub);

		ref1 += poyntingZ  (W, Wavelength * 2.5 - 100);
		tsm2 += -poyntingZ  (W, -h - h_SiO2 - dAg + 100);
		
	    updateE(W);
//		abs_Ag += objectAbsorption(W, Ag_sub);

		ref1 += poyntingZ  (W, Wavelength * 2.5 -100);
		tsm2 += -poyntingZ  (W, -h - h_SiO2 -dAg +100);
		
		if (!(n%(W->T))) { 	
			writeTxt(W, "/data", "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\r\n", W->dt*n/300, abs_Ag/W->T, ref1/W->T, tsm2/W->T, (abs_Ag+ref1+tsm2)/W->T);//simultion data
			abs_Ag = 0;
			ref1 = 0;
			tsm2 = 0;
		}
		
		//output images
		if (N-n <= 2*W->T) {
			sliceSnap(W, Ex, XZ, 20, png(dkbr, -1), "/XZ-Ex/");
			sliceSnap(W, Hy, XZ, 20, png(dkbr, -1), "/XZ-Hy/");
			
			sliceSnap(W, Ez, XZ, 20, png(dkbr, -1), "/XZ-Ez/");
			sliceSnap(W, Hy, XZ, 20, png(dkbr, -1), "/XZ-Hy/");

			sliceSnap(W, Ex, XY, 20, png(dkbr, -1), "/XY-Ex/");
			sliceSnap(W, Hy, XY, 20, png(dkbr, -1), "/XY-Hy/");
			
			sliceSnap(W, Ez, XY, 20, png(dkbr, -1), "/XY-Ez/");
			sliceSnap(W, Hy, XY, 20, png(dkbr, -1), "/XY-Hy/");

			sliceSnap(W, EE, XY, 20, png(jet, -1), "/XY-EE/");
			sliceSnap(W, HH, XY, 20, png(jet, -1), "/XY-HH/");
			sliceSnap(W, EE, XZ, 20, png(jet, -1), "/XZ-EE/");
			sliceSnap(W, HH, XZ, 20, png(jet, -1), "/XZ-HH/");
						
			}
		
		//for getting all electric field data for upper side of entire structure
		if (n > N-1)
		{
		for(float i = 0; i <=Wavelength * 2.5; i=i+W->dz)
			{
				writeRow(W, "/Total_EE", i, get(W, EE, 0, 0, i));
				writeRow(W, "/Total_Ex", i, get(W, Ex, 0, 0, i));
				writeRow(W, "/Total_Ez", i, get(W, Ez, 0, 0, i));
				if(i > Wavelength * 0.25 && i < Wavelength * 2.25) writeRow(W, "/Above_EE", i, get(W, EE, 0, 0, i));
				if(i > Wavelength * 0.25 && i < Wavelength * 2.25) writeRow(W, "/Above_Ex", i, get(W, Ex, 0, 0, i));
				if(i > Wavelength * 0.25 && i < Wavelength * 2.25) writeRow(W, "/Above_Ez", i, get(W, Ez, 0, 0, i));
			}		
		}
	}
}
