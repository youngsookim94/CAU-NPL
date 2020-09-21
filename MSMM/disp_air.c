#include <alis.h>


int main(int argc, char **argv)
{
	float wavelength = 1550;
	float k = atof(argv[1]);
	float h = atof(argv[2]); ; //input your structure height(nm)
	float t = atof(argv[3]);  //input your slot thickness(nm)


	dom Dom = {{0}, {0}, {1000}};
	// domain
	res Res = {2.5, {5}, {1240}};
	//set Resolution dt = 2.5, dx = 5
	sur Sur = {{PBC},{PBC, 2*PI*1000/k},{SYM, PML}, {24}};
	// set Boundary condition
	
	world W = createWorld(Dom, Res, Sur, "%s_%.0f_%.0f_%.0f", argv[0], k, h, t);
	
	object si = {Box, {{-INF, INF}, {-INF, INF}, {t/2, h/2}}};
	//define waveguide structure

	object sio2 = {Box, {{-INF, INF}, {-INF, INF}, {-5, t/2}}};
	//define waveguide slot

	putObjects(W, Air, sio2, Ag, si, Air);
	//input Object wg & slot

	pointDipole(W, Ez, 0, 0, 0, Pulse, 2500, 2000);
	//input Dipole Source

	slice XY = createSliceXY(W, 0);
	slice XZ = createSliceXZ(W, 0);
	slice YZ = createSliceYZ(W, 0);
	
	//slice test
	sliceSnap(W, LogRI, XY, png(dkbr,-1), "/%%");
	sliceSnap(W, LogRI, YZ, png(dkbr,-1), "/%%");
	sliceSnap(W, LogRI, XZ, png(dkbr,-1), "/%%");
	
	for (int n=1, N=15000; timer(n, N); n++) {
		updateH(W);
		updateE(W);
		
		writeRow(W, "/Totaltime", W->dt*n, get(W, Ez, 0, 0, 0));
		if (n > W->N) {
			writeRow(W,"/EzMode", W->dt*n, get(W, Ez, 0, 0, 0));
			//export data on txt file
			writeSpectrum(W, N, 100, 4500, "/Spectrum", get(W, Ez, 0, 0, 0));
			//export data on spectrum.txt 
		}
		if (N-n < 2*W->T) {
			//sliceSnap(W, Ez, XY, 25, png(dkbr,-1), "/%%/");
			sliceSnap(W, Ez, XZ, 25, png(dkbr,-1), "/%%/");
			sliceSnap(W, Ez, YZ, 25, png(dkbr,-1), "/%%/");
			//ploting image for each slices
		}
	}
}
