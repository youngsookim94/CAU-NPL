#include <alis.h>


int main(int argc, char **argv)
{
	float resolution = atof(argv[1]); 
	float R = atof(argv[2]);	
	float WC = atof(argv[3]);	
	float DW = atof(argv[4]);	
	float xc = atof(argv[5]);		
	float dom_x = R + 1000;

	dom Dom = {{dom_x}, {dom_x}, {0}};
	res Res = {0.5*resolution, {resolution}, {1240}};
	sur Sur = {{SYM,PML}, {SYM,PML}, {SYM}};
	world W = createWorld(Dom, Res, Sur, "%s_R%01.1f_res%03.0f_WC%03.0f_DW%03.0f_xc%03.0f", argv[0],R,resolution, WC, DW,xc);

	object disk  = {RodZ, {{ -R,  R}, {-R,  R}, {-INF, INF}}};
	object disk_air  = {RodZ, {{ -(R-300),  (R-300)}, {-(R-300),  (R-300)}, {-INF, INF}}};	
	putObjects(W,  n(3.3), disk, Air);//	putObjects(W, n(1.0), disk_air, n(3.3), disk, Air);

	pointDipole(W, Ex, R-xc, 0, 0, Pulse, WC,DW);
	slice XY = createSliceXY(W, 0);
	sliceSnap(W, LogRI, XY, png(gray,0), "/%%");	

	for (int n=1, N=32768*16; timer(n, W->N+N); n++) {
		updateE(W);
		updateH(W);

		writeRow(W, "/Totaltime", W->dt*n, get(W, Ex, R-xc, 0, 0));
		writeSpectrum(W, N, 500, 1500, "/JX", get(W, Jx, R-xc, 0, 0));
		writeSpectrum(W, N, 500, 1500, "/EX", get(W, Ex, R-xc, 0, 0));
		writeSpectrum(W, N, 500, 1500, "/JY", get(W, Jy, R-xc, 0, 0));
		writeSpectrum(W, N, 500, 1500, "/EY", get(W, Ey, R-xc, 0, 0));
		writeSpectrum(W, N, 500, 1500, "/JZ", get(W, Jz, R-xc, 0, 0));
		writeSpectrum(W, N, 500, 1500, "/EZ", get(W, Ez, R-xc, 0, 0));	
			
		if (n > W->N) {
			writeRow(W,"/ExMode", W->dt*n, get(W, Ex, R-80, 0, 0));
			//export data on txt file
			writeSpectrum(W, N, 500, 1500, "/Spectrum", get(W, Ex, R-80, 0, 0));
			
			
			//export data on spectrum.txt 
		}
		if (N-n < 2*W->T) {
			sliceSnap(W, EE, XY, 25, png(hot, -1), "/XY-EE/");
			//ploting image for each slices
		}
	}


}
