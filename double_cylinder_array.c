#include "alis.h"

int main(int argc, char **argv)
{
	//r = diameter 
	float r = 440;
	float a = 658.47;
	float H = 740;
	float H1 = H+440;
	float Wavelength = 1514;
	float delLamb = 50;

	res Res = {0.1646175, {0.329235, 0.329235, 20}, {1240}};
	dom TF  = {{-a/2-100, a/2+100}, {-a/2-100, a/2+100}, {-INF, 1500}};
	dom Dom = {{-a/2, a/2}, {-a/2, a/2}, {-2000, 2000}};
	sur Sur = {{SYM, PBC}, {SYM, PBC}, {PML, PML}, {24}};
	world W = createWorld(Dom, Res, Sur, "a%.2f_H%.0f_lamb%04.0f",a ,H ,Wavelength);		

	object disk1 = {RodZ, {{-r/2, r/2}, {-r/2, r/2}, {H/2, r+H/2}}};
	object disk2 = {RodZ, {{-r/2, r/2}, {-r/2, r/2}, {-r-H/2,-H/2}}};
	
	putObjects(W, n(3.479), disk1, n(3.479), disk2 , Air);
	planewave(W, TF, Ppol, 0, 0, Pulse, Wavelength, delLamb);
	
	slice XZ = createSliceXZ(W, 0);
	sliceSnap(W, LogRI, XZ, png(jet,2), "XZ");
	
	for (int n=0, N=Wavelength * delLamb / W->dx / W->dt; timer(n, N); n++) {
		updateH(W);
	    updateE(W);
		writeRow(W, "/TotalTime", W->dt*n, get(W, Ex, 0, 0, 0));
		
		if (n > W->N) {
//			writeRow(W, "/TimeH", W->dt*n, get(W, Ex, 0, 0, 1700));
//		    writeRow(W, "/TimeM", W->dt*n, get(W, Ex, 0, 0, 0));
//		    writeRow(W, "/TimeD", W->dt*n, get(W, Ex, 0, 0, -1700));
            writeSpectrum(W, N, 100, 4500, "/Spectrum", get(W, Ex, 0, 0, 0));
        }
		
		if (N-n <= 2*W->T) {
			sliceSnap(W, Ex, XZ, 50, png(dkbr, -1), "/XZ-Ex/");
			sliceSnap(W, EE, XZ, 50, png(jet, -1), "/XZ-EE/");
			sliceSnap(W, ScEx, XZ, 50, png(dkbr, -1), "/XZ-ScEx/");	
			writeRow(W, "/Time", W->dt*n, get(W, Ex, 0, 0, 0));
	
		}
	}
}
