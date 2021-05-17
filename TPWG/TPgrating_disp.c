#include <alis.h>


int main(int argc, char **argv)
{
	float k = atof(argv[1]);
	float length_Ag = 440; //atof(argv[2]);
	float w = 1250;
	float d_Si = 120;
	float d_SiO2 = 300;
	float d_DBR = (d_Si + d_SiO2);

	float nSi = 3.4684, nSiO2 = 1.4428;
	float aX = 520;

	dom Dom = {{-aX/2, aX/2}, {2500}, {-2000, 1700}};
	res Res = {2.5, {5, 25, 15}, {1240}};
	sur Sur = {{PBC, 2*PI*(1000-aX)/k}, {PML, PML}, {PML, PML}};
	world W = createWorld(Dom, Res, Sur, "%s_%03.1f_%.0f", argv[0], k, length_Ag);

	object grating = {Box, {{-length_Ag/2, length_Ag/2}, {-0.5 * w, 0.5 * w}, { 0,  50}}};
	object Substrate = {Box, {{-INF, INF}, {-0.5 * w, 0.5 * w}, {-INF, 0}}};
	object Substrate2 = {Box, {{-INF, INF}, {-INF, INF}, {-INF, -4*d_DBR}}};	
	object Si1 = {Box, {{-INF, INF}, {-0.5 * w, 0.5 * w}, {-d_Si, 0}}};	
	object Si2 = {Box, {{-INF, INF}, {-0.5 * w, 0.5 * w}, {-d_Si-d_DBR, -d_DBR}}};	
	object Si3 = {Box, {{-INF, INF}, {-0.5 * w, 0.5 * w}, {-d_Si-2*d_DBR, -2*d_DBR}}};	
	object Si4 = {Box, {{-INF, INF}, {-0.5 * w, 0.5 * w}, {-d_Si-3*d_DBR, -3*d_DBR}}};		
	putObjects(W, Ag, grating, n(nSi), Si4,n(nSi), Si3,n(nSi), Si2,n(nSi), Si1, n(nSi), Substrate2, n(nSiO2), Substrate, Air);

	pointDipole(W, Ey, 0, 0, -100, Band, 1500, 500);

	slice XY = createSliceXY(W, -100);
	slice XZ = createSliceXZ(W, 0);
	slice YZ = createSliceYZ(W, 0);
	sliceSnap(W, LogRI, XY, png(gray, 0), "/XY");
	sliceSnap(W, LogRI, YZ, png(gray, 0), "/YZ");
	sliceSnap(W, LogRI, XZ, png(gray, 0), "/XZ");	

	for (int n=1, N=1500*500/2.5/W->dt; timer(n, W->N+N); n++) {
		updateE(W);
		updateH(W);

		writeRow(W, "/Totaltime", W->dt*n, get(W, Ey, 0, 0, -100));

		if (n > W->N) {
			writeRow(W, "/EyMode", W->dt*n, get(W, Ey, 0, 0, -100));
			writeRow(W, "/EzMode", W->dt*n, get(W, Ez, 0, 0, -100));
			writeSpectrum(W, N, 500, 2500, "/SpectrumEy", get(W, Ey, 0, 0, -100));
			writeSpectrum(W, N, 500, 2500, "/SpectrumEz", get(W, Ez, 0, 0, -100));
		}
		if ( W->N+N-n < 2*W->T ) {
			sliceSnap(W, Ey, YZ, 25, png(dkbr, 0), "/YZ-Ey/");
			sliceSnap(W, EE, YZ, 25, png(hot, -1), "/YZ-EE/");
			sliceSnap(W, Ey, XY, 25, png(dkbr, 0), "/XY-Ey/");
			sliceSnap(W, EE, XY, 25, png(hot, -1), "/XY-EE/");
			sliceSnap(W, Ez, YZ, 25, png(dkbr, 0), "/YZ-Ez/");
			sliceSnap(W, Ey, XZ, 25, png(dkbr, 0), "/XZ-Ey/");
			sliceSnap(W, EE, XZ, 25, png(hot, -1), "/XZ-EE/");
/*//		sliceFreqDom(W, rawEx, XZ, N, 625, h5, "/%%/");
//			sliceFreqDom(W, rawEz, XZ, N, 625, h5, "/%%/");
			sliceFreqDom(W, Ex, XZ, N, 1253, png(dkbr,0), "/%%/");
			sliceFreqDom(W, Ey, XZ, N, 1253, png(dkbr,0), "/%%/");
			sliceFreqDom(W, Hz, XZ, N, 1253, png(dkbr,0), "/%%/"); */
		}
	}
}
