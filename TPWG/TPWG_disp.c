#include <alis.h>


int main(int argc, char **argv)
{
	float k = atof(argv[1]);
	float w = atof(argv[2]);
	float d_Si = 120;
	float d_SiO2 = 300;
	float d_DBR = (d_Si + d_SiO2);

	float nSi = 3.4684, nSiO2 = 1.4428;

	dom Dom = {{0}, {2500}, {-2000, 1700}};
	res Res = {7.5, {20, 25, 15}, {1240}};
	sur Sur = {{PBC, 2*PI*1000/k}, {SYM, PML}, {PML, PML}};
	world W = createWorld(Dom, Res, Sur, "%s_w%.0f_k%03.1f", argv[0], w, k);

	object grating = {Box, {{-INF, INF}, {-0.5 * w, 0.5 * w}, { 0,  50}}};
	object Substrate = {Box, {{-INF, INF}, {-0.5 * w, 0.5 * w}, {-INF, 0}}};
	object Substrate2 = {Box, {{-INF, INF}, {-INF, INF}, {-INF, -4*d_DBR}}};	
	object Si1 = {Box, {{-INF, INF}, {-0.5 * w, 0.5 * w}, {-d_Si, 0}}};	
	object Si2 = {Box, {{-INF, INF}, {-0.5 * w, 0.5 * w}, {-d_Si-d_DBR, -d_DBR}}};	
	object Si3 = {Box, {{-INF, INF}, {-0.5 * w, 0.5 * w}, {-d_Si-2*d_DBR, -2*d_DBR}}};	
	object Si4 = {Box, {{-INF, INF}, {-0.5 * w, 0.5 * w}, {-d_Si-3*d_DBR, -3*d_DBR}}};		
	putObjects(W, Ag, grating, n(nSi), Si4,n(nSi), Si3,n(nSi), Si2,n(nSi), Si1, n(nSi), Substrate2, n(nSiO2), Substrate, Air);

	pointDipole(W, Ey, 0, 0, -100, Band, 2000, 1500);

	slice XY = createSliceXY(W, 0);
	slice XZ = createSliceXZ(W, 0);
	slice YZ = createSliceYZ(W, 0);
	sliceSnap(W, LogRI, XY, png(gray, 0), "/XY");
	sliceSnap(W, LogRI, YZ, png(gray, 0), "/YZ");
	sliceSnap(W, LogRI, XZ, png(gray, 0), "/XZ");	

	for (int n=1, N=1600*600/7.5/W->dt; timer(n, W->N+N); n++) {
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
//			sliceSnap(W, Ey, YZ, 25, h5, "/YZ-Ey/");
			sliceSnap(W, EE, YZ, 25, png(dkbr, 0), "/YZ-EE/");
			sliceSnap(W, Ez, YZ, 25, png(dkbr, 0), "/YZ-Ez/");
//			sliceSnap(W, Ez, YZ, 25, h5, "/YZ-Ez/");
/*//		sliceFreqDom(W, rawEx, XZ, N, 625, h5, "/%%/");
//			sliceFreqDom(W, rawEz, XZ, N, 625, h5, "/%%/");
			sliceFreqDom(W, Ex, XZ, N, 1253, png(dkbr,0), "/%%/");
			sliceFreqDom(W, Ey, XZ, N, 1253, png(dkbr,0), "/%%/");
			sliceFreqDom(W, Hz, XZ, N, 1253, png(dkbr,0), "/%%/"); */
		}
	}
}
