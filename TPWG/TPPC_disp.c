#include <alis.h>


int main(int argc, char **argv)
{
	float k = atof(argv[1]);
	float w = 250.0;
	float d_Si = 108;
	float d_SiO2 = 256;
	float d_DBR = (d_Si + d_SiO2);

	dom Dom = {{-200,200}, {0}, {-1700, 1000}};
	res Res = {5, {10}, {1240}};
	sur Sur = {{BBC,1000/k}, {BBC}, {PML, PML}};
	world W = createWorld(Dom, Res, Sur, "%s-%03.1f", argv[0], k);

	object grating = {Box, {{ -0.5*w,  0.5*w}, {-INF, INF}, { 0,  50}}};
	object Substrate = {Box, {{-INF, INF}, {-INF, INF}, {-INF, 0}}};
	object Substrate2 = {Box, {{-INF, INF}, {-INF, INF}, {-INF, -4*d_DBR}}};	
	object Si1 = {Box, {{-INF, INF}, {-INF, INF}, {-d_Si, 0}}};	
	object Si2 = {Box, {{-INF, INF}, {-INF, INF}, {-d_Si-d_DBR, -d_DBR}}};	
	object Si3 = {Box, {{-INF, INF}, {-INF, INF}, {-d_Si-2*d_DBR, -2*d_DBR}}};	
	object Si4 = {Box, {{-INF, INF}, {-INF, INF}, {-d_Si-3*d_DBR, -3*d_DBR}}};		
	putObjects(W, Au, grating, n(3.482), Si4,n(3.482), Si3,n(3.482), Si2,n(3.482), Si1, n(3.482), Substrate2, n(1.5), Substrate, Air);

	pointDipole(W, Ey, 0, 0, -0.5*d_Si, Pulse, 4500, 4000);

	slice XY = createSliceXY(W, 0);
	slice XZ = createSliceXZ(W, 0);
	slice YZ = createSliceYZ(W, 0);
	sliceSnap(W, LogRI, XY, png(gray, 0), "/XY");
	sliceSnap(W, LogRI, YZ, png(gray, 0), "/YZ");
	sliceSnap(W, LogRI, XZ, png(gray, 0), "/XZ");	

	for (int n=1, N=1000*1500/2/W->dt; timer(n, W->N+N); n++) {
		updateE(W);
		updateH(W);

		writeRow(W, "/Totaltime", W->dt*n, get(W, Ey, 0, 0, -0.5*d_Si));

		if (n > W->N) {
			writeRow(W, "/EyMode", W->dt*n, get(W, Ey, 0, 0, -0.5*d_Si));
			writeRow(W, "/HzMode", W->dt*n, get(W, Ey, 0, 0, -0.5*d_Si));
			writeSpectrum(W, N, 500, 2500, "/SpectrumEy", get(W, Ey, 0, 0, -0.5*d_Si));
			writeSpectrum(W, N, 500, 2500, "/SpectrumHz", get(W, Ey, 0, 0, -0.5*d_Si));
		}
		if (N-n < 2*W->T) {
			sliceSnap(W, Ey, XZ, 50, png(dkbr, 0), "/XZ-Ey/");
			sliceSnap(W, Hz, XZ, 50, png(dkbr, 0), "/XZ-Hz/");
/*//		sliceFreqDom(W, rawEx, XZ, N, 625, h5, "/%%/");
//			sliceFreqDom(W, rawEz, XZ, N, 625, h5, "/%%/");
			sliceFreqDom(W, Ex, XZ, N, 1253, png(dkbr,0), "/%%/");
			sliceFreqDom(W, Ey, XZ, N, 1253, png(dkbr,0), "/%%/");
			sliceFreqDom(W, Hz, XZ, N, 1253, png(dkbr,0), "/%%/"); */
		}
	}
}
