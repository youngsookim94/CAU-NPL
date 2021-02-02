#include <alis.h>


int main(int argc, char **argv)
{
	dom Dom = {{2250}, {2250}, {1000}};
	res Res = {10, {20}};
	sur Sur = {{SYM, PML}, {SYM, PML}, {SYM, PML}};
	world W = createWorld(Dom, Res, Sur, "%s", argv[0]);

	object Slab = {Difference, {2}, objects {
						{Box, {{-INF, INF}, {-INF, INF}, {-100, 100}}},
						{Difference, {2}, objects {
							{DLattice, {{500, 9}, {500*sqrtf(3), 5}}, objects {
								{RodZ, {{-175, 175}, {-175, 175}, {-INF, INF}}},
							}},
							{RodZ, {{-175, 175}, {-175, 175}, {-INF, INF}}},
						}},
					}};
	putObjects(W, n(3.4), Slab, Air);
	pointDipole(W, Ey, 0, 0, 0, Band, 750, 2000);

	slice XY = createSliceXY(W, 0);
	slice XZ = createSliceXZ(W, 0);
	slice YZ = createSliceYZ(W, 0);
	sliceSnap(W, LogRI, XY, png(dkbr, -1), "/%%");
	sliceSnap(W, LogRI, XZ, png(dkbr, -1), "/%%");
	sliceSnap(W, LogRI, YZ, png(dkbr, -1), "/%%");

	for (int n=1, N=750*2000/10/W->dt; timer(n, W->N+N); n++) {
		updateE(W);
		updateH(W);

		writeSpectrum(W, N, 750, 2000, "/JX", get(W, Jx, 0, 0, 0));
		writeSpectrum(W, N, 750, 2000, "/EX", get(W, Ex, 0, 0, 0));
		writeSpectrum(W, N, 750, 2000, "/JY", get(W, Jy, 0, 0, 0));
		writeSpectrum(W, N, 750, 2000, "/EY", get(W, Ey, 0, 0, 0));
		writeSpectrum(W, N, 750, 2000, "/JZ", get(W, Jz, 0, 0, 0));
		writeSpectrum(W, N, 750, 2000, "/EZ", get(W, Ez, 0, 0, 0));	

		writeSpectrum(W, N, 750, 2000, "/JE", get(W, JE, 200, 0, 0));	

		if (n > W->N) {
			writeRow(W, "/Time", W->dt*n, get(W, Ey, 0, 0, 0));
			writeSpectrum(W, N, 750, 2000, "/Spectrum", get(W, Ey, 0, 0, 0));
		}
		if (n+2*W->T > W->N+N) {
			sliceSnap(W, EE, XY, 15, png(hot,-1), "/%%/");
			sliceSnap(W, Ey, XY, 15, png(dkbr, -1), "/%%/");
			sliceSnap(W, HH, XY, 15, png(hot, -1), "/%%/");
			sliceSnap(W, Hz, XY, 15, png(dkbr, -1), "/%%/");
		}
	}
}
