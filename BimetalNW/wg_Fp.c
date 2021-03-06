#include <alis.h>
#include <math.h>

int main(int argc, char **argv)
{
	float lambda = atof(argv[1]); // center wavelength

	float wx = atof(argv[2]); //width of dielectrics
    float wy = wx;
    float b0 = 150; //distance bet. surface of NW to PML (metal)

    float domLengX = (wx+b0*2)/2; //define domain scale parameter
    float domLengY = (wy+b0*2)/2;

    float sPosX = 0;
    float resParam = 2;

    float sPosZ = atof(argv[3]);

    res Res = {resParam/2, {resParam, resParam, 4*resParam}, {1240}, {{2, -wx/2, wx/2}, {2, -wy/2, wy/2}, {1, -INF, INF}}}; //variable grid by width of wg 
    dom Dom = {{domLengX}, {domLengY}, {-2000, 3000}};  
    sur Sur = {{SYM, PEC}, {SYM, PEC}, {PML}, {24}}; //su MyOldCode 1000/k
    world W = createWorld(Dom, Res, Sur, "%s_l%.0f_w%.0f_dx%.1f_X%.0f_Z%.0f", argv[0],lambda, wx, resParam, sPosX, sPosZ);
   
    //input object
    object Ag_Side = {Box, {{-INF, INF}, {-INF, INF}, {-INF, 0}}}; //metal
    object Au_Side = {Box, {{-INF, INF}, {-INF, INF}, {0, INF}}}; //metal
    object NW_wg = {Box, {{-wx/2, wx/2}, {-wy/2, wy/2}, {-INF, INF}}}; //dielectric NW

	object Ag_wire = {Difference, {2}, objects {Ag_Side, NW_wg}};
	object Au_wire = {Difference, {2}, objects {Au_Side, NW_wg}};

//    material def. for alis 1.0.2
    matter Drude_Ag = {{4.07666}, {9.2186, 0.02776}};
    matter Drude_Au = {{10.48449}, {9.0540, 0.07750}};

//    matter Drude_Ag = {{4.07666}, {9.2186, 0}};
//    matter Drude_Au = {{10.48449}, {9.0540, 0}};

    //input objects in world
 	putObjects(W, Drude_Ag, Ag_wire, Drude_Au, Au_wire, n(2.6));

    pointDipole(W, Ex, sPosX, 0, sPosZ, Band, lambda, 200, 0); //inducing fundamental mode 

    slice XZ = createSliceXZ(W, 0);
	slice XY = createSliceXY(W, 0);
    slice YZ = createSliceYZ(W, 0);
    sliceSnap(W, LogRI, XZ, png(jet, 2), "/%%");
    sliceSnap(W, LogRI, XZ, txt, "/%%");
	sliceSnap(W, LogRI, XY, png(jet, 2), "/%%");
	sliceSnap(W, LogRI, XY, txt, "/%%");
	sliceSnap(W, LogRI, YZ, png(jet, 2), "/%%");
    sliceSnap(W, LogRI, YZ, txt, "/%%");

	float totalOut = 0;
	float AgOut = 0, AgOut2 = 0, AgOut3 = 0, AgOut4 = 0;
	float AuOut = 0;

    for (int n = 1, N = lambda*100/resParam/2/W->dt; timer(n, W->N+N); n++) {//W->N+N
        updateH(W);
	    updateE(W);

        writeSpectrum(W, N, 600, 1200, "/JX", get(W, Jx, sPosX, 0, sPosZ));
		writeSpectrum(W, N, 600, 1200, "/EX", get(W, Ex, sPosX, 0, sPosZ));
		writeSpectrum(W, N, 600, 1200, "/JY", get(W, Jy, sPosX, 0, sPosZ));
		writeSpectrum(W, N, 600, 1200, "/EY", get(W, Ey, sPosX, 0, sPosZ));
		writeSpectrum(W, N, 600, 1200, "/JZ", get(W, Jz, sPosX, 0, sPosZ));
		writeSpectrum(W, N, 600, 1200, "/EZ", get(W, Ez, sPosX, 0, sPosZ));	

		writeSpectrum(W, N, 400, 2000, "/JE", get(W, JE, sPosX, 0, sPosZ));	

		if (n > W->N) {
			writeRow(W, "/Time", W->dt*n, get(W, Ex, sPosX, 0, sPosZ));
			writeSpectrum(W, N, 600, 1200, "/Spectrum", get(W, Ex, sPosX, 0, sPosZ));
        }
        if ( W->N+N-n < 2*W->T ) {//W->N+N
            sliceSnap(W, Ex, XZ, 25, png(dkbr, -1), "/XZ-Ex/");
            sliceSnap(W, Ey, XZ, 25, png(dkbr, -1), "/XZ-Ey/");
            sliceSnap(W, EE, XZ, 25, png(hot, -1), "/XZ-EE/");
            sliceSnap(W, Ex, XY, 25, png(dkbr, -1), "/XY-Ex/");
            sliceSnap(W, Ey, XY, 25, png(dkbr, -1), "/XY-Ey/");
            sliceSnap(W, EE, XY, 25, png(hot, -1), "/XY-EE/");
            //export mode image
        }
    }
}
