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

    float sPos = wx*0.5-5;
    float resParam = wx*0.03;

    float nLoop = atof(argv[3]);

    res Res = {resParam/2, {resParam, resParam, 4*resParam}, {1240}}; //variable grid by width of wg 
    dom Dom = {{domLengX}, {domLengY}, {-2000, 3000}};  
    sur Sur = {{SYM, PEC}, {SYM, PEC}, {PML}, {24}}; //su MyOldCode 1000/k
    world W = createWorld(Dom, Res, Sur, "%s_l%.0f_w%.0f_dx%.1f_N%.0f", argv[0],lambda, wx, resParam, nLoop);
   
    //input object
    object Ag_Side = {Box, {{-INF, INF}, {-INF, INF}, {-INF, INF}}}; //metal
//    object Au_Side = {Box, {{-INF, INF}, {-INF, INF}, {-INF, 0}}}; //metal
    object NW_wg = {Box, {{-wx/2, wx/2}, {-wy/2, wy/2}, {-INF, INF}}}; //dielectric NW

	object Ag_wire = {Difference, {2}, objects {Ag_Side, NW_wg}};
//	object Au_wire = {Difference, {2}, objects {Au_Side, NW_wg}};

    // material def. for alis 1.0.2
    matter Drude_Ag = {{4.07666}, {9.2186, 0}};
    matter Drude_Au = {{10.48449}, {9.0540, 0}};

    //input objects in world
 	putObjects(W, Drude_Ag, Ag_wire, n(2.6));

    guidedWaveZ(W, "Guided_Mode_Profiles/Ex2", Ey, 3000, Sine, lambda, 100, 1, 0);
    guidedWaveZ(W, "Guided_Mode_Profiles/Ez2", Ex, 3000, Sine, lambda, 100, 1, 0);

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
	float AuOut = 0, AuOut2 = 0, AuOut3 = 0, AuOut4 = 0;
	float AgOut = 0;

	writeTxt(W, "/Ratio", "Wavelength\tAgOut\tAuOut\tTotal\r\n");

    for (int n = 1, N = nLoop*lambda/W->dt; timer(n, W->N+N); n++) {//W->N+N
        updateH(W);
		AgOut = fabs(poyntingZ(W, 200, -wx/2, wx/2, -wx/2, wx/2));
		AuOut = fabs(poyntingZ(W, -200, -wx/2, wx/2, -wx/2, wx/2)); 
		AuOut2 = fabs(poyntingZ(W, -400, -wx/2, wx/2, -wx/2, wx/2)); 
		AuOut3 = fabs(poyntingZ(W, -800, -wx/2, wx/2, -wx/2, wx/2)); 
		AuOut4 = fabs(poyntingZ(W, -1600, -wx/2, wx/2, -wx/2, wx/2)); 
		//totalOut += poyntingOut(W, -wx/2, wx/2, -wx/2, wx/2, -700, -690);

        updateE(W);
		AgOut = fabs(poyntingZ(W, 200, -wx/2, wx/2, -wx/2, wx/2));
		AuOut = fabs(poyntingZ(W, -200, -wx/2, wx/2, -wx/2, wx/2));
        AuOut2 = fabs(poyntingZ(W, -400, -wx/2, wx/2, -wx/2, wx/2)); 
		AuOut3 = fabs(poyntingZ(W, -800, -wx/2, wx/2, -wx/2, wx/2)); 
		AuOut4 = fabs(poyntingZ(W, -1600, -wx/2, wx/2, -wx/2, wx/2)); 
		//totalOut += poyntingOut(W, -wx/2, wx/2, -wx/2, wx/2, -700, -690);

		if (!(n%(W->T))) {
			writeRow(W, "/Ratio", W->dt*n/nLoop, (AgOut), (AuOut), (AuOut2), (AuOut3), (AuOut4));
			AgOut = 0;
            AuOut = 0;
			totalOut = 0;
		}
		
        writeRow(W, "/Totaltime", W->dt*n, get(W, Ex, sPos, 3000 ,0));
        
        if (n > W->N) {
            writeRow(W, "/ExMode", W->dt*n, get(W, Ex, sPos, 3000, 0));
            writeSpectrum(W, N, 100, 4500, "/Spectrum", get(W, Ex, sPos, 3000, 0));
            //export data Ex to raw txt files
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
