#include <alis.h>

int main(int argc, char **argv)
{
    float k = atof(argv[1]) ; //wavevector, unit :2Pi/um
	float lambda = atof(argv[2]); // center wavelength
    float matlSel = atof(argv[3]); // 0 : Ag , 1 : Au

	float wx = atof(argv[4]); //width of dielectrics
    float wy = wx;
    float b0 = 150; //distance bet. surface of NW to PML (metal)

    float domLengX = (wx+b0*2)/2; //define domain scale parameter
    float domLengY = (wy+b0*2)/2;

    float sPos = wx*0.5-5;
    float resParam = wx*0.03;

    res Res = {resParam/2, {resParam}, {1240}}; //1 nm grid 
    dom Dom = {{domLengX}, {0}, {domLengY}};  
    sur Sur = {{SYM, PEC}, {PBC, 1000/k}, {SYM, PEC}}; //su MyOldCode 1000/k
    world W = createWorld(Dom, Res, Sur, "%s_lambdac%.0f_matl%.0f_k%.1f_w%.0f_dx%.1f", argv[0],lambda, matlSel, k, wx, resParam);
   
    //input object
    object metSide = {Box, {{-INF, INF}, {-INF, INF}, {-INF, INF}}}; //metal
    object NW_wg = {Box, {{-wx/2, wx/2}, {-INF, INF}, {-wy/2, wy/2}}}; //dielectric NW

    object Met_Slab = {Difference, {2}, objects {metSide, NW_wg}};

    // material def. for alis 1.0.2
    matter Drude_Ag = {{4.07666}, {9.2186, 0.02776}};
    matter Drude_Au = {{10.48449}, {9.0540, 0.07750}};

    //input objects in world
    if ( matlSel == 0 ) {
		putObjects(W, Drude_Ag, Met_Slab, n(2.6));
	} else if (matlSel == 1 ) {
		putObjects(W, Drude_Au, Met_Slab, n(2.6));
	}

    //pointDipole(W, Ex, sPos, 0, 0, Pulse, lambda, 100, 1); //input dipole source
    pointDipole(W, Ex, 0, 0, 0, Pulse, lambda, 100, 0); //inducing fundamental mode 

    slice XZ = createSliceXZ(W, 0);
    sliceSnap(W, LogRI, XZ, png(jet, 2), "/%%");
    sliceSnap(W, LogRI, XZ, txt, "/%%");

    for (int n = 1, N = 15000/W->dt; timer(n, W->N+N); n++) {//W->N+N
        updateH(W);
        updateE(W);

        writeRow(W, "/Totaltime", W->dt*n, get(W, Ex, sPos,0,0));
        
        if (n > W->N) {
            writeRow(W, "/ExMode", W->dt*n, get(W, Ex, sPos, 0, 0));
            writeSpectrum(W, N, 100, 4500, "/Spectrum", get(W, Ex, sPos, 0, 0));
            //export data Ex to raw txt files
        }

        if ( W->N+N-n < 2*W->T ) {//W->N+N
            sliceSnap(W, Ex, XZ, 25, png(dkbr, -1), "/XZ-Ex/");
//            sliceSnap(W, Ex, XZ, 25, h5, "/XZ-Ex/");
            sliceSnap(W, Hx, XZ, 25, png(dkbr, -1), "/XZ-Hx/");
//            sliceSnap(W, Hx, XZ, 25, h5, "/XZ-Hx/");
            sliceSnap(W, Ez, XZ, 25, png(dkbr, -1), "/XZ-Ez/");
//            sliceSnap(W, Ez, XZ, 25, h5, "/XZ-Ez/");
            sliceSnap(W, Hz, XZ, 25, png(dkbr, -1), "/XZ-Hz/");
//            sliceSnap(W, Hz, XZ, 25, h5, "/XZ-Hz/");
            sliceSnap(W, EE, XZ, 25, png(hot, -1), "/XZ-EE/");;
            //export mode image
        }
    }
     exec("cut -f2 %s/ExMode.txt | harminv -s amplitude -t %f %f-%f | tee %s/Hinv.txt", W->ID, W->dt, 1/(lambda+100), 1/(lambda-100), W->ID);

}
