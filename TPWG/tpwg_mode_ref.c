#include <alis.h>


int main(int argc, char **argv)
{
	float k = atof(argv[1]); //wavenumber
	float w = 2500; //width of waveguide
	float dom_x = w; // domain parameter(x)
	float sZ = atof(argv[2]);
	float iLambda = atof(argv[3]);

	//Material Index
	float n_GaAs = 3.7;
	float n_AlAs = 3.0;

	//Structure parameters
	float d_GaAs = 57.4324;
	float d_AlAs = 70.8334;
	float period = d_GaAs + d_AlAs;
	float height = period * 24 + 2000;

	dom Dom = {{-dom_x, dom_x}, {0}, {-4179.75, 2325}};
	// domain
	res Res = {5, {25, 25, 10}, {1240}};
	//set Resolution
	sur Sur = {{SYM, PML},{PBC, 2*PI*1000/k},{PML}, {24}};
	// set Boundary condition
	
	world W = createWorld(Dom, Res, Sur, "%s_k%01.1f_w%03.0f_sZ%.0f_L%0.2f", argv[0], k, w, sZ, iLambda);
	
	//under thisiline : objects declaration 
	object layer_GaAs = {Box, {{-dom_x, dom_x}, {-INF, INF}, {-4179.75, 0}}};
	object layer_AlAs = {Box, {{-dom_x, dom_x}, {-INF, INF}, {-d_AlAs, 0}}};	

	object AlAs_lattice = 
		{Translation, {0, 0, -12*period}, 
			objects{
				{Lattice, {{0, 1}, {0, 1}, {period, 24}}, objects{layer_AlAs}}
			}
		};

    object GaAs_lattice = {Difference, {2}, objects {layer_GaAs, AlAs_lattice}}; 

    object ag_Plate = {Box, {{-1250, 1250}, {-INF, INF}, {0, 50}}};

	//end objects declaration  

	putObjects(W, Ag, ag_Plate, n(n_GaAs), GaAs_lattice, n(n_AlAs), AlAs_lattice, Air);
	//input object in the World

	pointDipole(W, Hz, 0, 0, -sZ, Pulse, iLambda, 1);
	
	slice XY = createSliceXY(W, 0);
	slice XZ = createSliceXZ(W, 0);
	slice YZ = createSliceYZ(W, 0);
	
	//slice test
	sliceSnap(W, LogRI, XY, png(gray, 0), "/XY");
	sliceSnap(W, LogRI, YZ, png(gray, 0), "/YZ");
	sliceSnap(W, LogRI, XZ, png(gray, 0), "/XZ");
	
	for (int n=1, N=750000; timer(n, N); n++) {
		updateH(W);
		updateE(W);
		
		writeRow(W, "/Totaltime", W->dt*n, get(W, Ez, 0, 0, -sZ));
		if (n > W->N) {
			writeRow(W,"/EzMode", W->dt*n, get(W, Ez, 0, 0, -sZ));
			//export data on txt file
			writeSpectrum(W, N, 100, 4500, "/Spectrum", get(W, Ez, 0, 0, -sZ));
			//export data on spectrum.txt 
		}
		if (N-n < 2*W->T) {
			sliceSnap(W, Ez, XY, 25, png(dkbr,-1), "/XY-Ez/");
			sliceSnap(W, Ez, XZ, 25, png(dkbr,-1), "/XZ-Ez/");
			sliceSnap(W, Ez, YZ, 25, png(dkbr,-1), "/YZ-Ez/");
			sliceSnap(W, Ey, XY, 25, png(dkbr,-1), "/XZ-Ey/");
			sliceSnap(W, Ey, XZ, 25, png(dkbr,-1), "/XZ-Ey/");
			sliceSnap(W, Ey, YZ, 25, png(dkbr,-1), "/YZ-Ey/");
			sliceSnap(W, EE, XY, 25, png(hot, -1), "/XY-EE/");
			sliceSnap(W, EE, XZ, 25, png(hot, -1), "/XZ-EE/");
			sliceSnap(W, EE, YZ, 25, png(hot, -1), "/YZ-EE/");
			//ploting image for each slices
		}
	}
}
