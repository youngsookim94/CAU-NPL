#include <alis.h>


int main(int argc, char **argv)
{
	float k = atof(argv[1]); //wavenumber
	float w = 2500; //width of waveguide
	float dom_x = w; // domain parameter(x)
	float sZ = atof(argv[2]);

	//Material Index
	float n_Si = 3.4777;
	float n_SiO2 = 1.477;

	//Structure parameters
	float d_Si = 97.607;
	float d_SiO2 = 265.684;
	float period = d_Si + d_SiO2;
	float height = period * 6 + 2000;

	dom Dom = {{-dom_x, dom_x}, {0}, {-4179.75, 2325}};
	// domain
	res Res = {5, {25, 25, 10}, {1240}};
	//set Resolution
	sur Sur = {{SYM, PML},{PBC, 2*PI*1000/k},{PML}, {24}};
	// set Boundary condition
	
	world W = createWorld(Dom, Res, Sur, "%s_k%01.1f_w%03.0f_sZ%.0f", argv[0], k, w, sZ);
	
	//under thisiline : objects declaration 
	object layer_Si = {Box, {{-dom_x, dom_x}, {-INF, INF}, {-4179.75, 0}}};
	object layer_SiO2 = {Box, {{-dom_x, dom_x}, {-INF, INF}, {-d_SiO2, 0}}};	

	object SiO2_lattice = 
		{Translation, {0, 0, -3*period+d_Si}, 
			objects{
				{Lattice, {{0, 1}, {0, 1}, {period, 6}}, objects{layer_SiO2}}
			}
		};

    object Si_lattice = {Difference, {2}, objects {layer_Si, SiO2_lattice}}; 

    object ag_Plate = {Box, {{-1250, 1250}, {-INF, INF}, {0, 50}}};

	//end objects declaration  

	putObjects(W, Ag, ag_Plate, n(n_Si), Si_lattice, n(n_SiO2), SiO2_lattice, Air);
	//input object in the World

	pointDipole(W, Ez, 0, 0, -sZ, Pulse, 900, 100);
	
	slice XY = createSliceXY(W, 0);
	slice XZ = createSliceXZ(W, 0);
	slice YZ = createSliceYZ(W, 0);
	
	//slice test
	sliceSnap(W, LogRI, XY, png(gray, 0), "/XY");
	sliceSnap(W, LogRI, YZ, png(gray, 0), "/YZ");
	sliceSnap(W, LogRI, XZ, png(gray, 0), "/XZ");
	
	for (int n=1, N=500000; timer(n, N); n++) {
		updateH(W);
		updateE(W);
		
		writeRow(W, "/Totaltime", W->dt*n, get(W, Ez, 0, 0, -sZ));
		if (n > W->N) {
			writeRow(W,"/ExMode", W->dt*n, get(W, Ez, 0, 0, -sZ));
			//export data on txt file
			writeSpectrum(W, N, 100, 4500, "/Spectrum", get(W, Ez, 0, 0, -sZ));
			//export data on spectrum.txt 
		}
		if (N-n < 2*W->T) {
			sliceSnap(W, Ez, XY, 25, png(dkbr,-1), "/XY-Ex/");
			sliceSnap(W, Ez, XZ, 25, png(dkbr,-1), "/XZ-Ex/");
			sliceSnap(W, Ez, YZ, 25, png(dkbr,-1), "/YZ-Ex/");
			sliceSnap(W, EE, XY, 25, png(hot, -1), "/XY-EE/");
			sliceSnap(W, EE, XZ, 25, png(hot, -1), "/XZ-EE/");
			sliceSnap(W, EE, YZ, 25, png(hot, -1), "/YZ-EE/");
			//ploting image for each slices
		}
	}
}
