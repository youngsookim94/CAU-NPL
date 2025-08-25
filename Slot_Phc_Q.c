#include <alis.h>


int main(int argc, char **argv)
{
	float a = 500 ;
	float b =  atof(argv[1]) ; //radius chages by argument
	float r = b * a ;
	float r_0 = 0.25*a ; 
	
	//Parameters of Structure
	float h = 1000 ; //set height of cylinder
	float t =  10  ; //set slot thickness
	float d =  atof(argv[2]) ;

	dom Dom = {{10.5*a}, {10.5*a}, {h/2 + 700}};
	res Res = {2.5, {5}, {}, {{4, -9*a ,9*a}, {4, -9*a, 9*a}, {4, -400, 400}}};
	sur Sur = {{SYM, PML}, {SYM, PML}, {SYM, PML}, {96}} ;
	world W = createWorld(Dom, Res, Sur, "%s-r%.2f-d%3f", argv[0], b, d );
	//Set calculation domains & Resolution & Boundary conditition
	//Create World

	object cylinder = {RodZ, {{-r_0, r_0}, {-r_0, r_0}, {t/2, h/2}}};
	object slot = {RodZ, {{-r_0, r_0}, {-r_0, r_0}, {-t/2, t/2}}};
	//define object cylinder & slot

	object square_lattice =
	{Difference, {6}, objects 
		{
			{Lattice, {{a, 19}, {a, 19}, {a, 1}}, objects{cylinder}},
			{RodZ, {{-a/2, a/2}, {-a/2, a/2}, {-INF, INF}}},

			{RodZ, {{-3*a/2, -a/2}, {-a/2, a/2}, {-INF, INF}}},
			{RodZ, {{ 3*a/2, a/2}, {-a/2, a/2}, {-INF, INF}}},
			{RodZ, {{-a/2, a/2}, {-3*a/2, -a/2}, {-INF, INF}}},
			{RodZ, {{-a/2, a/2}, {3*a/2, a/2}, {-INF, INF}}}
		}
	};
	object slot_lattice =
	{Difference, {6}, objects 
		{
			{Lattice, {{a, 19}, {a, 19}, {a, 1}}, objects{slot}},
			{RodZ, {{-a/2, a/2}, {-a/2, a/2}, {-INF, INF}}},

			{RodZ, {{-3*a/2, -a/2}, {-a/2, a/2}, {-INF, INF}}},
			{RodZ, {{ 3*a/2, a/2}, {-a/2, a/2}, {-INF, INF}}},
			{RodZ, {{-a/2, a/2}, {-3*a/2, -a/2}, {-INF, INF}}},
			{RodZ, {{-a/2, a/2}, { 3*a/2, a/2}, {-INF, INF}}}
		}
	};
	//define Slot incerted square lattice

	object defect = 
	{Combination, {5}, objects
		{
			{RodZ, {{-r, r}, {-r, r}, {t/2, h/2}}},
			{RodZ, {{-r_0 -a +d, r_0 -a +d}, {-r_0, r_0}, {t/2, h/2}}},
			{RodZ, {{-r_0 +a -d, r_0 +a -d}, {-r_0, r_0}, {t/2, h/2}}},
			{RodZ, {{-r_0, r_0}, {-r_0 +a -d, r_0 +a -d}, {t/2, h/2}}},
			{RodZ, {{-r_0, r_0}, {-r_0 -a +d, r_0 -a +d}, {t/2, h/2}}}
		}
	};

	object defect_slot = 
	{Combination, {5}, objects
		{
			{RodZ, {{-r, r}, {-r, r}, {-t/2, t/2}}},
			{RodZ, {{-r_0 -a +d, r_0 -a +d}, {-r_0, r_0}, {-t/2, t/2}}},
			{RodZ, {{-r_0 +a -d, r_0 +a -d}, {-r_0, r_0}, {-t/2, t/2}}},
			{RodZ, {{-r_0, r_0}, {-r_0 +a -d, r_0 +a -d}, {-t/2, t/2}}},
			{RodZ, {{-r_0, r_0}, {-r_0 -a +d, r_0 -a +d}, {-t/2, t/2}}}
		}
	};

	putObjects(W, n(3.4777), square_lattice, n(1.444), slot_lattice, n(3.4777), defect, n(1.444), defect_slot, Air );
	pointDipole(W, Ez, 0, 0, 0, Pulse, 1550, 100);


	slice XY = createSliceXY(W, 0);
	slice XZ = createSliceXZ(W, 0);
	slice YZ = createSliceYZ(W, 0);
	sliceSnap(W, LogRI, XY, png(gray, 0), "/%%");
	sliceSnap(W, LogRI, YZ, png(gray, 0), "/%%");
	sliceSnap(W, LogRI, XZ, png(gray, 0), "/%%");

	float X_right = 0, X_left = 0, Y_right = 0, Y_left = 0, Z_right = 0 , Z_left = 0 ;
	float TE_poyntingout = 0, TE_poyntingbox = 0;
	float Plane_X = 0, Plane_Y = 0, Plane_Z = 0;
	float Slot_Area = 0, Slot_X = 0, Slot_Y = 0;

	writeTxt(W, "", "Time\tTE_poyntingout\tTE_poyntingbox\tHorizontal(slot)\n");

	for (int n=1, N=18000; timer(W->ID, n, W->N+N); n++) {
		updateH(W);
		updateE(W);

		writeRow(W, "/TotalTime", W->dt*n, get(W, Ez, 0, 0, 0));

		X_right += poyntingX(W, 10*a, - 10*a, 10*a, - h/2 - 400, h/2 + 400);
        X_left -= poyntingX(W, -10*a, - 10*a, 10*a, - h/2 - 400, h/2 + 400);
 
        Y_right += poyntingY(W, 10*a, -10*a, 10*a, - h/2 - 400, h/2 + 400);
        Y_left -= poyntingY(W, -10*a, -10*a, 10*a, - h/2 - 400, h/2 + 400);
  
        Z_right += poyntingZ(W, h/2 + 400, -10*a, 10*a, -10*a, 10*a);
        Z_left -= poyntingZ(W, -h/2 - 400, -10*a, 10*a, -10*a, 10*a);
  
        TE_poyntingout += poyntingOut (W, -10*a, 10*a, -10*a, 10*a, - h/2 - 400, h/2 + 400);
        TE_poyntingbox = X_right + X_left + Y_right + Y_left + Z_right + Z_left ;
  
        Plane_X = X_left + X_right;
        Plane_Y = Y_left + Y_right;
        Plane_Z = Z_left + Z_right;
  
        Slot_X += poyntingX(W, 10*a, -10*a, 10*a, -400, 400) - poyntingX(W, -10*a, -10*a, 10*a, -400, 400) ;
        Slot_Y += poyntingY(W, 10*a, -10*a, 10*a, -400, 400) - poyntingY(W, -10*a, -10*a, 10*a, -400, 400) ;


		Slot_Area = Slot_X + Slot_Y ;

		if (n > W->N) {
			writeRow(W, "/Time-Ez", W->dt*n, get(W, Ez, 0, 0, 0));
			writeSpectrum(W, N, 500, 2000, "/Spectrum-Ez", get(W, Ez, 0, 0, 0));
			
		}

		writeRow(W, "", W->dt*n/100, TE_poyntingout, TE_poyntingbox, Slot_Area);

		if (n+2*W->T > W->N+N) {
			sliceSnap(W, EE, XY, 15, png(hot, -1), "/%%/");
			sliceSnap(W, EE, XY, 15, txt, "/%%/");
			sliceSnap(W, EE, XZ, 15, png(hot, -1), "/%%/");
			sliceSnap(W, EE, XZ, 15, txt, "/%%/");	
			sliceSnap(W, Ez, XY, 15, png(dkbr, -1), "/%%/");
			sliceSnap(W, Ez, XY, 15, txt, "/%%/");
			sliceSnap(W, Ez, XZ, 15, png(dkbr, -1), "/%%/");
			sliceSnap(W, Ez, XZ, 15, txt, "/%%/");
		}
	}
}
