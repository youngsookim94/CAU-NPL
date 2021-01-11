#include <alis.h>
#include <math.h>

int main(int argc, char **argv) {
    float lambda = atof(argv[1]);
    float nSiO2 = 1.45;
    float nTiO2 = 2.48;
    float dTiO2 = 110;
    float dSiO2 = 190;
    float dUnit = dTiO2 + dSiO2;
    float numLayer = atof(argv[2]);
    float lengthY = dUnit*numLayer;

    res RES = {5, {20, 20, 10}, {1240}};
    dom DTF = {{INF}, {INF}, {lengthY/2 + 200}};
    dom DSF = {{1000}, {0}, {lengthY/2 + 1500}};
    sur SUR = {{SYM, PML}, {PML}, {PML}, {24}};
    world WLD = createWorld(DSF, RES, SUR, "%s_n%.0f_l%.0f", argv[0], numLayer, lambda);

    object boxSiO2 = {Box, {{-INF, INF}, {-INF, INF}, {-150, 40}}};
    object boxTiO2 = {Box, {{-INF, INF}, {-INF, INF}, {40, 150}}};

    object latticeTiO2 = {Lattice, {{0, 1}, {0, 1}, {dUnit, numLayer}}, objects{boxTiO2}};
    object latticeSiO2 = {Lattice, {{0, 1}, {0, 1}, {dUnit, numLayer}}, objects{boxSiO2}};

    putObjects (WLD, n(nTiO2), latticeTiO2, n(nSiO2), latticeSiO2, Air);
    planewave(WLD, DTF, Ppol, 0, 0, Sine, lambda, 10);

    phaser P = createPhaser(WLD, lambda);
    float mesR = 0, mesT = 0, mesA = 0, mesTotal = 0;
    
    slice XY = createSliceXY(WLD, 0);
    slice XZ = createSliceXZ(WLD, 0);
    slice YZ = createSliceYZ(WLD, 0);

    sliceSnap(WLD, LogRI, XY, png(gray, 1), "/%%");
    sliceSnap(WLD, LogRI, XZ, png(gray, 1), "/%%");
    sliceSnap(WLD, LogRI, YZ, png(gray, 1), "/%%");

    writeTxt(WLD, "/RTA", "Wavelength\tReflection\tTransmission\tAbsorption\tTotal\r\n");

    for (int n = 1, N = 50*lambda/WLD -> dt; timer(n, N); n++) {
        updateH(WLD);
        mesR += fabs(poyntingZ(WLD, lengthY/2 + 150));
        mesT += fabs(poyntingZ(WLD, lengthY/2 - 250));
        mesA += objectAbsorption(WLD, latticeSiO2) + objectAbsorption(WLD, latticeTiO2);
        mesTotal = mesR + mesT + mesA;

        updateE(WLD);
        mesR += fabs(poyntingZ(WLD, dUnit*numLayer/2 + 150));
        mesT += fabs(poyntingZ(WLD, dUnit*numLayer/2 - 250));
        mesA += objectAbsorption(WLD, latticeSiO2) + objectAbsorption(WLD, latticeTiO2);
        mesTotal = mesR + mesT + mesA;

        if (N-n < WLD -> T) updatePhaser(WLD, P);
        if (!(n%( WLD->T ))) {
            writeRow(WLD, "/RTA", WLD->dt*n/50, mesR, mesT, mesA, mesTotal);
            mesR = 0;
            mesT = 0;
            mesA = 0;
        }
        
        if (N-n <= 2*WLD -> T) {
            sliceSnap(WLD, Ex, XZ, 20, png(dkbr, -1), "/%%/");
        }
    }
}