#include <alis.h>

int main(int argc, char **argv) {
    float lambda = atof(argv[1]);
    float nSiO2 = 1.45;
    float nTiO2 = 2.48;
    float dTiO2 = 110;
    float dSiO2 = 190;
    float dUnit = dTiO2 + dSiO2;
    float numLayer = atof(argv[2]);
    float lengthY = dUnit*numLayer;

    res RES = {5, {10}, {1240}};
    dom DTF = {{INF}, {lengthY/2 + 200}, {INF}};
    dom DSF = {{1000}, {lengthY/2 + 1000}, {0}};
    sur SUR = {{SYM, PML}, {PML}, {PBC}};
    world WLD = createWorld(DSF, RES, SUR, "%s_n%.0f_l%.0f", argv[0], numLayer, lambda);

    object boxSiO2 = {Box, {{-INF, INF}, {-150, 40}, {-INF, INF}}};
    object boxTiO2 = {Box, {{-INF, INF}, {40, 150}, {-INF, INF}}};

    object latticeTiO2 = {Lattice, {{0, 1}, {dUnit, numLayer}, {0, 1}}, objects{boxTiO2}};
    object latticeSiO2 = {Lattice, {{0, 1}, {dUnit, numLayer}, {0, 1}}, objects{boxSiO2}};

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

    for (int n =1, N=100*lambda/WLD -> dt; timer(n, N); n++) {
        updateH(WLD);
        mesR += -poyntingZ(WLD, dUnit*numLayer/2 + 250);
        mesT += -poyntingZ(WLD, dUnit*numLayer/2 - 150);
        mesA += objectAbsorption(WLD, latticeSiO2) + objectAbsorption(WLD, latticeTiO2);
        mesTotal = mesR + mesT + mesA;

        updateE(WLD);
        mesR += -poyntingZ(WLD, dUnit*numLayer/2 + 250);
        mesT += -poyntingZ(WLD, dUnit*numLayer/2 - 150);
        mesA += objectAbsorption(WLD, latticeSiO2) + objectAbsorption(WLD, latticeTiO2);
        mesTotal = mesR + mesT + mesA;

        if (N-n < WLD -> T) updatePhaser(WLD, P);
        if (!(n%(WLD->T))) {
            writeRow(WLD, "/RTA", WLD->dt*n/100, mesR/mesTotal, mesT/mesTotal, mesA/mesTotal, mesTotal);
            mesR = 0;
            mesT = 0;
            mesA = 0;
        }
        
        if (N-n <= 2*WLD -> T) {
            sliceSnap(WLD, Ex, XY, 20, png(dkbr, -1), "/%%/");
        }
    }
}