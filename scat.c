#include <alis.h>

int main(int argc, char **argv) {
    float rNP = atof(argv[1]);

    dom DTF = {{rNP + 20}, {rNP + 20}, {rNP + 20}};
    dom DSF = {{rNP + 200}, {rNP + 200}, {rNP + 200}};
    res Res = {5, {10}, {1240}};
    sur Sur = {{SYM< PML}, {SYM, PML}, {PML, PML}};
    world W = createWorld(DSF, Res, Sur, "%s_r%.0f", argv[0], rNP);

    matter Si_Palik = {{4.4544},{{0.5359,5.6633,6.8761,2.6536},{3.6426,0.1318,7.4933,0.499}}};


    object SiNP = {Ball, {{-rNP, rNP}, {-rNP, rNP}, {-rNP, rNP}}};

    putObjects(W, Si_Palik, SiNP, Air);

    planewave(W, DTF, Spol, 0, 0, Band, 400, 2000);
    slice XZ = createSliceXZ(W, 0);
    phaser A = createPhasers(W, 400, 2000, 10, -rNP -10, rNP + 10, -rNP - 10, rNP + 10, -rNP - 10, rNP + 10);
    phaser S = createPhasers(W, 400, 2000, 10, -rNP - 30, rNP + 30, -rNP - 30, rNP + 30, -rNP - 30, rNP + 30);
    
    for (int n = 1; timer(n, S->N); n++) {
        updateE(W);
        updateH(W);
        updatePhaser(W, A);
        updatePhaser(W, S);
        sliceFreqDom(W, ScEy, XZ, S->N, 610, png(dkbr, 0), "/%%-610nm/");
        sliceFreqDom(W, ScEy, XZ, S->N, 770, png(dkbr, 0), "/%%_770nm/");
    }
   
    writePoyntingSpectrum(W, "", A, S);
}

