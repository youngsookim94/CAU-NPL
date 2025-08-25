#include <alis.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>

int readParticlePositions2D(const char* filename, float positions[][2], int maxParticles) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Cannot open file: %s\n", filename);
        return -1;
    }
    
    char line[256];
    int count = 0;
    
    while (fgets(line, sizeof(line), file) && count < maxParticles) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') {
            continue;
        }
        
        float x, z;
        if (sscanf(line, "%f %f", &x, &z) == 2) {
            positions[count][0] = x;
            positions[count][1] = z;
            count++;
        }
    }
    
    fclose(file);
    printf("\nRead %d of 2D particle positions from file.\n", count);
    return count;
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("Usage: %s <numParticles> [positionFile]\n", argv[0]);
        return 1;
    }
    
    float wL = atof(argv[1]);
    float dwL = 0.01;
    int numParticles = atoi(argv[2]);
    float particleRadius = 0.2;
    const char* positionFile = (argc >= 3) ? argv[3] : "positions_2d.txt";
    float Zpos = atof(argv[4]);
    float Xpos = atof(argv[5]);

    dom DSF = {{25}, {0}, {40}}; 
    res Res = {0.025, {0.05}, {1.240}};
    sur Sur = {{PML, PML}, {PBC, PBC}, {PML, PML}, {48}};
    dom DTF = {{-INF, INF}, {-INF, INF}, {-35, 35}};

    world W = createWorld(DSF, Res, Sur, "%s_l%.3f_n%d_Z%.0f_X%.0f", argv[0], wL, numParticles, Zpos, Xpos);

    object *particles = malloc(numParticles * sizeof(object));
    float (*positions)[2] = malloc(numParticles * sizeof(float[2]));

    object movingParticle = {RodY, {{Xpos - particleRadius, Xpos + particleRadius}, {-INF, INF}, {Zpos - particleRadius, Zpos + particleRadius}}};

    int actualParticles = readParticlePositions2D(positionFile, positions, numParticles);
    
    if (actualParticles <= 0) {
        printf("Failed to read particle positions. Using default random generation.\n");
        srand(time(NULL));
        
        float domainSize_x = 50/2 - particleRadius;
        float limitZ1 = -30;
        float limitZ2 = 30;
        float zRange = limitZ2 - limitZ1;

        for (int i = 0; i < numParticles; i++) {
            positions[i][0] = (float)rand() / RAND_MAX * 2 * domainSize_x - domainSize_x; // x 좌표
            positions[i][1] = limitZ1 + (float)rand() / RAND_MAX * zRange; // z 좌표
        }
        actualParticles = numParticles;
    }
     if (actualParticles < numParticles) {
        numParticles = actualParticles;
    }

    for (int i = 0; i < numParticles; i++) {
        float x = positions[i][0];
        float z = positions[i][1];
        particles[i] = (object){RodY, {{x - particleRadius, x + particleRadius}, {-INF, INF}, {z - particleRadius, z + particleRadius}}};


    }
    object randParticles[numParticles];
    for (int i = 0; i < numParticles; i++) {
        if (i < 10) {
            float x = positions[i][0];
            float z = positions[i][1];
        }
    }
    object particleGroup = {Combination, {numParticles}, particles};
    
    matter Si_Palik = {{4.4544},{{0.5359,5.6633,6.8761,2.6536},{3.6426,0.1318,7.4933,0.499}}};
    putObjects(W, Si_Palik, particleGroup, Si_Palik, movingParticle, Air);

    planewave(W, DTF, Spol, 0, 0, Sine, wL, dwL); 
    //phaser P = createPhaser(W, wL);

    slice XZ_domain = createSliceXZ(W, 0); 


    slice Line_measure1 = createSliceX(W, 0, -40);
    slice Line_measure2 = createSliceX(W, 0, 40);

    sliceSnap(W, LogRI, XZ_domain, png(gray, 1), "/%%");


    for (int n = 0, N = 2000 * W->dx / W->dt; timer(n, N); n++) {
        updateH(W);
        updateE(W);

        sliceSnap(W, ScEE, XZ_domain, 20, png(jet, -1), "/XZ_Field_Distribution/");
        //sliceSnap(W, ScEE, XZ_domain, 20, h5, "/XZ_Field_Distribution/");
        sliceSnap(W, Sz, XZ_domain, 20, png(gray, -1), "/XZ_Sz_Distribution/");

        sliceSnap(W, ScEE, Line_measure1, 20, png(jet, -1), "/Speckle_Pattern_Line_F/");
        //sliceSnap(W, ScEE, Line_measure1, 20, h5, "/Speckle_Pattern_Line_F/");
        sliceSnap(W, ScEE, Line_measure2, 20, png(jet, -1), "/Speckle_Pattern_Line_B/");
        //sliceSnap(W, ScEE, Line_measure2, 20, h5, "/Speckle_Pattern_Line_B/");
        sliceSnap(W, Sz, Line_measure1, 20, png(gray, -1), "/Speckle_Sz_Line_F/");
        //sliceSnap(W, Sz, Line_measure1, 20, h5, "/Speckle_Sz_Line_F/");
        sliceSnap(W, Sz, Line_measure2, 20, png(gray, -1), "/Speckle_Sz_Line_B/");
        //sliceSnap(W, Sz, Line_measure2, 20, h5, "/Speckle_Sz_Line_B/");
     }
    

    //farFieldTheta(W, P, 90, "");
    //farFieldProfile(W, P, AzimuthalMap, 360, 180, png(hot, 0), "/FF/");

    free(particles);
    free(positions);
    return 0;
}
