/********************************************************************/
/* Meta표면 계산을 위한 위상 측정 관련 코드 - 작성자: 김영수, 2025.05.09.  */ 
/* 이 코드는 ALIS 라이브러리를 사용하여 메타표면의 반사 및 투과 위상을 측정함 */
/* if문 내부에서 계산하는 결과는 직접필드측정방식으로 계산된 것임            */
/* Far-Field spectrum을 통해 계산한 것은 phaser 객체를 통해 계산된 것임   */
/* 두 방식은 서로 다른 위상기준점을 가짐으로, 결과값이 다름                 */
/* Far-field 위상의 경우 공간상에서 지연된 누적 phase 값을 출력함          */
/* 내부 모드 분석이나 필드의 동작을 분석할때는 직접필드방식이 더 적합함       */
/********************************************************************/


#include <alis.h>

int main(int argc, char **argv) {
    float H = 250; //atof(argv[1]);
    float Wid = 250; //atof(argv[2]);
    float tS = 10; //atof(argv[3]);
    float wL = 1550; //atof(argv[4]);
    float sDom = Wid + tS;

    res Res = {0.5, {5, 10, 10}, {1240}};
    dom Dom = {{sDom/2}, {0}, {1500}};
    sur Sur = {{PBC, PBC}, {PBC, PBC}, {PML, PML}};
    world W = createWorld(Dom, Res, Sur, "%s_H%.0f_W%.0f_t%.0f_L%.0f", argv[0], H, Wid, tS, wL);
    dom TF = {{-INF, INF}, {-INF, INF}, {-INF, H/2 + 100}};

    object Si_Block = {Lattice, {{Wid + tS, 2}, {0, 1}, {0, 1}}, objects {Box, {{-Wid/2, Wid/2},{-INF, INF}, {-H/2, H/2}}}};

    putObjects(W, n(3.47), Si_Block, Air);
    // 복소수 필드를 위해 복소 소스 활성화
    W->complexField = 1;
    planewave(W, TF, Ppol, 0, 0, Sine, wL, 10);

    phaser P = createPhaser(W, wL);

    float SiAbs = 0, Ref = 0, Tsm = 0, sOut = 0;
    float complex Ex_ref, Ex_trans; // 복소수 필드 값 저장용

    slice XZ = createSliceXZ(W, 0);
    slice XY = createSliceXY(W, 0);
    slice YZ = createSliceYZ(W, 0);
    sliceSnap(W, LogRI, XZ, png(gray, 1), "/%%");
    sliceSnap(W, LogRI, XY, png(gray, 1), "/%%");
    sliceSnap(W, LogRI, YZ, png(gray, 1), "/%%");

    writeTxt(W, "/RTA", "Wavelength\tTotal\tAbs\tRef\tTsm\r\n");
    writeTxt(W, "/Phase", "Wavelength\tRef_Phase\tTrans_Phase\r\n");

    for (int n = 0, N=100*wL/W->dt; timer(n, N, W->ID); n++) {
        updateH(W);
        Ref += abs(poyntingZ(W, H/2 + 150));
        Tsm += abs(poyntingZ(W, -H/2 -150));
        SiAbs += objectAbsorption(W, Si_Block);
        sOut = Ref + Tsm + SiAbs;

        updateE(W);
        Ref += abs(poyntingZ(W, H/2 + 150));
        Tsm += abs(poyntingZ(W, -H/2 -150));
        SiAbs += objectAbsorption(W, Si_Block);
        sOut = Ref + Tsm + SiAbs;

        if (N-n < W->T) {
            updatePhaser(W, P);
            
            // 반사파와 투과파의 위상 측정
            if (W->complexField) {
                // 복소수 필드 값 가져오기
                Ex_ref = get(W, Ex, 0, 0, H/2 + 150) + I * get(W, iEx, 0, 0, H/2 + 150);
                Ex_trans = get(W, Ex, 0, 0, -H/2 - 150) + I * get(W, iEx, 0, 0, -H/2 - 150);
                
                // 위상 계산 (라디안을 도로 변환)
                float phase_ref = cargf(Ex_ref) * 180.0 / PI;
                float phase_trans = cargf(Ex_trans) * 180.0 / PI;
                
                // 결과 기록
                writeRow(W, "/Phase", wL, phase_ref, phase_trans);
            }
        }
        
        if (!(n%(W->T))) {
            writeRow(W, "/RTA", W->dt*n/100, sOut/W->T, (SiAbs/sOut), (Ref/sOut), (Tsm/sOut));
            Ref = 0;
            Tsm = 0;
            SiAbs = 0;
        }
        
        if (N-n <= 2*W->T) {
            sliceSnap(W, Ex, XZ, 20, png(dkbr, -1), "/XZ-Ex/");
            sliceSnap(W, ScEx, XZ, 20, png(dkbr, -1), "/XZ_ScEx/");
            sliceSnap(W, EE, XZ, 20, png(hot, -1), "/XZ_EE/");
            sliceSnap(W, ScEE, XZ, 20, png(hot, -1), "/XZ_ScEE/");
        }
    }

    updatePhaser(W, P);
    // 원거리장 분석 - 위상 정보 포함
    farFieldTheta(W, P, 90, "");
    farFieldProfile(W, P, AzimuthalMap, 360, 180, png(hot,0), "/");

    // 추가: 반사파와 투과파 방향의 원거리장 스펙트럼 (위상 정보 포함)
    writeFarFieldSpectrum(W, Ppol, 180, 0, "/Reflection_Spectrum", P);
    writeFarFieldSpectrum(W, Ppol, 0, 0, "/Transmission_Spectrum", P);
}