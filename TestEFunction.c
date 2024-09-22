//#include <stdio.h>
//#include <math.h>
//#include <complex.h>
//
//#define M_PI 3.14159265358979323846
//// Define constants (if needed)
//#define NP 201  // Number of grid points
//#define WL 550e-9  // Wavelength (in meters)
//#define ZR 20e-6  // Rayleigh range
//#define m1 1  // Vortex charge
//
//// Define complex number type (if not already)
//typedef double complex cdouble;
//
//void calculate_electric_field(double *w0, double *w1, double *rho1, double *z1, double *xp, double *zp, double *phi1,
//                              double wL, double zR, double thetaX1, int np, cdouble *E);
//
//int main() {
//    // Define variables
//    double wL = WL;
//    double zR = ZR;
//    double w0[NP] = {20e-6};  // Beam waist
//    double xp[NP*NP], yp[NP*NP], zp[NP*NP];  // Coordinates
//    double w1[NP*NP], rho1[NP*NP], z1[NP*NP], phi1[NP*NP];  // Other variables
//    cdouble E[NP*NP];  // To store the electric field
//    double thetaX1 = 0.0;  // Incident angle
//
//    // Initialize variables (for simplicity, initializing all entries with some values)
//    for (int i = 0; i < NP; i++) {
//        for (int j = 0; j < NP; j++) {
//            int idx = i * NP + j;
//            xp[idx] = (i - NP / 2) * 80e-6 / NP;  // X coordinates
//            yp[idx] = (j - NP / 2) * 80e-6 / NP;  // Y coordinates
//            zp[idx] = 3e-3;  // Z distance to screen
//            z1[idx] = sqrt(zp[idx] * zp[idx] + xp[idx] * xp[idx]);  // z1 calculation
//            w1[idx] = w0[0] * sqrt(1 + pow(z1[idx] / zR, 2));  // w1 calculation
//            rho1[idx] = sqrt(pow(xp[idx], 2) + pow(yp[idx], 2));  // rho1 calculation
//            phi1[idx] = atan2(yp[idx], xp[idx]);  // phi1 calculation
//        }
//    }
//
//    // Call the function to calculate the electric field
//    calculate_electric_field(w0, w1, rho1, z1, xp, zp, phi1, wL, zR, thetaX1, NP, E);
//
//    // Print the electric field values (for testing)
//    for (int i = 0; i < NP; i++) {
//        for (int j = 0; j < NP; j++) {
//            int idx = i * NP + j;
//            printf("E[%d, %d] = %.5f + %.5fi\n", i, j, creal(E[idx]), cimag(E[idx]));
//        }
//    }
//
//    return 0;
//}
//
//
//void calculate_electric_field(double *w0, double *w1, double *rho1, double *z1, double *xp, double *zp, double *phi1,
//                              double wL, double zR, double thetaX1, int np, cdouble *E) {
//
//    int i, j;
//    double k = 2 * M_PI / wL;
//
//    // Step 1: Calculate the beam waist scaling
//    double beam_waist_scaling[np][np];
//    for (i = 0; i < np; i++) {
//        for (j = 0; j < np; j++) {
//            beam_waist_scaling[i][j] = w0[i] / w1[i * np + j];
//        }
//    }
//
//    // Step 2: Gaussian term for beam profile
//    cdouble gaussian_term[np][np];
//    for (i = 0; i < np; i++) {
//        for (j = 0; j < np; j++) {
//            gaussian_term[i][j] = cexp(-pow(rho1[i * np + j], 2) / pow(w1[i * np + j], 2));
//        }
//    }
//
//    // Step 3: Phase shift due to propagation and curvature
//    cdouble curvature_phase_shift[np][np];
//    for (i = 0; i < np; i++) {
//        for (j = 0; j < np; j++) {
//            curvature_phase_shift[i][j] = cexp(-I * 2 * M_PI / wL * z1[i * np + j] * pow(rho1[i * np + j], 2) /
//                                               (2 * (pow(z1[i * np + j], 2) + pow(zR, 2))));
//        }
//    }
//
//    // Step 4: Gouy phase shift
//    cdouble gouy_phase_shift[np][np];
//    for (i = 0; i < np; i++) {
//        for (j = 0; j < np; j++) {
//            gouy_phase_shift[i][j] = cexp(I * atan(z1[i * np + j] / zR));
//        }
//    }
//
//    // Step 5: Radial dependence term for vortex beam
//    cdouble radial_term[np][np];
//    for (i = 0; i < np; i++) {
//        for (j = 0; j < np; j++) {
//            radial_term[i][j] = cpow(rho1[i * np + j] * sqrt(2) / w1[i * np + j], m1);
//        }
//    }
//
//    // Step 6: Azimuthal phase shift due to vortex charge
//    cdouble azimuthal_phase_shift[np][np];
//    for (i = 0; i < np; i++) {
//        for (j = 0; j < np; j++) {
//            azimuthal_phase_shift[i][j] = cexp(I * m1 * phi1[i * np + j]);
//        }
//    }
//
//    // Step 7: Additional Gouy phase shift for vortex charge
//    cdouble vortex_gouy_phase_shift[np][np];
//    for (i = 0; i < np; i++) {
//        for (j = 0; j < np; j++) {
//            vortex_gouy_phase_shift[i][j] = cexp(I * m1 * atan(z1[i * np + j] / zR));
//        }
//    }
//
//    // Step 8: Plane wave tilt phase shift
//    cdouble tilt_phase_shift[np][np];
//    for (i = 0; i < np; i++) {
//        for (j = 0; j < np; j++) {
//            tilt_phase_shift[i][j] = cexp(I * 2 * M_PI / wL * xp[i * np + j] * sin(thetaX1));
//        }
//    }
//
//    // Step 9: Phase shift due to propagation along z-axis
//    cdouble propagation_phase_shift[np][np];
//    for (i = 0; i < np; i++) {
//        for (j = 0; j < np; j++) {
//            propagation_phase_shift[i][j] = cexp(-I * k * zp[i * np + j]);
//        }
//    }
//
//    // Step 10: Combine all the terms to calculate the electric field E
//    for (i = 0; i < np; i++) {
//        for (j = 0; j < np; j++) {
//            E[i * np + j] = beam_waist_scaling[i][j] * gaussian_term[i][j] * curvature_phase_shift[i][j] * gouy_phase_shift[i][j] *
//                            radial_term[i][j] * azimuthal_phase_shift[i][j] * vortex_gouy_phase_shift[i][j] *
//                            tilt_phase_shift[i][j] * propagation_phase_shift[i][j];
//        }
//    }
//}
//
//
