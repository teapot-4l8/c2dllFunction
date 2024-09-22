#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>  // For handling complex numbers

// Constants
#define wL 550e-9  // wavelength
#define np 201     // grid points
#define S1x 0e-3   // Source point x (non-zero for tilted incidence)
#define zSP 3e-3   // Axial distance from Source 1 to Detector Screen
#define pmax 80e-6 // maximum detector screen coordinate
#define pmin -pmax // minimum detector screen coordinate
#define w0 20e-6   // beam waist
#define m1 1       // vortex beam indicator (0 for Gaussian, non-zero for vortex)

// Helper function to initialize a 2D array dynamically
double** create_2d_array(int rows, int cols) {
    double** array = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        array[i] = (double*)malloc(cols * sizeof(double));
    }
    return array;
}

int main() {
    double k = 2 * M_PI / wL;
    double dp = (pmax - pmin) / (np - 1);

    // Create 2D arrays for xp, yp, and zp
    double **xp = create_2d_array(np, np);
    double **yp = create_2d_array(np, np);
    double **zp = create_2d_array(np, np);

    // Initialize xp, yp, zp
    for (int c = 0; c < np; c++) {
        double p = pmin + dp * c;
        for (int r = 0; r < np; r++) {
            xp[r][c] = p;
            yp[c][r] = -p;
            zp[r][c] = zSP;
        }
    }

    // Create z1, thetaX1, w1, phi1, rho1 arrays
    double z1[np][np];
    double w1[np][np];
    double phi1[np][np];
    double rho1[np][np];
    double complex E[np][np]; // To store the result of E (complex)

    double thetaX1 = atan(S1x / (2 * zSP));
    double zR = M_PI * w0 * w0 / wL;

    // Compute z1, w1, phi1, rho1, and E values
    for (int r = 0; r < np; r++) {
        for (int c = 0; c < np; c++) {
            z1[r][c] = sqrt(zp[r][c] * zp[r][c] + (xp[r][c] + S1x) * (xp[r][c] + S1x));
            w1[r][c] = w0 * sqrt(1 + (z1[r][c] / zR) * (z1[r][c] / zR));
            phi1[r][c] = atan2(yp[r][c], xp[r][c]);
            rho1[r][c] = sqrt((xp[r][c] * cos(thetaX1)) * (xp[r][c] * cos(thetaX1)) + yp[r][c] * yp[r][c]);

            // Now calculate E using the formula provided
            double complex term1 = w0 / w1[r][c];
            double complex term2 = exp(-rho1[r][c] * rho1[r][c] / (w1[r][c] * w1[r][c]));
            double complex term3 = exp(-I * 2 * M_PI / wL * z1[r][c] * rho1[r][c] * rho1[r][c] / 2 / (z1[r][c] * z1[r][c] + zR * zR));
            double complex term4 = exp(I * atan(z1[r][c] / zR));
            double complex term5 = cpow(rho1[r][c] * sqrt(2) / w1[r][c], m1);
            double complex term6 = exp(I * m1 * phi1[r][c]);
            double complex term7 = exp(I * m1 * atan(z1[r][c] / zR));
            double complex term8 = exp(I * 2 * M_PI / wL * xp[r][c] * sin(thetaX1));
            double complex term9 = exp(-I * k * zp[r][c]);

            // Combine all terms to calculate E
            E[r][c] = term1 * term2 * term3 * term4 * term5 * term6 * term7 * term8 * term9;


        }
        // Print the result (real and imaginary parts)
        printf("E[%d][%d] = %e + %ei\n", r, 0, creal(E[r][0]), cimag(E[r][0]));
//        printf("E[%d][%d] = %e + %ei\n", r, 0, creal(E[r][0]), cimag(E[r][0]));
    }

    // Remember to free the allocated memory when done
    for (int i = 0; i < np; i++) {
        free(xp[i]);
        free(yp[i]);
        free(zp[i]);
    }
    free(xp);
    free(yp);
    free(zp);

    return 0;
}
