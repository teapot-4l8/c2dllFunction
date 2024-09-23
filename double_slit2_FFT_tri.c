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
#define TRIANGLE_HEIGHT 50  // Height of the triangle

double** create_2d_array(int rows, int cols);

// Function to calculate and apply the triangular slit mask
void calculate_triangular_slit(double** slit , double complex E[np][np]);


int main() {
    double k = 2 * M_PI / wL;
    double dp = (pmax - pmin) / (np - 1);

    // Create 2D arrays for xp, yp, and zp
    double **xp = create_2d_array(np, np);
    double **yp = create_2d_array(np, np);
    double **zp = create_2d_array(np, np);
    double **slit = create_2d_array(np, np);
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

            // term1: w0 ./ w1  Calculate the beam waist scaling
            double complex term1 = w0 / w1[r][c];
            // term2: exp(-rho1.^2 ./ (w1).^2)  Gaussian term for beam profile
            double complex term2 = exp(-rho1[r][c] * rho1[r][c] / (w1[r][c] * w1[r][c]));
            // term3: exp(-1i * 2 * pi / wL * z1 * rho1.^2 / (2 * (z1.^2 + zR.^2))) Phase shift due to propagation and curvature
            double complex term3 = cexp(-I * 2 * M_PI / wL * z1[r][c] * (rho1[r][c] * rho1[r][c]) / (2 * (z1[r][c] * z1[r][c] + zR * zR)));
            // term4: exp(1i * atan(z1/zR))
            double complex term4 = cexp(I * atan(z1[r][c] / zR));
            // term5: (rho1 * sqrt(2) ./ w1).^m1
            double complex term5 = cpow(rho1[r][c] * sqrt(2) / w1[r][c], m1);
            // term6: exp(1i * m1 * phi1)
            double complex term6 = cexp(I * m1 * phi1[r][c]);
            // term7: exp(1i * m1 * atan(z1/zR))
            double complex term7 = cexp(I * m1 * atan(z1[r][c] / zR));
            // term8: exp(1i * 2 * pi / wL * xp * sin(thetaX1))
            double complex term8 = cexp(I * 2 * M_PI / wL * xp[r][c] * sin(thetaX1));
            // term9: exp(-1i * k * zp)
            double complex term9 = cexp(-I * k * zp[r][c]);
            // Combine all terms to calculate E
            E[r][c] = term1 * term2 * term3 * term4 * term5 * term6 * term7 * term8 * term9;
        }
    }

    calculate_triangular_slit(slit, E);




    free(xp);
    free(yp);
    free(zp);
    free(slit);
    return 0;
}


// Helper function to initialize a 2D array dynamically
double** create_2d_array(int rows, int cols) {
    double** array = (double**)malloc(rows * sizeof(double*));
    for (int i = 0; i < rows; i++) {
        array[i] = (double*)malloc(cols * sizeof(double));
    }
    return array;
}

// Function to calculate and apply the triangular slit mask
void calculate_triangular_slit(double** slit , double complex E[np][np]) {
    double center_x = np / 2.0;
    double center_y = np / 2.0;

    double vertex1[2] = {center_x, center_y - 2.0 / 3.0 * TRIANGLE_HEIGHT};
    double vertex2[2] = {center_x + TRIANGLE_HEIGHT / sqrt(3), center_y + 1.0 / 3.0 * TRIANGLE_HEIGHT};
    double vertex3[2] = {center_x - TRIANGLE_HEIGHT / sqrt(3), center_y + 1.0 / 3.0 * TRIANGLE_HEIGHT};

    // Initialize slit array with zeros
    for (int i = 0; i < np; i++) {
        for (int j = 0; j < np; j++) {
            slit[i][j] = 0;
        }
    }

    // Iterate over the grid points
    for (int nn1 = 0; nn1 < np; nn1++) {
        for (int nn2 = 0; nn2 < np; nn2++) {
            // Compute vectors v0, v1, and v2
            double v0[2] = {vertex3[0] - vertex1[0], vertex3[1] - vertex1[1]};
            double v1[2] = {vertex2[0] - vertex1[0], vertex2[1] - vertex1[1]};
            double v2[2] = {nn1 - vertex1[0] + 1, nn2 - vertex1[1] + 1};

            // Dot products
            double dot00 = v0[0] * v0[0] + v0[1] * v0[1];
            double dot01 = v0[0] * v1[0] + v0[1] * v1[1];
            double dot02 = v0[0] * v2[0] + v0[1] * v2[1];
            double dot11 = v1[0] * v1[0] + v1[1] * v1[1];
            double dot12 = v1[0] * v2[0] + v1[1] * v2[1];


            // Compute barycentric coordinates u and v
            double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
            double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
            double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

//             Check if point is inside the triangle
            if (u >= 0 && v >= 0 && (u + v) < 1) {
                slit[nn1][nn2] = 1;  // Mark as inside
            }
        }
    }

//    // Apply the slit mask to the electric field (element-wise multiplication)
//    for (int i = 0; i < np; i++) {
//        for (int j = 0; j < np; j++) {
//            E[i][j] *= slit[i][j];
//        }
//    }
}