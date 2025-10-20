#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

// number of grids
const int x = 129; 
const int y = 129; 

const double lx = 1.0; // Cavity length
const double ly = 1.0; // Cavity height

const double dx = lx / (x - 1);
const double dy = ly / (y - 1);
const double max_tol = 1e-8;            // Convergence tolerance

// --- Fluid Properties ---
const double rho = 10;    // Density
const double mu = 0.01;  // Dynamic viscosity
const double nu = mu / rho; // Kinematic viscosity

// --- Boundary Conditions ---
const double u_lid = 1.0; // Top wall velocity

// --- Relaxation Factors ---
const double alpha_p = 0.3;
const double alpha_u = 0.5;
const double alpha_v = 0.5;
const double alpha_k = 0.5;
const double alpha_omega = 0.5;

// solver constants
const double beta_star = 0.09;

const double beta_1 = 0.075;
const double sigma_k_1 = 0.85;
const double sigma_omega_1 = 0.5;
const double kappa = 0.41;

const double beta_2 = 0.0828;
const double sigma_k_2 = 1.0;
const double sigma_omega_2 = 0.856;

const double gamma_1 = (beta_1 / beta_star) - (sigma_omega_1 * kappa * kappa / sqrt(beta_star));
const double gamma_2 = (beta_2 / beta_star) - (sigma_omega_2 * kappa * kappa / sqrt(beta_star));
const double a_1_const = 0.31;

// Function to get the normal distance from the closest wall
double getWallDistance(int i, int j) {
    double x_pos = (i - 0.5) * dx;
    double y_pos = (j - 0.5) * dy;
    double dist_x = min(x_pos, lx - x_pos);
    double dist_y = min(y_pos, ly - y_pos);
    return min(dist_x, dist_y);
}

// Van Leer flux limiter function for TVD scheme
double fluxLimiter(double r) {
    if (abs(1.0 + abs(r)) < 1e-9) {
        return 0.0;
    }
    return (r + abs(r)) / (1.0 + abs(r));
}

void setBoundaryConditions(vector<vector<double>>& u, vector<vector<double>>& v, vector<vector<double>>& k, vector<vector<double>>& omega) {
    // Wall omega value based on distance to the first cell center
    const double omega_wall_y = (6.0 * nu) / (beta_1 * pow(0.5 * dy, 2));
    const double omega_wall_x = (6.0 * nu) / (beta_1 * pow(0.5 * dx, 2));

    // top and bottom u velocity
    for (int i = 1; i < x; ++i) {
        u[i-1][0] = -u[i-1][1];         
        u[i-1][y] = 2.0 * u_lid - u[i-1][y - 1]; 
    }
    for (int i = 0; i <= x; ++i) {
        // v on top and bottom wall
        v[i][0] = 0;
        v[i][y-1] = 0;  
        
        // k = 0 at walls
        k[i][0] = 0;
        k[i][y] = 0;    
        k[i][1] = 0;
        k[i][y-1] = 0;

        // setting omega wall value
        omega[i][0] = omega_wall_y;
        omega[i][y] = omega_wall_y;
        omega[i][1] = omega_wall_y;
        omega[i][y-1] = omega_wall_y;
    }

    // v velocity at left and right wall
    for (int j = 1; j < y; ++j) {
        v[0][j-1] = -v[1][j-1];
        v[x][j-1] = -v[x-1][j-1];
    }
    for (int j = 0; j <= y; ++j) {
        // u-velocity at left and right wall
        u[0][j] = 0;
        u[x-1][j] = 0;

        // k = 0 at walls
        k[0][j] = 0;
        k[x][j] = 0;
        k[1][j] = 0;
        k[x-1][j] = 0;

        // setting omega wall value
        omega[0][j] = omega_wall_x;
        omega[x][j] = omega_wall_x;
        omega[1][j] = omega_wall_x;
        omega[x-1][j] = omega_wall_x;
    }
}

void calculateTurbulentViscosity(const vector<vector<double>>& u, const vector<vector<double>>& v,const vector<vector<double>>& k, const vector<vector<double>>& omega, vector<vector<double>>& tensor_dot, vector<vector<double>>& mu_t) {
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y; j++) {
            // All derivatives are calculated at the cell center (i, j)
            double dudx = (u[i][j] - u[i-1][j]) / dx;
            double dvdy = (v[i][j] - v[i][j-1]) / dy;
            
            double u_n = 0.25 * (u[i][j+1] + u[i-1][j+1] + u[i][j] + u[i-1][j]);
            double u_s = 0.25 * (u[i][j-1] + u[i-1][j-1] + u[i][j] + u[i-1][j]);
            double dudy = (u_n - u_s) / dy;

            double v_e = 0.25 * (v[i][j] + v[i+1][j] + v[i][j-1] + v[i+1][j-1]);
            double v_w = 0.25 * (v[i][j] + v[i-1][j] + v[i][j-1] + v[i-1][j-1]);
            double dvdx = (v_e - v_w) / dx;

            tensor_dot[i][j] = pow(dudx, 2) + pow(dvdy, 2) + 0.5 * pow(dudy + dvdx, 2);
            double S = sqrt(2 * tensor_dot[i][j]);

            double wall_dist = getWallDistance(i, j);
            double arg2 = max((2 * sqrt(k[i][j])) / (beta_star * omega[i][j] * wall_dist), (500 * nu) / (omega[i][j] * pow(wall_dist, 2)));
            double F2 = tanh(pow(arg2, 2));

            mu_t[i][j] = (rho * a_1_const * k[i][j]) / max((a_1_const * omega[i][j]), (S * F2));
        }
    }
}

void solveMomentum(const vector<vector<double>>& u, const vector<vector<double>>& v, vector<vector<double>>& u_star, vector<vector<double>>& v_star, vector<vector<double>>& ap_u, vector<vector<double>>& ap_v, const vector<vector<double>>& p, const vector<vector<double>>& mu_t) {
    // --- Solve u-momentum equation --- //
    for (int i = 1; i < x - 1; i++) {
        for (int j = 1; j < y; j++) {
            // Convective fluxes
            double Fe = rho * 0.5 * (u[i+1][j] + u[i][j]) * dy;
            double Fw = rho * 0.5 * (u[i-1][j] + u[i][j]) * dy;
            double Fn = rho * 0.5 * (v[i][j] + v[i+1][j]) * dx;
            double Fs = rho * 0.5 * (v[i][j-1] + v[i+1][j-1]) * dx;

            // Diffusive fluxes 
            double mu_t_e = mu_t[i+1][j];
            double mu_t_w = mu_t[i][j];
            double mu_t_n = 0.25 * (mu_t[i][j] + mu_t[i+1][j] + mu_t[i][j+1] + mu_t[i+1][j+1]);
            double mu_t_s = 0.25 * (mu_t[i][j] + mu_t[i+1][j] + mu_t[i][j-1] + mu_t[i+1][j-1]);
            
            double De = (mu + mu_t_e) * dy / dx;
            double Dw = (mu + mu_t_w) * dy / dx;
            double Dn = (mu + mu_t_n) * dx / dy;
            double Ds = (mu + mu_t_s) * dx / dy;

            // tvd scheme
            double r_e = 0.0, r_w = 0.0, r_n = 0.0, r_s = 0.0;

            // east face
            if (Fe > 0) {
                r_e = (u[i][j] - u[i-1][j]) / (u[i+1][j] - u[i][j] + 1e-9);
            } else {
                if (i < x - 2) {
                    r_e = (u[i+1][j] - u[i+2][j]) / (u[i][j] - u[i+1][j] + 1e-9);
                }
            }

            // west face
            if (Fw > 0) {
                if (i > 1) {
                    r_w = (u[i-1][j] - u[i-2][j]) / (u[i][j] - u[i-1][j] + 1e-9);
                }
            } else {
                if (i < x - 1) {
                    r_w = (u[i][j] - u[i+1][j]) / (u[i-1][j] - u[i][j] + 1e-9);
                }
            }

            // north face
            if (Fn > 0) {
                r_n = (u[i][j] - u[i][j-1]) / (u[i][j+1] - u[i][j] + 1e-9);
            } else {
                if (j < y - 1) {
                    r_n = (u[i][j+1] - u[i][j+2]) / (u[i][j] - u[i][j+1] + 1e-9);
                }
            }

            // south face
            if (Fs > 0) {
                if (j > 1) {
                    r_s = (u[i][j-1] - u[i][j-2]) / (u[i][j] - u[i][j-1] + 1e-9);
                }
            } else { 
                if (j < y) {
                    r_s = (u[i][j] - u[i][j+1]) / (u[i][j-1] - u[i][j] + 1e-9);
                }
            }

            double alpha_e = (Fe > 0) ? 1 : 0;
            double alpha_w = (Fw > 0) ? 1 : 0;
            double alpha_n = (Fn > 0) ? 1 : 0;
            double alpha_s = (Fs > 0) ? 1 : 0;

            // Deferred correction fluxes
            double F_corr_e = 0.5 * Fe * ((1 - alpha_e) * fluxLimiter(r_e) - fluxLimiter(r_e) * alpha_e) * (u[i+1][j] - u[i][j]);
            double F_corr_w = 0.5 * Fw * (fluxLimiter(r_w) * alpha_w - (1 - alpha_w) * fluxLimiter(r_w)) * (u[i][j] - u[i-1][j]);
            double F_corr_n = 0.5 * Fn * ((1 - alpha_n) * fluxLimiter(r_n) - fluxLimiter(r_n) * alpha_n) * (u[i][j+1] - u[i][j]);
            double F_corr_s = 0.5 * Fs * (fluxLimiter(r_s) * alpha_s - (1 - alpha_s) * fluxLimiter(r_s)) * (u[i][j] - u[i][j-1]);

            double S_u_DC = F_corr_w + F_corr_e + F_corr_s + F_corr_n;

            double aw = Dw + max(Fw, 0.0);
            double ae = De + max(-Fe, 0.0);
            double as = Ds + max(Fs, 0.0);
            double an = Dn + max(-Fn, 0.0);
            double ap_u_unrelaxed = ae + aw + an + as + (Fe - Fw + Fn - Fs);
            
            double pres_u = (p[i][j] - p[i + 1][j]) * dy;
            
            if (abs(ap_u_unrelaxed) > 1e-9) {
                ap_u[i][j] = ap_u_unrelaxed / alpha_u;
                double source = ae * u[i+1][j] + aw * u[i-1][j] + an * u[i][j+1] + as * u[i][j-1] + S_u_DC + pres_u + (1 - alpha_u) * ap_u[i][j] * u[i][j];
                u_star[i][j] = source / ap_u[i][j];
            } else {
                ap_u[i][j] = 0;
                u_star[i][j] = u[i][j];
            }
        }
    }

    // --- Solve v-momentum equation --- //
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y - 1; j++) {
            // Convective fluxes
            double Fe = rho * 0.5 * (u[i][j + 1] + u[i][j]) * dy;
            double Fw = rho * 0.5 * (u[i - 1][j + 1] + u[i - 1][j]) * dy;
            double Fn = rho * 0.5 * (v[i][j + 1] + v[i][j]) * dx;
            double Fs = rho * 0.5 * (v[i][j - 1] + v[i][j]) * dx;

            // Diffusive fluxes
            double mu_t_n = mu_t[i][j+1];
            double mu_t_s = mu_t[i][j];
            double mu_t_e = 0.25 * (mu_t[i][j] + mu_t[i][j+1] + mu_t[i+1][j] + mu_t[i+1][j+1]);
            double mu_t_w = 0.25 * (mu_t[i][j] + mu_t[i][j+1] + mu_t[i-1][j] + mu_t[i-1][j+1]);
            
            double De = (mu + mu_t_e) * dy / dx;
            double Dw = (mu + mu_t_w) * dy / dx;
            double Dn = (mu + mu_t_n) * dx / dy;
            double Ds = (mu + mu_t_s) * dx / dy;

            // tvd scheme
            double r_e = 0.0, r_w = 0.0, r_n = 0.0, r_s = 0.0;

            // east face
            if (Fe > 0) {
                r_e = (v[i][j] - v[i-1][j]) / (v[i+1][j] - v[i][j] + 1e-9);
            } else {
                if (i < x - 1) { 
                    r_e = (v[i+1][j] - v[i+2][j]) / (v[i][j] - v[i+1][j] + 1e-9);
                }
            }

            // west face
            if (Fw > 0) {
                if (i > 1) {
                    r_w = (v[i-1][j] - v[i-2][j]) / (v[i][j] - v[i-1][j] + 1e-9);
                }
            } else {
                if (i < x) {
                    r_w = (v[i][j] - v[i+1][j]) / (v[i-1][j] - v[i][j] + 1e-9);
                }
            }

            // north face
            if (Fn > 0) {
                r_n = (v[i][j] - v[i][j-1]) / (v[i][j+1] - v[i][j] + 1e-9);
            } else {
                if (j < y - 2) {
                    r_n = (v[i][j+1] - v[i][j+2]) / (v[i][j] - v[i][j+1] + 1e-9);
                }
            }

            // south face
            if (Fs > 0) {
                if (j > 1) {
                    r_s = (v[i][j-1] - v[i][j-2]) / (v[i][j] - v[i][j-1] + 1e-9);
                }
            } else {
                if (j < y - 1) {
                    r_s = (v[i][j] - v[i][j+1]) / (v[i][j-1] - v[i][j] + 1e-9);
                }
            }

            double alpha_e = (Fe > 0) ? 1 : 0;
            double alpha_w = (Fw > 0) ? 1 : 0;
            double alpha_n = (Fn > 0) ? 1 : 0;
            double alpha_s = (Fs > 0) ? 1 : 0;

            double F_corr_e = 0.5 * Fe * ((1 - alpha_e) * fluxLimiter(r_e) - fluxLimiter(r_e) * alpha_e) * (v[i+1][j] - v[i][j]);
            double F_corr_w = 0.5 * Fw * (fluxLimiter(r_w) * alpha_w - (1 - alpha_w) * fluxLimiter(r_w)) * (v[i][j] - v[i-1][j]);
            double F_corr_n = 0.5 * Fn * ((1 - alpha_n) * fluxLimiter(r_n) - fluxLimiter(r_n) * alpha_n) * (v[i][j+1] - v[i][j]);
            double F_corr_s = 0.5 * Fs * (fluxLimiter(r_s) * alpha_s - (1 - alpha_s) * fluxLimiter(r_s)) * (v[i][j] - v[i][j-1]);

            double S_v_DC = F_corr_w + F_corr_e + F_corr_s + F_corr_n;

            double aw = Dw + max(Fw, 0.0);
            double ae = De + max(-Fe, 0.0);
            double as = Ds + max(Fs, 0.0);
            double an = Dn + max(-Fn, 0.0);
            double ap_v_unrelaxed = ae + aw + an + as + (Fe - Fw + Fn - Fs);

            double pres_v = (p[i][j] - p[i][j + 1]) * dx;
            
            if (abs(ap_v_unrelaxed) > 1e-9) {
                 ap_v[i][j] = ap_v_unrelaxed / alpha_v;
                 double source = ae * v[i+1][j] + aw * v[i-1][j] + an * v[i][j+1] + as * v[i][j-1] + S_v_DC + pres_v + (1.0 - alpha_v) * ap_v[i][j] * v[i][j];
                 v_star[i][j] = source / ap_v[i][j];
            } else {
                ap_v[i][j] = 0.0;
                v_star[i][j] = v[i][j];
            }
        }
    }
}

void pressureCorrection(const vector<vector<double>>& u_star, const vector<vector<double>>& v_star, vector<vector<double>>& p_prime, const vector<vector<double>>& ap_u, const vector<vector<double>>& ap_v) {
    // Setting all values of p_prime to zero
    for (auto& row : p_prime) {
        fill(row.begin(), row.end(), 0.0);
    }

    double p_tol = 1e-5;
    int p_max_iter = 100;

    for (int iter = 0; iter < p_max_iter; iter++) {
        double max_p_err = 0;
        for (int i = 1; i < x; i++) {
            for (int j = 1; j < y; j++) {
                double p_prime_old = p_prime[i][j];
                
                // d coefficients for pressure correction equation
                double de = (i < x - 1 && abs(ap_u[i][j]) > 1e-9) ? dy * dy / ap_u[i][j] : 0;
                double dw = (i > 1 && abs(ap_u[i-1][j]) > 1e-9) ? dy * dy / ap_u[i-1][j] : 0;
                double dn = (j < y - 1 && abs(ap_v[i][j]) > 1e-9) ? dx * dx / ap_v[i][j] : 0;
                double ds = (j > 1 && abs(ap_v[i][j-1]) > 1e-9) ? dx * dx / ap_v[i][j-1] : 0;

                double ae = rho * de;
                double aw = rho * dw;
                double an = rho * dn;
                double as = rho * ds;
                
                double ap = ae + aw + an + as;

                double b = (rho * u_star[i-1][j] * dy - rho * u_star[i][j] * dy) + (rho * v_star[i][j-1] * dx - rho * v_star[i][j] * dx);
                
                if (abs(ap) > 1e-9) {
                    p_prime[i][j] = (ae * p_prime[i+1][j] + aw * p_prime[i-1][j] + an * p_prime[i][j+1] + as * p_prime[i][j-1] + b) / ap;
                }
                max_p_err = max(max_p_err, abs(p_prime[i][j] - p_prime_old));
            }
        }

        // Applying pressure boundary conditions (zero gradient)
        for (int i = 0; i <= x; i++) {
            p_prime[i][0] = p_prime[i][1];
            p_prime[i][y] = p_prime[i][y-1];
        }

        for (int j = 0; j <= y; j++) {
            p_prime[0][j] = p_prime[1][j];
            p_prime[x][j] = p_prime[x-1][j];
        }

        if (max_p_err < p_tol) {
            break;
        }
    }
}

void correctFields(vector<vector<double>>& u, vector<vector<double>>& v, vector<vector<double>>& p, const vector<vector<double>>& u_star, const vector<vector<double>>& v_star, const vector<vector<double>>& p_prime, const vector<vector<double>>& ap_u, const vector<vector<double>>& ap_v) {
    // u velocity
    for (int i = 1; i < x - 1; i++) {
        for (int j = 1; j < y; j++) {
            if (abs(ap_u[i][j]) > 1e-9) {
                double d = dy / ap_u[i][j];
                u[i][j] = u_star[i][j] + d * (p_prime[i][j] - p_prime[i + 1][j]);
            }
        }
    }

    // v velocity
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y - 1; j++) {
            if (abs(ap_v[i][j]) > 1e-9) {
                double d = dx / ap_v[i][j];
                v[i][j] = v_star[i][j] + d * (p_prime[i][j] - p_prime[i][j + 1]);
            }
        }
    }

    // pressure
    for (int i = 0; i <= x; i++) {
        for (int j = 0; j <= y; j++) {
            p[i][j] += alpha_p * p_prime[i][j];
        }
    }
}

void solveKOmegaSST(const vector<vector<double>>& u, const vector<vector<double>>& v, const vector<vector<double>>& tensor_dot, const vector<vector<double>>& mu_t, vector<vector<double>>& k, vector<vector<double>>& omega) {
    
    vector<vector<double>> Pk(x + 1, vector<double>(y + 1, 0));
    vector<vector<double>> F1(x + 1, vector<double>(y + 1, 0));
    
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y; j++) {
            // Gradient terms for cross-diffusion
            double dkdx = (k[i+1][j] - k[i-1][j]) / (2 * dx);
            double dkdy = (k[i][j+1] - k[i][j-1]) / (2 * dy);
            double domegadx = (omega[i+1][j] - omega[i-1][j]) / (2 * dx);
            double domegady = (omega[i][j+1] - omega[i][j-1]) / (2 * dy);
            double grad_k_dot_grad_omega = dkdx * domegadx + dkdy * domegady;

            double CD_k_omega = max(2 * rho * sigma_omega_2 * grad_k_dot_grad_omega / (omega[i][j] + 1e-9), 1e-10);

            // calculating F1 production term
            double wall_dist = getWallDistance(i, j);
            double arg1 = min(max(sqrt(k[i][j]) / (beta_star * omega[i][j] * wall_dist), (500 * nu) / (pow(wall_dist, 2) * omega[i][j])), (4 * rho * sigma_k_2 * k[i][j]) / (CD_k_omega * pow(wall_dist, 2)));
            F1[i][j] = tanh(pow(arg1, 4));

            // calculating production term
            Pk[i][j] = min(mu_t[i][j] * 2.0 * tensor_dot[i][j], 10.0 * beta_star * rho * k[i][j] * omega[i][j]);
        }
    }

    // --- Solve k equation --- //
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y; j++) {

            // blending sigma_k
            double sigma_k = F1[i][j] * sigma_k_1 + (1 - F1[i][j]) * sigma_k_2;
            
            // convective fluxes
            double Fe = rho * u[i][j] * dy;
            double Fw = rho * u[i-1][j] * dy;
            double Fn = rho * v[i][j] * dx;
            double Fs = rho * v[i][j-1] * dx;

            // diffusive fluxes
            double mu_t_e = 0.5 * (mu_t[i][j] + mu_t[i+1][j]);
            double mu_t_w = 0.5 * (mu_t[i][j] + mu_t[i-1][j]);
            double mu_t_n = 0.5 * (mu_t[i][j] + mu_t[i][j+1]);
            double mu_t_s = 0.5 * (mu_t[i][j] + mu_t[i][j-1]);

            double De = (mu + sigma_k * mu_t_e) * dy / dx;
            double Dw = (mu + sigma_k * mu_t_w) * dy / dx;
            double Dn = (mu + sigma_k * mu_t_n) * dx / dy;
            double Ds = (mu + sigma_k * mu_t_s) * dx / dy;

            // tvd scheme
            double r_e = 0.0, r_w = 0.0, r_n = 0.0, r_s = 0.0;

            // east face
            if (Fe > 0) {
                r_e = (k[i][j] - k[i-1][j]) / (k[i+1][j] - k[i][j] + 1e-9);
            } else {
                if (i < x - 1) { 
                    r_e = (k[i+1][j] - k[i+2][j]) / (k[i][j] - k[i+1][j] + 1e-9);
                }
            }

            // west face
            if (Fw > 0) {
                if (i > 1) {
                    r_w = (k[i-1][j] - k[i-2][j]) / (k[i][j] - k[i-1][j] + 1e-9);
                }
            } else {
                r_w = (k[i][j] - k[i+1][j]) / (k[i-1][j] - k[i][j] + 1e-9);
            }

            // north face
            if (Fn > 0) {
                r_n = (k[i][j] - k[i][j-1]) / (k[i][j+1] - k[i][j] + 1e-9);
            } else {
                if (j < y - 1) {
                    r_n = (k[i][j+1] - k[i][j+2]) / (k[i][j] - k[i][j+1] + 1e-9);
                }
            }

            // south face
            if (Fs > 0) {
                if (j > 1) {
                    r_s = (k[i][j-1] - k[i][j-2]) / (k[i][j] - k[i][j-1] + 1e-9);
                }
            } else {
                r_s = (k[i][j] - k[i][j+1]) / (k[i][j-1] - k[i][j] + 1e-9);
            }

            double alpha_e = (Fe > 0) ? 1 : 0;
            double alpha_w = (Fw > 0) ? 1 : 0;
            double alpha_n = (Fn > 0) ? 1 : 0;
            double alpha_s = (Fs > 0) ? 1 : 0;

            double F_corr_e = 0.5 * Fe * ((1 - alpha_e) * fluxLimiter(r_e) - fluxLimiter(r_e) * alpha_e) * (k[i+1][j] - k[i][j]);
            double F_corr_w = 0.5 * Fw * (fluxLimiter(r_w) * alpha_w - (1 - alpha_w) * fluxLimiter(r_w)) * (k[i][j] - k[i-1][j]);
            double F_corr_n = 0.5 * Fn * ((1 - alpha_n) * fluxLimiter(r_n) - fluxLimiter(r_n) * alpha_n) * (k[i][j+1] - k[i][j]);
            double F_corr_s = 0.5 * Fs * (fluxLimiter(r_s) * alpha_s - (1 - alpha_s) * fluxLimiter(r_s)) * (k[i][j] - k[i][j-1]);

            double S_k_DC = F_corr_w + F_corr_e + F_corr_s + F_corr_n;
            
            double aw = Dw + max(Fw, 0.0);
            double ae = De + max(-Fe, 0.0);
            double as = Ds + max(Fs, 0.0);
            double an = Dn + max(-Fn, 0.0);
            
            double su = Pk[i][j] * dx * dy;
            double sp = -beta_star * rho * omega[i][j] * dx * dy;
            double ap_k_unrelaxed = ae + aw + an + as + (Fe - Fw + Fn - Fs) - sp;

            if (abs(ap_k_unrelaxed) > 1e-9) {
                double ap_k = ap_k_unrelaxed / alpha_k;
                double source = ae * k[i+1][j] + aw * k[i-1][j] + an * k[i][j+1] + as * k[i][j-1] + S_k_DC + su + ((1-alpha_k) * ap_k * k[i][j]);
                k[i][j] = source / ap_k;
            }
        }
    }

    // --- Solve omega equation --- //
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y; j++) {

            // blended constants with F1
            double sigma_omega = F1[i][j] * sigma_omega_1 + (1 - F1[i][j]) * sigma_omega_2;
            double beta = F1[i][j] * beta_1 + (1 - F1[i][j]) * beta_2;
            double gamma = F1[i][j] * gamma_1 + (1 - F1[i][j]) * gamma_2;
            
            // convective fluxes
            double Fe = rho * u[i][j] * dy;
            double Fw = rho * u[i-1][j] * dy;
            double Fn = rho * v[i][j] * dx;
            double Fs = rho * v[i][j-1] * dx;

            // diffusive fluxes
            double mu_t_e = 0.5 * (mu_t[i][j] + mu_t[i+1][j]);
            double mu_t_w = 0.5 * (mu_t[i][j] + mu_t[i-1][j]);
            double mu_t_n = 0.5 * (mu_t[i][j] + mu_t[i][j+1]);
            double mu_t_s = 0.5 * (mu_t[i][j] + mu_t[i][j-1]);

            double De = (mu + sigma_omega * mu_t_e) * dy / dx;
            double Dw = (mu + sigma_omega * mu_t_w) * dy / dx;
            double Dn = (mu + sigma_omega * mu_t_n) * dx / dy;
            double Ds = (mu + sigma_omega * mu_t_s) * dx / dy;

            // tvd scheme
            double r_e = 0.0, r_w = 0.0, r_n = 0.0, r_s = 0.0;

            // east face
            if (Fe > 0) {
                r_e = (omega[i][j] - omega[i-1][j]) / (omega[i+1][j] - omega[i][j] + 1e-9);
            } else {
                if (i < x - 1) {
                    r_e = (omega[i+1][j] - omega[i+2][j]) / (omega[i][j] - omega[i+1][j] + 1e-9);
                }
            }

            // west face
            if (Fw > 0) {
                if (i > 1) {
                    r_w = (omega[i-1][j] - omega[i-2][j]) / (omega[i][j] - omega[i-1][j] + 1e-9);
                }
            } else {
                r_w = (omega[i][j] - omega[i+1][j]) / (omega[i-1][j] - omega[i][j] + 1e-9);
            }

            // north face
            if (Fn > 0) {
                r_n = (omega[i][j] - omega[i][j-1]) / (omega[i][j+1] - omega[i][j] + 1e-9);
            } else {
                if (j < y - 1) {
                    r_n = (omega[i][j+1] - omega[i][j+2]) / (omega[i][j] - omega[i][j+1] + 1e-9);
                }
            }

            // south face
            if (Fs > 0) {
                if (j > 1) {
                    r_s = (omega[i][j-1] - omega[i][j-2]) / (omega[i][j] - omega[i][j-1] + 1e-9);
                }
            } else {
                r_s = (omega[i][j] - omega[i][j+1]) / (omega[i][j-1] - omega[i][j] + 1e-9);
            }

            double alpha_e = (Fe > 0) ? 1 : 0;
            double alpha_w = (Fw > 0) ? 1 : 0;
            double alpha_n = (Fn > 0) ? 1 : 0;
            double alpha_s = (Fs > 0) ? 1 : 0;

            double F_corr_e = 0.5 * Fe * ((1 - alpha_e) * fluxLimiter(r_e) - fluxLimiter(r_e) * alpha_e) * (omega[i+1][j] - omega[i][j]);
            double F_corr_w = 0.5 * Fw * (fluxLimiter(r_w) * alpha_w - (1 - alpha_w) * fluxLimiter(r_w)) * (omega[i][j] - omega[i-1][j]);
            double F_corr_n = 0.5 * Fn * ((1 - alpha_n) * fluxLimiter(r_n) - fluxLimiter(r_n) * alpha_n) * (omega[i][j+1] - omega[i][j]);
            double F_corr_s = 0.5 * Fs * (fluxLimiter(r_s) * alpha_s - (1 - alpha_s) * fluxLimiter(r_s)) * (omega[i][j] - omega[i][j-1]);

            double S_omega_DC = F_corr_w + F_corr_e + F_corr_s + F_corr_n;

            double aw = Dw + max(Fw, 0.0);
            double ae = De + max(-Fe, 0.0);
            double as = Ds + max(Fs, 0.0);
            double an = Dn + max(-Fn, 0.0);
            
            // Cross-diffusion term
            double dkdx = (k[i+1][j] - k[i-1][j]) / (2 * dx);
            double dkdy = (k[i][j+1] - k[i][j-1]) / (2 * dy);
            double domegadx = (omega[i+1][j] - omega[i-1][j]) / (2 * dx);
            double domegady = (omega[i][j+1] - omega[i][j-1]) / (2 * dy);
            double grad_k_dot_grad_omega = dkdx * domegadx + dkdy * domegady;
            double cross_diffusion_term = 2 * (1 - F1[i][j]) * rho * sigma_omega_2 * (grad_k_dot_grad_omega / (omega[i][j] + 1e-9));

            double su = (gamma / (mu_t[i][j]/rho + 1e-9) * Pk[i][j] + cross_diffusion_term) * dx * dy;
            double sp = -beta * rho * omega[i][j] * dx * dy;
            
            double ap_omega_unrelaxed = ae + aw + an + as + (Fe - Fw + Fn - Fs) - sp;

            if (abs(ap_omega_unrelaxed) > 1e-9) {
                double ap_o = ap_omega_unrelaxed / alpha_omega;
                double source = ae * omega[i+1][j] + aw * omega[i-1][j] + an * omega[i][j+1] + as * omega[i][j-1] + S_omega_DC + su + ((1-alpha_omega) * ap_o * omega[i][j]);
                omega[i][j] = source / ap_o;
            }
        }
    }
}

void writeResults(const vector<vector<double>>& u, const vector<vector<double>>& v, const vector<vector<double>>& p, const vector<vector<double>>& k, const vector<vector<double>>& omega, const vector<vector<double>>& mu_t, int iter = 1) {
    string name = "results-" + to_string(iter) + ".csv";
    ofstream outfile(name);
    if (!outfile.is_open()) {
        cerr << "Error: Could not open results.csv for writing." << endl;
        return;
    }

    // Write CSV header
    outfile << "x,y,X,Y,U,V,P,k,omega,mu_t\n";

    // Output data at cell centers
    for (int j = 0; j < y; ++j) {
        for (int i = 0; i < x; ++i) {
            double x_pos = i * dx;
            double y_pos = j * dy;

            // Interpolate velocities and pressure from faces to cell center
            double u_center = 0.5 * (u[i][j] + u[i][j+1]);
            double v_center = 0.5 * (v[i][j] + v[i+1][j]);
            double p_center = 0.25 * (p[i][j] + p[i+1][j] + p[i][j+1] + p[i+1][j+1]);
            double k_center = 0.25 * (k[i][j] + k[i+1][j] + k[i][j+1] + k[i+1][j+1]);
            double omega_center = 0.25 * (omega[i][j] + omega[i+1][j] + omega[i][j+1] + omega[i+1][j+1]);
            double mu_t_center = 0.25 * (mu_t[i][j] + mu_t[i+1][j] + mu_t[i][j+1] + mu_t[i+1][j+1]);

            outfile << i+1 << "," << j+1 << "," << x_pos << "," << y_pos << "," << u_center << "," << v_center << "," << p_center << "," << k_center << "," << omega_center << "," << mu_t_center << "\n";
        }
    }

    outfile.close();
    cout << "\nResults written to results.csv" << endl;
}

int main() {
    // --- Initialize Flow Field Variables ---
    // Staggered grid arrangement
    vector<vector<double>> u(x, vector<double>(y + 1, 0));      // u-velocity on vertical faces
    vector<vector<double>> u_old(x, vector<double>(y + 1, 0));
    vector<vector<double>> u_star(x, vector<double>(y + 1, 0));
    vector<vector<double>> ap_u(x, vector<double>(y + 1, 0));

    vector<vector<double>> v(x + 1, vector<double>(y, 0));      // v-velocity on horizontal faces
    vector<vector<double>> v_star(x + 1, vector<double>(y, 0));
    vector<vector<double>> ap_v(x + 1, vector<double>(y, 0));

    vector<vector<double>> p(x + 1, vector<double>(y + 1, 0));  // Pressure at cell centers
    vector<vector<double>> p_prime(x + 1, vector<double>(y + 1, 0));

    // Turbulence quantities at cell centers
    vector<vector<double>> k(x + 1, vector<double>(y + 1, 1e-6));
    vector<vector<double>> omega(x + 1, vector<double>(y + 1, 1e-6));
    vector<vector<double>> mu_t(x + 1, vector<double>(y + 1, 0));
    vector<vector<double>> tensor_dot(x + 1, vector<double>(y + 1, 0));
    
    int iter = 0;
    double max_err = 1.0;
    int centerline_index = (x + 1) / 2; // For convergence check

    while (max_err > max_tol) {
        setBoundaryConditions(u, v, k, omega);

        calculateTurbulentViscosity(u, v, k, omega, tensor_dot, mu_t);

        solveMomentum(u, v, u_star, v_star, ap_u, ap_v, p, mu_t);

        pressureCorrection(u_star, v_star, p_prime, ap_u, ap_v);

        correctFields(u, v, p, u_star, v_star, p_prime, ap_u, ap_v);
        
        solveKOmegaSST(u, v, tensor_dot, mu_t, k, omega);

        // --- Check for convergence ---
        max_err = 0;
        for (int j = 0; j <= y; j++) {
            max_err = max(max_err, abs(u_old[centerline_index][j] - u[centerline_index][j]));
        }
        
        cout << "Iteration: " << iter << ", Residual (u-centerline): " << max_err << "\r" << flush;

        if (iter > 0 && iter % 1000 == 0) {
            writeResults(u, v, p, k, omega, mu_t, iter);
        }
        
        iter++;
        u_old = u;
    }

    cout << "\n\nConvergence reached after " << iter << " iterations." << endl;
    writeResults(u, v, p, k, omega, mu_t, iter);

    return 0;
}

