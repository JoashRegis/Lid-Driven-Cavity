#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

// number of grids
const int x = 400; 
const int y = 200;

const double lx = 2.0; // Cavity length
const double ly = 1.0; // Cavity height

const double dx = lx / (x - 1);
const double dy = ly / (y - 1);

// --- Relaxation Factors ---
const double alpha_p = 0.1;
const double alpha_u = 0.5;
const double alpha_v = 0.5;
const double alpha_k = 0.5;
const double alpha_omega = 0.5;

double u_inlet = 1.0;
const double u_lid = 1; // velocity of top squeeze plate in downward direction
const double temperature = 299.904; // assuming isothermal process
const double R = 287.05; // gas constant (air)
const double atm = 101325; // atomspheric pressure
double mu = 1.8e-5; // dynamic viscosity (air)

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

// semicircle placement
const double cx = 0.5; // Center x-coordinate (placed 1/3 downstream)
const double cy = 0.5; // Center y-coordinate (centered vertically)
// const double cy2 = ly;
const double cr = 0.05;      // Cylinder radius
const double cr_sq = cr * cr; // Radius squared (for distance check)

const double T = 10; // total time of solver
const double dt = 0.005; // time step
const double t_loop = 5; // iterations per time step

struct linearSystem {
    vector<vector<double>> aP;
    vector<vector<double>> aE;
    vector<vector<double>> aW;
    vector<vector<double>> aN;
    vector<vector<double>> aS;
    vector<vector<double >> b;

    void resize(int row, int col) {
        aP.assign(row, vector<double>(col, 0.0));
        aE.assign(row, vector<double>(col, 0.0));
        aW.assign(row, vector<double>(col, 0.0));
        aN.assign(row, vector<double>(col, 0.0));
        aS.assign(row, vector<double>(col, 0.0));
        b.assign(row, vector<double>(col, 0.0));
    }
};

void successiveOverRelaxation(vector<vector<double>>& x_sol, const linearSystem sys, int i_start, int i_end, int j_start, int j_end, double omega=1.0) {
    double tol = 1e-6;
    int max_iter = 50;

    for (int iter = 0; iter < max_iter; iter++) {
        double max_residual = 0.0;

        for (int i = i_start; i < i_end; i++) {
            for (int j = j_start; j < j_end; j++) {
                if (sys.aP[i][j] < 1e-9) continue;

                double surr = sys.aE[i][j] * x_sol[i+1][j] + sys.aW[i][j] * x_sol[i-1][j] + sys.aN[i][j] * x_sol[i][j+1] + sys.aS[i][j] * x_sol[i][j-1];
                double x_new = (surr + sys.b[i][j]) / sys.aP[i][j];
                double x_curr = x_sol[i][j];
                x_sol[i][j] = x_curr + omega * (x_new - x_curr);

                double diff = abs(x_sol[i][j] - x_curr);
                if (diff > max_residual) {
                    max_residual = diff;
                }
            }
        }
        if (max_residual < tol) break;
    }
}

// Van Leer flux limiter function for TVD scheme
double fluxLimiter(double r) {
    if (abs(1.0 + abs(r)) < 1e-9) {
        return 0.0;
    }
    return (r + abs(r)) / (1.0 + abs(r));
}

// double getWallDistance(int i, int j) {
//     double x_dis = min(abs((i - 0.5) * dx), lx - abs((i - 0.5) * dx));
//     double y_dis = min(abs((j - 0.5) * dy), ly - abs((j - 0.5) * dy));
//     return min(x_dis, y_dis);
// }

double getWallDistance(int i, int j) {
    double x_pos = abs(i - 0.5) * dx;
    double y_pos = abs(j - 0.5) * dy;

    double dist_x = min(x_pos, lx - x_pos);
    double dist_y = min(y_pos, ly - y_pos);
    double min_dist = min(dist_x, dist_y);

    double dist_to_center = sqrt(pow(x_pos - cx, 2) + pow(y_pos - cy, 2));
    double dist_cylinder = abs(dist_to_center - cr);

    return min(min_dist, dist_cylinder);
}

void solveDensity(vector<vector<double>>& rho, const vector<vector<double>>& p) {
    for (int i = 0; i <= x; i++) {
        for (int j = 0; j <= y; j++) {
            rho[i][j] = p[i][j] / (R * temperature);
        }
    }
}

double massFlowRateRatio(const vector<vector<double>>& u, const vector<vector<double>>& rho) {
    double massFlowIn = 0.0;
    double massFlowOut = 0.0;
    for (int j = 1; j < y; j++) {
        massFlowIn += u[0][j] + 0.5 * (rho[0][j] + rho[1][j]);
        massFlowOut += u[x-1][j] + 0.5 * (rho[x-1][j] + rho[x][j]);
    }

    return massFlowIn / massFlowOut;
}

void setBoundaryConditons(vector<vector<double>>& u, vector<vector<double>>& v, vector<vector<double>>& k, vector<vector<double>>& omega, const vector<vector<double>>& rho) {
    // u for top and bottom wall
    for (int i = 0; i < x; i++) {
        u[i][0] = -u[i][1]; // bottom wall
        u[i][y] = - u[i][y-1]; // top wall
    }

    // u for left and right wall
    for (int j = 0; j <= y; j++) {
        u[0][j] = u_inlet; // left wall
        u[x-1][j] = u[x-2][j]; // right wall
    }

    // v for top and bottom wall
    for (int i = 0; i <= x; i++) {
        v[i][0] = 0.0; // bottom wall
        v[i][y-1] = 0.0; // top wall
    }

    // v for left and right wall
    for (int j = 0; j < y; j++) {
        v[0][j] = -v[1][j]; // left wall
        v[x][j] = v[x-1][j]; // right wall
    }

    // k and omega for top and bottom wall
    for (int i = 1; i < x; i++) {
        k[i][0] = -k[i][1]; // bottom wall
        k[i][y] = -k[i][y-1]; // top wall

        omega[i][0] = 2 * ((60 * mu) / (rho[i][1] * beta_1 * pow(getWallDistance(i, 1), 2.0))) - omega[i][1]; // bottom wall
        omega[i][y] = 2 * ((60 * mu) / (rho[i][y-1] * beta_1 * pow(getWallDistance(i, y-1), 2.0))) - omega[i][y-1]; // top wall
    }

    // k and omega for left and right wall
    // for (int j = 1; j < y; j++) {
    //     k[0][j] = -k[1][j]; // left wall
    //     k[x][j] = -k[x-1][j]; // right wall

    //     omega[0][j] = 2 * ((60 * mu) / (rho[1][j] * beta_1 * pow(getWallDistance(1, j), 2.0))) - omega[1][j]; // left wall
    //     omega[x][j] = 2 * ((60 * mu) / (rho[x-1][j] * beta_1 * pow(getWallDistance(x-1, j), 2.0))) - omega[x-1][j]; // right wall
    // }
}

void solveEddyViscosity(const vector<vector<double>>& u, const vector<vector<double>>& v, const vector<vector<double>>& rho, const vector<vector<double>>& k, const vector<vector<double>>& omega, vector<vector<double>>& mu_t, vector<vector<double>>& tensor_dot) {
    // calculate tensor dot
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y; j++) {
            double dudx = (u[i][j] - u[i-1][j]) / dx;
            double dvdy = (v[i][j] - v[i][j-1]) / dy;

            double dudy = 0.25 * (u[i][j+1] + u[i-1][j+1] - u[i][j-1] - u[i-1][j-1]) / dy;
            double dvdx = 0.25 * (v[i+1][j] + v[i+1][j-1] - v[i-1][j] - v[i-1][j-1]) / dx;

            tensor_dot[i][j] = pow(dudx, 2.0) + pow(dvdy, 2.0) + 0.5 * pow(dudy + dvdx, 2.0);
            double S = sqrt(2.0 * tensor_dot[i][j]);
            double arg_2 = max((2 * sqrt(k[i][j]) / (beta_star * omega[i][j] * getWallDistance(i, j))), (500 * mu) / (rho[i][j] * omega[i][j] * pow(getWallDistance(i, j), 2.0)));
            double F_2 = tanh(arg_2);

            mu_t[i][j] = (a_1_const * rho[i][j] * k[i][j]) / max(a_1_const * omega[i][j], S * F_2);
        }
    }
}

void solveMomentum(vector<vector<double>>& u, vector<vector<double>>& v, const vector<vector<double>>& p, const vector<vector<double>>& rho, const vector<vector<double>>& mu_t, vector<vector<double>>& u_star, vector<vector<double>>& v_star, vector<vector<double>>& ap_u, vector<vector<double>>& ap_v, const vector<vector<double>>& u_old, const vector<vector<double>>& v_old, const vector<vector<double>>& rho_old) {
    static linearSystem U;
    static linearSystem V;

    U.resize(x, y+1);
    V.resize(x+1, y);

    // solve u momentum
    for (int i = 1; i < x-1; i++) {
        for (int j = 1; j < y; j++) {
            double Fe = rho[i+1][j] * 0.5 * (u[i+1][j] + u[i][j]) * dy;
            double Fw = rho[i][j] * 0.5 * (u[i][j] + u[i-1][j]) * dy;
            double Fn = 0.25 * (rho[i+1][j] + rho[i][j] + rho[i+1][j+1] + rho[i][j+1]) * 0.5 * (v[i+1][j] + v[i][j]) * dx;
            double Fs = 0.25 * (rho[i+1][j] + rho[i][j] + rho[i+1][j-1] + rho[i][j-1]) * 0.5 * (v[i+1][j-1] + v[i][j-1]) * dx;

            double De = (mu + mu_t[i+1][j]) * dy / dx;
            double Dw = (mu + mu_t[i][j]) * dy / dx;
            double Dn = (mu + 0.25 * (mu_t[i+1][j] + mu_t[i][j] + mu_t[i][j+1] + mu_t[i+1][j+1])) * dx / dy;
            double Ds = (mu + 0.25 * (mu_t[i+1][j] + mu_t[i][j] + mu_t[i+1][j-1] + mu_t[i][j-1])) * dx / dy;

            // TVD neighbour coefficients
            double aw = Dw + max(Fw, 0.0);
            double ae = De + max(-Fe, 0.0);
            double as = Ds + max(Fs, 0.0);
            double an = Dn + max(-Fn, 0.0);

            double ap0 = 0.5 * (rho_old[i][j] + rho_old[i+1][j])  * dx * dy / dt;
            double delta_f = Fe - Fw + Fn - Fs;
            double ap_u_unrelaxed = ae + aw + an + as + ap0 + delta_f;
            double pres_u = (p[i][j] - p[i + 1][j]) * dy;

            double re = 0.0, rw = 0.0, rn = 0.0, rs = 0.0;

            if (Fe > 0.0) {
                re = (u[i][j] - u[i-1][j]) / (u[i+1][j] - u[i][j] + 1e-9);
            } else {
                if (i+2 <= x-1) {
                    re = (u[i+1][j] - u[i+2][j]) / (u[i][j] - u[i+1][j] + 1e-9);
                }
            }

            if (Fw > 0.0) {
                if (i-2 >= 0) {
                    rw = (u[i-1][j] - u[i-2][j]) / (u[i][j] - u[i-1][j] + 1e-9);
                }
            } else {
                rw = (u[i][j] - u[i+1][j]) / (u[i-1][j] - u[i][j] + 1e-9);
            }

            if (Fn > 0.0) {
                rn = (u[i][j] - u[i][j-1]) / (u[i][j+1] - u[i][j] + 1e-9);
            } else {
                if (j+2 <= y) {
                    rn = (u[i][j+1] - u[i][j+2]) / (u[i][j] - u[i][j+1] + 1e-9);
                }
            }

            if (Fs > 0.0) {
                if (j-2 >= 0) {
                    rs = (u[i][j-1] - u[i][j-2]) / (u[i][j] - u[i][j-1] + 1e-9);
                }
            } else {
                rs = (u[i][j] - u[i][j+1]) / (u[i][j-1] - u[i][j] + 1e-9);
            }

            double alpha_e = (Fe > 0) ? 1 : 0;
            double alpha_w = (Fw > 0) ? 1 : 0;
            double alpha_n = (Fn > 0) ? 1 : 0;
            double alpha_s = (Fs > 0) ? 1 : 0;

            double S_u_DC = 0.5 * Fe * ((1 - alpha_e) * fluxLimiter(re) - alpha_e * fluxLimiter(re)) * (u[i+1][j] - u[i][j]) + 0.5 * Fw * (alpha_w * fluxLimiter(rw) - (1 - alpha_w) * fluxLimiter(rw)) * (u[i][j] - u[i-1][j]) + 0.5 * Fn * ((1 - alpha_n) * fluxLimiter(rn) - alpha_n * fluxLimiter(rn)) * (u[i][j+1] - u[i][j]) + 0.5 * Fs * (alpha_s * fluxLimiter(rs) - (1 - alpha_s) * fluxLimiter(rs)) * (u[i][j] - u[i][j-1]);

            if (i == x-2) {
                u_star[i+1][j] = u_star[i][j] * massFlowRateRatio(u_star, rho);
            }

            double x_pos = abs(i - 0.5) * dx;
            double y_pos = abs(j - 0.5) * dy;

            double dist_sq1 = (x_pos - cx)*(x_pos - cx) + (y_pos - cy)*(y_pos - cy);
            // double dist_sq2 = (x_pos - cx)*(x_pos - cx) + (y_pos - cy2)*(y_pos - cy2);

            if (dist_sq1 <= cr_sq) {
                U.aE[i][j] = 0.0;
                U.aW[i][j] = 0.0;
                U.aN[i][j] = 0.0;
                U.aS[i][j] = 0.0;
                U.aP[i][j] = 1e50;
                U.b[i][j] = 0.0;
            } else {
                U.aE[i][j] = ae;
                U.aW[i][j] = aw;
                U.aN[i][j] = an;
                U.aS[i][j] = as;
                U.aP[i][j] = ap_u_unrelaxed / alpha_u;
                U.b[i][j] = pres_u + S_u_DC + ap0 * u_old[i][j] + (1.0 - alpha_u) * U.aP[i][j] * u[i][j];
            }
            ap_u[i][j] = U.aP[i][j];
        }
    }
    successiveOverRelaxation(u_star, U, 1, x-1, 1, y, alpha_u);

    // solve v momentum equation
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y-1; j++) {
            double Fe = 0.25 * (rho[i][j] + rho[i+1][j] + rho[i][j+1] + rho[i+1][j+1]) * 0.5 * (u[i][j] + u[i][j+1]) * dy;
            double Fw = 0.25 * (rho[i][j] + rho[i-1][j] + rho[i][j+1] + rho[i-1][j+1]) * 0.5 * (u[i-1][j] + u[i-1][j+1]) * dy;
            double Fn = rho[i][j+1] * 0.5 * (v[i][j] + v[i][j+1]) * dx;
            double Fs = rho[i][j] * 0.5 * (v[i][j] + v[i][j-1]) * dx;

            double De = (mu + 0.25 * (mu_t[i][j] + mu_t[i+1][j] + mu_t[i][j+1] + mu_t[i+1][j+1])) * dy / dx;
            double Dw = (mu + 0.25 * (mu_t[i][j] + mu_t[i-1][j] + mu_t[i][j+1] + mu_t[i-1][j+1])) * dy / dx;
            double Dn = (mu + mu_t[i][j+1]) * dx / dy;
            double Ds = (mu + mu_t[i][j]) * dx / dy;

            // TVD neighbour coefficients
            double aw = Dw + max(Fw, 0.0);
            double ae = De + max(-Fe, 0.0);
            double as = Ds + max(Fs, 0.0);
            double an = Dn + max(-Fn, 0.0);

            double ap0 = 0.5 * (rho_old[i][j] + rho_old[i+1][j])  * dx * dy / dt;
            double delta_f = Fe - Fw + Fn - Fs;
            double ap_v_unrelaxed = ae + aw + an + as + ap0 + delta_f;
            double pres_v = (p[i][j] - p[i][j + 1]) * dx;

            double re = 0.0, rw = 0.0, rn = 0.0, rs = 0.0;

            if (Fe > 0.0) {
                re = (v[i][j] - v[i-1][j]) / (v[i+1][j] - v[i][j] + 1e-9);
            } else {
                if (i+2 <= x) {
                    re = (v[i+1][j] - v[i+2][j]) / (v[i][j] - v[i+1][j] + 1e-9);
                }
            }

            if (Fw > 0.0) {
                if (i-2 >= 0) {
                    rw = (v[i-1][j] - v[i-2][j]) / (v[i][j] - v[i-1][j] + 1e-9);
                }
            } else {
                rw = (v[i][j] - v[i+1][j]) / (v[i-1][j] - v[i][j] + 1e-9);
            }

            if (Fn > 0.0) {
                rn = (v[i][j] - v[i][j-1]) / (v[i][j+1] - v[i][j] + 1e-9);
            } else {
                if (j+2 <= y-1) {
                    rn = (v[i][j+1] - v[i][j+2]) / (v[i][j] - v[i][j+1] + 1e-9);
                }
            }

            if (Fs > 0.0) {
                if (j-2 >= 0) {
                    rs = (v[i][j-1] - v[i][j-2]) / (v[i][j] - v[i][j-1] + 1e-9);
                }
            } else {
                rs = (v[i][j] - v[i][j+1]) / (v[i][j-1] - v[i][j] + 1e-9);
            }

            double alpha_e = (Fe > 0) ? 1 : 0;
            double alpha_w = (Fw > 0) ? 1 : 0;
            double alpha_n = (Fn > 0) ? 1 : 0;
            double alpha_s = (Fs > 0) ? 1 : 0;

            double S_v_DC = 0.5 * Fe * ((1 - alpha_e) * fluxLimiter(re) - alpha_e * fluxLimiter(re)) * (v[i+1][j] - v[i][j]) + 0.5 * Fw * (alpha_w * fluxLimiter(rw) - (1 - alpha_w) * fluxLimiter(rw)) * (v[i][j] - v[i-1][j]) + 0.5 * Fn * ((1 - alpha_n) * fluxLimiter(rn) - alpha_n * fluxLimiter(rn)) * (v[i][j+1] - v[i][j]) + 0.5 * Fs * (alpha_s * fluxLimiter(rs) - (1 - alpha_s) * fluxLimiter(rs)) * (v[i][j] - v[i][j-1]);

            double x_pos = abs(i - 0.5) * dx;
            double y_pos = abs(j - 0.5) * dy;

            double dist_sq1 = (x_pos - cx)*(x_pos - cx) + (y_pos - cy)*(y_pos - cy);
            // double dist_sq2 = (x_pos - cx)*(x_pos - cx) + (y_pos - cy2)*(y_pos - cy2);

            if (dist_sq1 <= cr_sq) {
                V.aE[i][j] = 0.0;
                V.aW[i][j] = 0.0;
                V.aN[i][j] = 0.0;
                V.aS[i][j] = 0.0;
                V.aP[i][j] = 1e50;
                V.b[i][j] = 0.0;
            } else {
                V.aE[i][j] = ae;
                V.aW[i][j] = aw;
                V.aN[i][j] = an;
                V.aS[i][j] = as;
                V.aP[i][j] = ap_v_unrelaxed / alpha_v;
                V.b[i][j] = pres_v + S_v_DC + ap0 * v_old[i][j] + (1.0 - alpha_v) * V.aP[i][j] * v[i][j];
            }
            ap_v[i][j] = U.aP[i][j];
        }
    }
    successiveOverRelaxation(v_star, V, 1, x, 1, y-1, alpha_v);
}

void pressureCorrection(const vector<vector<double>>& u_star, const vector<vector<double>>& v_star, vector<vector<double>>& p_prime, const vector<vector<double>>& ap_u, const vector<vector<double>>& ap_v, const vector<vector<double>>& rho, const vector<vector<double>>& rho_old) {
    // setting all values of  p_prime to zero
    for (auto& row : p_prime) {
        fill(row.begin(), row.end(), 0);
    }

    static linearSystem P;
    P.resize(x+1, y+1);
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y; j++) {
                // d coefficients for pressure correction equation
            double de = (i < x - 1 && abs(ap_u[i][j]) > 1e-9) ? dy / ap_u[i][j] : 0;
            double dw = (i > 0 && abs(ap_u[i - 1][j]) > 1e-9) ? dy / ap_u[i - 1][j] : 0;
            double dn = (j < y - 1 && abs(ap_v[i][j]) > 1e-9) ? dx / ap_v[i][j] : 0;
            double ds = (j > 0 && abs(ap_v[i][j - 1]) > 1e-9) ? dx / ap_v[i][j - 1] : 0;

            double rho_e = 0.5 * (rho[i][j] + rho[i+1][j]);
            double rho_w = 0.5 * (rho[i][j] + rho[i-1][j]);
            double rho_n = 0.5 * (rho[i][j] + rho[i][j+1]);
            double rho_s = 0.5 * (rho[i][j] + rho[i][j-1]);

            double ae = rho_e * de * dy;
            double aw = rho_w * dw * dy;
            double an = rho_n * dn * dx;
            double as = rho_s * ds * dx;
                
            double ap = ae + aw + an + as;

            double b = (rho_w * u_star[i-1][j] * dy) - (rho_e * u_star[i][j] * dy) + (rho_s * v_star[i][j - 1] * dx) - (rho_n * v_star[i][j] * dx) + ((rho_old[i][j] - rho[i][j]) * dx * dy / dt);
                
            if (abs(ap) > 1e-12) {
                P.aE[i][j] = ae;
                P.aW[i][j] = aw;
                P.aN[i][j] = an;
                P.aS[i][j] = as;
                P.aP[i][j] = ap;
                P.b[i][j] = b;
            } else {
                P.aP[i][j] = 1;
                P.b[i][j] = 0;
            }
        }
    }

    successiveOverRelaxation(p_prime, P, 1, x, 1, y, 1.5);

    // p_prime[1][1] = 0;

    // applying pressure boundary conditions (zerogradient)
    for (int i = 0; i < x + 1; i++) {
        p_prime[i][0] = p_prime[i][1];
        p_prime[i][y] = p_prime[i][y - 1];
    }

    for (int j = 0; j < y + 1; j++) {
        p_prime[0][j] = p_prime[1][j];
        p_prime[x][j] = 0.0;
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
    for (int i = 0; i < x + 1; i++) {
        for (int j = 0; j < y + 1; j++) {
            p[i][j] += alpha_p * p_prime[i][j];
        }
    }
}

void solveKOmegaSST(const vector<vector<double>>& u, const vector<vector<double>>& v, const vector<vector<double>>& rho, const vector<vector<double>>& rho_old, const vector<vector<double>>& mu_t, const vector<vector<double>>& tensor_dot, vector<vector<double>>& k, const vector<vector<double>>& k_old, vector<vector<double>>& omega, const vector<vector<double>>& omega_old){
    vector<vector<double>> P(x+1, vector<double>(y+1, 0.0));
    vector<vector<double>> F1 = P;
    vector<vector<double>> grad_k_dot_grad_omega = P;

    static linearSystem K;
    static linearSystem OMEGA;

    K.resize(x+1, y+1);
    OMEGA.resize(x+1, y+1);

    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y; j++) {
            // production term
            P[i][j] = min(2 * mu_t[i][j] * tensor_dot[i][j], 10 * beta_star * k[i][j] * omega[i][j]);

            // solve blending function F1
            double dkdx = 0.5 * (k[i+1][j] - k[i-1][j]) / dx;
            double dkdy = 0.5 * (k[i][j+1] - k[i][j-1]) / dy;
            double domegadx = 0.5 * (omega[i+1][j] - omega[i-1][j]) / dx;
            double domegady = 0.5 * (omega[i][j+1] - omega[i][j-1]) / dy;
            grad_k_dot_grad_omega[i][j] = dkdx * domegadx + dkdy * domegady;

            double cross_diffusion = max(2 * rho[i][j] * sigma_omega_2 * grad_k_dot_grad_omega[i][j] / (omega[i][j] + 1e-9), 1.0e-10);
            double arg_1 = min(max(sqrt(k[i][j]) / (beta_star * omega[i][j] * getWallDistance(i, j)), (500 * mu) / (pow(getWallDistance(i, j), 2.0) * omega[i][j])), (4 * rho[i][j] * sigma_omega_2 * k[i][j]) / (cross_diffusion * pow(getWallDistance(i, j), 2.0)));
            F1[i][j] = tanh(pow(arg_1, 4.0));
        }
    }

    // solve k
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y; j++) {
            double sigma_k = F1[i][j] * sigma_k_1 + (1 - F1[i][j]) * sigma_k_2;

            double Fe = 0.5 * (rho[i+1][j] + rho[i][j]) * u[i][j] * dy;
            double Fw = 0.5 * (rho[i-1][j] + rho[i][j]) * u[i-1][j] * dy;
            double Fn = 0.5 * (rho[i][j+1] + rho[i][j]) * v[i][j] * dx;
            double Fs = 0.5 * (rho[i][j-1] + rho[i][j]) * v[i][j-1] * dx;

            double De = (mu + sigma_k * 0.5 * (mu_t[i+1][j] + mu_t[i][j])) * dy / dx;
            double Dw = (mu + sigma_k * 0.5 * (mu_t[i-1][j] + mu_t[i][j])) * dy / dx;
            double Dn = (mu + sigma_k * 0.5 * (mu_t[i][j+1] + mu_t[i][j])) * dx / dy;
            double Ds = (mu + sigma_k * 0.5 * (mu_t[i][j-1] + mu_t[i][j])) * dx / dy;

            // TVD neighbour coefficients
            double aw = Dw + max(Fw, 0.0);
            double ae = De + max(-Fe, 0.0);
            double as = Ds + max(Fs, 0.0);
            double an = Dn + max(-Fn, 0.0);

            double re = 0.0, rw = 0.0, rn = 0.0, rs = 0.0;

            if (Fe > 0.0) {
                re = (k[i][j] - k[i-1][j]) / (k[i+1][j] - k[i][j] + 1e-9);
            } else {
                if (i+2 <= x) {
                    re = (k[i+1][j] - k[i+2][j]) / (k[i][j] - k[i+1][j] + 1e-9);
                }
            }

            if (Fw > 0.0) {
                if (i-2 >= 0) {
                    rw = (k[i-1][j] - k[i-2][j]) / (k[i][j] - k[i-1][j] + 1e-9);
                }
            } else {
                rw = (k[i][j] - k[i+1][j]) / (k[i-1][j] - k[i][j] + 1e-9);
            }

            if (Fn > 0.0) {
                rn = (k[i][j] - k[i][j-1]) / (k[i][j+1] - k[i][j] + 1e-9);
            } else {
                if (j+2 <= y) {
                    rn = (k[i][j+1] - k[i][j+2]) / (k[i][j] - k[i][j+1] + 1e-9);
                }
            }

            if (Fs > 0.0) {
                if (j-2 >= 0) {
                    rs = (k[i][j-1] - k[i][j-2]) / (k[i][j] - k[i][j-1] + 1e-9);
                }
            } else {
                rs = (k[i][j] - k[i][j+1]) / (k[i][j-1] - k[i][j] + 1e-9);
            }

            double alpha_e = (Fe > 0) ? 1 : 0;
            double alpha_w = (Fw > 0) ? 1 : 0;
            double alpha_n = (Fn > 0) ? 1 : 0;
            double alpha_s = (Fs > 0) ? 1 : 0;

            double S_k_DC = 0.5 * Fe * ((1 - alpha_e) * fluxLimiter(re) - alpha_e * fluxLimiter(re)) * (k[i+1][j] - k[i][j]) + 0.5 * Fw * (alpha_w * fluxLimiter(rw) - (1 - alpha_w) * fluxLimiter(rw)) * (k[i][j] - k[i-1][j]) + 0.5 * Fn * ((1 - alpha_n) * fluxLimiter(rn) - alpha_n * fluxLimiter(rn)) * (k[i][j+1] - k[i][j]) + 0.5 * Fs * (alpha_s * fluxLimiter(rs) - (1 - alpha_s) * fluxLimiter(rs)) * (k[i][j] - k[i][j-1]);
            
            
            double su = P[i][j] * dx * dy;
            double sp = -beta_star * rho[i][j] * omega[i][j] * dx * dy;
            double ap0 = 0.5 * (rho_old[i][j] + rho_old[i+1][j])  * dx * dy / dt;
            double ap_k_unrelaxed = ae + aw + an + as + ap0 + (Fe - Fw + Fn - Fs) - sp;

            double x_pos = abs(i - 0.5) * dx;
            double y_pos = abs(j - 0.5) * dy;

            double dist_sq1 = (x_pos - cx)*(x_pos - cx) + (y_pos - cy)*(y_pos - cy);
            if (dist_sq1 <= cr_sq) {
                K.aE[i][j] = 0;
                K.aW[i][j] = 0;
                K.aN[i][j] = 0;
                K.aS[i][j] = 0;
                K.aP[i][j] = 1.0e+50;
                K.b[i][j] = 0;
            } else {
                K.aE[i][j] = ae;
                K.aW[i][j] = aw;
                K.aN[i][j] = an;
                K.aS[i][j] = as;
                K.aP[i][j] = ap_k_unrelaxed / alpha_k;
                K.b[i][j] = sp + S_k_DC + ap0 * k_old[i][j] + (1.0 - alpha_k) * K.aP[i][j] * k[i][j];
            }
        }
    }
    successiveOverRelaxation(k, K, 1, x, 1, y, alpha_k);
    for(int i=0; i<=x; i++) for(int j=0; j<=y; j++) k[i][j] = max(k[i][j], 1e-12);

    // solve omega
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y; j++) {
            double sigma_omega = F1[i][j] * sigma_omega_1 + (1 - F1[i][j]) * sigma_omega_2;
            double gamma = F1[i][j] * gamma_1 + (1 - F1[i][j]) * gamma_2;
            double beta = F1[i][j] * beta_1 + (1 - F1[i][j]) * beta_2;

            double Fe = 0.5 * (rho[i+1][j] + rho[i][j]) * u[i][j] * dy;
            double Fw = 0.5 * (rho[i-1][j] + rho[i][j]) * u[i-1][j] * dy;
            double Fn = 0.5 * (rho[i][j+1] + rho[i][j]) * v[i][j] * dx;
            double Fs = 0.5 * (rho[i][j-1] + rho[i][j]) * v[i][j-1] * dx;

            double De = (mu + sigma_omega * 0.5 * (mu_t[i+1][j] + mu_t[i][j])) * dy / dx;
            double Dw = (mu + sigma_omega * 0.5 * (mu_t[i-1][j] + mu_t[i][j])) * dy / dx;
            double Dn = (mu + sigma_omega * 0.5 * (mu_t[i][j+1] + mu_t[i][j])) * dx / dy;
            double Ds = (mu + sigma_omega * 0.5 * (mu_t[i][j-1] + mu_t[i][j])) * dx / dy;

            // TVD neighbour coefficients
            double aw = Dw + max(Fw, 0.0);
            double ae = De + max(-Fe, 0.0);
            double as = Ds + max(Fs, 0.0);
            double an = Dn + max(-Fn, 0.0);

            double re = 0.0, rw = 0.0, rn = 0.0, rs = 0.0;

            if (Fe > 0.0) {
                re = (omega[i][j] - omega[i-1][j]) / (omega[i+1][j] - omega[i][j] + 1e-9);
            } else {
                if (i+2 <= x) {
                    re = (omega[i+1][j] - omega[i+2][j]) / (omega[i][j] - omega[i+1][j] + 1e-9);
                }
            }

            if (Fw > 0.0) {
                if (i-2 >= 0) {
                    rw = (omega[i-1][j] - omega[i-2][j]) / (omega[i][j] - omega[i-1][j] + 1e-9);
                }
            } else {
                rw = (omega[i][j] - omega[i+1][j]) / (omega[i-1][j] - omega[i][j] + 1e-9);
            }

            if (Fn > 0.0) {
                rn = (omega[i][j] - omega[i][j-1]) / (omega[i][j+1] - omega[i][j] + 1e-9);
            } else {
                if (j+2 <= y) {
                    rn = (omega[i][j+1] - omega[i][j+2]) / (omega[i][j] - omega[i][j+1] + 1e-9);
                }
            }

            if (Fs > 0.0) {
                if (j-2 >= 0) {
                    rs = (omega[i][j-1] - omega[i][j-2]) / (omega[i][j] - omega[i][j-1] + 1e-9);
                }
            } else {
                rs = (omega[i][j] - omega[i][j+1]) / (omega[i][j-1] - omega[i][j] + 1e-9);
            }

            double alpha_e = (Fe > 0) ? 1 : 0;
            double alpha_w = (Fw > 0) ? 1 : 0;
            double alpha_n = (Fn > 0) ? 1 : 0;
            double alpha_s = (Fs > 0) ? 1 : 0;

            double S_omega_DC = 0.5 * Fe * ((1 - alpha_e) * fluxLimiter(re) - alpha_e * fluxLimiter(re)) * (omega[i+1][j] - omega[i][j]) + 0.5 * Fw * (alpha_w * fluxLimiter(rw) - (1 - alpha_w) * fluxLimiter(rw)) * (omega[i][j] - omega[i-1][j]) + 0.5 * Fn * ((1 - alpha_n) * fluxLimiter(rn) - alpha_n * fluxLimiter(rn)) * (omega[i][j+1] - omega[i][j]) + 0.5 * Fs * (alpha_s * fluxLimiter(rs) - (1 - alpha_s) * fluxLimiter(rs)) * (omega[i][j] - omega[i][j-1]);
            
            double cross_diffusion_term = 2 * (1 - F1[i][j]) * rho[i][j] * sigma_omega_2 * (grad_k_dot_grad_omega[i][j] / (omega[i][j] + 1e-9));

            double su = (gamma / (mu_t[i][j]/rho[i][j] + 1e-9) * P[i][j] + cross_diffusion_term) * dx * dy;
            double sp = -beta * rho[i][j] * omega[i][j] * dx * dy;
            double ap0 = 0.5 * (rho_old[i][j] + rho_old[i+1][j])  * dx * dy / dt;
            double ap_omega_unrelaxed = ae + aw + an + as + ap0 + (Fe - Fw + Fn - Fs) - sp;

            double x_pos = abs(i - 0.5) * dx;
            double y_pos = abs(j - 0.5) * dy;

            double dist_sq1 = (x_pos - cx)*(x_pos - cx) + (y_pos - cy)*(y_pos - cy);
            if (dist_sq1 <= cr_sq) {
                OMEGA.aE[i][j] = 0;
                OMEGA.aW[i][j] = 0;
                OMEGA.aN[i][j] = 0;
                OMEGA.aS[i][j] = 0;
                OMEGA.aP[i][j] = 1.0e+50;
                OMEGA.b[i][j] = 0;
            } else {      
                OMEGA.aE[i][j] = ae;
                OMEGA.aW[i][j] = aw;
                OMEGA.aN[i][j] = an;
                OMEGA.aS[i][j] = as;
                OMEGA.aP[i][j] = ap_omega_unrelaxed / alpha_k;
                OMEGA.b[i][j] = sp + S_omega_DC + ap0 * omega_old[i][j] + (1.0 - alpha_omega) * OMEGA.aP[i][j] * omega[i][j];
            }
        }
    }
    successiveOverRelaxation(omega, OMEGA, 1, x, 1, y, alpha_omega);
    for(int i=0; i<=x; i++) for(int j=0; j<=y; j++) omega[i][j] = max(omega[i][j], 1e-12);
}

void writeResults(const vector<vector<double>>& u, const vector<vector<double>>& v, const vector<vector<double>>& p, const vector<vector<double>>& rho, const vector<vector<double>>& k, const vector<vector<double>>& omega, const vector<vector<double>>& mu_t, int iter=1) {
    string name = "results-" + to_string(iter) + ".csv";
    ofstream outfile(name);
    if (!outfile.is_open()) {
        cerr << "Error: Could not open results.csv for writing." << endl;
        return;
    }

    // Write CSV header
    outfile << "x,y,X,Y,U,V,P,rho,k,omega,mu_t\n";

    // Output data at cell centers
    for (int j = 0; j < y; ++j) {
        for (int i = 0; i < x; ++i) {
            double x_pos = i * dx;
            double y_pos = j * dy;

            // Interpolate velocities and pressure from faces/nodes to cell center
            double u_center = 0.5 * (u[i][j] + u[i][j + 1]);
            double v_center = 0.5 * (v[i][j] + v[i + 1][j]);
            double p_center = 0.25 * (p[i][j] + p[i + 1][j] + p[i][j + 1] + p[i + 1][j + 1]);
            double rho_center = 0.25 * (rho[i][j] + rho[i + 1][j] + rho[i][j + 1] + rho[i + 1][j + 1]);
            double k_center = 0.25 * (k[i][j] + k[i+1][j] + k[i][j+1] + k[i+1][j+1]);
            double omega_center = 0.25 * (omega[i][j] + omega[i+1][j] + omega[i][j+1] + omega[i+1][j+1]);
            double mu_t_center = 0.25 * (mu_t[i][j] + mu_t[i+1][j] + mu_t[i][j+1] + mu_t[i+1][j+1]);

            outfile << i + 1 << "," << j + 1 << "," << x_pos << "," << y_pos << "," << u_center << "," << v_center << "," << p_center << "," << rho_center << "," << k_center << "," << omega_center << "," << mu_t_center << "\n";
        }
    }

    outfile.close();
    cout << "\nResults written to results.csv" << endl;
}

int main() {
    // --- Initialize Flow Field Variables ---
    // Staggered grid arrangement
    vector<vector<double>> u(x, vector<double>(y + 1, 0.0));      // u-velocity on vertical faces
    vector<vector<double>> u_old = u;
    vector<vector<double>> u_history = u;
    vector<vector<double>> u_star = u;
    vector<vector<double>> ap_u = u;

    vector<vector<double>> v(x + 1, vector<double>(y, 0.0));      // v-velocity on horizontal faces
    vector<vector<double>> v_star = v;
    vector<vector<double>> v_old = v;
    vector<vector<double>> ap_v = v;

    vector<vector<double>> p(x + 1, vector<double>(y + 1, 101325));  // Pressure at cell centers
    vector<vector<double>> p_prime = p;

    vector<vector<double>> rho(x + 1, vector<double>(y + 1, 1.177));
    vector<vector<double>> rho_old = rho;

    // Turbulence quantities at cell centers
    vector<vector<double>> k(x + 1, vector<double>(y + 1, 1e-6));
    vector<vector<double>> k_old = k;
    vector<vector<double>> omega(x + 1, vector<double>(y + 1, 1e-6));
    vector<vector<double>> omega_old = omega;
    vector<vector<double>> mu_t(x + 1, vector<double>(y + 1, 0));
    vector<vector<double>> tensor_dot(x + 1, vector<double>(y + 1, 0));
    
    int iteration = 0;
    double max_err = 1.0;
    int centerline_index = (x + 1) / 2; // For convergence check

    setBoundaryConditons(u, v, k, omega, rho);
    
    // outer loop
    for (double time = 0; time < T; time += dt) {
        u_old = u;
        v_old = v;
        rho_old = rho;
        k_old = k;
        omega_old = omega;
        iteration++;

        cout << time << endl;
        // inner loop inside a time step
        for (int iter = 0; iter < t_loop; iter++) {
            max_err = 0;
            u_history = u;

            u_star = u;
            v_star = v;

            solveDensity(rho, p);
            
            solveMomentum(u, v, p, rho, mu_t, u_star, v_star, ap_u, ap_v, u_old, v_old, rho_old);

            pressureCorrection(u_star, v_star, p_prime, ap_u, ap_v, rho, rho_old);

            correctFields(u, v, p, u_star, v_star, p_prime, ap_u, ap_v);

            solveKOmegaSST(u, v, rho, rho_old, mu_t, tensor_dot, k, k_old, omega, omega_old);

            setBoundaryConditons(u, v, k, omega, rho);

            for (int j = 0; j <= y; j++) {
                max_err = max(max_err, abs(u_history[51][j] - u[51][j]));
            }
            
            cout << "Iteration: " << iter << " " << max_err << '\r' << flush;
        }
        if (iteration % 40 == 0) {
            writeResults(u, v, p, rho, k, omega, mu_t, iteration);
        }
    }

    //writeResults(u, v, p, rho, dy);

    return 0;
}