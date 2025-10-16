#include <iostream>
#include <vector>
#include <fstream>
#include <string>

using namespace std;

// number of grids
const int x = 129;
const int y = 129;

// dimension of cavity
const double lx = 1;
const double ly = 1;

const double dx = lx / (x - 1); //  x spacing between grids
const double dy = ly / (y - 1); // y spacing between grids

// const double dt = 0.0005; // time step
// const int max_iter = 20000;
const double max_tol = 1e-8;

// relaxation factors
const double alpha_p = 0.3;
const double alpha_u = 0.5;
const double alpha_v = 0.5;
const double alpha_k = 0.5;
const double alpha_epsilon = 0.5; 

// constants for k-epsilon model
const double kappa = 0.4187; // Von Karman's factor
// const double E = 9.793;
const double C_mu = 0.09;
const double sigma_k = 1.00;
const double sigma_epsilon = 1.30;
const double C_1epsilon = 1.44;
const double C_2epsilon = 1.92;

const double rho = 50; // density
const double mu = 0.01; // dynamic viscosity
const double u_lid = 1; // top wall velocity
const double nu = mu / rho;
double Re = rho * u_lid * lx / mu; // Re = 100 for rho = 1

double getWallDistance(int i, int j) {
    double y_dist = min((j - 0.5) * dy, ly - (j - 0.5) * dy);
    double x_dist = min((i - 0.5) * dx, lx - (i - 0.5) * dx);
    return min(x_dist, y_dist);
}

void setBoundaryCondition(vector<vector<double>>& u, vector<vector<double>>& v, vector<vector<double>>& k, vector<vector<double>>& epsilon) {
    // u for bottom and top wall
    for (int i = 0; i < x; i++) {
        u[i][0] = -u[i][1];
        u[i][y] = 2.0 * u_lid - u[i][y - 1];
    }

    // v for left and right wall
    for (int j = 0; j < y; j++) {
        v[0][j] = -v[1][j];
        v[x][j] = -v[x - 1][j];
    }

    // u for left and right wall
    for (int j = 0; j < y + 1; j++) {
        u[0][j] = 0;
        u[x - 1][j] = 0;
    }

    // v for top and bottom wall
    for (int i = 0; i < x + 1; i++) {
        v[i][0] = 0;
        v[i][y - 1] = 0;
    }

    // Top wall
    double y_p = dy / 2;
    for (int i = 1; i < x; i++) {
        k[i][y] = k[i][y-1]; // Zero gradient for k
        double u_tau = pow(C_mu, 0.25) * pow(k[i][y-1], 0.5);
        double y_plus = rho * u_tau * y_p / mu;
        epsilon[i][y-1] = y_plus > 11.63 ? (2 * k[i][y-1] * mu / rho * pow(y_p, 2)) : pow(C_mu, 0.75) * pow(k[i][y-1], 0.5) / (kappa * y_p);
    }

    for (int i = 1; i < x; i++) {
        k[i][0] = k[i][1]; // Zero gradient for k
        double u_tau = pow(C_mu, 0.25) * pow(k[i][1], 0.5);
        double y_plus = rho * u_tau * y_p / mu;
        epsilon[i][1] = y_plus > 11.63 ? (2 * k[i][1] * mu / rho * pow(y_p, 2)) : pow(C_mu, 0.75) * pow(k[i][1], 1.5) / (kappa * y_p);
    }

    // Left wall
    double x_p = dx / 2;
    for (int j = 1; j < y; j++) {
        k[0][j] = k[1][j]; // Zero gradient for k
        double u_tau = pow(C_mu, 0.25) * pow(k[1][j], 0.5);
        double x_plus = rho * u_tau * x_p / mu;
        epsilon[1][j] = x_plus > 11.63 ? (2 * k[1][j] * mu / rho * pow(x_p, 2)) : pow(C_mu, 0.75) * pow(k[1][j], 0.5) / (kappa * x_p);
    }
    
    // Right wall
    for (int j = 1; j < y; j++) {
        k[x][j] = k[x-1][j]; // Zero gradient for k
        double u_tau = pow(C_mu, 0.25) * pow(k[1][j], 0.5);
        double x_plus = rho * u_tau * x_p / mu;
        epsilon[x-1][j] = x_plus > 11.63 ? (2 * k[x-1][j] * mu / rho * pow(x_p, 2)) : pow(C_mu, 0.75) * pow(k[x-1][j], 0.5) / (kappa * x_p);
    }
}

void calculateTurbulentViscosity(const vector<vector<double>>& k, const vector<vector<double>>& epsilon, vector<vector<double>>& mu_t) {
    for (int i = 0; i <=x; i++) {
        for (int j = 0; j <= y; j++) {
            double f_mu = 1;

            if (Re <= 3200) {
            double y_wall = getWallDistance(i, j);
            double Re_y = sqrt(k[i][j]) * y_wall / nu;
            double Re_t = k[i][j] * k[i][j] / ((epsilon[i][j] + 1e-10) * nu);
            f_mu = pow(1 - exp(-0.0165 * Re_y), 2) * (1 + (20.5 / (Re_t + 1e-10)));
            }

            mu_t[i][j] = rho * C_mu  * f_mu * k[i][j] * k[i][j] / (epsilon[i][j] + 1e-10);
        }
    }
}

void solveMomentum(const vector<vector<double>>& u, const vector<vector<double>>& v, vector<vector<double>>& u_star, vector<vector<double>>& v_star, vector<vector<double>>& ap_u, vector<vector<double>>& ap_v, const vector<vector<double>>& p, const vector<vector<double>>& mu_t) {
    // solve u momentum equation
    for (int i = 1; i < x - 1; i++) {
        for (int j = 1; j < y; j++) {
            // convective coefficients
            double Fe = rho * 0.5 * (u[i + 1][j] + u[i][j]) * dy;
            double Fw = rho * 0.5 * (u[i - 1][j] + u[i][j]) * dy;
            double Fn = rho * 0.5 * (v[i][j] + v[i + 1][j]) * dx;
            double Fs = rho * 0.5 * (v[i][j - 1] + v[i + 1][j - 1]) * dx;

            // diffusion coefficients
            double De = (mu + mu_t[i][j]) * dy / dx;
            double Dw = (mu + mu_t[i][j]) * dy / dx;
            double Dn = (mu + mu_t[i][j]) * dx / dy;
            double Ds = (mu + mu_t[i][j]) * dx / dy;

            // hybrid scheme
            double ae = max({-Fe, De - 0.5 * Fe, 0.0});
            double aw = max({Fw, Dw + 0.5 * Fw, 0.0});
            double an = max({-Fn, Dn - 0.5 * Fn, 0.0});
            double as = max({Fs, Ds + 0.5 * Fs, 0.0});
            double delta_f = Fe - Fw + Fn - Fs;
            double ap_u_unrelaxed = ae + aw + an + as + delta_f;

            double pres_u = (p[i][j] - p[i + 1][j]) * dy;
            
            if (abs(ap_u_unrelaxed) > 1e-9) {
                ap_u[i][j] = ap_u_unrelaxed / alpha_u;
                double source = ae * u[i + 1][j] + aw * u[i - 1][j] + an * u[i][j + 1] + as * u[i][j - 1] + pres_u + (1.0 - alpha_u) * ap_u[i][j] * u[i][j];
                u_star[i][j] = source / ap_u[i][j];
            } else {
                ap_u[i][j] = 0.0;
                u_star[i][j] = 0.0;
            }
        }
    }

    // v momentum
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y - 1; j++) {
            // convective coefficients
            double Fe = rho * 0.5 * (u[i][j + 1] + u[i][j]) * dy;
            double Fw = rho * 0.5 * (u[i - 1][j + 1] + u[i - 1][j]) * dy;
            double Fn = rho * 0.5 * (v[i][j + 1] + v[i][j]) * dx;
            double Fs = rho * 0.5 * (v[i][j - 1] + v[i][j]) * dx;

            // diffusion coefficients
            double De = (mu + mu_t[i][j]) * dy / dx;
            double Dw = (mu + mu_t[i][j]) * dy / dx;
            double Dn = (mu + mu_t[i][j]) * dx / dy;
            double Ds = (mu + mu_t[i][j]) * dx / dy;

            // hybrid scheme
            double ae = max({-Fe, De - 0.5 * Fe, 0.0});
            double aw = max({Fw, Dw + 0.5 * Fw, 0.0});
            double an = max({-Fn, Dn - 0.5 * Fn, 0.0});
            double as = max({Fs, Ds + 0.5 * Fs, 0.0});
            double delta_f = Fe - Fw + Fn - Fs;
            double ap_v_unrelaxed = ae + aw + an + as + delta_f;

            double pres_v = (p[i][j] - p[i][j + 1]) * dx;
            
            if (abs(ap_v_unrelaxed) > 1e-9) {
                 ap_v[i][j] = ap_v_unrelaxed / alpha_v;
                 double source = ae * v[i + 1][j] + aw * v[i - 1][j] + an * v[i][j + 1] + as * v[i][j - 1] + pres_v + (1.0 - alpha_v) * ap_v[i][j] * v[i][j];
                 v_star[i][j] = source / ap_v[i][j];
            } else {
                ap_v[i][j] = 0.0;
                v_star[i][j] = 0.0;
            }
        }
    }
}

void pressureCorrection(const vector<vector<double>>& u_star, const vector<vector<double>>& v_star, vector<vector<double>>& p_prime, const vector<vector<double>>& ap_u, const vector<vector<double>>& ap_v) {
    // setting all values of  p_prime to zero
    for (auto& row : p_prime) {
        fill(row.begin(), row.end(), 0);
    }

    double p_tol = 1e-5;
    int p_iter = 100;

    for (int iter = 0; iter < p_iter; iter++) {
        double max_p_err = 0;
        for (int i = 1; i < x; i++) {
            for (int j = 1; j < y; j++) {
                double p_prime_old = p_prime[i][j];
                // d coefficients for pressure correction equation
                double de = (i < x - 1 && abs(ap_u[i][j]) > 1e-9) ? dy / ap_u[i][j] : 0;
                double dw = (i > 0 && abs(ap_u[i - 1][j]) > 1e-9) ? dy / ap_u[i - 1][j] : 0;
                double dn = (j < y - 1 && abs(ap_v[i][j]) > 1e-9) ? dx / ap_v[i][j] : 0;
                double ds = (j > 0 && abs(ap_v[i][j - 1]) > 1e-9) ? dx / ap_v[i][j - 1] : 0;

                double ae = rho * de * dy;
                double aw = rho * dw * dy;
                double an = rho * dn * dx;
                double as = rho * ds * dx;
                
                double ap = ae + aw + an + as;

                double b = (rho * u_star[i - 1][j] * dy - rho * u_star[i][j] * dy) + (rho * v_star[i][j - 1] * dx - rho * v_star[i][j] * dx);
                
                if (abs(ap) > 1e-9) {
                    p_prime[i][j] = (ae * p_prime[i + 1][j] + aw * p_prime[i - 1][j] + an * p_prime[i][j + 1] + as * p_prime[i][j - 1] + b) / ap;
                }
                max_p_err = max(max_p_err, abs(p_prime[i][j] - p_prime_old));
            }
        }

        p_prime[1][1] = 0;

        // applying pressure boundary conditions (zerogradient)
        for (int i = 0; i < x + 1; i++) {
            p_prime[i][0] = p_prime[i][1];
            p_prime[i][y] = p_prime[i][y - 1];
        }

        for (int j = 0; j < y + 1; j++) {
            p_prime[0][j] = p_prime[1][j];
            p_prime[x][j] = p_prime[x - 1][j];
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
    for (int i = 0; i < x + 1; i++) {
        for (int j = 0; j < y + 1; j++) {
            p[i][j] += alpha_p * p_prime[i][j];
        }
    }
}

void solveKEpsilon(const vector<vector<double>>& u, const vector<vector<double>>& v, vector<vector<double>>& k, vector<vector<double>>& epsilon, const vector<vector<double>>& mu_t) {
    vector<vector<double>> Pk(x + 1, vector<double>(y + 1, 0));
    vector<vector<double>> k_star = k;
    vector<vector<double>> epsilon_star = epsilon;

    // production term Pk
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y; j++) {
            double dudx = (u[i][j] - u[i-1][j]) / dx;
            double dvdy = (v[i][j] - v[i][j-1]) / dy;
            
            double u_at_nj = 0.5 * (u[i][j] + u[i-1][j]);
            double u_at_sj = 0.5 * (u[i][j-1] + u[i-1][j-1]);
            double dudy = (u_at_nj - u_at_sj) / dy;

            double v_at_ei = 0.5 * (v[i][j] + v[i][j-1]);
            double v_at_wi = 0.5 * (v[i-1][j] + v[i-1][j-1]);
            double dvdx = (v_at_ei - v_at_wi) / dx;

            Pk[i][j] = mu_t[i][j] * (2 * (dudx * dudx + dvdy * dvdy) + (dudy + dvdx)*(dudy + dvdx));
        }
    }

    // solving for k
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y; j++) {
            // convection coefficients
            double Fe = rho * u[i][j] * dy;
            double Fw = rho * u[i-1][j] * dy;
            double Fn = rho * v[i][j] * dx;
            double Fs = rho * v[i][j-1] * dx;

            // diffusion coefficients
            double De = (mu + (0.5 * (mu_t[i][j] + mu_t[i+1][j]) / sigma_k)) * dy / dx;
            double Dw = (mu + (0.5 * (mu_t[i][j] + mu_t[i-1][j]) / sigma_k)) * dy / dx;
            double Dn = (mu + (0.5 * (mu_t[i][j] + mu_t[i][j+1]) / sigma_k)) * dx / dy;
            double Ds = (mu + (0.5 * (mu_t[i][j] + mu_t[i][j-1]) / sigma_k)) * dx / dy;

            // hybrid scheme
            double ae = max({-Fe, De - 0.5 * Fe, 0.0});
            double aw = max({Fw, Dw + 0.5 * Fw, 0.0});
            double an = max({-Fn, Dn - 0.5 * Fn, 0.0});
            double as = max({Fs, Ds + 0.5 * Fs, 0.0});
            double delta_f = Fe - Fw + Fn - Fs;

            double su = Pk[i][j] * dx * dy;
            double sp = - rho * epsilon[i][j] * dx * dy / k[i][j];
            double ap_k = ae + aw + an + as + delta_f - sp;

            if (abs(ap_k) > 1e-9) {
                double source = ae * k[i+1][j] + aw * k[i-1][j] + an * k[i][j+1] + as * k[i][j-1] + su;
                k_star[i][j] = source / ap_k;
            } else {
                k_star[i][j] = 0.0;
            }
            k[i][j] = (1 - alpha_k) * k[i][j] + alpha_k * k_star[i][j];
            k[i][j] = max(k[i][j], 1e-10);
        }
    }

    // solving epsilon
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y; j++) {

            double f_1 = 1;
            double f_2 = 1;

            if (Re <= 3200) {
            double y_wall = getWallDistance(i,j);
            double Re_y = sqrt(k[i][j]) * y_wall / nu;
            double Re_t = k[i][j]*k[i][j] / ((epsilon[i][j] + 1e-10) * nu);
            double f_mu = pow(1 - exp(-0.0165 * Re_y), 2) * (1+ (20.5 / (Re_t + 1e-10)));
            f_1 = pow(1 + (0.05 / f_mu), 3);
            f_2 = 1 - exp(-pow(Re_t, 2));
            }

            // convection coefficients
            double Fe = rho * u[i][j] * dy;
            double Fw = rho * u[i-1][j] * dy;
            double Fn = rho * v[i][j] * dx;
            double Fs = rho * v[i][j-1] * dx;

            // diffusion coefficients
            double De = (mu + (0.5 * (mu_t[i][j] + mu_t[i+1][j]) / sigma_epsilon)) * dy / dx;
            double Dw = (mu + (0.5 * (mu_t[i][j] + mu_t[i-1][j]) / sigma_epsilon)) * dy / dx;
            double Dn = (mu + (0.5 * (mu_t[i][j] + mu_t[i][j+1]) / sigma_epsilon)) * dx / dy;
            double Ds = (mu + (0.5 * (mu_t[i][j] + mu_t[i][j-1]) / sigma_epsilon)) * dx / dy;

            // hybrid scheme
            double ae = max({-Fe, De - 0.5 * Fe, 0.0});
            double aw = max({Fw, Dw + 0.5 * Fw, 0.0});
            double an = max({-Fn, Dn - 0.5 * Fn, 0.0});
            double as = max({Fs, Ds + 0.5 * Fs, 0.0});
            double delta_f = Fe - Fw + Fn - Fs;

            double su = C_1epsilon * f_1 * epsilon[i][j] * Pk[i][j] * dx * dy / k[i][j];
            double sp = - C_2epsilon * f_2 * rho * epsilon[i][j] * dx * dy / k[i][j];
            double ap_e = ae + aw + an + as + delta_f - sp;

            if (abs(ap_e) > 1e-9) {
                double source = ae * epsilon[i+1][j] + aw * epsilon[i-1][j] + an * epsilon[i][j+1] + as * epsilon[i][j-1] + su;
                epsilon_star[i][j] = source / ap_e;
            } else {
                epsilon_star[i][j] = 0.0;
            }
            epsilon[i][j] = (1 - alpha_epsilon) * epsilon[i][j] + alpha_epsilon * epsilon_star[i][j];
            epsilon[i][j] = max(epsilon[i][j], 1e-10);
        }
    }
}

void writeResults(const vector<vector<double>>& u, const vector<vector<double>>& v, const vector<vector<double>>& p, const vector<vector<double>>& k, const vector<vector<double>>& epsilon, const vector<vector<double>>& mu_t, int iter = 1) {
    string name = "results-" + to_string(iter) + ".csv";
    ofstream outfile(name);
    if (!outfile.is_open()) {
        cerr << "Error: Could not open results.csv for writing." << endl;
        return;
    }

    // Write CSV header
    outfile << "x,y,X,Y,U,V,P,k,epsilon,mu_t\n";

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
            double epsilon_center = 0.25 * (epsilon[i][j] + epsilon[i+1][j] + epsilon[i][j+1] + epsilon[i+1][j+1]);
            double mu_t_center = 0.25 * (mu_t[i][j] + mu_t[i+1][j] + mu_t[i][j+1] + mu_t[i+1][j+1]);

            outfile << i+1 << "," << j+1 << "," << x_pos << "," << y_pos << "," << u_center << "," << v_center << "," << p_center << "," << k_center << "," << epsilon_center << "," << mu_t_center << "\n";
        }
    }

    outfile.close();
    cout << "\nResults written to results.csv" << endl;
}

int main() {
    
    // Initializing grids for different flow fields
    vector<vector<double>> u(x, vector<double>(y + 1, 0));
    vector<vector<double>> u_old(x, vector<double>(y + 1, 0));
    vector<vector<double>> u_star(x, vector<double>(y + 1, 0));
    vector<vector<double>> ap_u(x, vector<double>(y + 1, 0));

    vector<vector<double>> v(x + 1, vector<double>(y, 0));
    vector<vector<double>> v_star(x + 1, vector<double>(y, 0));
    vector<vector<double>> ap_v(x + 1, vector<double>(y, 0));

    vector<vector<double>> p(x + 1, vector<double>(y + 1, 0));
    vector<vector<double>> p_prime(x + 1, vector<double>(y + 1, 0));

    vector<vector<double>> k(x + 1, vector<double>(y + 1, 0));
    vector<vector<double>> epsilon(x + 1, vector<double>(y + 1, 0));
    vector<vector<double>> mu_t(x + 1, vector<double>(y + 1, 0));

    setBoundaryCondition(u, v, k, epsilon);
    int iter = 0;
    double max_err = 1;

    while (max_err > max_tol) {
        max_err = 0;

        calculateTurbulentViscosity(k, epsilon, mu_t);

        solveMomentum(u, v, u_star, v_star, ap_u, ap_v, p, mu_t);

        pressureCorrection(u_star, v_star, p_prime, ap_u, ap_v);

        correctFields(u, v, p, u_star, v_star, p_prime, ap_u, ap_v);

        solveKEpsilon(u, v, k, epsilon, mu_t);

        setBoundaryCondition(u, v, k, epsilon);

        for (int j = 0; j <= y; j++) {
            max_err = max(max_err, abs(u_old[65][j] - u[65][j]));
        }
        
        // if (iter % 100 == 0) {
        //     cout << "Iteration: " <<iter << " " << max_err << "\r" << flush;
        // }
        cout << "Iteration: " <<iter << " " << max_err << "\r" << flush;

        if (iter % 2000 == 0) {
            writeResults(u, v, p, k, epsilon, mu_t, iter);
        }
        
        iter += 1;
        u_old = u;
    }

    writeResults(u, v, p, k, epsilon, mu_t);

    return 0;

}