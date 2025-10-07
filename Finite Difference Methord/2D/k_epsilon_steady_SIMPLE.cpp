#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm> // For std::max

using namespace std;

// number of grids
const int x = 129;
const int y = 129;

// dimension of cavity
const double lx = 1;
const double ly = 1;

const double dx = lx / (x - 1); //  x spacing between grids
const double dy = ly / (y - 1); // y spacing between grids

const double dt = 0.0005; // time step
const int max_iter = 20000;

// relaxation factors
const double alpha_p = 0.8;
const double alpha_u = 0.7;
const double alpha_v = 0.7;
const double alpha_k = 0.7;
const double alpha_epsilon = 0.7; 

// constants for k-epsilon model
const double kappa = 0.4187; // More precise value from textbook
const double E = 9.793;
const double C_mu = 0.09;
const double sigma_k = 1.00;
const double sigma_epsilon = 1.30;
const double C_1epsilon = 1.44;
const double C_2epsilon = 1.92;

const double rho = 1; // density
const double mu = 0.01; // dynamic viscosity
const double u_lid = 1; // top wall velocity

void setBoundaryConditions(vector<vector<double>>& u, vector<vector<double>>& v, vector<vector<double>>& k, vector<vector<double>>& epsilon) {
    // normal velocity to left and right wall
    for (int j = 0; j <= y; j++) {
        u[0][j] = 0;
        u[x-1][j] = 0;
    }

    // normal velocity to top and bottom wall
    for (int i = 0; i <= x; i++) {
        v[i][0] = 0;
        v[i][y-1] = 0;
    }

    // tangential velocity of right and left wall
    for (int j = 1; j < y - 1; j++) {
        v[0][j] = -v[1][j];
        v[x][j] = -v[x-1][j];
    }

    // tangential velocity for top and bottom wall
    for (int i = 1; i < x - 1; i++) {
        u[i][0] = -u[i][1];
        u[i][y] = 2 * u_lid - u[i][y-1];
    }

    // Boundary conditions k-epsilon
    // Bottom wall
    double y_p_bottom = dy / 2.0;
    for (int i = 1; i < x; i++) {
        k[i][0] = k[i][1]; // Zero gradient for k
        epsilon[i][1] = pow(C_mu, 0.75) * pow(k[i][1], 1.5) / (kappa * y_p_bottom);
    }
    
    // Top wall
    double y_p_top = dy / 2.0;
    for (int i = 1; i < x; i++) {
        k[i][y] = k[i][y-1]; // Zero gradient for k
        epsilon[i][y-1] = pow(C_mu, 0.75) * pow(k[i][y-1], 1.5) / (kappa * y_p_top);
    }


    // Left wall
    double x_p_left = dx / 2.0;
    for (int j = 1; j < y; j++) {
        k[0][j] = k[1][j]; // Zero gradient for k
        epsilon[1][j] = pow(C_mu, 0.75) * pow(k[1][j], 1.5) / (kappa * x_p_left);
    }
    
    // Right wall
    double x_p_right = dx / 2.0;
    for (int j = 1; j < y; j++) {
        k[x][j] = k[x-1][j]; // Zero gradient for k
        epsilon[x-1][j] = pow(C_mu, 0.75) * pow(k[x-1][j], 1.5) / (kappa * x_p_right);
    }
}

void calculateTurbulentViscosity(const vector<vector<double>>& k, const vector<vector<double>>& epsilon, vector<vector<double>>& mu_t) {
    for (int i = 0; i <=x; i++) {
        for (int j = 0; j <= y; j++) {
            mu_t[i][j] = rho * C_mu * k[i][j] * k[i][j] / (epsilon[i][j] + 1e-10);
        }
    }
}

void solveMomentum(const vector<vector<double>>& u, const vector<vector<double>>& v, const vector<vector<double>>& p, vector<vector<double>>& u_star, vector<vector<double>>& v_star, const vector<vector<double>>& mu_t) {
    // solving for x-momentum equation
    for (int i = 1; i < x - 1; i++) {
        for (int j = 1; j < y; j++) {
            // average values of u, v in u-grid cell
            double u_avg_n = 0.5 * (u[i][j] + u[i][j+1]);
            double u_avg_s = 0.5 * (u[i][j] + u[i][j-1]);
            double u_avg_e = 0.5 * (u[i][j] + u[i+1][j]); 
            double u_avg_w = 0.5 * (u[i][j] + u[i-1][j]); 
            double v_avg_n = 0.5 * (v[i][j] + v[i+1][j]);
            double v_avg_s = 0.5 * (v[i][j-1] + v[i+1][j-1]);

            double mu_eff = mu + 0.5 * (mu_t[i][j] + mu_t[i+1][j]);

            // u convective term
            double u_conv = (u_avg_e * u_avg_e - u_avg_w * u_avg_w) / dx + (u_avg_n * v_avg_n - u_avg_s * v_avg_s) / dy;

            // u diffusion term
            double u_diff = (mu_eff / rho) * ((u[i+1][j] - 2 * u[i][j] + u[i-1][j]) / (dx * dx) + (u[i][j+1] - 2 * u[i][j] + u[i][j-1]) / ( dy * dy));

            // u pressure gradient
            double u_pres_grad = (1 / rho) * (p[i+1][j] - p[i][j]) / dx;

            u_star[i][j] = (u_diff - u_conv - u_pres_grad) * dt + u[i][j];
        }
    }

    // solving for y-momentum equation
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y - 1; j++) {
            // average values of u, v in v-grid cell
            double v_avg_n = 0.5 * (v[i][j] + v[i][j+1]);
            double v_avg_s = 0.5 * (v[i][j] + v[i][j-1]);
            double v_avg_e = 0.5 * (v[i][j] + v[i+1][j]);
            double v_avg_w = 0.5 * (v[i][j] + v[i-1][j]);
            double u_avg_w = 0.5 * (u[i-1][j] + u[i-1][j+1]);
            double u_avg_e = 0.5 * (u[i][j] + u[i][j+1]);

            double mu_eff = mu + 0.5 * (mu_t[i][j] + mu_t[i][j+1]);

            // v convective term
            double v_conv = (u_avg_e * v_avg_e - u_avg_w * v_avg_w) / dx + (v_avg_n * v_avg_n - v_avg_s * v_avg_s) / dy;

            // v diffusion term
            double v_diff = (mu_eff / rho) * ((v[i+1][j] - 2 * v[i][j] + v[i-1][j]) / (dx * dx) + (v[i][j+1] - 2 * v[i][j] + v[i][j-1]) / (dy * dy));

            // v pressure gradient
            double v_pres_grad = (1 / rho) * (p[i][j+1] - p[i][j]) / dy;

            v_star[i][j] = (v_diff - v_conv - v_pres_grad) * dt + v[i][j];
        }
    }
}

void pressureCorrection(const vector<vector<double>>& u_star, const vector<vector<double>>& v_star, vector<vector<double>>& p_prime) {
    for (auto& row : p_prime) {
        fill(row.begin(), row.end(), 0);
    }

    double p_tol = 1e-5;
    int p_iter = 100;

    double a = 2 * ((dt / (dx * dx)) + (dt / (dy * dy)));
    double b = -1 * (dt / (dx * dx));
    double c = -1 * (dt / (dy * dy));

    vector<vector<double>> d(x + 1, vector<double>(y + 1, 0));
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y; j++) {
            d[i][j] = ((1 / dx) * (rho * u_star[i][j] - rho * u_star[i-1][j])) + ((1 / dy) * (rho * v_star[i][j] - rho * v_star[i][j-1]));
        }
    }

    for (int iter = 0; iter < p_iter; iter++) {
        double max_err = 0;
        for (int i = 1; i < x; i++) {
            for (int j = 1; j < y; j++) {
                double p_prime_old = p_prime[i][j];
                p_prime[i][j] = (-d[i][j] - (b * p_prime[i+1][j]) - (b * p_prime[i-1][j]) - (c * p_prime[i][j+1]) - (c * p_prime[i][j-1])) / a;
                max_err = max(max_err, abs(p_prime[i][j] - p_prime_old));
            }
        }

        for (int i = 1; i < x; i++) {
            p_prime[i][0] = p_prime[i][1];
            p_prime[i][y] = p_prime[i][y-1];
        }

        for (int j = 1; j < y; j++) {
            p_prime[0][j] = p_prime[1][j];
            p_prime[x][j] = p_prime[x-1][j];
        }

        if (max_err < p_tol) {
            break;
        }
    } 
}

void correctFields(vector<vector<double>>& u, vector<vector<double>>& v, vector<vector<double>>& p, const vector<vector<double>>& u_star, const vector<vector<double>>& v_star, const vector<vector<double>>& p_prime) {
    for (int i = 1; i < x - 1; i++) {
        for (int j = 1; j < y; j++) {
            u[i][j] = (1 - alpha_u) * u[i][j] + alpha_u * (u_star[i][j] - (dt / rho) * (p_prime[i+1][j] - p_prime[i][j]) / dx);
        }
    }

    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y - 1; j++) {
            v[i][j] = (1 - alpha_v) * v[i][j] + alpha_v * (v_star[i][j] - (dt / rho) * (p_prime[i][j+1] - p_prime[i][j]) / dy);
        }
    }

    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y; j++) {
            p[i][j] += alpha_p * p_prime[i][j];
        }
    }
}

void solveKEpsilon(const vector<vector<double>>& u, const vector<vector<double>>& v, vector<vector<double>>& k, vector<vector<double>>& epsilon, const vector<vector<double>>& mu_t) {
    vector<vector<double>> Pk(x + 1, vector<double>(y + 1, 0));
    vector<vector<double>> k_star = k;
    vector<vector<double>> epsilon_star = epsilon;

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
            double u_e = 0.5 * (u[i][j] + u[i][j+1]);     
            double u_w = 0.5 * (u[i-1][j] + u[i-1][j+1]); 
            double v_n = 0.5 * (v[i][j+1] + v[i+1][j+1]);     
            double v_s = 0.5 * (v[i][j] + v[i+1][j]);     

            double k_conv = (u_e * (u_e > 0 ? k[i][j] : k[i + 1][j]) - u_w * (u_w > 0 ? k[i - 1][j] : k[i][j])) / dx + (v_n * (v_n > 0 ? k[i][j] : k[i][j + 1]) - v_s * (v_s > 0 ?k[i][j - 1] : k[i][j])) / dy;

            double mu_eff_k = mu + mu_t[i][j] / sigma_k;
            double k_diff = (mu_eff_k / rho) * ((k[i+1][j] - 2 * k[i][j] + k[i-1][j]) / (dx * dx) + (k[i][j+1] - 2 * k[i][j] + k[i][j-1]) / (dy * dy));

            double k_source = (Pk[i][j] - rho * epsilon[i][j]) / rho;

            k_star[i][j] = k[i][j] + dt * (k_diff - k_conv + k_source);
        }
    }

    // solving for epsilon
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y; j++) {
            double u_e = 0.5 * (u[i][j] + u[i][j+1]);     
            double u_w = 0.5 * (u[i-1][j] + u[i-1][j+1]); 
            double v_n = 0.5 * (v[i][j+1] + v[i+1][j+1]);     
            double v_s = 0.5 * (v[i][j] + v[i+1][j]);   

            double epsilon_conv = (u_e * (u_e > 0 ? epsilon[i][j] : epsilon[i + 1][j]) - u_w * (u_w > 0 ? epsilon[i - 1][j] : epsilon[i][j])) / dx + (v_n * (v_n > 0 ? epsilon[i][j] : epsilon[i][j + 1]) - v_s * (v_s > 0 ? epsilon[i][j - 1] : epsilon[i][j])) / dy;

            double mu_eff_eps = mu + mu_t[i][j] / sigma_epsilon;
            double epsilon_diff = (mu_eff_eps / rho) * ((epsilon[i+1][j] - 2 * epsilon[i][j] + epsilon[i-1][j]) / (dx * dx) + (epsilon[i][j+1] - 2 * epsilon[i][j] + epsilon[i][j-1]) / (dy * dy));

            double k_val = max(k[i][j], 1e-10);
            double epsilon_source = (1 / rho) * ((C_1epsilon * epsilon[i][j] * Pk[i][j] / k_val) -  (C_2epsilon * rho * pow(epsilon[i][j], 2) / k_val));

            epsilon_star[i][j] = epsilon[i][j] + dt * (epsilon_diff - epsilon_conv + epsilon_source);
        }
    }

    // updating k and epsilon with under-relaxation
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y; j++) {
            k[i][j] = (1 - alpha_k) * k[i][j] + alpha_k * k_star[i][j];
            epsilon[i][j] = (1 - alpha_epsilon) * epsilon[i][j] + alpha_epsilon * epsilon_star[i][j];

            // ensuring it is not zero
            k[i][j] = max(k[i][j], 1e-10);
            epsilon[i][j] = max(epsilon[i][j], 1e-10);
        }
    }
}

void writeResults(const vector<vector<double>>& u, const vector<vector<double>>& v, const vector<vector<double>>& p, const vector<vector<double>>& k, const vector<vector<double>>& epsilon, const vector<vector<double>>& mu_t) {
    ofstream outfile("results.csv");
    if (!outfile.is_open()) {
        cerr << "Error: Could not open results.csv for writing." << endl;
        return;
    }

    outfile << "x_coord,y_coord,U,V,P,k,epsilon,mu_t\n";

    for (int j = 0; j < y; ++j) {
        for (int i = 0; i < x; ++i) {
            double x_pos = i * dx;
            double y_pos = j * dy;
            
            // Interpolate velocities to cell center
            double u_center = 0.5 * (u[i][j] + u[i][j + 1]);
            double v_center = 0.5 * (v[i][j] + v[i + 1][j]);
            double p_center = 0.25 * (p[i][j] + p[i+1][j] + p[i][j+1] + p[i+1][j+1]);
            double k_center = 0.25 * (k[i][j] + k[i+1][j] + k[i][j+1] + k[i+1][j+1]);
            double epsilon_center = 0.25 * (epsilon[i][j] + epsilon[i+1][j] + epsilon[i][j+1] + epsilon[i+1][j+1]);
            double mu_t_center = 0.25 * (mu_t[i][j] + mu_t[i+1][j] + mu_t[i][j+1] + mu_t[i+1][j+1]);

            outfile << x_pos << "," << y_pos << "," << u_center << "," << v_center << "," << p_center << "," << k_center << "," << epsilon_center << "," << mu_t_center << "\n";
        }
    }

    outfile.close();
    cout << "\nResults written to results.csv" << endl;
}

int main() {
    // Initializing grids for different flow fields
    vector<vector<double>> u(x, vector<double>(y + 1, 0));
    vector<vector<double>> v(x + 1, vector<double>(y, 0));
    vector<vector<double>> p(x + 1, vector<double>(y + 1, 0));
    vector<vector<double>> u_star(x, vector<double>(y + 1, 0));
    vector<vector<double>> v_star(x + 1, vector<double>(y, 0));
    vector<vector<double>> p_prime(x + 1, vector<double>(y + 1, 0));
    vector<vector<double>> k(x + 1, vector<double>(y + 1, 1e-6));
    vector<vector<double>> epsilon(x + 1, vector<double>(y + 1, 1e-6));
    vector<vector<double>> mu_t(x + 1, vector<double>(y + 1, 0));

    setBoundaryConditions(u, v, k, epsilon);

    for (int iter = 0; iter < max_iter; iter++) {

        calculateTurbulentViscosity(k, epsilon, mu_t);

        solveMomentum(u, v, p, u_star, v_star, mu_t);

        pressureCorrection(u_star, v_star, p_prime);

        correctFields(u, v, p, u_star, v_star, p_prime);

        solveKEpsilon(u, v, k, epsilon, mu_t);

        setBoundaryConditions(u, v, k, epsilon);
        
        if (iter % 100 == 0) {
            cout << "Iteration: " << iter <<  "\r" << flush;
        }
    }
    
    writeResults(u, v, p, k, epsilon, mu_t);

    return 0;
}