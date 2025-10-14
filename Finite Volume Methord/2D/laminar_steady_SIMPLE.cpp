#include <algorithm>
#include <vector>
#include <iostream>
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
const double alpha_p = 0.8; // relaxation factor for pressure
const double alpha_u = 0.7; // relaxation factor for u-velocity
const double alpha_v = 0.7; // relaxation factor for v-velocity
// const double max_tol = 1e-8;
const int max_iter = 100000;

// double max_err = 1;

const double rho = 1; // density
const double mu = 0.01; // dynamic viscosity
const double u_lid = 1; // top wall velocity

void setBoundaryCondition(vector<vector<double>>& u, vector<vector<double>>& v) {
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
}


void solveMomentum(const vector<vector<double>>& u, const vector<vector<double>>& v, vector<vector<double>>& u_star, vector<vector<double>>& v_star, vector<vector<double>>& ap_u, vector<vector<double>>& ap_v, const vector<vector<double>>& p) {
    // solve u momentum equation
    for (int i = 1; i < x - 1; i++) {
        for (int j = 1; j < y; j++) {
            // convective coefficients
            double Fe = rho * 0.5 * (u[i + 1][j] + u[i][j]) * dy;
            double Fw = rho * 0.5 * (u[i - 1][j] + u[i][j]) * dy;
            double Fn = rho * 0.5 * (v[i][j] + v[i + 1][j]) * dx;
            double Fs = rho * 0.5 * (v[i][j - 1] + v[i + 1][j - 1]) * dx;

            // diffusion coefficients
            double De = mu * dy / dx;
            double Dw = mu * dy / dx;
            double Dn = mu * dx / dy;
            double Ds = mu * dx / dy;

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
            double De = mu * dy / dx;
            double Dw = mu * dy / dx;
            double Dn = mu * dx / dy;
            double Ds = mu * dx / dy;

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

void writeResults(const vector<vector<double>>& u, const vector<vector<double>>& v, const vector<vector<double>>& p, int iter=1) {
    string name = "results-" + to_string(iter) + ".csv";
    ofstream outfile(name);
    if (!outfile.is_open()) {
        cerr << "Error: Could not open results.csv for writing." << endl;
        return;
    }

    // Write CSV header
    outfile << "x,y,X,Y,U,V,P\n";

    // Output data at cell centers
    for (int j = 0; j < y; ++j) {
        for (int i = 0; i < x; ++i) {
            double x_pos = i * dx;
            double y_pos = j * dy;

            // Interpolate velocities and pressure from faces/nodes to cell center
            double u_center = 0.5 * (u[i][j] + u[i][j + 1]);
            double v_center = 0.5 * (v[i][j] + v[i + 1][j]);
            double p_center = 0.25 * (p[i][j] + p[i + 1][j] + p[i][j + 1] + p[i + 1][j + 1]);

            outfile << i + 1 << "," << j + 1 << "," << x_pos << "," << y_pos << "," << u_center << "," << v_center << "," << p_center << "\n";
        }
    }

    outfile.close();
    cout << "\nResults written to results.csv" << endl;
}

int main() {
    
    // Initializing grids for different flow fields
    vector<vector<double>> u(x, vector<double>(y + 1, 0));
    vector<vector<double>> u_old(x, vector<double>(y + 1, 0));
    vector<vector<double>> v(x + 1, vector<double>(y, 0));
    vector<vector<double>> p(x + 1, vector<double>(y + 1, 0));
    vector<vector<double>> u_star(x, vector<double>(y + 1, 0));
    vector<vector<double>> v_star(x + 1, vector<double>(y, 0));
    vector<vector<double>> p_prime(x + 1, vector<double>(y + 1, 0));
    vector<vector<double>> ap_u(x, vector<double>(y + 1, 0));
    vector<vector<double>> ap_v(x + 1, vector<double>(y, 0));

    double max_err = 1e-7;
    int iter = 0;

    setBoundaryCondition(u, v);
    
    while (max_err > 1e-8) {
    // for (int iter = 0; iter < max_iter; iter++) {
        max_err = 0;
        u_old = u;
        
        solveMomentum(u, v, u_star, v_star, ap_u, ap_v, p);

        pressureCorrection(u_star, v_star, p_prime, ap_u, ap_v);

        correctFields(u, v, p, u_star, v_star, p_prime, ap_u, ap_v);

        setBoundaryCondition(u, v);

        for (int j = 0; j <= y; j++) {
            max_err = max(max_err, abs(u_old[64][j] - u[64][j]));
        }

        if (iter % 100 == 0) {
            cout << "Iteration: " << iter << " " << max_err << "\r" << flush;
            writeResults(u, v, p, iter);
        }
        iter++;
    }

    writeResults(u, v, p);

    return 0;
}