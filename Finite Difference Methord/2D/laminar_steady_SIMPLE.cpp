#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

// number of grids
const int x = 129;
const int y = 129;

// dimension of cavity
const double lx = 1;
const double ly = 1;

const double dx = lx / (x - 1); //  x spaceing between grids
const double dy = ly / (y - 1); // y spacing between grids
const double dt = 0.0005; // time step
const double alpha_p = 0.3; // relaxation factor for pressure
const double alpha_u = 0.5;
const double alpha_v = 0.5;
const int max_iter = 20000;

// double max_err = 1;

const double rho = 1; // density 
const double mu = 0.01; // dynamic viscosity
const double u_lid = 1; // top wall velocity

void setBoundaryConditions(vector<vector<double>>& u, vector<vector<double>>& v) {
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
}

void solveMomentum(const vector<vector<double>>& u, const vector<vector<double>>& v, const vector<vector<double>>& p, vector<vector<double>>& u_star, vector<vector<double>>& v_star) {

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

            // u convective term
            double u_conv = (u_avg_e * u_avg_e - u_avg_w * u_avg_w) / dx + (u_avg_n * v_avg_n - u_avg_s * v_avg_s) / dy;

            // u diffusion term
            double u_diff = (mu / rho) * ((u[i+1][j] - 2 * u[i][j] + u[i-1][j]) / (dx * dx) + (u[i][j+1] - 2 * u[i][j] + u[i][j-1]) / ( dy * dy));

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

            // v convective term
            double v_conv = (u_avg_e * v_avg_e - u_avg_w * v_avg_w) / dx + (v_avg_n * v_avg_n - v_avg_s * v_avg_s) / dy;

            // v diffusion term
            double v_diff = (mu / rho) * ((v[i+1][j] - 2 * v[i][j] + v[i-1][j]) / (dx * dx) + (v[i][j+1] - 2 * v[i][j] + v[i][j-1]) / (dy * dy));

            // v pressure gradient
            double v_pres_grad = (1 / rho) * (p[i][j+1] - p[i][j]) / dy;

            v_star[i][j] = (v_diff - v_conv - v_pres_grad) * dt + v[i][j];
        }
    }
}

void pressureCorrection(const vector<vector<double>>& u_star, const vector<vector<double>>& v_star, vector<vector<double>>& p_prime) {
    // setting all values of  p_prime to zero
    for (auto& row : p_prime) {
        fill(row.begin(), row.end(), 0);
    }

    double p_tol = 1e-5;
    int p_iter = 100;

    double a = 2 * ((dt / (dx * dx)) + (dt / (dy * dy)));
    double b = -1 * (dt / (dx * dx));
    double c = -1 * (dt / (dy * dy));

    // calculating mass source for every point
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

        // applying boundary conditions
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
    
    // u field
    for (int i = 1; i < x - 1; i++) {
        for (int j = 1; j < y; j++) {
            u[i][j] = (u_star[i][j] - (dt / rho) * (p_prime[i+1][j] - p_prime[i][j]) / dx) * alpha_u + (1 - alpha_u) * u[i][j];
        }
    }

    // v field
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y - 1; j++) {
            v[i][j] = (v_star[i][j] - (dt / rho) * (p_prime[i][j+1] - p_prime[i][j]) / dy) * alpha_v + (1 - alpha_v) * v[i][j];
        }
    }

    // p field
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y; j++) {
            p[i][j] += alpha_p * p_prime[i][j];
        }
    }
}

void writeResults(const vector<vector<double>>& u, const vector<vector<double>>& v, const vector<vector<double>>& p) {
    ofstream outfile("results.csv");
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

            // Interpolate velocities and pressure from faces to cell center
            double u_center = 0.5 * (u[i][j] + u[i][j+1]);
            double v_center = 0.5 * (v[i][j] + v[i+1][j]);
            double p_center = 0.25 * (p[i][j] + p[i+1][j] + p[i][j+1] + p[i+1][j+1]);

            outfile << i+1 << "," << j+1 << "," << x_pos << "," << y_pos << "," << u_center << "," << v_center << "," << p_center << "\n";
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

    setBoundaryConditions(u, v);

    double max_err = 1e-7;
    int iter = 0;
    while (max_err > 1e-8) {
    //for (int iter = 0; iter < max_iter; iter++) {
        max_err = 0;
        u_old = u;
        solveMomentum(u, v, p, u_star, v_star);

        pressureCorrection(u_star, v_star, p_prime);

        correctFields(u, v, p, u_star, v_star, p_prime);

        setBoundaryConditions(u, v);
        for (int j = 0; j <= y; j++) {
            max_err = max(max_err, abs(u_old[64][j] - u[64][j]));
        }
        if (iter % 100 == 0) {
            cout << "Iteration: " << iter << " " << max_err <<  "\r" << flush;
        }
        iter++;    
    }

    writeResults(u, v, p);

    return 0;
}