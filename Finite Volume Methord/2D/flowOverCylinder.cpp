#include <vector>
#include <iostream>
#include <fstream>
#include <string>


using namespace std;

int x = 160; // number of grid points in radial direction (across diameter)
int y = 60; // number of grid points in film thickness

const double lx = 8.0; // length
const double ly = 3.0; // height

double dx = lx / (x - 1); // distance between each grid point in x direction
double dy = ly / (y - 1); // distance between each grid point in y direction

// relaxation factors
const double alpha_p = 0.3;
const double alpha_u = 0.5;
const double alpha_v = 0.5;

// double v_vel = 0.05;
// const double u_lid = 0; // velocity of top squeeze plate in downward direction
const double u_in = 0.5;
const double temperature = 299.904; // assuming isothermal flow (K)
const double R = 287.05; // gas constant (air)
const double atm = 101325; // atomspheric pressure
double mu = 0.00392; // dynamic viscosity (air)

const double T = 100; // total time of solver
const double dt = 0.1; // time step
const double t_loop = 50; // iterations per time step

// Cylinder Geometry
const double cx = lx / 4.5; // Center x-coordinate 
const double cy = ly / 2.0; // Center y-coordinate 
const double cr = 0.3;      // Cylinder radius
const double cr_sq = cr * cr;

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
                if (sys.aP[i][j] < 1e-12) continue;

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

void setBoundaryCondition(vector<vector<double>>& u, vector<vector<double>>& v) {
    // u for left and right wall
    for (int j = 0; j < y + 1; j++) {
        // inlet
        u[0][j] = u_in;
        u[x-1][j] = u[x-2][j];
    }

    // v for left and right wall
    for (int j = 0; j < y; j++) {
        v[0][j] = -v[1][j];
        v[x][j] = v[x - 1][j];
    }
    for (int j = 0; j < y; j++) v[0][j] = 0.0;

    for (int i = 0; i < x; i++) {
        u[i][0] = u[i][1];
        u[i][y] = u[i][y - 1]; 
    }
    
    // v for top and bottom wall
    for (int i = 0; i < x + 1; i++) {
        v[i][0] = 0;
        v[i][y - 1] = 0;
    }
}

void solveMomentum(vector<vector<double>>& u, vector<vector<double>>& v, const vector<vector<double>>& p, const vector<vector<double>>& rho, vector<vector<double>>& u_star, vector<vector<double>>& v_star, vector<vector<double>>& ap_u, vector<vector<double>>& ap_v, const vector<vector<double>>& u_old, const vector<vector<double>>& v_old, const vector<vector<double>>& rho_old) {

    static linearSystem U;
    static linearSystem V;

    U.resize(x, y+1);
    V.resize(x+1, y);

    // solve u momentum
    for (int i = 1; i < x-1; i++) {
        for (int j = 1; j < y; j++) {
            // convective coefficients
            double Fw = 0.5 * ((rho[i+1][j] + rho[i][j]) * 0.5 * u[i][j] + (rho[i-1][j] + rho[i][j]) * 0.5 * u[i-1][j]) * dy;
            double Fe = 0.5 * ((rho[i+2][j] + rho[i+1][j]) * 0.5 * u[i+1][j] + (rho[i][j] + rho[i+1][j]) * 0.5 * u[i][j]) * dy;
            double Fs = 0.5 * ((rho[i][j] + rho[i][j-1]) * 0.5 * v[i][j-1] + (rho[i+1][j] + rho[i+1][j-1]) * 0.5 * v[i+1][j-1]) * dx;
            double Fn = 0.5 * ((rho[i][j] + rho[i][j+1]) * 0.5 * v[i][j] + (rho[i+1][j] + rho[i+1][j+1]) * 0.5 * v[i+1][j]) * dx;

            // diffusion coefficients
            double De = mu * dy / dx;
            double Dw = mu * dy / dx;
            double Dn = mu * dx / dy;
            double Ds = mu * dx / dy;

            // hybrid scheme
            double aw = max({Fw, Dw + Fw * 0.5, 0.0});
            double ae = max({-Fe, De - Fe * 0.5, 0.0});
            double as = max({Fs, Ds + Fs * 0.5, 0.0});
            double an = max({-Fn, Dn - Fn * 0.5, 0.0});
            double ap0 = 0.5 * (rho_old[i][j] + rho_old[i+1][j])  * dx * dy / dt;
            double delta_f = Fe - Fw + Fn - Fs;
            double ap_u_unrelaxed = ae + aw + an + as + ap0 + delta_f;

            double pres_u = (p[i][j] - p[i + 1][j]) * dy;

            

            double x_pos = i * dx; 
            double y_pos = j * dy; 

            // Check distance to cylinder center
            double dist_sq = (x_pos - cx)*(x_pos - cx) + (y_pos - cy)*(y_pos - cy);

            if (dist_sq <= cr_sq) {
                // grid point inside obstacle
                U.aE[i][j] = 0.0;
                U.aW[i][j] = 0.0;
                U.aN[i][j] = 0.0;
                U.aS[i][j] = 0.0;
                U.aP[i][j] = 1.0e30;
                U.b[i][j] = 0.0;
                ap_u[i][j] = 1.0e30;
            } else {
                if (i == x-2) {
                    u_star[i+1][j] = u_star[i][j] * massFlowRateRatio(u_star, rho);
                }
                U.aE[i][j] = ae;
                U.aW[i][j] = aw;
                U.aN[i][j] = an;
                U.aS[i][j] = as;
                U.aP[i][j] = ap_u_unrelaxed / alpha_u;
                U.b[i][j] = pres_u + ap0 * u_old[i][j] + (1.0 - alpha_u) * U.aP[i][j] * u[i][j];
                ap_u[i][j] = U.aP[i][j];
            }       
        }
    }
    successiveOverRelaxation(u_star, U, 1, x-1, 1, y, alpha_u);

    // solve v momentum
    for (int i = 1; i < x; i++) {
        for (int j = 1; j < y-1; j++) {
            // convective coefficients
            double Fw = 0.5 * ((rho[i][j] + rho[i-1][j]) * 0.5 * u[i-1][j] + (rho[i][j+1] + rho[i-1][j+1]) * 0.5 * u[i-1][j+1]) * dy;
            double Fe = 0.5 * ((rho[i+1][j] + rho[i][j]) * 0.5 * u[i][j] + (rho[i+1][j+1] + rho[i][j+1]) * 0.5 * u[i][j+1]) * dy;
            double Fs = 0.5 * ((rho[i][j-1] + rho[i][j]) * 0.5 * v[i][j-1] + (rho[i][j+1] + rho[i][j]) * 0.5 * v[i][j]) * dx;
            double Fn = 0.5 * ((rho[i][j+2] + rho[i][j+1]) * 0.5 * v[i][j+1] + (rho[i+1][j] + rho[i][j]) * 0.5 * v[i][j]) * dx;

            // diffusion coefficients
            double De = mu * dy / dx;
            double Dw = mu * dy / dx;
            double Dn = mu * dx / dy;
            double Ds = mu * dx / dy;

            // hybrid scheme
            double aw = max({Fw, Dw + Fw * 0.5, 0.0});
            double ae = max({-Fe, De - Fe * 0.5, 0.0});
            double as = max({Fs, Ds + Fs * 0.5, 0.0});
            double an = max({-Fn, Dn - Fn * 0.5, 0.0});

            double delta_f = Fe - Fw + Fn - Fs;
            double ap0 = 0.5 * (rho_old[i][j] + rho_old[i][j+1])  * dx * dy / dt;
            double ap_v_unrelaxed = ae + aw + an + as + ap0 + delta_f;

            double pres_v = (p[i][j] - p[i][j + 1]) * dx;

            double x_pos = i * dx;
            double y_pos = j * dy; 

            double dist_sq = (x_pos - cx)*(x_pos - cx) + (y_pos - cy)*(y_pos - cy);

            if (dist_sq <= cr_sq) {
                V.aE[i][j] = 0.0;
                V.aW[i][j] = 0.0;
                V.aN[i][j] = 0.0;
                V.aS[i][j] = 0.0;
                V.aP[i][j] = 1.0e30;
                V.b[i][j] = 0.0;
                ap_v[i][j] = 1.0e30;
            } else {
                V.aE[i][j] = ae;
                V.aW[i][j] = aw;
                V.aN[i][j] = an;
                V.aS[i][j] = as;
                V.aP[i][j] = ap_v_unrelaxed / alpha_v;
                V.b[i][j] = pres_v + ap0 * v_old[i][j] + (1.0 - alpha_v) * V.aP[i][j] * v[i][j];
                ap_v[i][j] = V.aP[i][j];
            }
        }
        successiveOverRelaxation(v_star, V, 1, x, 1, y-1, alpha_v);
    }
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
    for (int i = 1; i < x; i++) {
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

void writeResults(const vector<vector<double>>& u, const vector<vector<double>>& v, const vector<vector<double>>& p, const vector<vector<double>>& rho, int iter=1) {
    string name = "results-" + to_string(iter) + ".csv";
    ofstream outfile(name);
    if (!outfile.is_open()) {
        cerr << "Error: Could not open results.csv for writing." << endl;
        return;
    }

    // Write CSV header
    outfile << "x,y,X,Y,Z,U,V,P,rho\n";

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

            outfile << i + 1 << "," << j + 1 << "," << x_pos << "," << y_pos << "," << 0.0 << "," << u_center << "," << v_center << "," << p_center << "," << rho_center << "\n";
        }
    }

    outfile.close();
    cout << "\nResults written to results.csv" << endl;
}

int main() {

    // initializing grids
    vector<vector<double>> u(x, vector<double>(y+1, 0.0));
    vector<vector<double>> u_old = u;
    vector<vector<double>> u_history = u;
    vector<vector<double>> u_star = u;
    vector<vector<double>> ap_u = u;

    vector<vector<double>> v(x+1, vector<double>(y, 0.0));
    vector<vector<double>> v_star = v;
    vector<vector<double>> v_old = v;
    vector<vector<double>> ap_v = v;

    vector<vector<double>> p(x+1, vector<double>(y+1, atm));
    vector<vector<double>> p_prime(x+1, vector<double>(y+1, 0.0));

    vector<vector<double>> rho(x+1, vector<double>(y+1, 1.177));
    vector<vector<double>> rho_old = rho;

    double max_err = 1e-7;

    int iter = 0;

    setBoundaryCondition(u, v);
    
    // outer loop
    for (double time = 0; time < T; time += dt) {
        u_old = u;
        v_old = v;
        rho_old = rho;

        cout << time << endl;
        // inner loop inside a time step
        for (int iter = 0; iter < t_loop; iter++) {
            max_err = 0;
            u_history = u;

            u_star = u;
            v_star = v;

            solveDensity(rho, p);
            
            solveMomentum(u, v, p, rho, u_star, v_star, ap_u, ap_v, u_old, v_old, rho_old);

            pressureCorrection(u_star, v_star, p_prime, ap_u, ap_v, rho, rho_old);

            correctFields(u, v, p, u_star, v_star, p_prime, ap_u, ap_v);

            setBoundaryCondition(u, v);

            for (int j = 0; j <= y; j++) {
                max_err = max(max_err, abs(u_history[70][j] - u[70][j]));
            }
            
            cout << "Iteration: " << iter << " " << max_err << '\r' << flush;
        }
        writeResults(u, v, p, rho, time);
    }

    //writeResults(u, v, p, rho, dy);

    return 0;
}
