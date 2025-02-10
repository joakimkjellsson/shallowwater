/* 
   2D shallow water model 
   in C++ 
*/   

// iostream for print statement etc
#include <iostream>
// vector for arrays
#include <vector>
// standard math
#include <cmath>
// file streams
#include <fstream>

/* 
With this we can write vector instead of std::vector
However, it is bad practice to load in the entire 
namespace of libraries since they may have conflicting 
names, or I might write a function with a conflicting
name. If cmath includes a sin function, and I write 
a function called sin as well, it would confuse the
compiler. 
Better to leave this line commented out and always
write std::vector etc. It also helps a reader to 
understand where functions come from.
*/
// using namespace std;

// Define constants
const double g = 9.81;  // acceleration due to gravity (m/s^2)
const double dx = 1.0;  // grid spacing in x direction
const double dy = 1.0;  // grid spacing in y direction
const double dt = 0.1;  // time step (s)
const int nx = 100;     // number of grid points in x direction
const int ny = 100;     // number of grid points in y direction
const int nt = 200;     // number of time steps to simulate

/* 
2D array for water depth (h), and velocity components (u_x, u_y)
C++ does not have matrices, so we have arrays of arrays, e.g. 
1st element of h (i=1) is an array of size ny. 
We also initialise to 0.  
*/
vector<vector<double>> h(nx, vector<double>(ny, 0.0));  // water height
vector<vector<double>> u_x(nx, vector<double>(ny, 0.0)); // x-velocity
vector<vector<double>> u_y(nx, vector<double>(ny, 0.0)); // y-velocity

// Initialize the water depth (h) and velocity (u_x, u_y)
void initialize() {
    for (int i = 30; i < 70; i++) {
        for (int j = 30; j < 70; j++) {
            h[i][j] = 10.0;  // Initial wave in the middle of the grid
        }
    }
}

// Update the boundary conditions (simple reflective boundaries)
void applyBoundaryConditions() {
    for (int i = 0; i < nx; i++) {
        h[i][0] = h[i][1]; h[i][ny - 1] = h[i][ny - 2];
        u_x[i][0] = -u_x[i][1]; u_x[i][ny - 1] = -u_x[i][ny - 2];
        u_y[i][0] = -u_y[i][1]; u_y[i][ny - 1] = -u_y[i][ny - 2];
    }
    for (int j = 0; j < ny; j++) {
        h[0][j] = h[1][j]; h[nx - 1][j] = h[nx - 2][j];
        u_x[0][j] = -u_x[1][j]; u_x[nx - 1][j] = -u_x[nx - 2][j];
        u_y[0][j] = -u_y[1][j]; u_y[nx - 1][j] = -u_y[nx - 2][j];
    }
}

// Step function to update water depth and velocity based on shallow water equations
void step() {
    vector<vector<double>> h_new = h;
    vector<vector<double>> u_x_new = u_x;
    vector<vector<double>> u_y_new = u_y;

    // Loop over all grid points
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            // Update the water height (continuity equation)
            double dhdx = (h[i+1][j] - h[i-1][j]) / (2.0 * dx);
            double dhdy = (h[i][j+1] - h[i][j-1]) / (2.0 * dy);
            h_new[i][j] = h[i][j] - dt * (u_x[i][j] * dhdx + u_y[i][j] * dhdy);

            // Update the x-velocity (momentum equation)
            double dudx = (u_x[i+1][j] - u_x[i-1][j]) / (2.0 * dx);
            double dudy = (u_x[i][j+1] - u_x[i][j-1]) / (2.0 * dy);
            u_x_new[i][j] = u_x[i][j] - dt * (u_x[i][j] * dudx + u_y[i][j] * dudy + g * dhdx / h[i][j]);

            // Update the y-velocity (momentum equation)
            double dvdx = (u_y[i+1][j] - u_y[i-1][j]) / (2.0 * dx);
            double dvdy = (u_y[i][j+1] - u_y[i][j-1]) / (2.0 * dy);
            u_y_new[i][j] = u_y[i][j] - dt * (u_x[i][j] * dvdx + u_y[i][j] * dvdy + g * dhdy / h[i][j]);
        }
    }

    // Apply the boundary conditions again
    applyBoundaryConditions();

    // Copy updated values back into the original arrays
    h = h_new;
    u_x = u_x_new;
    u_y = u_y_new;
}

// Function to write the results to a file for visualization
void writeOutput(int step) {
    ofstream file("output_" + to_string(step) + ".dat");
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            file << h[i][j] << " ";
        }
        file << endl;
    }
    file.close();
}

int main() {
    // Initialize the system
    initialize();

    // Time stepping loop
    for (int t = 0; t < nt; t++) {
        cout << "Time step: " << t << endl;
        step();
        if (t % 10 == 0) {
            writeOutput(t);
        }
    }

    return 0;
}
