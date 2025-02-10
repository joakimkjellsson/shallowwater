/* 
*/ 
#include <iostream>
#include <vector>
#include <netcdf.h>

// Define parameters
const double g = 9.81; // Gravitational acceleration
const double H = 10.0; // Mean water depth
const double dx = 1.0; // Grid spacing in x-direction
const double dy = 1.0; // Grid spacing in y-direction
const double dt = 0.01; // Time step
const int nx = 100;    // Number of grid points in x-direction
const int ny = 100;    // Number of grid points in y-direction
const int nt = 100;    // Number of time steps

// Error handler for NetCDF
void check_nc(int status) {
    if (status != NC_NOERR) {
        std::cerr << "NetCDF Error: " << nc_strerror(status) << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Initialize the height field with a small perturbation
void initialize(std::vector<std::vector<double>> &h) {
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            h[j][i] = H + 0.1 * std::exp(-0.01 * ((i - nx / 2) * (i - nx / 2) + (j - ny / 2) * (j - ny / 2)));
        }
    }
}

// Write data to a NetCDF file
void write_netcdf(const std::string &filename, const std::vector<std::vector<double>> &h, int timestep) {
    int ncid, varid, dimids[3];
    int x_dim, y_dim, t_dim;
    size_t start[3] = {timestep, 0, 0};
    size_t count[3] = {1, ny, nx};

    if (timestep == 0) {
        // Create the file
        check_nc(nc_create(filename.c_str(), NC_CLOBBER, &ncid));

        // Define dimensions
        check_nc(nc_def_dim(ncid, "x", nx, &x_dim));
        check_nc(nc_def_dim(ncid, "y", ny, &y_dim));
        check_nc(nc_def_dim(ncid, "time", NC_UNLIMITED, &t_dim));

        dimids[0] = t_dim;
        dimids[1] = y_dim;
        dimids[2] = x_dim;

        // Define variable
        check_nc(nc_def_var(ncid, "h", NC_DOUBLE, 3, dimids, &varid));

        // End define mode
        check_nc(nc_enddef(ncid));
    } else {
        // Open the file
        check_nc(nc_open(filename.c_str(), NC_WRITE, &ncid));
        check_nc(nc_inq_varid(ncid, "h", &varid));
    }

    // Write data
    check_nc(nc_put_vara_double(ncid, varid, start, count, &h[0][0]));

    // Close the file
    check_nc(nc_close(ncid));
}

int main() {
    
    /*
    Grid is from index 0 to nx, i.e. not 1 as in Fortran
    */
    
    // Allocate fields
    std::vector<std::vector<double>> h(ny, std::vector<double>(nx, H));
    std::vector<std::vector<double>> u(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> v(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> h_new = h;
    std::vector<std::vector<double>> u_new = u;
    std::vector<std::vector<double>> v_new = v;

    initialize(h);

    std::string filename = "shallow_water.nc";

    // Time integration loop
    for (int t = 0; t < nt; ++t) {
        // Update u, v, h 
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                u_new[j][i] = u[j][i] - dt * g * (h[j][i + 1] - h[j][i]) / dx;
                v_new[j][i] = v[j][i] - dt * g * (h[j + 1][i] - h[j][i]) / dy;
                h_new[j][i] = h[j][i] - dt * H * ((u[j][i] - u[j][i - 1]) / dx + (v[j][i] - v[j - 1][i]) / dy);
            }
        }

        // Swap arrays
        u.swap(u_new);
        v.swap(v_new);
        h.swap(h_new);

        // Output to NetCDF
        write_netcdf(filename, h, t);
    }

    std::cout << "Simulation complete. Data written to " << filename << std::endl;
    return 0;
}
 
