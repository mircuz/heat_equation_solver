//
//  support_heat.cpp
//  heat_transfer_eq
//
//  Created by Mirco Meazzo on 14/10/2019.
//  Copyright © 2019 Mirco Meazzo. All rights reserved.
//

#include "support_heat.hpp"
#include <iostream>
#include <fstream>

void grid::init_grid() {
    x_pos[0] = 0;
    for (int i=1; i<XDIM; i++) {
        x_pos[i]= x_pos[i-1] + stepx;
    }
    
    y_pos[0] = 0;
    for (int j=1; j<YDIM; j++) {
        y_pos[j] = y_pos[j-1] + stepy;
    }
}


grid::grid( int len_x_, int len_y_, float step_x_=1, float step_y_=1) :
lenx (len_x_), leny (len_y_), stepx (step_x_), stepy (step_y_)
{
    init_grid();
    printf( "Grid initialized\n");
};


grid::~grid() {
    delete [] xn;
    delete [] x;
}

double grid::compute_dx(int i) {
    return (this->x_pos[i+1] - this->x_pos[i-1]);
}

double grid::compute_dy(int i) {
    return (this->y_pos[i+1] - this->y_pos[i-1]);
}


void grid::visualize_values() {
    using namespace std;
    for (int j=0; j<YDIM; j++) {
        for (int i=0; i<XDIM; i++) {
           printf("%g\t", x[i+XDIM*j]);
        }
        printf("\n");
    }
}




thermal_domain::thermal_domain(int x_, int y_) : grid(x_,y_,x_/XDIM,y_/YDIM)
{
};


void thermal_domain::set_conditions(double temp) {
    
    for (int j=0; j<YDIM; j++) {
         for (int i=0; i<XDIM; i++) {
             x [j*XDIM+i]=temp;
             xn[j*XDIM+i]=temp;
         }
     }
}

double _temp, _time_of_contact;
int _x_coord, _y_coord;
void thermal_domain::set_source(int x_coord, int y_coord, double temp, double time_of_contact) {
    _temp = temp;
    _time_of_contact = time_of_contact;
    _x_coord = x_coord;
    _y_coord = y_coord;
    set_source_backend();
}

void thermal_domain::set_source_backend() {
    x[_y_coord*YDIM+_x_coord] = _temp;
}

void thermal_domain::set_alpha() {
    std::cout << "Insert Thermal Conductivity: " << std::flush;
    double alpha, k, d, cp;
    std::cin >> k;
    std::cout << "Insert Mass Density: " << std::flush;
    std::cin >> d;
    std::cout << "Insert Specific Heat Capacity: " << std::flush;
    std::cin >> cp;
    alpha = k/(d*cp);
    std::cout << "alpha="<< alpha << std::endl;
    this->a = alpha;
}


void thermal_domain::solve_pde(double time) {
    double t=0;
    double dx,dy,dt, C,Clim=0.5;
    dx = this->compute_dx(int(XDIM/2));     // Can be useful to define dx for every
    dy = this->compute_dy(int(YDIM/2));     // grid point in case of varying mesh
    dt = (Clim/this->a) * 1/((1/(dx))+(1/(dy)));
    C  = this->a*dt/dx + this->a*dt/dy;

    std::ofstream results;
    results.open("thermal_results.dat",std::ios::trunc);
    if (results.is_open()) {
        results << "##### DOMAIN DIMENSION #####\n";
        results << XDIM << "\n";
        results << YDIM << "\n";
        results << "##### END OF HEADER #####\n";
        
        while (t < time) {
            std::cout << "Time " << t << ", dt " << dt << ", Courant: " << C << std::endl;
            
    //      First row
    //      i=0,j=0
            int i=0,j=0;
            xn[i*YDIM+j] =   x[i*YDIM+j] + this->a*dt * (
                            (x[i*YDIM+j+2] -2*x[i*YDIM+j+1] + x[i*YDIM+j])/(dx*dx) +
                            (-2*x[(i+1)*YDIM+j] +x[i*YDIM+j] + x[(i+2)*YDIM+j])/(dy*dy) );
    //      i=0
            i=0;
#pragma omp parallel for
            for (int j=1; j<YDIM-1; j++) {
                xn[i*YDIM+j] =   x[i*YDIM+j] + this->a*dt * (
                                (x[i*YDIM+j+1] -2*x[i*YDIM+j] + x[i*YDIM+j-1])/(dx*dx) +
                                (-2*x[(i+1)*YDIM+j] +x[i*YDIM+j] + x[(i+2)*YDIM+j])/(dy*dy) );
            }
    //      i=0,j=YDIM-1
            i=0; j=YDIM-1;
            xn[i*YDIM+j] =   x[i*YDIM+j] + this->a*dt * (
                            (x[i*YDIM+j-2] -2*x[i*YDIM+j-1] + x[i*YDIM+j])/(dx*dx) +
                            (-2*x[(i+1)*YDIM+j] +x[i*YDIM+j] + x[(i+2)*YDIM+j])/(dy*dy) );
          
    //      i=XDIM-1,j=0
            i=XDIM-1; j=0;
            xn[i*YDIM+j] =   x[i*YDIM+j] + this->a*dt * (
                            (x[i*YDIM+j+2] -2*x[i*YDIM+j+1] + x[i*YDIM+j])/(dx*dx) +
                            (-2*x[(i-1)*YDIM+j] +x[i*YDIM+j] + x[(i-2)*YDIM+j])/(dy*dy) );
    //      i=XDIM-1
            i=XDIM-1;
#pragma omp parallel for
            for (int j=1; j<YDIM-1; j++) {
                xn[i*YDIM+j] =   x[i*YDIM+j] + this->a*dt * (
                                (x[i*YDIM+j+1] -2*x[i*YDIM+j] + x[i*YDIM+j-1])/(dx*dx) +
                                (-2*x[(i-1)*YDIM+j] +x[i*YDIM+j] + x[(i-2)*YDIM+j])/(dy*dy) );
            }
    //      i=XDIM-1,j=YDIM-1
            i=XDIM-1; j=YDIM-1;
            xn[i*YDIM+j] =   x[i*YDIM+j] + this->a*dt * (
                            (x[i*YDIM+j] -2*x[i*YDIM+j-1] + x[i*YDIM+j-2])/(dx*dx) +
                            (-2*x[(i-1)*YDIM+j] +x[i*YDIM+j] + x[(i-2)*YDIM+j])/(dy*dy) );
    //      Column 0
            j=0;
            for (int i=1; i<XDIM-1; i++) {
                xn[i*YDIM+j] =   x[i*YDIM+j] + this->a*dt * (
                                (x[i*YDIM+j+2] -2*x[i*YDIM+j+1] + x[i*YDIM+j])/(dx*dx) +
                                (-2*x[(i)*YDIM+j] +x[(i+1)*YDIM+j] + x[(i-1)*YDIM+j])/(dy*dy) );
            }
   
     //      Column YDIM-1
             j=YDIM-1;
             for (int i=1; i<XDIM-1; i++) {
                 xn[i*YDIM+j] =   x[i*YDIM+j] + this->a*dt * (
                                 (x[i*YDIM+j] -2*x[i*YDIM+j-1] + x[i*YDIM+j-2])/(dx*dx) +
                                 (-2*x[(i)*YDIM+j] +x[(i+1)*YDIM+j] + x[(i-1)*YDIM+j])/(dy*dy) );
             }
    
            
    //      Inner points
#pragma omp parallel for
            for (int i=1; i<XDIM-1; i++) {
                for (int j=1; j<YDIM-1; j++) {
                    xn[i*YDIM+j] =   x[i*YDIM+j] + this->a*dt * (
                                    (x[i*YDIM+j+1] -2*x[i*YDIM+j] + x[i*YDIM+j-1])/(dx*dx) +
                                    (x[(i+1)*YDIM+j] -2*x[i*YDIM+j] + x[(i-1)*YDIM+j])/(dy*dy) );
                }
            }
#pragma omp parallel for
            for (unsigned int i=0; i<XDIM; i++) {
                for (unsigned int j=0; j<YDIM; j++) {
                    x[i*YDIM+j] = xn[i*YDIM+j];
                }
            }
            if (t < _time_of_contact ) set_source_backend();
            t = t+dt;
            
//          Write data
            for (unsigned int i=0; i<XDIM; i++) {
                for (unsigned int j=0; j<YDIM; j++) {
                    results << x[i*YDIM+j] << "\t";
                }
                results << "\n";
            }
            results << "\n";
        }
        results.close();
    }
    else {
        std::cout << "Error opening file\nSimulation Aborted..." << std::endl;
    }
}
