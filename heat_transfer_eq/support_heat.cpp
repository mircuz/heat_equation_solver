//
//  support_heat.cpp
//  heat_transfer_eq
//
//  Created by Mirco Meazzo on 14/10/2019.
//  Copyright © 2019 Mirco Meazzo. All rights reserved.
//

#include "support_heat.hpp"
#include <iostream>

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


void thermal_domain::set_source(int x_coord, int y_coord, double temp) {
    x[y_coord*YDIM+x_coord] = temp;
}

double thermal_domain::set_alpha() {
    std::cout << "Insert Thermal Conductivity: " << std::flush;
    double alpha, k, d, cp;
    std::cin >> k;
    std::cout << "Insert Mass Density: " << std::flush;
    std::cin >> d;
    std::cout << "Insert Specific Heat Capacity: " << std::flush;
    std::cin >> cp;
    alpha = k/(d*cp);
    std::cout << "alpha="<< alpha << std::endl;
    return alpha;
}


void thermal_domain::solve_pde(double time, double alpha) {
    double t=0;
    double dx,dy,dt;
    
    dx = this->compute_dx(XDIM/2);     // Can be useful to define dx for every
    dy = this->compute_dy(YDIM/2);     // grid point in case of varying mesh
    dt = (1/alpha) * 1/((1/(dx*dx))+(1/(dy*dy))) *0.55;     //Safety factor
    
    while (t < time) {
        std::cout << "Time " << t << " dt " << dt << std::endl;
        
//      First row
        for (int j=0; j<YDIM; j++) {
            int i=0;
            xn[i*YDIM+j] =   x[i*YDIM+j] + alpha*dt * (
                            (x[i*YDIM+j+1] -2*x[i*YDIM+j] + x[i*YDIM+j-1])/(dx*dx) +
                            (-2*x[(i+1)*YDIM+j] +x[i*YDIM+j] + x[(i+2)*YDIM+j])/(dy*dy) );
        }
      
//      Last row
        for (int j=0; j<YDIM; j++) {
            int i=XDIM-1;
            xn[i*YDIM+j] =   x[i*YDIM+j] + alpha*dt * (
                            (x[i*YDIM+j+1] -2*x[i*YDIM+j] + x[i*YDIM+j-1])/(dx*dx) +
                            (-2*x[(i-1)*YDIM+j] +x[i*YDIM+j] + x[(i-2)*YDIM+j])/(dy*dy) );
        
//      Inner points
        for (int i=1; i<XDIM-1; i++) {
            for (int j=0; j<YDIM; j++) {

                xn[i*YDIM+j] =   x[i*YDIM+j] + alpha*dt * (
                                (x[i*YDIM+j+1] -2*x[i*YDIM+j] + x[i*YDIM+j-1])/(dx*dx) +
                                (x[(i+1)*YDIM+j] -2*x[i*YDIM+j] + x[(i-1)*YDIM+j])/(dy*dy) );

            }
        }
        t = t+dt;
        for (unsigned int i=0; i<XDIM; i++) {
            for (unsigned int j=0; j<YDIM; j++) {
                x[i*YDIM+j] = xn[i*YDIM+j];
            }
        }
        visualize_values();
    }
}
