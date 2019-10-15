//
//  support_heat.cpp
//  heat_transfer_eq
//
//  Created by Mirco Meazzo on 14/10/2019.
//  Copyright © 2019 Mirco Meazzo. All rights reserved.
//

#include "support_heat.hpp"

void grid::init_grid() {
    x_pos[0] = 0;
    for (int i=1; i<XDIM; i++) {
        x_pos[i]= x_pos[i-1] + stepx;
    }
    
    y_pos[0] = 0;
    for (int j=1; j<YDIM; j++) {
        y_pos[j] = y_pos[j-1] +stepy;
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
