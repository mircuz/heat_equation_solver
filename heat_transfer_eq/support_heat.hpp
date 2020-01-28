//
//  support_heat.hpp
//  heat_transfer_eq
//
//  Created by Mirco Meazzo on 14/10/2019.
//  Copyright © 2019 Mirco Meazzo. All rights reserved.
//

#ifndef support_heat_hpp
#define support_heat_hpp

#ifndef XDIM
    #define XDIM 10
#endif
#ifndef YDIM
    #define YDIM 10
#endif

#include <stdio.h>


class grid {
public:
    int lenx, leny;
    float stepx, stepy;
    void init_grid();
    void visualize_values();
    double *xn = new double[XDIM*YDIM];
    double *x = new double[XDIM*YDIM];
    double *x_pos = new double[XDIM];
    double *y_pos = new double[YDIM];
    double compute_dx(int i);
    double compute_dy(int i);
    grid( int len_x_, int len_y_, float step_x_, float step_y_);
    ~grid();
};


class thermal_domain: public grid {
public:
    int nx, ny;
    thermal_domain(int x_, int y_);
    void set_conditions(double temp);
    void set_source(int x_coord, int y_coord, double temp);
    double set_alpha();
    void solve_pde(double time, double alpha);
};


#endif /* support_heat_hpp */
