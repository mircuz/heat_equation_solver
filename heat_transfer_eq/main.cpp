//
//  main.cpp
//  heat_transfer_eq
//
//  Created by Mirco Meazzo on 14/10/2019.
//  Copyright © 2019 Mirco Meazzo. All rights reserved.
//

#include <iostream>
#include "support_heat.hpp"


int main(int argc, const char * argv[]) {
    
    thermal_domain domain(10,10);  // Bug here, remove 100,100 and use X/YDIM instead
    domain.set_conditions(273);
    domain.set_source(5, 5, 350, 1);
    domain.set_alpha();
    domain.solve_pde(3);
}
