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
    
    thermal_domain domain(10,10);
    domain.set_conditions(273);
    domain.set_source(5, 5, 350);
    domain.visualize_values();

}
