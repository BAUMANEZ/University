//
//  norms.hpp
//  Parallels
//
//  Created by Арсений Токарев on 12.10.2022.
//

#ifndef norms_hpp
#define norms_hpp

#include <iostream>
#include "../core/matrix.hpp"
#include "../algorithm/algorithm.hpp"

struct norms {
    static double frobenius(const matrix& lhs, const matrix& rhs);
    static double max_error(const matrix& lhs, const matrix& rhs);
};

#endif /* norms_hpp */
