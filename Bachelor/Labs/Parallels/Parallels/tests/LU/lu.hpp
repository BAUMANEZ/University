//
//  lu.hpp
//  Parallels
//
//  Created by Арсений Токарев on 12.10.2022.
//

#ifndef lu_hpp
#define lu_hpp

#include <iostream>
#include <time.h>
#include "../../matrix/core/matrix.hpp"
#include "../../matrix/algorithm/algorithm.hpp"
#include "../../matrix/norms/norms.hpp"

struct lu_test {
    enum method { block, non_block };
    
    static void run(method method, bool paralleled);
};

#endif /* lu_hpp */
