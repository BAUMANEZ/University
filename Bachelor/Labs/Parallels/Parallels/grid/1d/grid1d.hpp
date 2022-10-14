//
//  grid1d.hpp
//  Parallels
//
//  Created by Арсений Токарев on 14.10.2022.
//

#ifndef grid1d_hpp
#define grid1d_hpp

#include <iostream>

struct grid1d final {
private:
    double a;
    double b;
    double h;
    size_t n;
    
//MARK: - Getters
public:
    inline double start() const { return a; }
    inline double end()   const { return b; }
    inline size_t step()  const { return h; }
    inline size_t steps() const { return n; }
    
//MARK: - Init
public:
    grid1d(double start, double end, double step);
    
//MARK: - Setters
    void set_step(double h);
    void set_nodes(size_t n);
    void resize(double a, double b);
    
//MARK: - Operators
    double operator[](size_t i) const;
};

#endif /* grid1d_hpp */
