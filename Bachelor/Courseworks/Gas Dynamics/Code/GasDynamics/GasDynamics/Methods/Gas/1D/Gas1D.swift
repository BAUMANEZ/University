//
//  Gas1D.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 31.10.2022.
//

import Foundation

final class Gas1D: JSONConvertableAlgorithm {
    
    // MARK: - Nested types
    
    typealias Mesh = [BoundaryValue: V]
    typealias Solution = [Time: Mesh]
    
    struct V: Hashable {
        
        // MARK: - Type properties
        
        static let zero = V(density: .zero, speed: .zero, pressure: .zero)
        
        // MARK: - Type methods
        
        /// e_i stands for eigen vector value (either left eigen vector or right eigen vector)
        static func makeFrom(alpha: Double, e1: Double, e2: Double, e3: Double) -> V {
            return Gas1D.V(density: alpha * e1, speed: alpha * e2, pressure: alpha * e3)
        }
        
        // MARK: - Internal properties
        
        let density : Double
        let speed   : Double
        let pressure: Double
    }
    
    struct U: Hashable {
        
        // MARK: - Type properties
        
        static let zero = U(physical: .zero)
        
        // MARK: - Internal properties
        
        /// - rho
        let u1: Double
        
        /// - rho * u
        let u2: Double
        
        /// - E
        let u3 : Double
        
        public init(physical v: V, gamma: Double = Constants.gamma) {
            self.init(density: v.density, speed: v.speed, pressure: v.pressure, gamma: gamma)
        }
        
        public init(density: Double, speed: Double, pressure: Double, gamma: Double = Constants.gamma) {
            u1 = density
            u2 = density * speed
            u3 = Formulas.energy(density: density, speed: speed, pressure: pressure, gamma: gamma)
        }
    }
    
    struct F {
        
        // MARK: - Internal properties
        
        /// - rho * u
        let f1: Double
        
        /// - rho * u^2 + p
        let f2: Double
        
        /// - (E + p) * u
        let f3: Double
        
        public init(density: Double, speed: Double, pressure: Double, gamma: Double = Constants.gamma) {
            f1 = density*speed
            f2 = density*pow(speed, 2) + pressure
            f3 = (Formulas.energy(density: density, speed: speed, pressure: pressure, gamma: gamma) + pressure) * speed
        }
    }
    
    // MARK: - Internal properties
    
    private(set) var solution: Solution = [:]
    
    let space: Grid
    let T: Time
    let v0: (Double) -> V
    
    let profile: String?
    let gamma = Constants.gamma
    
    // MARK: - Initializers
    
    init(test: Test) {
        self.profile = test.description
        self.space = Grid(start: test.a, end: test.b, steps: test.N)
        self.T = test.T
        self.v0 = test.initial
        
        solve()
    }
    
    // MARK: - Internal methods
    
    func solve() {
        var t: Double = 0
        
        solution[t] = space.nodes().reduce(into: Mesh()) { mesh, x in
            mesh[BoundaryValue(value: x, side: .middle)] = v0(x)
        }
        
        while t <= T {
            guard let meshP = solution[t] else { assertionFailure("Time iteration error"); return }
            
            let tau = self.tau(average: Set(meshP.values))
            
        }
    }
    
    // MARK: - Private methods
    
    private func tau(average vSet: Set<V>) -> Double {
        let maxLamba = vSet.map { abs($0.speed + Constants.speedOfSound(physical: $0, gamma: gamma)) }.max()!
        
        return Constants.sigmaGas * space.step / maxLamba
    }
}

extension Gas1D {
    
    // MARK: - Nested types
    
    enum Formulas {
        
        // MARK: - Type methods
        
        static func energy(
            density: Double,
            speed: Double,
            pressure: Double,
            gamma: Double = Constants.gamma
        ) -> Double {
            return pressure / (gamma - 1) + 0.5 * density * pow(speed, 2)
        }
        
        static func energy(physical v: V, gamma: Double = Constants.gamma) -> Double {
            return energy(density: v.density, speed: v.speed, pressure: v.pressure, gamma: gamma)
        }
    }
    
    enum Eigen {
        
        // MARK: - Type methods
        
        static func Lambda(
            density: Double,
            speed: Double,
            pressure: Double,
            gamma: Double = Constants.gamma
        ) -> [Double] {
            let c = Constants.speedOfSound(density: density, pressure: pressure, gamma: gamma)
            
            return [speed - c, speed, speed + c]
        }
        
        static func Lambda(physical v: V, gamma: Double = Constants.gamma) -> [Double] {
            return Lambda(density: v.density, speed: v.speed, pressure: v.pressure, gamma: gamma)
        }
        
        static func R(density: Double, speed: Double, pressure: Double, gamma: Double = Constants.gamma) -> Matrix {
            let c = Constants.speedOfSound(density: density, pressure: pressure, gamma: gamma)
            var result = Matrix(n: 3, initial: 0.0)
            
            result[0, 0] = 1
            result[0, 1] = 1
            result[0, 2] = 1
            
            result[1, 0] = -c / density
            result[1, 1] = 0
            result[1, 2] = c / density
            
            result[2, 0] = c * c
            result[2, 1] = 0
            result[2, 2] = c * c
            
            return result
        }
        
        static func R(physical v: V, gamma: Double = Constants.gamma) -> Matrix {
            return R(density: v.density, speed: v.speed, pressure: v.pressure, gamma: gamma)
        }
        
        static func L(density: Double, speed: Double, pressure: Double, gamma: Double = Constants.gamma) -> Matrix {
            let c = Constants.speedOfSound(density: density, pressure: pressure, gamma: gamma)
            var result = Matrix(n: 3, initial: 0.0)
            
            result[0, 0] = 0
            result[0, 1] = -density / (2 * c)
            result[0, 2] = 1 / (2 * c * c)
            
            result[1, 0] = 1
            result[1, 1] = 0
            result[1, 2] = -1 / (c * c)
            
            result[2, 0] = 0
            result[2, 1] = density / (2 * c)
            result[2, 2] = 1 / (2 * c * c)
            
            return result
        }
        
        static func L(physical v: V, gamma: Double = Constants.gamma) -> Matrix {
            return L(density: v.density, speed: v.speed, pressure: v.pressure, gamma: gamma)
        }
    }
}
