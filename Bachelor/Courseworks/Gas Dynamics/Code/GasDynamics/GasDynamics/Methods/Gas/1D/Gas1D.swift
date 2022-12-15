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
        
        static func *(l: [Double], v: V) -> Double {
            return l[0] * v.density + l[1] * v.speed + l[2] * v.pressure
        }
        
        static func *(v: V, r: [Double]) -> Double {
            return r * v
        }
        
        static func +(v1: V, v2: V) -> V {
            return V(density: v1.density + v2.density, speed: v1.speed + v2.speed, pressure: v1.pressure + v2.pressure)
        }
        
        static func /(v: V, x: Double) -> V {
            return V(density: v.density / x, speed: v.speed / x, pressure: v.pressure / x)
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
        
        public init(physical v: V, gamma: Double = Constants.gamma) {
            self.init(density: v.density, speed: v.speed, pressure: v.pressure, gamma: gamma)
        }
        
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
        
        solution[t] = space.nodes().reduce(into: Mesh()) { mesh, _x in
            let xL = BoundaryValue(value: _x - space.halfed, side: .right)
            let x  = BoundaryValue(value: _x, side: .middle)
            let xR = BoundaryValue(value: _x + space.halfed, side: .left)
            let vL = v0(xL.value)
            let vR = v0(xR.value)
            let v  = (vL + vR) / 2.0

            mesh[xL] = vL
            mesh[x]  = v
            mesh[xR] = vR
        }
        
        while t <= T {
            guard let meshP = solution[t] else { assertionFailure("Time iteration error"); return }
            
            let tau = self.tau(average: Set(meshP.filter{ $0.key.side == .middle }.values))
            let tP = t
            t += tau
            solution[t] = [:]
            
            space.nodes().forEach { node in
                let x = BoundaryValue(value: node, side: .middle)
                let xL = BoundaryValue(value: node - space.halfed, side: .right)
                let xR = BoundaryValue(value: node + space.halfed, side: .left)
                
                
            }
        }
    }
    
    // MARK: - Private methods
    
    private func averageFlow(FL: F, FR: F, vStar: V, vR: V, vL: V) {
        
    }
    
    private func tau(average vSet: Set<V>) -> Double {
        let maxLamba = vSet.map { abs($0.speed + Constants.speedOfSound(physical: $0, gamma: gamma)) }.max()!
        
        return Constants.sigmaGas * space.step / maxLamba
    }
        
    private func Nv(
        t: Double,
        x: Double,
        xL: BoundaryValue,
        xM: BoundaryValue,
        xR: BoundaryValue
    ) -> V {
        guard
            let vL = solution[t]?[xL],
            let vM  = solution[t]?[xM],
            let vR = solution[t]?[xR]
        else { return .zero }
        
        return V(
            density: Nv(x: x, xL: xL.value, vL: vL.density, vM: vM.density, vR: vR.density),
            speed: Nv(x: x, xL: xL.value, vL: vL.speed, vM: vM.speed, vR: vR.speed),
            pressure: Nv(x: x, xL: xL.value, vL: vL.pressure, vM: vM.pressure, vR: vR.pressure)
        )
    }
    
    private func Nv(
        x: Double,
        xL: Double,
        vL: Double,
        vM: Double,
        vR: Double
    ) -> Double {
        let xi = (x - xL) / space.step
        let delta = vR - vL
        let sixth = 6.0 * (vM - 0.5 * (vR + vL))
        
        return vL + xi * (delta + sixth * (1.0 - xi))
    }
    
    private func shifted(x: Double, lambda: Double, tau: Double) -> Double {
        return x - lambda * tau
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
        
        static func R(
            density: Double,
            speed: Double,
            pressure: Double,
            gamma: Double = Constants.gamma
        ) -> [[Double]] {
            let c = Constants.speedOfSound(density: density, pressure: pressure, gamma: gamma)
            
            return [
                [1, -c / density, c * c],
                [1, 0, 0],
                [1, c / density, c * c]
            ]
        }
        
        static func R(physical v: V, gamma: Double = Constants.gamma) -> [[Double]] {
            return R(density: v.density, speed: v.speed, pressure: v.pressure, gamma: gamma)
        }
        
        static func L(
            density: Double,
            speed: Double,
            pressure: Double,
            gamma: Double = Constants.gamma
        ) -> [[Double]] {
            let c = Constants.speedOfSound(density: density, pressure: pressure, gamma: gamma)
            
            return [
                [0, -density / (2 * c), 1 / (2 * c * c)],
                [1, 0, -1 / (c * c)],
                [0, density / (2 * c), 1 / (2 * c * c)]
            ]
        }
        
        static func L(physical v: V, gamma: Double = Constants.gamma) -> [[Double]] {
            return L(density: v.density, speed: v.speed, pressure: v.pressure, gamma: gamma)
        }
    }
}
