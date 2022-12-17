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
        
        static func makeFrom(alpha: Double, eigen vector: [Double]) -> V {
            return Gas1D.V(density: alpha * vector[0], speed: alpha * vector[1], pressure: alpha * vector[2])
        }
        
        static func *(l: [Double], v: V) -> Double {
            return l[0] * v.density + l[1] * v.speed + l[2] * v.pressure
        }
        
        static func *(v: V, r: [Double]) -> Double {
            return r * v
        }
        
        static func *(v: V, x: Double) -> V {
            return V(density: v.density * x, speed: v.speed * x, pressure: v.pressure * x)
        }
        
        static func *(x: Double, v: V) -> V {
            return v * x
        }
        
        static func *=(v: inout V, x: Double) {
            v.density *= x
            v.speed *= x
            v.pressure *= x
        }
        
        static func *=(x: Double, v: inout V) {
            v *= x
        }
        
        static func +(v1: V, v2: V) -> V {
            return V(density: v1.density + v2.density, speed: v1.speed + v2.speed, pressure: v1.pressure + v2.pressure)
        }
        
        static func +=(v1: inout V, v2: V) {
            v1.density += v2.density
            v1.speed += v2.speed
            v1.pressure += v2.pressure
        }
        
        static func -(v1: V, v2: V) -> V {
            return V(density: v1.density - v2.density, speed: v1.speed - v2.speed, pressure: v1.pressure - v2.pressure)
        }
        
        // MARK: - Internal properties
        
        var density : Double
        var speed   : Double
        var pressure: Double
        
        public init(density: Double, speed: Double, pressure: Double) {
            self.density = density
            self.speed = speed
            self.pressure = pressure
        }
        
        public init(vector: [Double]) {
            self.density = vector[0]
            self.speed = vector[1]
            self.pressure = vector[2]
        }
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
    
    // MARK: - Type methods
    
    static func solve(for test: Test) {
        let problem = Gas1D(test: test)
        problem.solve()
        problem.export()
    }
    
    // MARK: - Internal properties
    
    private(set) var solution: Solution = [:]
    
    let space: Grid
    let T: Time
    let v0: (Double) -> V
    
    let profile: String?
    let gamma = Constants.gamma
    
    // MARK: - Private properties
    
    private lazy var metadata: Data? = {
        let json: [String: Any] = [
            "test": profile as Any,
            "grid": [
                "time" : T,
                "space": space.json
            ]
        ]
        
        return try? JSONSerialization.data(withJSONObject: json, options: [.sortedKeys, .prettyPrinted])
    }()
    
    // MARK: - Initializers
    
    private init(test: Test) {
        self.profile = test.description
        self.space = Grid(start: test.a, end: test.b, steps: test.N)
        self.T = test.T
        self.v0 = test.initial
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
            
            let v: V
            if abs(vR.density) < .ulpOfOne {
                v = .zero
            } else {
                let system = Eigen.system(physical: vR, gamma: gamma)
                v = system.reduce(into: V.zero) { v, system in

                    v += V(vector: 0.5 * (system.l * vL + system.l * vR) * system.r)
                }
            }
            
            mesh[xL] = vL
            mesh[x]  = v
            mesh[xR] = vR
        }
        while t <= T {
            guard let meshP = solution[t] else { assertionFailure("Time iteration error"); return }
            
            let tau = tau(average: Set(meshP.filter{ $0.key.side == .middle }.values))
            let tP = t
            t += tau
            solution[t] = [:]
            
            space.nodes().forEach { node in
                let xPL = BoundaryValue(value: node - space.step - space.halfed, side: .right)
                let xP = BoundaryValue(value: node - space.step, side: .middle)
                let xPR = BoundaryValue(value: node - space.step + space.halfed, side: .left)
                
                let xL = BoundaryValue(value: node - space.halfed, side: .right)
                let x = BoundaryValue(value: node, side: .middle)
                let xR = BoundaryValue(value: node + space.halfed, side: .left)
                
                let xNL = BoundaryValue(value: node + space.step - space.halfed, side: .right)
                let xN = BoundaryValue(value: node + space.step, side: .middle)
                let xNR = BoundaryValue(value: node + space.step + space.halfed, side: .left)
                
                /// left border
                let vStarL = vStar(t: tP, x: xL.value)
                let systemL = Eigen.system(physical: vStarL, gamma: gamma)
                let vL = boundaryV(
                    t: tP,
                    tau: tau,
                    vStar: vStarL,
                    system: systemL,
                    xLL: xPL,
                    xLM: xP,
                    xLR: xPR,
                    xRL: xL,
                    xRM: x,
                    xRR: xR
                )
                
                solution[t]?[xL] = vL
                
                /// right border
                let vStarR = vStar(t: tP, x: xR.value)
                let systemR = Eigen.system(physical: vStarR, gamma: gamma)
                let vR = boundaryV(
                    t: tP,
                    tau: tau,
                    vStar: vStarR,
                    system: systemR,
                    xLL: xL,
                    xLM: x,
                    xLR: xPR,
                    xRL: xNL,
                    xRM: xN,
                    xRR: xNR
                )
                
                solution[t]?[xR] = vR
                
                /// middle through flow F
                
            }
        }
    }
    
    // MARK: - Private methods
    
    private func tau(average vSet: Set<V>) -> Double {
        let maxLamba = vSet.map { abs($0.speed + Constants.speedOfSound(physical: $0, gamma: gamma)) }.max()!
        
        return Constants.sigmaGas * space.step / maxLamba
    }
    
    private func boundaryV(
        t: Double,
        tau: Double,
        vStar: V,
        system: [Eigen.System],
        xLL: BoundaryValue,
        xLM: BoundaryValue,
        xLR: BoundaryValue,
        xRL: BoundaryValue,
        xRM: BoundaryValue,
        xRR: BoundaryValue
    ) -> V {
        return system.reduce(into: V.zero) { result, system in
            let shiftedX = shifted(x: xLR.value, lambda: system.lambda, tau: tau)
            
            if system.lambda > 0 {
                let vL = Nv(
                    t: t,
                    x: shifted(x: xLR.value, lambda: system.lambda, tau: tau),
                    xL: xLL,
                    xM: xLM,
                    xR: xLR
                )
                
                result += V(vector: ((system.l * (0.5 * (vL + vStar))) * system.r))
            } else {
                let vR = Nv(
                    t: t,
                    x: shifted(x: xLR.value, lambda: system.lambda, tau: tau),
                    xL: xRL,
                    xM: xRM,
                    xR: xRR
                )
                
                result += V(vector: ((system.l * (0.5 * (vR - vStar))) * system.r))
            }
        }
    }
    
    private func vStar(t: Double, x: Double) -> V {
        let xL = BoundaryValue(value: x, side: .left)
        let xR = BoundaryValue(value: x, side: .right)
        
        guard
            let vL = solution[t]?[xL],
            let vR = solution[t]?[xR]
        else { return .zero }
        
        let sqrtDensityL = sqrt(vL.density)
        let sqrtDensityR = sqrt(vR.density)
        let sqrtSum = sqrtDensityL + sqrtDensityR
        
        let density = sqrtDensityL * sqrtDensityR
        let speed = (sqrtDensityL * vL.speed + sqrtDensityR * vR.speed) / sqrtSum
        
        let enthalpyL = Formulas.enthalpy(physical: vL, gamma: gamma)
        let enthalpyR = Formulas.enthalpy(physical: vR, gamma: gamma)
        let enthalpy = (sqrtDensityL * enthalpyL + sqrtDensityR * enthalpyR) / sqrtSum
        let pressure = Formulas.pressure(enthalpy: enthalpy, density: density, speed: speed, gamma: gamma)
        
        return V(density: density, speed: speed, pressure: pressure)
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
        
        return Nv(x: x, xL: xL.value, vL: vL, vM: vM, vR: vR)
    }
    
    private func Nv(
        x: Double,
        xL: Double,
        vL: V,
        vM: V,
        vR: V
    ) -> V {
        let xi = xi(x: x, xL: xL)
        let delta = delta(vL: vL, vR: vR)
        let sixth = sixth(vL: vL, vM: vM, vR: vR)
        
        return Nv(x: x, vL: vL, xi: xi, delta: delta, sixth: sixth)
    }
    
    private func Nv(
        x: Double,
        vL: V,
        xi: Double,
        delta: V,
        sixth: V
    ) -> V {
        return V(
            density: Nv(x: x, vL: vL.density, xi: xi, delta: delta.density, sixth: sixth.density),
            speed: Nv(x: x, vL: vL.speed, xi: xi, delta: delta.speed, sixth: sixth.speed),
            pressure: Nv(x: x, vL: vL.pressure, xi: xi, delta: delta.pressure, sixth: sixth.pressure)
        )
    }
    
    private func Nv(
        x: Double,
        xL: Double,
        vL: Double,
        vM: Double,
        vR: Double
    ) -> Double {
        let xi = xi(x: x, xL: xL)
        let delta = delta(vL: vL, vR: vR)
        let sixth = sixth(vL: vL, vM: vM, vR: vR)
        
        return Nv(x: x, vL: vL, xi: xi, delta: delta, sixth: sixth)
    }
    
    private func Nv(
        x: Double,
        vL: Double,
        xi: Double,
        delta: Double,
        sixth: Double
    ) -> Double {
        return vL +  xi * (delta + sixth * (1.0 - xi))
    }
    
    private func xi(x: Double, xL: Double) -> Double {
        return (x - xL) / space.step
    }
    
    private func delta(vL: Double, vR: Double) -> Double {
        return vR - vL
    }
    
    private func delta(vL: V, vR: V) -> V {
        return vR - vL
    }
    
    private func sixth(vL: Double, vM: Double, vR: Double) -> Double {
        return 6.0 * (vM - 0.5 * (vR + vL))
    }
    
    private func sixth(vL: V, vM: V, vR: V) -> V {
        return V(
            density: sixth(vL: vL.density, vM: vM.density, vR: vR.density),
            speed: sixth(vL: vL.speed, vM: vM.speed, vR: vR.speed),
            pressure: sixth(vL: vL.pressure, vM: vM.pressure, vR: vR.pressure)
        )
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
        
        static func enthalpy(
            density: Double,
            speed: Double,
            pressure: Double,
            gamma: Double = Constants.gamma
        ) -> Double {
            return (gamma / (gamma - 1) + density * speed * speed / 2) / density
        }
        
        static func enthalpy(physical v: V, gamma: Double = Constants.gamma) -> Double {
            return enthalpy(density: v.density, speed: v.speed, pressure: v.pressure, gamma: gamma)
        }
        
        static func pressure(
            enthalpy: Double,
            density: Double,
            speed: Double,
            gamma: Double = Constants.gamma
        ) -> Double {
            return (gamma - 1) / gamma * density * (enthalpy - speed * speed / 2)
        }
    }
    
    enum Eigen {
        
        // MARK: - Nested types
        
        struct System {
            
            let lambda: Double
            let r: [Double]
            let l: [Double]
        }
        
        // MARK: - Type methods
        
        static func system(
            density: Double,
            speed: Double,
            pressure: Double,
            gamma: Double = Constants.gamma
        ) -> [System] {
            let lambdas = Lambda(density: density, speed: speed, pressure: pressure, gamma: gamma)
            let R = R(density: density, speed: speed, pressure: pressure, gamma: gamma)
            let L = L(density: density, speed: speed, pressure: pressure, gamma: gamma)
            
            return [
                System(lambda: lambdas[0], r: R[0], l: L[0]),
                System(lambda: lambdas[1], r: R[1], l: L[1]),
                System(lambda: lambdas[2], r: R[2], l: L[2])
            ]
        }
        
        static func system(physical v: V, gamma: Double = Constants.gamma) -> [System] {
            return system(density: v.density, speed: v.speed, pressure: v.pressure, gamma: gamma)
        }
        
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

extension Gas1D {
    
    private func export() {
        let json = solution.keys.reduce(into: [String: Any]()) { json, t in
            json[String(t)] = space.nodes().reduce(into: [String: String]()) { mesh, node in
                let xL = BoundaryValue(value: node - space.halfed, side: .right)
                let x  = BoundaryValue(value: node, side: .middle)
                let xR = BoundaryValue(value: node + space.halfed, side: .left)

                guard
                    let vL = solution[t]?[xL],
                    let v  = solution[t]?[x],
                    let vR = solution[t]?[xR]
                else { return }

                mesh[String(xL.value)] = String(vL.density)
                mesh[String(x.value)] = String(v.density)
                mesh[String(xR.value)] = String(vR.density)
                
                let delta = delta(vL: vL, vR: vR)
                let sixth = sixth(vL: vL, vM: v, vR: vR)
                let step = space.step / 250
                for k in stride(from: xL.value + step, to: xR.value, by: step) {
                    let xi = xi(x: k, xL: xL.value)
                    let vK = Nv(x: k, vL: vL, xi: xi, delta: delta, sixth: sixth)
                    
                    mesh[String(k)] = String(vK.density)
                }
            }
        }
        
        let data = try? JSONSerialization.data(withJSONObject: json, options: [.sortedKeys, .prettyPrinted])
        save(file: "density", data: data, meta: metadata)
    }
}
