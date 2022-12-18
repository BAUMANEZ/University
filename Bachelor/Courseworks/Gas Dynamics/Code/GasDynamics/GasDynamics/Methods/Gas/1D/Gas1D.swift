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
        
        var plain: [Double] { [density, speed, pressure] }
        
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
        var u1: Double
        
        /// - rho * u
        var u2: Double
        
        /// - E
        var u3 : Double
        
        var plain: [Double] { [u1, u2, u3] }
        
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
        
        // MARK: - Type properties
        
        static let zero = F(f1: .zero, f2: .zero, f3: .zero)
        
        // MARK: - Operators
        
        static func +(lhs: F, rhs: F) -> F {
            return F(f1: lhs.f1 + rhs.f1, f2: lhs.f2 + rhs.f2, f3: lhs.f3 + rhs.f3)
        }
        
        static func -(lhs: F, rhs: F) -> F {
            return F(f1: lhs.f1 - rhs.f1, f2: lhs.f2 - rhs.f2, f3: lhs.f3 - rhs.f3)
        }
        
        static func *(x: Double, f: F) -> F {
            return F(f1: x * f.f1, f2: x * f.f2, f3: x * f.f3)
        }
        
        static func *(f: F, x: Double) -> F {
            return x * f
        }
        
        // MARK: - Internal properties
        
        /// - rho * u
        let f1: Double
        
        /// - rho * u^2 + p
        let f2: Double
        
        /// - (E + p) * u
        let f3: Double
        
        var plain: [Double] { [f1, f2, f3] }
        
        public init(physical v: V, gamma: Double = Constants.gamma) {
            self.init(density: v.density, speed: v.speed, pressure: v.pressure, gamma: gamma)
        }
        
        public init(density: Double, speed: Double, pressure: Double, gamma: Double = Constants.gamma) {
            self.init(
                f1: density * speed,
                f2: density * speed * speed + pressure,
                f3: (
                    Formulas.energy(
                        density: density,
                        speed: speed,
                        pressure: pressure,
                        gamma: gamma
                    ) + pressure
                ) * speed
            )
        }
        
        public init(vector: [Double]) {
            self.init(f1: vector[0], f2: vector[1], f3: vector[2])
        }
        
        private init(f1: Double, f2: Double, f3: Double) {
            self.f1 = f1
            self.f2 = f2
            self.f3 = f3
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
                
                guard let v = solution[tP]?[x] else { solution[t]?[x] = .zero; return }
                
                /// left border
                let leftFlow = flow(
                    x: xL, t: t, tP: tP, tau: tau,
                    xLL: xPL, xLM: xP, xLR: xPR,
                    xRL: xL, xRM: x, xRR: xR
                )
                
                /// right border
                let rightFlow = flow(
                    x: xR, t: t, tP: tP, tau: tau,
                    xLL: xL, xLM: x, xLR: xR,
                    xRL: xNL, xRM: xN, xRR: xNR
                )
                
                solution[t]?[x] = v - tau / space.step * (rightFlow - leftFlow)
            }
        }
    }
    
    // MARK: - Private methods
    
    private func tau(average vSet: Set<V>) -> Double {
        let maxLamba = vSet.map { abs($0.speed + Constants.speedOfSound(physical: $0, gamma: gamma)) }.max()!
        
        return Constants.sigmaGas * space.step / maxLamba
    }
    
    private func flow(
        x: BoundaryValue,
        t: Double,
        tP: Double,
        tau: Double,
        xLL: BoundaryValue,
        xLM: BoundaryValue,
        xLR: BoundaryValue,
        xRL: BoundaryValue,
        xRM: BoundaryValue,
        xRR: BoundaryValue
    ) -> V {
        let vStar = vStar(t: tP, x: x.value)
        
        guard vStar.density > .ulpOfOne else { return .zero }
        
        let system = Eigen.system(physical: vStar, gamma: gamma)
        
        let vBoundary = boundaryV(
            t: tP,
            tau: tau,
            vStar: vStar,
            system: system,
            xLL: xLL,
            xLM: xLM,
            xLR: xLR,
            xRL: xRL,
            xRM: xRM,
            xRR: xRR
        )
        
        solution[t]?[x] = vBoundary
        
        let vL: V = {
            guard
                let maxSystem = system.max(by: { $0.lambda < $1.lambda }),
                maxSystem.lambda >= 0
            else { return .zero }
            
            let maxV = averageBoundaryVp(
                t: tP, tau: tau, system: maxSystem,
                xLL: xLL, xLM: xLM, xLR: xLR,
                xRL: xRL, xRM: xRM, xRR: xRR
            )
            
            let sumVL = system.filter { $0.lambda >= 0 }.reduce(into: V.zero) { v, system in
                let averageVp = averageBoundaryVp(
                    t: tP, tau: tau, system: system,
                    xLL: xLL, xLM: xLM, xLR: xLR,
                    xRL: xRL, xRM: xRM, xRR: xRR
                )
                v += V(vector: system.r * (system.l * (averageVp - maxV)))
            }
            
            return sumVL + maxV
        }()
        
        let vR: V = {
            guard
                let minSystem = system.min(by: { $0.lambda < $1.lambda }),
                minSystem.lambda < 0
            else { return .zero }
            
            let minV = averageBoundaryVp(
                t: tP, tau: tau, system: minSystem,
                xLL: xLL, xLM: xLM, xLR: xLR,
                xRL: xRL, xRM: xRM, xRR: xRR
            )
            
            let sumVR = system.filter { $0.lambda < 0 }.reduce(into: V.zero) { v, system in
                let averageVp = averageBoundaryVp(
                    t: tP, tau: tau, system: system,
                    xLL: xLL, xLM: xLM, xLR: xLR,
                    xRL: xRL, xRM: xRM, xRR: xRR
                )
                v += V(vector: system.r * (system.l * (averageVp - minV)))
            }
            
            return minV - sumVR
        }()
        
        let fL = F(physical: vL, gamma: gamma)
        let fR = F(physical: vR, gamma: gamma)
        
        let M = Eigen.M(physical: vStar, gamma: gamma)
        
        let deltaV = delta(vL: vL, vR: vR)
        let sum = system.reduce(into: [Double](repeating: 0, count: 3)) { sum, system in
            sum += 0.5 * system.lambda.magnitude * (system.l * deltaV) * (M * system.r)
        }
        
        let flow = 0.5 * (fL + fR) - F(vector: sum)
        let invM = Eigen.InvM(physical: vStar, gamma: gamma)
        let _flow = invM * flow.plain
        
        return V(vector: _flow)
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
            guard system.lambda.magnitude > .ulpOfOne else { return () }
            
            let shiftedX = shifted(x: xLR.value, lambda: system.lambda, tau: tau)
            
            if system.lambda > 0 {
                let vL = Nv(
                    t: t,
                    x: shiftedX,
                    xL: xLL,
                    xM: xLM,
                    xR: xLR
                )
                
                result += V(vector: ((system.l * (0.5 * (vL + vStar))) * system.r))
            } else {
                let vR = Nv(
                    t: t,
                    x: shiftedX,
                    xL: xRL,
                    xM: xRM,
                    xR: xRR
                )
                
                result += V(vector: ((system.l * (0.5 * (vR - vStar))) * system.r))
            }
        }
    }
    
    private func averageBoundaryVp(
        t: Double,
        tau: Double,
        system: Eigen.System,
        xLL: BoundaryValue,
        xLM: BoundaryValue,
        xLR: BoundaryValue,
        xRL: BoundaryValue,
        xRM: BoundaryValue,
        xRR: BoundaryValue
    ) -> V {
        
        guard system.lambda.magnitude > .ulpOfOne else { return .zero }
        
        let delta = system.lambda.magnitude * tau
        let shiftedX = shifted(x: xLR.value, lambda: system.lambda, tau: tau)
        
        let shiftedLeft: V
        let shiftedMiddle: V
        let shiftedRight: V
        
        if system.lambda > 0 {
            shiftedLeft = Nv(
                t: t,
                x: shiftedX,
                xL: xLL,
                xM: xLM,
                xR: xLR
            )
            
            shiftedMiddle = Nv(
                t: t,
                x: shiftedX + 0.5 * delta,
                xL: xLL,
                xM: xLM,
                xR: xLR
            )
            
            shiftedRight = Nv(
                t: t,
                x: shiftedX + delta,
                xL: xLL,
                xM: xLM,
                xR: xLR
            )
        } else {
            shiftedLeft = Nv(
                t: t,
                x: shiftedX - delta,
                xL: xLL,
                xM: xLM,
                xR: xLR
            )
            
            shiftedMiddle = Nv(
                t: t,
                x: shiftedX - 0.5 * delta,
                xL: xLL,
                xM: xLM,
                xR: xLR
            )
            
            shiftedRight = Nv(
                t: t,
                x: shiftedX,
                xL: xLL,
                xM: xLM,
                xR: xLR
            )
        }
        
        return V(vector: ((system.l * (1 / 6 * (shiftedLeft + 4 * shiftedMiddle + shiftedRight))) * system.r))
    }
    
    /// Note: Nv is not applicable for shifting vector of variables V, it is supposed to initialize V in a known parabolic equation on a previous time stamp,
    /// so in order to initialize V on a current time stamp you have to define a basis in a cell and calculate a new V like (lp * Vprev) * rp
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
            return pressure / (gamma - 1) + 0.5 * density * speed * speed
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
            return (pressure * gamma / (gamma - 1) + density * speed * speed / 2) / density
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
        
        static func M(
            density: Double,
            speed: Double,
            pressure: Double,
            gamma: Double = Constants.gamma
        ) -> Matrix {
            var result = Matrix(n: 3, initial: 0)

            result[0,0] = 1
            result[0,1] = 0
            result[0,2] = 0

            result[1,0] = speed
            result[1,1] = density
            result[1,2] = 0

            result[2,0] = speed * speed / 2
            result[2,1] = density * speed
            result[2,2] = 1 / (gamma - 1)
            
            return result
        }
        
        static func M(physical v: V, gamma: Double = Constants.gamma) -> Matrix {
            return M(density: v.density, speed: v.speed, pressure: v.pressure, gamma: gamma)
        }
        
        static func InvM(
            density: Double,
            speed: Double,
            pressure: Double,
            gamma: Double = Constants.gamma
        ) -> Matrix {
            var result = Matrix(n: 3, initial: 0)

            result[0,0] = 1
            result[0,1] = 0
            result[0,2] = 0

            result[1,0] = -speed / density
            result[1,1] = 1 / density
            result[1,2] = 0

            result[2,0] = speed * speed / 2 * (gamma - 1)
            result[2,1] = -speed * (gamma - 1)
            result[2,2] = gamma - 1
            
            return result
        }
        
        static func InvM(physical v: V, gamma: Double = Constants.gamma) -> Matrix {
            return InvM(density: v.density, speed: v.speed, pressure: v.pressure, gamma: gamma)
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
                let step = space.step / 5
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
