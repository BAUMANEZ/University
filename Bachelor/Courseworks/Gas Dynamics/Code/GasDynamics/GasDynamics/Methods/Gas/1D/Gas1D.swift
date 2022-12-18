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
    typealias EigenMesh = [BoundaryValue: [Eigen.System]]
    
    typealias Solution = [Time: Mesh]
    typealias EigenCell = [Time: EigenMesh]
    
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
        
        static func -=(v1: inout V, v2: V) {
            v1.density -= v2.density
            v1.speed -= v2.speed
            v1.pressure -= v2.pressure
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
        
        static func +=(lhs: inout F, rhs: F) {
            lhs.f1 += rhs.f1
            lhs.f2 += rhs.f2
            lhs.f3 += rhs.f3
        }
        
        // MARK: - Internal properties
        
        /// - rho * u
        var f1: Double
        
        /// - rho * u^2 + p
        var f2: Double
        
        /// - (E + p) * u
        var f3: Double
        
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
    private(set) var eigenCell: EigenCell = [:]
    
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
        eigenCell[t] = [:]
        
        solution[t] = space.nodes().reduce(into: Mesh()) { mesh, node in
            let xL = BoundaryValue(value: node - space.halfed, side: .right)
            let x  = BoundaryValue(value: node, side: .middle)
            let xR = BoundaryValue(value: node + space.halfed, side: .left)
            let vL = v0(xL.value)
            let vR = v0(xR.value)
            
            let system = Eigen.system(physical: vR.density > .ulpOfOne ? vR : vL, gamma: gamma)
            eigenCell[t]?[x] = system
            
            let v = system.reduce(into: V.zero) { v, system in
                v += V(vector: 0.5 * (system.l * vL + system.l * vR) * system.r)
            }
            
            mesh[xL] = vL
            mesh[x]  = v
            mesh[xR] = vR
        }
        
        while t <= T {
            guard let meshP = solution[t] else { assertionFailure("Time iteration error"); return }
            
            let tau = tau(average: Set(meshP.filter{ $0.key.side == .middle}.values))
            let tP = t
            t += tau
            solution[t] = [:]
            eigenCell[t] = [:]
            
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
                
                guard
                    let v = solution[tP]?[safe: x],
                    let systemP = eigenCell[tP]?[safe: x]
                else {
                    assertionFailure("Eigen system must be defined")
                    return
                }
                
                /// right border
                let (averageVL, VL, Vl): (V, V, V) = {
                    guard
                        let maxSystemP = systemP.last,
                        maxSystemP.lambda > 0,
                        let vL = solution[tP]?[safe: xL],
                        let vR = solution[tP]?[safe: xR]
                    else { return (.zero, .zero, .zero) }
                    
                    let averageVL = averageBoundaryVp(
                        t: t, tP: tP, tau: tau, lambda: maxSystemP.lambda,
                        xLL: xL, xLM: x, xLR: xR,
                        xRL: .zero, xRM: .zero, xRR: .zero,
                        vLL: vL, vLM: v, vLR: vR,
                        vRL: .zero, vRM: .zero, vRR: .zero
                    )
                    
                    let averageVLSystem = Eigen.system(physical: averageVL, gamma: gamma).filter { $0.lambda > 0 }
                    
                    let VL = averageVLSystem.reduce(into: V.zero) { v, system in
                        let x = shifted(x: xR.value, lambda: system.lambda, tau: tau)
                        v += Nvp(x: x, system: system, xL: xL.value, vL: vL, vM: v, vR: vR)
                    }
                    
                    let Vl = averageVL + averageVLSystem.reduce(into: V.zero) { v, system in
                        let averageVp = averageBoundaryVp(
                            t: t, tP: tP, tau: tau, lambda: system.lambda,
                            xLL: xL, xLM: x, xLR: xR,
                            xRL: .zero, xRM: .zero, xRR: .zero,
                            vLL: vL, vLM: v, vLR: vR,
                            vRL: .zero, vRM: .zero, vRR: .zero
                        )
                        
                        v += V(vector: system.r * (system.l * (averageVp - averageVL)))
                    }
                    
                    return (averageVL, VL, Vl)
                }()
                                
                let (averageVR, VR, Vr): (V, V, V) = {
                    guard
                        let minSystemP = systemP.first,
                        minSystemP.lambda < 0,
                        let vNL = solution[tP]?[safe: xNL],
                        let vN = solution[tP]?[safe: xN],
                        let vNR = solution[tP]?[safe: xNR]
                    else {
                        return (.zero, .zero, .zero)
                    }

                    let averageVR = averageBoundaryVp(
                        t: t, tP: tP, tau: tau, lambda: minSystemP.lambda,
                        xLL: .zero, xLM: .zero, xLR: .zero,
                        xRL: xNL, xRM: xN, xRR: xNR,
                        vLL: .zero, vLM: .zero, vLR: .zero,
                        vRL: vNL, vRM: vN, vRR: vNR
                    )
                    
                    let averageVRSystem = Eigen.system(physical: averageVR, gamma: gamma).filter {
                        $0.lambda < 0
                    }

                    let VR = averageVRSystem.reduce(into: V.zero) { v, system in
                        let x = shifted(x: xR.value, lambda: system.lambda, tau: tau)
                        v += Nvp(x: x, system: system, xL: xNL.value, vL: vNL, vM: vN, vR: vNR)
                    }

                    let Vr = averageVRSystem.reduce(into: V.zero) { v, system in
                        let averageVp = averageBoundaryVp(
                            t: t, tP: tP, tau: tau, lambda: system.lambda,
                            xLL: .zero, xLM: .zero, xLR: .zero,
                            xRL: xNL, xRM: xN, xRR: xNR,
                            vLL: .zero, vLM: .zero, vLR: .zero,
                            vRL: vNL, vRM: vN, vRR: vNR
                        )

                        v += V(vector: system.r * (system.l * (averageVp - averageVL)))
                    } - averageVR
                    
                    return (averageVR, VR, Vr)
                }()
                
                let vStar = vStar(vL: averageVL, vR: averageVR)
                let newSystem = Eigen.system(physical: vStar, gamma: gamma)
                let proportionalVStar = 0.5 * newSystem.reduce(into: V.zero) { v, system in
                    guard system.lambda.magnitude > .ulpOfOne else { return () }

                    if system.lambda > 0 {
                        v += V(vector: (system.l * vStar) * system.r)
                    } else {
                        v -= V(vector: (system.l * vStar) * system.r)
                    }
                }
                let newVR = 0.5 * (VL + VR) + proportionalVStar

                solution[t]?[xR] = newVR
                eigenCell[t]?[x] = newSystem

                let flowVL = F(physical: averageVL, gamma: gamma)
                let flowVR = F(physical: averageVR, gamma: gamma)
                let flowVRVLHalved = 0.5 * (flowVL + flowVR)

                let M = Eigen.M(physical: vStar, gamma: gamma)

                let flowUR = flowVRVLHalved - 0.5 * newSystem.reduce(into: F.zero) { flow, system in
                    flow += F(vector: system.lambda.magnitude * (system.l * (VR - VL)) * (M * system.r) )
                }
                
                let invM = Eigen.InvM(physical: vStar, gamma: gamma)
                let flowR = invM * flowUR.plain
                var newV = v - tau / space.step * V(vector: flowR)
                
                /// left border
//                guard
//                    let vNL = solution[tP]?[xNL],
//                    let vN = solution[tP]?[xN],
//                    let vNR = solution[tP]?[xNR]
//                else {
//                    solution[t]?[xL] = .zero
//                    solution[t]?[x] = newV
//                    return
//                }
                
                solution[t]?[x] = newV
                solution[t]?[xL] = newV
            }
        }
    }
    
    // MARK: - Private methods
    
    private func tau(average vSet: Set<V>) -> Double {
        let maxLamba = vSet.map {
            let c = Constants.speedOfSound(physical: $0, gamma: gamma)
            return max(abs($0.speed + c), abs($0.speed - c))
        }.max()!
        
        return Constants.sigmaGas * space.step / maxLamba
    }
    
    private func alpha(v: V, l: [Double]) -> Double {
        return l * v
    }
    
    private func vStar(vL: V, vR: V) -> V {
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

    private func averageBoundaryVp(
        t: Double,
        tP: Double,
        tau: Double,
        lambda: Double,
        xLL: BoundaryValue,
        xLM: BoundaryValue,
        xLR: BoundaryValue,
        xRL: BoundaryValue,
        xRM: BoundaryValue,
        xRR: BoundaryValue,
        vLL: V,
        vLM: V,
        vLR: V,
        vRL: V,
        vRM: V,
        vRR: V
    ) -> V {
        guard lambda.magnitude > .ulpOfOne  else { return .zero }

        let delta = lambda.magnitude * tau
        let shiftedX = shifted(x: xLR.value, lambda: lambda, tau: tau)

        if lambda > 0 {
            guard let systemP = eigenCell[tP]?[safe: xLM] else { return .zero }
            
            return systemP.reduce(into: V.zero) { v, system in
                let alphaLp = alpha(v: vLL, l: system.l)
                let alphaMp = alpha(v: vLM, l: system.l)
                let alphaRp = alpha(v: vLR, l: system.l)
                
                let alpha1 = NalphaP(x: shiftedX, xL: xLL.value, alphaLp: alphaLp, alphaMp: alphaMp, alphaRp: alphaRp)
                let alpha2 = NalphaP(x: shiftedX + 0.5 * delta, xL: xLL.value, alphaLp: alphaLp, alphaMp: alphaMp, alphaRp: alphaRp)
                let alpha3 = NalphaP(x: shiftedX + delta, xL: xLL.value, alphaLp: alphaLp, alphaMp: alphaMp, alphaRp: alphaRp)
                
                v += V(vector: (1 / 6 * (alpha1 + 4 * alpha2 + alpha3)) * system.r)
            }
        } else {
            guard let systemP = eigenCell[tP]?[safe: xRM] else { return .zero }
            
            return systemP.reduce(into: V.zero) { v, system in
                let alphaLp = alpha(v: vRL, l: system.l)
                let alphaMp = alpha(v: vRM, l: system.l)
                let alphaRp = alpha(v: vRR, l: system.l)
                
                let alpha1 = NalphaP(x: shiftedX - delta, xL: xRL.value, alphaLp: alphaLp, alphaMp: alphaMp, alphaRp: alphaRp)
                let alpha2 = NalphaP(x: shiftedX - 0.5 * delta, xL: xRL.value, alphaLp: alphaLp, alphaMp: alphaMp, alphaRp: alphaRp)
                let alpha3 = NalphaP(x: shiftedX, xL: xRL.value, alphaLp: alphaLp, alphaMp: alphaMp, alphaRp: alphaRp)
                
                v += V(vector: (1 / 6 * (alpha1 + 4 * alpha2 + alpha3)) * system.r)
            }
        }
    }
    
    private func Nv(
        t: Double,
        x: Double,
        system: [Eigen.System],
        xL: BoundaryValue,
        xM: BoundaryValue,
        xR: BoundaryValue
    ) -> V {
        guard
            let vL = solution[t]?[xL],
            let vM  = solution[t]?[xM],
            let vR = solution[t]?[xR]
        else { return .zero }
        
        return Nv(x: x, system: system, xL: xL.value, vL: vL, vM: vM, vR: vR)
    }
    
    private func Nv(
        x: Double,
        system: [Eigen.System],
        xL: Double,
        vL: V,
        vM: V,
        vR: V
    ) -> V {
        return system.reduce(into: V.zero) { v, system in
            v += Nvp(x: x, system: system, xL: xL, vL: vL, vM: vM, vR: vR)
        }
    }
    
    private func Nvp(
        x: Double,
        system: Eigen.System,
        xL: Double,
        vL: V,
        vM: V,
        vR: V
    ) -> V {
        let alphaLp = alpha(v: vL, l: system.l)
        let alphaMp = alpha(v: vM, l: system.l)
        let alphaRp = alpha(v: vR, l: system.l)
        
        return Nvp(x: x, system: system, xL: xL, alphaLp: alphaLp, alphaMp: alphaMp, alphaRp: alphaRp)
    }
    
    private func Nvp(
        x: Double,
        system: Eigen.System,
        xL: Double,
        alphaLp: Double,
        alphaMp: Double,
        alphaRp: Double
    ) -> V {
        return V(vector: NalphaP(x: x, xL: xL, alphaLp: alphaLp, alphaMp: alphaMp, alphaRp: alphaRp) * system.r)
    }
    
    private func NalphaP(
        x: Double,
        xL: Double,
        alphaLp: Double,
        alphaMp: Double,
        alphaRp: Double
    ) -> Double {
        let xi = xi(x: x, xL: xL)
        let delta = delta(alphaL: alphaLp, alphaR: alphaRp)
        let sixth = sixth(alphaL: alphaLp, alphaM: alphaMp, alphaR: alphaRp)
        
        return NalphaP(x: x, alphaLp: alphaLp, xi: xi, delta: delta, sixth: sixth)
    }
    
    private func NalphaP(
        x: Double,
        alphaLp: Double,
        xi: Double,
        delta: Double,
        sixth: Double
    ) -> Double {
        return alphaLp +  xi * (delta + sixth * (1.0 - xi))
    }
    
    private func xi(x: Double, xL: Double) -> Double {
        return (x - xL) / space.step
    }
    
    private func delta(alphaL: Double, alphaR: Double) -> Double {
        return alphaR - alphaL
    }
    
    private func sixth(alphaL: Double, alphaM: Double, alphaR: Double) -> Double {
        return 6.0 * (alphaM - 0.5 * (alphaR + alphaL))
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
//
//                let delta = delta(vL: vL, vR: vR)
//                let sixth = sixth(vL: vL, vM: v, vR: vR)
//                let step = space.step / 5
//                for k in stride(from: xL.value + step, to: xR.value, by: step) {
//                    let xi = xi(x: k, xL: xL.value)
//                    let vK = Nv(x: k, vL: vL, xi: xi, delta: delta, sixth: sixth)
//
//                    mesh[String(k)] = String(vK.density)
//                }
            }
        }
        
        let data = try? JSONSerialization.data(withJSONObject: json, options: [.sortedKeys, .prettyPrinted])
        save(file: "density", data: data, meta: metadata)
    }
}
