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
            self.init(
                u1: density,
                u2: density * speed,
                u3: Formulas.energy(density: density, speed: speed, pressure: pressure, gamma: gamma)
            )
        }
        
        public init(vector: [Double]) {
            self.init(u1: vector[0], u2: vector[1], u3: vector[2])
        }
        
        public init(u1: Double, u2: Double, u3: Double) {
            self.u1 = u1
            self.u2 = u2
            self.u3 = u3
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
        
        static func -=(lhs: inout F, rhs: F) {
            lhs.f1 -= rhs.f1
            lhs.f2 -= rhs.f2
            lhs.f3 -= rhs.f3
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
            let v = v0(x.value)
            let vR = v0(xR.value)
            
            mesh[xL] = vL
            mesh[x]  = v
            mesh[xR] = vR
            
            eigenCell[t]?[x] = Eigen.system(physical: v, gamma: gamma)
        }
                
        while t < T {
            guard let meshP = solution[t] else { assertionFailure("Time iteration error"); return }
            
            let tau = tau(average: Set(meshP.filter{ $0.key.side == .middle}.values))
            let tP = t
            t += tau
            solution[t] = [:]
            eigenCell[t] = [:]
            
            print("-------\(t)-------")
            
            space.nodes().forEach { node in
                let x = BoundaryValue(value: node, side: .middle)
                
                guard let v = solution[tP]?[safe: x] else { fatalError() }
                
                let xPL = BoundaryValue(value: node - space.step - space.halfed, side: .right)
                let xP = BoundaryValue(value: node - space.step, side: .middle)
                let xPR = BoundaryValue(value: node - space.step + space.halfed, side: .left)
                let xL = BoundaryValue(value: node - space.halfed, side: .right)
                let xR = BoundaryValue(value: node + space.halfed, side: .left)
                let xNL = BoundaryValue(value: node + space.step - space.halfed, side: .right)
                let xN = BoundaryValue(value: node + space.step, side: .middle)
                let xNR = BoundaryValue(value: node + space.step + space.halfed, side: .left)
                
                let (vR, flowR) = calculateFlowAndNewBoundaryV(
                    tP: tP, tau: tau,
                    xPL: xL, xP: x, xPR: xR,
                    xNL: xNL, xN: xN, xNR: xNR
                )
                
                let (vL, flowL) = calculateFlowAndNewBoundaryV(
                    tP: tP, tau: tau,
                    xPL: xPL, xP: xP, xPR: xPR,
                    xNL: xL, xN: x, xNR: xR
                )
                                
                solution[t]?[xL] = vL
                solution[t]?[xR] = vR
                
                let u = U(physical: v, gamma: gamma)
                let newUVector = u.plain - tau / space.step * (flowR - flowL).plain
                let invM = Eigen.InvM(conservative: U(vector: newUVector), gamma: gamma)
                let newV = V(vector: invM * newUVector)
                
                eigenCell[t]?[x] = Eigen.system(physical: newV, gamma: gamma)
                solution[t]?[x] = newV
            }

            return
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
    
    private func calculateFlowAndNewBoundaryV(
        tP: Double,
        tau: Double,
        xPL: BoundaryValue,
        xP : BoundaryValue,
        xPR: BoundaryValue,
        xNL: BoundaryValue,
        xN : BoundaryValue,
        xNR: BoundaryValue
    ) -> (V, F) {
        
        let (VL, Vl): (V, V) = {
            guard
                let systemL = eigenCell[tP]?[safe: xP],
                let vL = solution[tP]?[safe: xPL],
                let v = solution[tP]?[safe: xP],
                let vR = solution[tP]?[safe: xPR]
            else {
                let vl = v0(space.start)
                
                return (vl, vl)
            }
            
            guard !systemL.contains(where: { $0.lambda.isZero }) else {
                let vl = solution[tP]?[safe: xPR] ?? .zero
                
                return (vl, vl)
            }
            
            return leftSidePair(system: systemL, tau: tau, xL: xPL, x: xP, xR: xPR, vL: vL, v: v, vR: vR)
        }()
        
        let (VR, Vr): (V, V) = {
            guard
                let systemR = eigenCell[tP]?[safe: xN],
                let vL = solution[tP]?[safe: xNL],
                let v = solution[tP]?[safe: xN],
                let vR = solution[tP]?[safe: xNR]
            else {
                let vr = v0(space.end)
                
                return (vr, vr)
            }
            
            guard !systemR.contains(where: { $0.lambda.isZero }) else {
                let vl = solution[tP]?[safe: xNR] ?? .zero
                
                return (vl, vl)
            }

            return rightSidePair(system: systemR, tau: tau, xL: xNL, x: xN, xR: xNR, vL: vL, v: v, vR: vR)
        }()
        
        let vStar = vStar(vL: VL, vR: VR)
        let systemVStar = Eigen.system(physical: vStar, gamma: gamma)
        
        var newBoundaryV = 0.5 * (VL + VR)
        var flow = 0.5 * (F(physical: Vl, gamma: gamma) + F(physical: Vr, gamma: gamma))
        
        guard !systemVStar.contains(where: { $0.lambda.isZero }) else {
            return (newBoundaryV, flow)
        }
        
        let M = Eigen.M(physical: vStar, gamma: gamma)
        
        systemVStar.forEach { system in
            let vp = 0.5 * V(vector: (system.l * vStar) * system.r)
    
            if system.lambda > 0 {
                newBoundaryV += vp
            } else {
                newBoundaryV -= vp
            }
            
            flow -= 0.5 * system.lambda.magnitude * F(vector: (system.l * (Vr - Vl)) * (M * system.r))
        }
        
        return (newBoundaryV, flow)
    }
    
    private func leftSidePair(
        system systemP: [Eigen.System],
        tau: Double,
        xL: BoundaryValue,
        x: BoundaryValue,
        xR: BoundaryValue,
        vL: V,
        v: V,
        vR: V
    ) -> (V, V) {
        guard
            let maxSystemP = systemP.last,
            maxSystemP.lambda > 0
        else { return (.zero, .zero) }
        
        let averageVL = averageBoundaryVp(
            tau: tau, system: maxSystemP, boundary: xR,
            xL: xL, x: x, xR: xR,
            vL: vL, v: v, vR: vR
        )
                
        /// for boundary
        let VL = systemP.reduce(into: V.zero) { result, system in
            guard system.lambda > 0 else { return () }
            
            let x = shifted(x: xR.value, lambda: system.lambda, tau: tau)
            let vp = Nvp(x: x, system: system, xL: xL.value, vL: vL, vM: v, vR: vR)
            
            result += vp
        }
        
        /// for flow
        let Vl = averageVL + systemP.reduce(into: V.zero) { result, system in
            guard system.lambda > 0 else { return () }
            
            let averageVp = averageBoundaryVp(
                tau: tau, system: system, boundary: xR,
                xL: xL, x: x, xR: xR,
                vL: vL, v: v, vR: vR
            )
                        
            result += V(vector: system.r * (system.l * (averageVp - averageVL)))
        }
        
        return (VL, Vl)
    }
    
    private func rightSidePair(
        system systemP: [Eigen.System],
        tau: Double,
        xL: BoundaryValue,
        x: BoundaryValue,
        xR: BoundaryValue,
        vL: V,
        v: V,
        vR: V
    ) -> (V, V) {
        guard
            let minSystemP = systemP.first,
            minSystemP.lambda < 0
        else { return (.zero, .zero) }

        let averageVR = averageBoundaryVp(
            tau: tau, system: minSystemP, boundary: xL,
            xL: xL, x: x, xR: xR,
            vL: vL, v: v, vR: vR
        )
        
        /// for boundary
        let VR = systemP.reduce(into: V.zero) { result, system in
            guard system.lambda < 0 else { return () }
                        
            let x = shifted(x: xL.value, lambda: system.lambda, tau: tau)
            let vp = Nvp(x: x, system: system, xL: xL.value, vL: vL, vM: v, vR: vR)
            
            result += vp
        }

        /// for flow
        let Vr = averageVR - systemP.reduce(into: V.zero) { result, system in
            guard system.lambda < 0 else { return () }
            
            let averageVp = averageBoundaryVp(
                tau: tau, system: system, boundary: xL,
                xL: xL, x: x, xR: xR,
                vL: vL, v: v, vR: vR
            )

            result += V(vector: system.r * ((averageVp - averageVR) * system.l))
        }
        
        return (VR, Vr)
    }
    
    private func vStar(vL: V, vR: V) -> V {
        guard !vL.density.isZero else {
            return vR
        }
        
        guard !vR.density.isZero else {
            return vL
        }
        
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
        tau: Double,
        system: Eigen.System,
        boundary: BoundaryValue,
        xL: BoundaryValue,
        x: BoundaryValue,
        xR: BoundaryValue,
        vL: V,
        v: V,
        vR: V
    ) -> V {
        guard !system.lambda.isZero else { return .zero }
        
        let delta = system.lambda * tau
        let shiftedX = shifted(x: boundary.value, lambda: system.lambda, tau: tau)
        
        let alphaLp = alpha(v: vL, l: system.l)
        let alphaMp = alpha(v: v, l: system.l)
        let alphaRp = alpha(v: vR, l: system.l)
        
        let alpha1 = NalphaP(x: shiftedX, xL: xL.value, alphaLp: alphaLp, alphaMp: alphaMp, alphaRp: alphaRp)
        let alpha2 = NalphaP(x: shiftedX + 0.5 * delta, xL: xL.value, alphaLp: alphaLp, alphaMp: alphaMp, alphaRp: alphaRp)
        let alpha3 = NalphaP(x: shiftedX + delta, xL: xL.value, alphaLp: alphaLp, alphaMp: alphaMp, alphaRp: alphaRp)
        
        return V(vector: (1 / 6 * (alpha1 + 4 * alpha2 + alpha3)) * system.r)
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
        return system.reduce(into: V.zero) { result, system in
            result += Nvp(x: x, system: system, xL: xL, vL: vL, vM: vM, vR: vR)
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
            guard !density.isZero else { fatalError() }
            
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
            guard !density.isZero else { fatalError() }
            
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
            guard !density.isZero else { fatalError() }
            
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
            gamma: Double = Constants.gamma
        ) -> Matrix {
            guard !density.isZero else { fatalError() }
            
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
            return M(density: v.density, speed: v.speed, gamma: gamma)
        }
        
        static func InvM(
            density: Double,
            speed: Double,
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
            return InvM(density: v.density, speed: v.speed, gamma: gamma)
        }
        
        static func InvM(conservative u: U, gamma:Double = Constants.gamma) -> Matrix {
            return InvM(density: u.u1, speed: u.u2 / u.u1, gamma: gamma)
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
                    let vR = solution[t]?[xR],
                    let system = eigenCell[t]?[safe: x]
                else { return }

                mesh[String(xL.value)] = String(vL.density)
                mesh[String(x.value)] = String(v.density)
                mesh[String(xR.value)] = String(vR.density)

//                let step = space.step / 20
//                for k in stride(from: xL.value + step, to: xR.value, by: step) {
//                    let vK = Nv(t: t, x: k, system: system, xL: xL, xM: x, xR: xR)
//
//                    mesh[String(k)] = String(vK.density)
//                }
            }
        }
        
        let data = try? JSONSerialization.data(withJSONObject: json, options: [.sortedKeys, .prettyPrinted])
        save(file: "density", data: data, meta: metadata)
    }
}
