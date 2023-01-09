//
//  Gas1D.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 31.10.2022.
//

import Foundation

protocol ArrayToVariableConverable {
    
    var values: [Double] { get set }
    
    init(vector: [Double])
}

extension ArrayToVariableConverable {
    
    // MARK: - Type methods
    
    static func *(l: [Double], v: Self) -> Double { l * v.values }
    
    static func *(v: Self, r: [Double]) -> Double { v.values * r }
    
    static func *(v: Self, x: Double) -> Self { Self(vector: v.values * x) }
    
    static func *(x: Double, v: Self) -> Self { Self(vector: x * v.values) }

    static func *=(v: inout Self, x: Double) { v.values *= x }

    static func *=(x: Double, v: inout Self) { v.values *= x }

    static func +(v1: Self, v2: Self) -> Self { Self(vector: v1.values + v2.values) }

    static func +=(v1: inout Self, v2: Self) { v1.values += v2.values }

    static func -=(v1: inout Self, v2: Self) { v1.values -= v2.values }

    static func -(v1: Self, v2: Self) -> Self { Self(vector: v1.values - v2.values) }
}

final class Gas1D: JSONConvertableAlgorithm {
    
    // MARK: - Nested types
    
    typealias Mesh = [Double: V]
    typealias EigenMesh = [GridCell: [Eigen.System]]
    
    typealias Solution = [Time: Mesh]
    typealias EigenCell = [Time: EigenMesh]
    
    struct V: ArrayToVariableConverable, Hashable {
        
        // MARK: - Type properties
        
        static let zero = V(density: .zero, speed: .zero, pressure: .zero)
        
        // MARK: - Internal properties
        
        var values: [Double]
        
        var density: Double { values[0] }
        var speed: Double { values[1] }
        var pressure: Double { values[2] }
        
        // MARK: - Initializers
                
        public init(conservative u: U, gamma: Double = Constants.gamma) {
            let density = u.u1
            let speed = u.u2 / u.u1
            let _energy = u.u3 - 0.5 * density * speed * speed
            
            self.init(density: density, speed: speed, pressure: _energy * (gamma - 1))
        }
        
        public init(density: Double, speed: Double, pressure: Double) {
            self.init(vector: [density, speed, pressure])
        }
        
        init(vector: [Double]) {
            self.values = vector
        }
    }
    
    struct U: ArrayToVariableConverable, Hashable {
        
        // MARK: - Type properties
        
        static let zero = U(physical: .zero)
        
        // MARK: - Internal properties
        
        var values: [Double]
        
        /// - rho
        var u1: Double { values[0] }
        
        /// - rho * u
        var u2: Double { values[1] }
        
        /// - E
        var u3 : Double { values[2] }
                
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
        
        public init(u1: Double, u2: Double, u3: Double) {
            self.init(vector: [u1, u2, u3])
        }
        
        init(vector: [Double]) {
            self.values = vector
        }
    }
    
    struct F: ArrayToVariableConverable, Hashable {
        
        // MARK: - Type properties
        
        static let zero = F(f1: .zero, f2: .zero, f3: .zero)
        
        // MARK: - Internal properties
        
        var values: [Double]
        
        /// - rho * u
        var f1: Double { values[0] }
        
        /// - rho * u^2 + p
        var f2: Double { values[1] }
        
        /// - (E + p) * u
        var f3: Double { values[2] }
                
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
        
        private init(f1: Double, f2: Double, f3: Double) {
            self.init(vector: [f1, f2, f3])
        }
        
        init(vector: [Double]) {
            self.values = vector
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
            let meshCell = GridCell(value: node, grid: space)
            
            let xL = meshCell.left
            let vL = v0(xL)
            mesh[xL] = vL
            
            let xR = meshCell.right
            let vR = v0(xR)
            mesh[xR] = vR
            
            let xM = meshCell.value
            let v = 0.5 * (vL + vR)
            mesh[xM] = v
            eigenCell[t]?[meshCell] = Eigen.system(physical: v, gamma: gamma)
        }
        
        var iterations = 0
        
        while t < T {
            guard let meshP = solution[t] else { assertionFailure("Time iteration error"); return }
            
            let tau = tau(average: Set(space.nodes().compactMap{ meshP[safe: $0] }))
            
            guard tau.isNormal && tau > 0 else { fatalError() }
            
            let tP = t
            t += tau
            solution[t] = [:]
            eigenCell[t] = [:]
            
            print("-------\(t)-------")
            
            space.nodes().enumerated().forEach { i, node in
                let cellP = GridCell(value: space.node(for: i - 1), grid: space)
                let cell = GridCell(value: node, grid: space)
                let cellN = GridCell(value: space.node(for: i + 1), grid: space)
                
                guard let v = solution[tP]![safe: cell.value] else { fatalError() }
                
                let (vL, flowL) = calculateFlowAndNewBoundaryV(
                    tP: tP, tau: tau,
                    cellL: cellP,
                    cellR: cell
                )
                
                let (vR, flowR) = calculateFlowAndNewBoundaryV(
                    tP: tP, tau: tau,
                    cellL: cell,
                    cellR: cellN
                )
                
                solution[t]?[cell.left] = vL
                solution[t]?[cell.right] = vR
                
                let u = U(physical: v, gamma: gamma)
                let newU = U(vector: u.values - tau / space.step * (flowR - flowL).values)
                let newV = V(conservative: newU, gamma: gamma)
                                
                guard newV.values.allSatisfy({ $0.isFinite }) else { fatalError() }
                
                eigenCell[t]?[cell] = Eigen.system(physical: newV, gamma: gamma)
                solution[t]?[cell.value] = newV
            }
                        
            iterations += 1
            
            if iterations == 2 {
                return
            }
        }
    }
    
    // MARK: - Private methods
    
    private func tau(average vSet: Set<V>) -> Double {
        let maxLamba = vSet.map {
            let c = Constants.speedOfSound(physical: $0, gamma: gamma)
            return abs($0.speed + c)
        }.max()!
        
        return Constants.sigmaGas * space.step / maxLamba
    }
    
    private func calculateFlowAndNewBoundaryV(
        tP: Double,
        tau: Double,
        cellL: GridCell,
        cellR: GridCell
    ) -> (V, F) {
        guard let (VL, Vl) = leftSidePair(tP: tP, tau: tau, cell: cellL) else {
            let (VR, Vr) = rightSidePair(tP: tP, tau: tau, cell: cellR)!
            
            return (VR, F(physical: Vr, gamma: gamma))
        }
        
        guard let (VR, Vr) = rightSidePair(tP: tP, tau: tau, cell: cellR) else {
            return (VL, F(physical: Vl, gamma: gamma))
        }
        
        let vStar = vStar(vL: VL, vR: VR)
        let systemVStar = Eigen.system(physical: vStar, gamma: gamma)
                
        var newBoundaryV = 0.5 * (VL + VR)
        var flow = 0.5 * (F(physical: Vl, gamma: gamma) + F(physical: Vr, gamma: gamma))
        
        let M = Eigen.M(physical: vStar, gamma: gamma)
        
        systemVStar.forEach { system in
            guard !system.lambda.isZero else { return }
            
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
        tP: Double,
        tau: Double,
        cell: GridCell
    ) -> (forBoundary: V, forFlow: V)? {
        guard
            let systemP = eigenCell[tP]?[cell],
            systemP.contains(where: { $0.lambda > 0 }),
            let maxSystemP = systemP.last,
            maxSystemP.lambda > 0
        else {
            return nil
        }
        
        guard
            let vL = solution[tP]?[safe: cell.left],
            let vM = solution[tP]?[safe: cell.value],
            let vR = solution[tP]?[safe: cell.right]
        else {
            let vl = v0(space.start)
            
            return (vl, vl)
        }
        
        /// for boundary
        let VL = systemP.reduce(into: V.zero) { result, system in
            guard system.lambda > 0 else { return () }
            
            let x = shifted(x: cell.right, lambda: system.lambda, tau: tau)
                                    
            result += Nv(x: x, xL: cell.left, vL: vL, vM: vM, vR: vR)
        }
        
        /// for flow
        let averageVL = averageBoundaryVp(
            tau: tau, system: maxSystemP, cell: cell,
            vL: vL, vM: vM, vR: vR,
            boundaryX: cell.right, boundaryV: vR
        )
        
        let Vl = averageVL + systemP.reduce(into: V.zero) { result, system in
            guard system.lambda > 0 else { return () }
            
            let averageVp = averageBoundaryVp(
                tau: tau, system: system, cell: cell,
                vL: vL, vM: vM, vR: vR,
                boundaryX: cell.right, boundaryV: vR
            )

            result += V(vector: system.r * (system.l * (averageVp - averageVL)))
        }
                        
        return (VL, Vl)
    }
    
    private func rightSidePair(
        tP: Double,
        tau: Double,
        cell: GridCell
    ) -> (forBoundary: V, forFlow: V)? {
        guard
            let systemP = eigenCell[tP]?[cell],
            systemP.contains(where: { $0.lambda < 0 }),
            let minSystemP = systemP.first,
            minSystemP.lambda < 0
        else {
            return nil
        }
        
        guard
            let vL = solution[tP]?[safe: cell.left],
            let vM = solution[tP]?[safe: cell.value],
            let vR = solution[tP]?[safe: cell.right]
        else {
            let vr = v0(space.end)
            
            return (vr, vr)
        }
        
        /// for boundary
        let VR = systemP.reduce(into: V.zero) { result, system in
            guard system.lambda < 0 else { return () }
            
            let x = shifted(x: cell.left, lambda: system.lambda, tau: tau)
            result += Nv(x: x, xL: cell.left, vL: vL, vM: vM, vR: vR)
        }
    
        /// for flow
        let averageVR = averageBoundaryVp(
            tau: tau, system: minSystemP, cell: cell,
            vL: vL, vM: vM, vR: vR,
            boundaryX: cell.left, boundaryV: vL
        )

        let Vr = averageVR - systemP.reduce(into: V.zero) { result, system in
            guard system.lambda < 0 else { return () }
            
            let averageVp = averageBoundaryVp(
                tau: tau, system: system, cell: cell,
                vL: vL, vM: vM, vR: vR,
                boundaryX: cell.left, boundaryV: vL
            )
            
            result += V(vector: system.r * ((averageVp - averageVR) * system.l))
        }
                                
        return (VR, Vr)
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
        tau: Double,
        system: Eigen.System,
        cell: GridCell,
        vL: V,
        vM: V,
        vR: V,
        boundaryX: Double,
        boundaryV: V
    ) -> V {
        guard !system.lambda.isZero else {
            return boundaryV
        }

        let delta = system.lambda * tau
        let shiftedX = shifted(x: boundaryX, lambda: system.lambda, tau: tau)
        
        let v1 = Nv(x: shiftedX, xL: cell.left, vL: vL, vM: vM, vR: vR)
        let v2 = Nv(x: shiftedX + 0.5 * delta, xL: cell.left, vL: vL, vM: vM, vR: vR)
        let v3 = Nv(x: shiftedX + delta, xL: cell.left, vL: vL, vM: vM, vR: vR)
        
        return 1 / 6 * (v1 + 4 * v2 + v3)
    }
    
    private func Nv(
        t: Double,
        x: Double,
        cell: GridCell
    ) -> V {
        guard
            let vL = solution[t]?[safe: cell.left],
            let vM  = solution[t]?[safe: cell.value],
            let vR = solution[t]?[safe: cell.right]
        else { return .zero }
        
        return Nv(x: x, xL: cell.left, vL: vL, vM: vM, vR: vR)
    }
    
    private func Nv(
        x: Double,
        xL: Double,
        vL: V,
        vM: V,
        vR: V
    ) -> V {
        let xi = xi(x: x, xL: xL)
        
        return V(
            density: NVar(x: x, xL: xL, xi: xi, varL: vL.density, varM: vM.density, varR: vR.density),
            speed: NVar(x: x, xL: xL, xi: xi, varL: vL.speed, varM: vM.speed, varR: vR.speed),
            pressure: NVar(x: x, xL: xL, xi: xi, varL: vL.pressure, varM: vM.pressure, varR: vR.pressure)
        )
    }
    
    private func NVar(
        x: Double,
        xL: Double,
        varL: Double,
        varM: Double,
        varR: Double
    ) -> Double {
        let xi = xi(x: x, xL: xL)
        
        return NVar(x: x, xL: xL, xi: xi, varL: varL, varM: varM, varR: varR)
    }
    
    private func NVar(
        x: Double,
        xL: Double,
        xi: Double,
        varL: Double,
        varM: Double,
        varR: Double
    ) -> Double {
        let delta = delta(varL: varL, varR: varR)
        let sixth = sixth(varL: varL, varM: varM, varR: varR)
        
        return NVar(x: x, varL: varL, xi: xi, delta: delta, sixth: sixth)
    }
    
    private func NVar(
        x: Double,
        varL: Double,
        xi: Double,
        delta: Double,
        sixth: Double
    ) -> Double {
        return varL +  xi * (delta + sixth * (1.0 - xi))
    }
    
    private func xi(x: Double, xL: Double) -> Double {
        return (x - xL) / space.step
    }
    
    private func delta(varL: Double, varR: Double) -> Double {
        return varR - varL
    }
    
    private func sixth(varL: Double, varM: Double, varR: Double) -> Double {
        return 6.0 * (varM - 0.5 * (varR + varL))
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
            return pressure * gamma / (density * (gamma - 1)) + 0.5 * speed * speed
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
        var json = "{\n"
        
        solution.keys.sorted().forEach { t in
            var mesh = "\t\"\(t)\" : {\n"
            
            space.nodes().reduce(into: [Double: Double]()) { mesh, node in
               let cell = GridCell(value: node, grid: space)
                
                guard
                    let vL = solution[t]?[safe: cell.left],
                    let v  = solution[t]?[safe: cell.value],
                    let vR = solution[t]?[safe: cell.right]
                else { return }
                
                mesh[cell.left] = vL.density
                mesh[cell.value] = v.density
                mesh[cell.right] = vR.density
                
                let step = space.step / 4
                for k in stride(from: cell.left + step, to: cell.right, by: step) {
                    let vK = Nv(x: k, xL: cell.left, vL: vL, vM: v, vR: vR)
                    mesh[k] = vK.density
                }
            }.sorted(by: { $0.key < $1.key }).forEach { x, y in
                mesh.append("\t\t\"\(x)\": \"\(y)\",\n")
            }
            mesh.removeLast(2)
            mesh.append("\n\t},\n")
            json.append("\(mesh)")
        }
        json.removeLast(2)
        json.append("\n}")
        
        let data = json.data(using: .utf8)
        save(file: "density", data: data, meta: metadata)
    }
}
