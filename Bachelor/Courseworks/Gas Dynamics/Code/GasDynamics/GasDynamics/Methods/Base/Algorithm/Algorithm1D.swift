//
//  Algorithm.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

 class Algorithm1D: Algorithm {
     
     typealias Mesh = [BoundaryValue: Double]
     typealias Solution = [Time: Mesh]
     
     let space: Grid
     var solution: Solution = [:]
    
     init(a: Double, b: Double, h: Double, tau: Double, deadline: Double) {
        self.space = Grid(start: a, end: b, step: h)
        super.init(tau: tau, deadline: deadline)
    }
    
     final func data(for solution: Solution) -> Data? {
        let json = solution.reduce(into: [String: Any]()) { json, solution in
            let mesh = solution.value.reduce(into: [String: String]()) { mesh, pair in
                let x = pair.key.value
                let y = pair.value
                mesh[String(x)] = String(y)
            }
            if !mesh.isEmpty {
                json[String(solution.key)] = mesh
            }
        }
        return try? JSONSerialization.data(withJSONObject: json, options: [.prettyPrinted])
    }
    
     func f(x: Double, t: Double) -> Double? {
        return nil
    }
}
