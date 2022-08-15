//
//  Algorithm.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

public class Algorithm1D: Algorithm {
    public typealias Time = Double
    public typealias Mesh = Dictionary<Double, Double>
    
    public let space: Grid
    
    public init(a: Double, b: Double, h: Double, tau: Double, deadline: Double) {
        self.space = Grid(start: a, end: b, step: h)
        super.init(tau: tau, deadline: deadline)
    }
    
    public final func data(for solutions: [Time: Mesh]) -> Data? {
        let json = solutions.reduce(into: [String: Any]()) { json, solution in
            let mesh = solution.value.reduce(into: [String: String]()) { mesh, pair in
                let x = pair.key
                let y = pair.value
                mesh[String(x)] = String(y)
            }
            if !mesh.isEmpty {
                json[String(solution.key)] = mesh
            }
        }
        return try? JSONSerialization.data(withJSONObject: json, options: [.sortedKeys, .prettyPrinted])
    }
    
    public func f(x: Double, t: Double) -> Double? {
        return nil
    }
}
