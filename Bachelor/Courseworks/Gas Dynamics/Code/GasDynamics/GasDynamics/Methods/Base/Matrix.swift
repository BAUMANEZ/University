//
//  Matrix.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 09.11.2022.
//

import Foundation


struct Matrix: CustomStringConvertible,
               CustomDebugStringConvertible {
    // MARK: - Operator
    
    static func *(left: Matrix, right: [Double]) -> [Double] {
        return (0 ..< left.m).reduce(into: [Double]()) { vector, i in
            var result: Double = 0
            
            for j in stride(from: 0, to: right.count, by: 1) {
                result += left[i, j] * right[j]
            }
            
            vector.append(result)
        }
    }
    
    // MARK: - Internal properties
    
    fileprivate(set) var n: Int
    fileprivate(set) var m: Int
    
    var description: String {
        return """
        MATRIX, dim=(\(n), \(m))
        \(values.description)
        """
    }
    
    var debugDescription: String {
        return description
    }
    
    // MARK: - Private properties
    
    fileprivate var values: [Double]
    
    // MARK: - Initializers
    
    init(n: Int = 0, initial: Double = 0.0) {
        self.init(n: n, m: n, initial: initial)
    }
    
    init(n: Int, m: Int, initial: Double = 0.0) {
        self.n = n
        self.m = m
        self.values = [Double](repeating: initial, count: n*m)
    }
    
    // MARK: - Operators
    
    subscript(i: Int, j: Int) -> Double {
        get {
            assert(i < n, "Wrong row")
            assert(j < m, "Wrong col")
            
            return values[i*m+j]
        }
        mutating set {
            assert(i < n, "Wrong row")
            assert(j < m, "Wrong col")
            
            return values[i*m+j] = newValue
        }
    }
    
    // MARK: - Internal methods
    
    func row(at i: Int) -> [Double] {
        let row = i * m
        return Array(values[row ..< row + m ])
    }
    
    func col(at j: Int) -> [Double] {
        var result = [Double](repeating: 0, count: n)
        
        for i in 0 ..< n {
            result[i] = self[i, j]
        }
        
        return result
    }
}

extension Array where Element == Double {
    static func *(left: Self, right: Self) -> Double {
        var result: Double = 0
        
        for i in stride(from: 0, to: left.count, by: 1) {
            result += left[i] * right[i]
        }
        
        return result
    }
    
    static func *(left: Self, right: Double) -> [Double] {
        return left.map { $0 * right }
    }
    
    static func *(left: Double, right: Self) -> [Double] {
        return right * left
    }
    
    static func +=(left: inout Self, right: Self) {
        right.enumerated().forEach { left[$0.offset] += $0.element }
    }
    
    static func +(left: Self, right: Self) -> Self {
        var result = Array(repeating: 0, count: left.count)
        
        left.indices.forEach { i in
            result[i] = left[i] + right[i]
        }
        
        return result
    }
    
    static func -(left: Self, right: Self) -> Self {
        var result = Array(repeating: 0, count: left.count)
        
        left.indices.forEach { i in
            result[i] = left[i] - right[i]
        }
        
        return result
    }
}
