//
//  Node.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 07.06.2022.
//

import Foundation

struct BoundaryValue: Hashable, Comparable {
    
    // MARK: - Type properties
    
    static let zero = BoundaryValue(value: .zero, side: .middle)
    
    //MARK: - Nested types
    
    enum Side: Int, Hashable, Comparable {
        case left
        case middle
        case right
        
        static func < (lhs: BoundaryValue.Side, rhs: BoundaryValue.Side) -> Bool {
            switch lhs {
            case .left:
                return rhs != .left
            case .middle :
                return rhs == .right
            case .right:
                return false
            }
        }
    }
    
    //MARK: - Static functions
    
    static func < (lhs: BoundaryValue, rhs: BoundaryValue) -> Bool {
        guard Swift.abs(lhs.value-rhs.value) > Double.leastNormalMagnitude else {
            return lhs.side < rhs.side
        }
        
        return lhs.value < rhs.value
    }
    static func == (lhs: BoundaryValue, rhs: BoundaryValue) -> Bool {
        return lhs.side == rhs.side && Swift.abs(lhs.value-rhs.value) < .ulpOfOne
    }
    
    //MARK: - Properties
    
    let value: Double
    let side: Side
}

extension Dictionary where Key == BoundaryValue {
    
    subscript(safe key: BoundaryValue) -> Value? {
        get {
            guard let key = self.first(where: { $0.key == key })?.key else { return nil }
            
            return self[key]
        }
        set {
            guard let key = self.first(where: { $0.key == key })?.key else { return }
            
            self[key] = newValue
        }
    }
}
