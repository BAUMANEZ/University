//
//  Node.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 07.06.2022.
//

import Foundation

public struct Node: Hashable, Comparable {
    public let value: Double
    public let side: Side
    
    public static func < (lhs: Node, rhs: Node) -> Bool {
        guard Swift.abs(lhs.value-rhs.value) > Double.leastNormalMagnitude else {
            return lhs.side < rhs.side
        }
        return lhs.value < rhs.value
    }
    
    public enum Side: Int, Hashable, Comparable {
        case left
        case middle
        case right
        
        public static func < (lhs: Node.Side, rhs: Node.Side) -> Bool {
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
}
