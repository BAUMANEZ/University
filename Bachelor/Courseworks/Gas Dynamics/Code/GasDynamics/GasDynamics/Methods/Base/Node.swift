//
//  Node.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 07.06.2022.
//

import Foundation

 struct Node: Hashable, Comparable {
    
    //MARK: - Nested types
    
     enum Side: Int, Hashable, Comparable {
        case left
        case middle
        case right
        
         static func < (lhs: Node.Side, rhs: Node.Side) -> Bool {
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
    
     static func < (lhs: Node, rhs: Node) -> Bool {
        guard Swift.abs(lhs.value-rhs.value) > Double.leastNormalMagnitude else {
            return lhs.side < rhs.side
        }
        return lhs.value < rhs.value
    }
     static func == (lhs: Node, rhs: Node) -> Bool {
        return lhs.side == rhs.side && Swift.abs(lhs.value-rhs.value) < Double.leastNormalMagnitude
    }
    
    //MARK: - Properties
    
     let value: Double
     let side: Side
}
