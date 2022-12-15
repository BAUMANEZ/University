//
//  Tests.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 15.12.2022.
//

import Foundation

extension Gas1D {
    
    struct Test {
        
        let description: String?
        let a: Double
        let b: Double
        let T: Time
        let N: Int
        let initial: (Double) -> V
    }
    
    enum Tests {
        
        static let sod = Test(
            description: "Sod",
            a: -1,
            b: 1,
            T: 0.4,
            N: 100,
            initial: { x in
                switch x {
                case -1 ... 0:
                    return V(density: 1, speed: 0, pressure: 1)
                    
                case 0 ... 1:
                    return V(density: 0.125, speed: 0, pressure: 0.1)
                    
                default:
                    return .zero
                }
            }
        )
        
        static let lux = Test(
            description: "Lux",
            a: -1,
            b: 1,
            T: 0.32,
            N: 100,
            initial: { x in
                switch x {
                case -1 ... 0:
                    return V(density: 0.445, speed: 0.698, pressure: 3.528)
                    
                case 0 ... 1:
                    return V(density: 0.5, speed: 0, pressure: 0.571)
                    
                default:
                    return .zero
                }
            }
        )
        
        static let shu = Test(
            description: "Shu",
            a: -1,
            b: 1,
            T: 0.36,
            N: 100,
            initial: { x in
                switch x {
                case -1 ... -0.8:
                    return V(density: 3.857143, speed: 2.629369, pressure: 10.3333)
                    
                case -0.8 ... 1:
                    return V(density: 1 + 0.2 * sin(5 * .pi * x), speed: 0, pressure: 1)
                    
                default:
                    return .zero
                }
            }
        )
        
        static let shockWave = Test(
            description: "Shock Wave",
            a: 0,
            b: 1,
            T: 0.038,
            N: 400,
            initial: { x in
                switch x {
                case 0 ..< 0.1:
                    return V(density: 1, speed: 0, pressure: 10e3)
                    
                case 0.1 ..< 0.9:
                    return V(density: 1, speed: 0, pressure: 10e-2)
                    
                case 0.9 ... 1:
                    return V(density: 1, speed: 0, pressure: 10e2)
                    
                default:
                    return .zero
                }
            }
        )
    }
}
