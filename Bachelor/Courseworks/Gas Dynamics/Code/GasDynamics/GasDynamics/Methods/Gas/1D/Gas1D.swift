//
//  Gas1D.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 31.10.2022.
//

import Foundation

final class Gas1D: Algorithm1D {
    
    // MARK: - Nested types
    
    struct V {
        
        let density : Double
        let speed   : Double
        let pressure: Double
    }
    
    struct U {
        
        /// - rho
        let u1: Double
        
        /// - rho * u
        let u2: Double
        
        /// - E
        let u3 : Double
        
        public init(conservative: V, energyDensity: Double) {
            self.init(
                density: conservative.density,
                speed: conservative.speed,
                pressure: conservative.pressure,
                energyDensity: energyDensity
            )
        }
        
        public init(density: Double, speed: Double, pressure: Double, energyDensity: Double) {
            u1 = density
            u2 = density * speed
            u3 = density*energyDensity + 0.5*density*pow(speed, 2)
        }
    }
    
    struct F {
        
        /// - rho * u
        let f1: Double
        
        /// - rho * u^2 + p
        let f2: Double
        
        /// - (E + p) * u
        let f3: Double
        
        public init(density: Double, speed: Double, energy: Double, pressure: Double) {
            f1 = density*speed
            f2 = density*pow(speed, 2) + pressure
            f3 = (energy+pressure)*speed
        }
    }
    
    // MARK: - Internal properties
    
//    let flow: (Double, Time) -> F
//    let u   : (Double, Time) -> U
    
    // MARK: - Internal methods
    
}
