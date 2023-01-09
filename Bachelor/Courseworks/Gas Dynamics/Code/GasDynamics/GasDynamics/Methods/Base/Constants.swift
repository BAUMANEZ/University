import Foundation

typealias Time = Double

enum Constants {
    
    // MARK: - Advection
    
    static let c: Double = 1.0
    static let h: Double = 1.0/2.0
    static let sigmaAdvection: Double = 1.0/1.0
    static let profile: Advection1D.Profile = .cos
    
    // MARK: - Gas
    
    static let gamma: Double = 1.4
    static let sigmaGas: Double = 0.5
    
    static func speedOfSound(density: Double, pressure: Double, gamma: Double = Constants.gamma) -> Double {
        return sqrt(gamma * pressure / density)
    }
    
    static func speedOfSound(physical v: Gas1D.V, gamma: Double = Constants.gamma) -> Double {
        return speedOfSound(density: v.density, pressure: v.pressure, gamma: gamma)
    }
}
