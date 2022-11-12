import Foundation

func speedOfSound(density: Double, pressure: Double, gamma: Double) -> Double {
    return sqrt(gamma*pressure/density)
}
