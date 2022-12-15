//
//  PPM.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 02.06.2022.
//

import Foundation

extension Advection1D {
    final class PPM: InterpolatedAdvection1D {
        override var identifier: String? {
            return "advectionPPM"
        }
        
        override func solve(for j: Int) {
            super.solve(for: j)
            
            let t = time.node(for: j)
            
            guard solution[t] != nil else { return }
            
            space.nodes().forEach { node in
                let x   = BoundaryValue(value: node, side: .middle)
                let xMM = BoundaryValue(value: node-2.0*space.step, side: .middle)
                let xM  = BoundaryValue(value: node-space.step, side: .middle)
                let xP  = BoundaryValue(value: node+space.step, side: .middle)
                let xPP = BoundaryValue(value: node+2.0*space.step, side: .middle)
                let yMM = solution[t]?[xMM] ?? 0
                let yM  = solution[t]?[xM]  ?? 0
                let y   = solution[t]?[x]   ?? 0
                let yP  = solution[t]?[xP]  ?? 0
                let yPP = solution[t]?[xPP] ?? 0
                let xL  = BoundaryValue(value: node-space.halfed, side: .right)
                let xR  = BoundaryValue(value: node+space.halfed, side: .left)
                let yL  = 0.5*(yM+y)-(1.0/6.0)*(deltaM(yL: yM, y: y, yR: yP)-deltaM(yL: yMM, y: yM, yR: y))
                let yR  = 0.5*(y+yP)-(1.0/6.0)*(deltaM(yL: y, y: yP, yR: yPP)-deltaM(yL: yM, y: y, yR: yP))
                
                guard (yR-y)*(y-yL) > 0 else {
                    solution[t]?[xL] = y
                    solution[t]?[xR] = y
                    return
                }
                
                let delta = yR-yL
                let sixth = 6.0*(y-0.5*(yL+yR))
                let deltaSix = delta*sixth
                let deltaSq = delta*delta
                solution[t]?[xL] = deltaSix > deltaSq ? (3.0*y-2.0*yR) : yL
                solution[t]?[xR] = deltaSix < -deltaSq ? (3.0*y-2.0*yL) : yR
            }
        }
        
        private func deltaM(yL: Double, y: Double, yR: Double) -> Double {
            let deltaM = 0.5*(yL+yR)
            guard ((yR-y)*(y-yL)).sign > 0 else { return .zero }
            let min1 = 2.0*(yR-y)
            let min2 = 2.0*(y-yL)
            return deltaM.sign*min(Swift.abs(deltaM), Swift.abs(min1), Swift.abs(min2))
        }
    }
}

extension Double {
    var sign: Double {
        return self >= .zero ? 1.0 : -1.0
    }
}
