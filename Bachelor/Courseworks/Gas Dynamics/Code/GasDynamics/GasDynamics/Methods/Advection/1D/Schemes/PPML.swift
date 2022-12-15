//
//  PPML.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 02.06.2022.
//

import Foundation

extension Advection1D {
    final class PPML: InterpolatedAdvection1D {
        override var identifier: String? {
            return "advectionPPML"
        }
        
        override func solve(for j: Int) {
            super.solve(for: j)
            
            let t = time.node(for: j)
            let jP = j - 1
            let tP = time.node(for: jP)
            
            guard solution[t] != nil, solution[tP] != nil else { return }
            
            space.nodes().forEach { node in
                let x  = BoundaryValue(value: node, side: .middle)
                let xL = BoundaryValue(value: node-space.halfed, side: .right)
                let xR = BoundaryValue(value: node+space.halfed, side: .left)
                let xP  = BoundaryValue(value: node-space.step, side: .middle)
                let xLP = BoundaryValue(value: node-space.step-space.halfed, side: .right)
                let xRP = BoundaryValue(value: node-space.halfed, side: .left)
                
                guard let yL = Nf(t: tP, x: xL.value-c*time.step, xL: xLP, xM: xP, xR: xRP),
                      let yR = Nf(t: tP, x: xR.value-c*time.step, xL: xL, xM: x, xR: xR),
                      let y  = solution[t]?[x]
                else { return }
                
                guard (yR-y)*(y-yL) > 0 else {
                    solution[t]?[xL] = y
                    solution[t]?[xR] = y
                    return
                }
                
                let delta = delta(yL: yL, yR: yR)
                let sixth = sixth(yL: yL, y: y, yR: yR)
                let deltaSix = delta*sixth
                let deltaSq = delta*delta
                solution[t]?[xL] = deltaSix > deltaSq ? (3.0*y-2.0*yR) : yL
                solution[t]?[xR] = deltaSix < -deltaSq ? (3.0*y-2.0*yL) : yR
            }
        }
    }
}
