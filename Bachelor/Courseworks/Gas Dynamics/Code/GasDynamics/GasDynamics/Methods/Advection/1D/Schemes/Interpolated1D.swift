//
//  Interpolated1D.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 07.06.2022.
//

import Foundation

class InterpolatedAdvection1D: Advection1D {
            
    var identifier: String? {
        return nil
    }
    
    private final var substep: Double {
        return 100.0
    }
    
    private final var plotStep: Int {
        return time.steps/5
    }
    
    private final var metadata: Data? {
        let json: [String: Any] = [
            "type"   : "Advection",
            "scheme" : identifier ?? "none",
            "profile": profile ?? "unknown",
            "grid"   : [
                "time" : time.json,
                "space": space.json
            ],
            "speed": c,
            "norms"  : [
                "C"  : String(normC),
                "L"  : nil,//String(normL),
                "L2" : String(normL2),
                "W2" : nil
            ]
        ]
        return try? JSONSerialization.data(withJSONObject: json, options: [.prettyPrinted, .sortedKeys])
    }
    
    final var normC: Double {
        var result = -Double.greatestFiniteMagnitude
        time.nodes(starting: 1).forEach { t in
            guard let mesh = solution[t] else { return }
            
            space.nodes().forEach { node in
                let xL = BoundaryValue(value: node-space.halfed, side: .right)
                let xM = BoundaryValue(value: node, side: .middle)
                let xR = BoundaryValue(value: node+space.halfed, side: .left)
                guard let y = mesh[xM],
                      let yL = mesh[xL],
                      let yR = mesh[xR]
                else { return }
                let delta = yR-yL
                let sixth = 6.0*(y-0.5*(yR+yL))
                for k in stride(from: xL.value, to: xR.value, by: (xR.value-xL.value)/substep) {
                    let xi = (k-xL.value)/space.step
                    result = Swift.max(result, (u(k, t)-(yL+xi*(delta+sixth*(1.0-xi)))).magnitude)
                }
            }
        }
        return result
    }
    
    final var normL: Double {
        return time.nodes(starting: 1).enumerated().reduce(into: Double()) { global, t in
            guard
                let mesh = solution[t.element],
                let meshP = solution[t.element - 1]
            else { return }
            
            global += space.nodes().reduce(into: Double()) { local, node in
                let xL = BoundaryValue(value: node-space.halfed, side: .right)
                let x  = BoundaryValue(value: node, side: .middle)
                let xR = BoundaryValue(value: node+space.halfed, side: .left)
                
                guard
                    let y = mesh[x],
                    let yL = mesh[xL],
                    let yR = mesh[xR],
                    let yP = mesh[x],
                    let yPL = meshP[xL],
                    let yPR = meshP[xR]
                else { return }
                
                let delta = yR-yL
                let deltaP = yPR-yPL
                let sixth = 6.0*(y-0.5*(yR+yL))
                let sixthP = 6.0*(yP-0.5*(yPR+yPL))
                
                for k in stride(from: xL.value, through: xR.value, by: (xR.value-xL.value)/substep) {
                    let xi = (k-xL.value)/space.step
                    let z = Swift.abs(u(k, t.element)-(yL+xi*(delta+sixth*(1.0-xi))))
                    let zP = Swift.abs(u(k, t.element - time.step)-(yPL+xi*(deltaP+sixthP*(1.0-xi))))
                    local += space.step*time.step*(zP-z)
                }
            }
        }
    }
    
    final var normL2: Double {
        return sqrt(time.nodes(starting: 1).enumerated().reduce(into: Double()) { global, t in
            guard
                let mesh = solution[t.element],
                let meshP = solution[t.element - 1]
            else { return }
            
            global += space.nodes().reduce(into: Double()) { local, node in
                let xL = BoundaryValue(value: node-space.halfed, side: .right)
                let x  = BoundaryValue(value: node, side: .middle)
                let xR = BoundaryValue(value: node+space.halfed, side: .left)
                
                guard
                    let y = mesh[x],
                    let yL = mesh[xL],
                    let yR = mesh[xR],
                    let yP = mesh[x],
                    let yPL = meshP[xL],
                    let yPR = meshP[xR]
                else { return }
                
                let delta = yR-yL
                let deltaP = yPR-yPL
                let sixth = 6.0*(y-0.5*(yR+yL))
                let sixthP = 6.0*(yP-0.5*(yPR+yPL))
                
                for k in stride(from: xL.value, through: xR.value, by: (xR.value-xL.value)/substep) {
                    let xi = (k-xL.value)/space.step
                    let z = pow((u(k, t.element)-(yL+xi*(delta+sixth*(1.0-xi)))), 2)
                    let zP = pow(Swift.abs(u(k, t.element - time.step)-(yPL+xi*(deltaP+sixthP*(1.0-xi)))), 2)
                    local += space.step*time.step*(zP-z)
                }
            }
        })
    }
    
    final var normW2: Double {
        return .zero
    }
    
    final override func solve() {        
        solution[time.start] = space.nodes().reduce(into: Mesh()) {
            let xL = BoundaryValue(value: $1-space.halfed, side: .right)
            let x  = BoundaryValue(value: $1, side: .middle)
            let xR = BoundaryValue(value: $1+space.halfed, side: .left)
            let yL = u0(xL.value)
            let yR = u0(xR.value)
            let y  = (yL+yR)/2.0
            $0[xL] = yL
            $0[x]  = y
            $0[xR] = yR
        }
        
        guard time.steps > 1 else { return }
        
        time.range(starting: 1).forEach{ solve(for: $0) }
        
        guard let identifier else { return }
        
        var solutions: Solution = [:]
        
        for j in Swift.stride(from: 0, through: time.steps, by: plotStep) {
            let t = time.node(for: j)
            
            guard let mesh = solution[t] else { continue }
            
            solutions[t] = space.nodes().reduce(into: Mesh()) { _mesh, node in
                let xL = BoundaryValue(value: node-space.halfed, side: .right)
                let x  = BoundaryValue(value: node, side: .middle)
                let xR = BoundaryValue(value: node+space.halfed, side: .left)
                
                guard let yL = mesh[xL],
                      let y  = mesh[x],
                      let yR = mesh[xR]
                else { return }
                
                _mesh[x] = y
                _mesh[BoundaryValue(value: xL.value + 10e-11, side: .right)] = yL
                _mesh[BoundaryValue(value: xR.value - 10e-11, side: .right)] = yR
                
                let delta = delta(yL: yL, yR: yR)
                let sixth = sixth(yL: yL, y: y, yR: yR)
                let substep = (xR.value-xL.value)/substep
                
                Swift.stride(
                    from: xL.value + substep,
                    through: xR.value - substep,
                    by: substep
                ).forEach { x_k in
                    let xi = (x_k - xL.value) / space.step
                    _mesh[BoundaryValue(value: x_k, side: .middle)] = yL + xi * (delta + sixth * (1.0 - xi))
                }
            }
        }
        
        save(file: identifier, data: data(for: solutions), meta: metadata)
    }
        
    func solve(for j: Int) {
        let jP = j-1
        let t = time.node(for: j)
        let tP = time.node(for: jP)
        solution[t] = [:]
        
        space.nodes().forEach { node in
            let x  = BoundaryValue(value: node, side: .middle)
            let xL = BoundaryValue(value: node - space.halfed, side: .right)
            let xR = BoundaryValue(value: node + space.halfed, side: .left)
            
            guard let leftFlow = flow(t: tP, xB: xL),
                  let rightFlow = flow(t: tP, xB: xR),
                  let y = solution[tP]?[x]
            else { return }
            
            let f = y-gamma*(rightFlow-leftFlow)
            solution[t]?[x] = f
        }
    }
    
    private func flow(t: Double, xB: BoundaryValue) -> Double? {
        let xL: BoundaryValue
        let x : BoundaryValue
        let xR: BoundaryValue
        
        if c > 0 {
            xL = BoundaryValue(value: xB.value-space.step, side: .right)
            x  = BoundaryValue(value: xB.value-space.halfed, side: .middle)
            xR = BoundaryValue(value: xB.value, side: .left)
        } else {
            xL = BoundaryValue(value: xB.value, side: .right)
            x  = BoundaryValue(value: xB.value+space.halfed, side: .middle)
            xR = BoundaryValue(value: xB.value+space.step, side: .left)
        }
        
        guard let y  = solution[t]?[x],
              let yL = solution[t]?[xL],
              let yR = solution[t]?[xR]
        else { return .zero }
        
        let boundary = c > 0 ? xR.value : xL.value
        let f1 = c*(Nf(x: boundary-c*time.step, xL: xL, yL: yL, y: y, yR: yR))
        let f2 = c*(Nf(x: boundary-0.5*c*time.step, xL: xL, yL: yL, y: y, yR: yR))
        let f3 = c*(Nf(x: boundary, xL: xL, yL: yL, y: y, yR: yR))
        let flow = 1.0/6.0*(f1+4.0*f2+f3)
        return flow
    }
}

extension InterpolatedAdvection1D {
    
    final func delta(t: Double, xL: BoundaryValue, xR: BoundaryValue) -> Double? {
        guard let yL = solution[t]?[xL],
              let yR = solution[t]?[xR]
        else { return nil }
        
        return delta(yL: yL, yR: yR)
    }
    
    final func delta(yL: Double, yR: Double) -> Double {
        return yR-yL
    }
    
    final func sixth(t: Double, xL: BoundaryValue, x: BoundaryValue, xR: BoundaryValue) -> Double? {
        guard let y  = solution[t]?[x],
              let yL = solution[t]?[xL],
              let yR = solution[t]?[xR]
        else { return nil }
        return sixth(yL: yL, y: y, yR: yR)
    }
    
    final func sixth(yL: Double, y: Double, yR: Double) -> Double {
        return 6.0*(y-0.5*(yR+yL))
    }
    
    final func xi(x: Double, xL: BoundaryValue) -> Double {
        return (x-xL.value)/space.step
    }
    
    final func Nf(t: Double, x: Double, xL: BoundaryValue, xM: BoundaryValue, xR: BoundaryValue) -> Double? {
        guard let y  = solution[t]?[xM],
              let yL = solution[t]?[xL],
              let yR = solution[t]?[xR]
        else { return .zero }
        
        let xi = (x-xL.value)/space.step
        let delta = yR-yL
        let sixth = 6.0*(y-0.5*(yR+yL))
        
        return yL+xi*(delta+sixth*(1.0-xi))
    }
    
    final func Nf(x: Double, xL: BoundaryValue, yL: Double, y: Double, yR: Double) -> Double {
        let xi = xi(x: x, xL: xL)
        let delta = delta(yL: yL, yR: yR)
        let sixth = sixth(yL: yL, y: y, yR: yR)
        
        return yL+xi*(delta+sixth*(1-xi))
    }
}
