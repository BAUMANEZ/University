//
//  Grid.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 02.06.2022.
//

import Foundation

final class Grid: Equatable {
    
    static let empty = Grid(start: 0.0, end: 0.0, step: 0.0)
    
    static func == (lhs: Grid, rhs: Grid) -> Bool {
        return lhs.start == rhs.start && lhs.end == rhs.end && lhs.steps == rhs.steps
    }
    
    let start: Double
    let end  : Double
    let step : Double
    let steps: Int
    
    final func range(starting: Int? = nil, ending: Int? = nil) -> ClosedRange<Int> {
        let starting = starting ?? 0
        var ending = ending ?? steps
        if ending < 0 { ending = steps-ending }
        
        guard ending >= starting else { return Int.max...Int.max }
        
        return starting...ending
    }
    
    final func nodes(starting: Int? = nil, ending: Int? = nil) -> [Double] {
        return range(starting: starting, ending: ending).compactMap{ self.node(for: $0) }
    }
    
    final var halfed: Double {
        return step/2
    }
    
    final var json: [String: String] {
        return [
            "a": String(start),
            "b": String(end),
            "h": String(step),
            "n": String(steps),
        ]
    }
    
    init(start: Double, end: Double, step: Double) {
        self.start = start
        self.end = end
        self.step = step
        self.steps = max(0, Int((end-start)/step) - 1)
    }
    
    init(start: Double, end: Double, steps: Int) {
        self.start = start
        self.end = end
        self.step = (end - start) / Double(steps - 1)
        self.steps = steps
    }
    
    final func node(for i: Int) -> Double {
        return start+Double(i)*step
    }
}
