//
//  main.swift
//  GasDynamics
//
//  Created by Арсений Токарев on 29.05.2022.
//

import Foundation

let c: Double = 1.0
let h: Double = 1.0/2.0
let sigma: Double = 1.0/1.0
let profile: Advection1D.Profile = .cos

Advection1D.PPM(c: c, h: h, sigma: sigma, profile: profile)
Advection1D.PPML(c: c, h: h, sigma: sigma, profile: profile)
