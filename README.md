# Rocket Flight Simulator

A C++ program that simulates a model rocket's vertical flight, including realistic 
physics and an autonomous control system that corrects for crosswind drift.

## Purpose

This program predicts how a model rocket will behave during flight; altitude, 
velocity, and stability, using real physics calculations. It's designed to test 
rocket flight behavior computationally before ever building or launching a physical 
rocket, the same way real aerospace engineers use simulation to validate designs.

## What it does

- Simulates the full flight from launch to landing using RK4 numerical integration
- Models thrust, gravity, and altitude dependent air drag
- Accounts for changing rocket mass as propellant burns during the engine's burn phase
- Simulates a crosswind disturbance pushing the rocket off course
- Implements a PID controller that acts as the rocket's automatic fin correction, 
  steering it back toward a stable trajectory
- Accepts user input for rocket mass, thrust, burn time, drag coefficient, 
  cross-sectional area, wind strength, and PID gains (with sensible defaults)
- Logs flight data at 1-millisecond resolution to a CSV file
- Automatically generates a 4 panel graph (altitude, velocity, mass, and 
  drift/correction over time) using the included Python script

## Files

- `rocket_sim.cpp` - the main simulation program (C++17)
- `plot_results.py` - generates graphs from the simulation's CSV output

## How to run it

**Requirements:** a C++17 compiler (e.g. g++) and Python 3 with `pandas` and 
`matplotlib` installed (`pip install pandas matplotlib`).
g++ -std=c++17 -o rocket_sim.exe rocket_sim.cpp
./rocket_sim.exe


The program will prompt for input values (press Enter to accept the shown default), 
then run the simulation, print a flight summary, and automatically generate the graph.

## Example results

With default values: max altitude ~145.6 m, max velocity ~53.7 m/s, max horizontal 
drift ~1.76 m, total flight time ~12.1 s.

## Known limitation

The PID controller's integral term causes the fin correction to settle at a small 
nonzero value rather than returning fully to zero, even after drift has mostly 
subsided. This is a common real-world PID tuning issue, reducing the integral 
gain (Ki) would reduce this lingering offset.
