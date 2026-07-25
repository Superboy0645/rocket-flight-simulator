#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
// Constants and model-rocket parameters
// TIME_STEP is shared by the RK4 integrator and CSV logger. A value of
// 0.001 seconds gives one state calculation and one data row per millisecond.
constexpr double TIME_STEP = 0.001;            // Time increment (s)
// Standard gravitational acceleration near Earth's surface. It is treated as
// constant because a model rocket's altitude is tiny relative to Earth.
constexpr double GRAVITY = 9.80665;            // Acceleration (m/s^2)
// Reference air density at zero altitude in the exponential atmosphere model.
constexpr double SEA_LEVEL_DENSITY = 1.225;    // Air density (kg/m^3)
// Controls how quickly density falls with altitude according to
// rho = rho_0 * exp(-altitude / scale height).
constexpr double ATMOSPHERE_SCALE_HEIGHT = 8500.0; // Scale height (m)
// The body diameter is used to calculate the default circular frontal area.
// The user may replace that area directly during program setup.
constexpr double ROCKET_DIAMETER = 0.08;        // Body diameter (m)
constexpr double PI = 3.14159265358979323846;   // Dimensionless constant
// Frontal area A = pi*d^2/4, which appears in the aerodynamic drag equation.
constexpr double DEFAULT_CROSS_SECTIONAL_AREA =
    PI * ROCKET_DIAMETER * ROCKET_DIAMETER / 4.0; // Frontal area (m^2)
// Converts commanded fin deflection into an equivalent lateral force. This is
// a simplified control-effectiveness model rather than a full fin flow model.
constexpr double FIN_FORCE_PER_RADIAN = 4.0;      // N/rad
// Physical fin travel limit. PID commands are clipped to this magnitude.
constexpr double MAX_FIN_DEFLECTION_DEG = 12.0;   // Mechanical fin limit
constexpr double MAX_FIN_DEFLECTION_RAD =
    MAX_FIN_DEFLECTION_DEG * PI / 180.0;           // Limit (rad)
// Caps accumulated position error to prevent integral "windup" while a fin
// command is limited at its maximum physical deflection.
constexpr double MAX_INTEGRAL_ERROR = 20.0; // Error integral limit (m*s)
// These example values describe a small model rocket. They appear in the input
// prompts; pressing Enter keeps a default for the current simulation run.
struct SimulationParameters {
    double initialMass = 0.75;        // Liftoff mass including fuel (kg)
    double propellantMass = 0.25;     // Fuel consumed during burn (kg)
    double thrust = 25.0;             // Constant motor force in burn (N)
    double burnTime = 2.0;            // Powered-flight duration (s)
    double dragCoefficient = 0.75;    // Dimensionless drag coefficient
    double crossSectionalArea = DEFAULT_CROSS_SECTIONAL_AREA; // Area (m^2)
    double crosswindForce = 0.30;      // Lateral disturbance force (N)
    double pidKp = 0.040;              // Proportional gain (rad/m)
    double pidKi = 0.005;              // Integral gain (rad/(m*s))
    double pidKd = 0.080;              // Derivative gain (rad*s/m)
    // Dry mass is the structure, casing, electronics, and hardware that remain
    // after the modeled propellant has been consumed.
    double dryMass() const {
        return initialMass - propellantMass;
    }
};
// State contains the four variables advanced by RK4. Upward and the selected
// crosswind direction are positive; zero horizontal position is the target.
struct State {
    double altitude;           // m
    double velocity;           // m/s; upward is positive
    double horizontalPosition; // m; intended flight line is position zero
    double horizontalVelocity; // m/s
};
// A Derivative stores the instantaneous rates of change returned by the
// equations of motion. RK4 combines four such estimates per timestep.
struct Derivative {
    double altitudeRate;           // m/s
    double velocityRate;           // m/s^2
    double horizontalPositionRate; // m/s
    double horizontalVelocityRate; // m/s^2
};
// The PID controller must retain accumulated error and the previous error
// between samples to calculate its integral and derivative actions.
struct PIDState {
    double integralError;
    double previousError;
};
// This structure keeps each controller contribution for later CSV analysis,
// along with the limited fin command and its modeled lateral force.
struct ControlOutput {
    double error;             // Current horizontal position error (m)
    double proportionalTerm; // P contribution to commanded angle (rad)
    double integralTerm;     // I contribution to commanded angle (rad)
    double derivativeTerm;   // D contribution to commanded angle (rad)
    double finAngle;          // Limited fin command (rad)
    double finForce;          // Resulting horizontal correction force (N)
};
// Interactive input handling
// This input layer lets the same executable represent different rockets and
// controller tunings without recompilation. Every prompt displays the default
// stored in SimulationParameters. Pressing Enter supplies a blank line and
// accepts that default; otherwise the text must contain one finite number.
// The Validator template argument lets each parameter apply its own physical
// rule, such as positive initial mass or nonnegative thrust. Invalid text or an
// out-of-range value prints the supplied explanation and repeats the prompt.
template <typename Validator>
double promptForValue(const std::string& label,
                      double defaultValue,
                      Validator isValid,
                      const std::string& validationMessage) {
    // Continue until this one parameter receives a valid value or its default.
    while (true) {
        std::cout << label << " [default " << std::setprecision(8)
                  << defaultValue << "]: ";
        std::string line;
        // End-of-file can occur when input is redirected; using the default
        // keeps that behavior equivalent to pressing Enter.
        if (!std::getline(std::cin, line)) {
            std::cout << "\nInput ended; using the default value.\n";
            return defaultValue;
        }
        // A line containing no visible characters means "accept the default."
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) {
            return defaultValue;
        }
        // Parse exactly one number. Any trailing non-whitespace text, NaN,
        // infinity, or value rejected by isValid causes a new prompt.
        std::istringstream input(line);
        double value = 0.0;
        if (input >> value) {
            input >> std::ws;
            if (input.eof() && std::isfinite(value) && isValid(value)) {
                return value;
            }
        }
        std::cerr << "Invalid input: " << validationMessage << '\n';
    }
}
SimulationParameters readSimulationParameters() {
    // Start with all documented defaults, then replace only the values for
    // which the user enters an explicit number.
    SimulationParameters parameters;
    std::cout << "Rocket Flight Simulation Setup\n"
              << "Press Enter at any prompt to accept its default.\n\n";
    // Mass validation ensures the dry mass remains positive, preventing a
    // division by zero in F = m*a after burnout.
    parameters.initialMass = promptForValue(
        "Initial rocket mass (kg)",
        parameters.initialMass,
        [](double value) { return value > 0.0; },
        "initial mass must be greater than zero.");
    parameters.propellantMass = promptForValue(
        "Propellant mass (kg)",
        parameters.propellantMass,
        [&parameters](double value) {
            return value >= 0.0 && value < parameters.initialMass;
        },
        "propellant mass must be nonnegative and less than initial mass.");
    // Motor inputs may be zero to represent an unpowered or zero-duration case,
    // but negative thrust or time has no physical meaning in this model.
    parameters.thrust = promptForValue(
        "Thrust (N)",
        parameters.thrust,
        [](double value) { return value >= 0.0; },
        "thrust cannot be negative.");
    parameters.burnTime = promptForValue(
        "Burn time (s)",
        parameters.burnTime,
        [](double value) { return value >= 0.0; },
        "burn time cannot be negative.");
    // Aerodynamic area must remain positive and Cd must not be negative because
    // both multiply the standard drag equation.
    parameters.dragCoefficient = promptForValue(
        "Drag coefficient",
        parameters.dragCoefficient,
        [](double value) { return value >= 0.0; },
        "drag coefficient cannot be negative.");
    parameters.crossSectionalArea = promptForValue(
        "Cross-sectional area (m^2)",
        parameters.crossSectionalArea,
        [](double value) { return value > 0.0; },
        "cross-sectional area must be greater than zero.");
    // A signed wind force selects the disturbance direction. PID gains remain
    // unrestricted so the user can deliberately explore different tunings.
    parameters.crosswindForce = promptForValue(
        "Crosswind force (N; positive or negative)",
        parameters.crosswindForce,
        [](double) { return true; },
        "crosswind force must be a finite number.");
    parameters.pidKp = promptForValue(
        "PID gain Kp",
        parameters.pidKp,
        [](double) { return true; },
        "Kp must be a finite number.");
    parameters.pidKi = promptForValue(
        "PID gain Ki",
        parameters.pidKi,
        [](double) { return true; },
        "Ki must be a finite number.");
    parameters.pidKd = promptForValue(
        "PID gain Kd",
        parameters.pidKd,
        [](double) { return true; },
        "Kd must be a finite number.");
    std::cout << "\nDry mass after burnout: "
              << parameters.dryMass() << " kg\n\n";
    return parameters;
}
// Rocket mass falls because hot propellant leaves the vehicle during powered
// flight. With no detailed motor mass-flow curve available, the simulator
// assumes the entered propellant mass is consumed at a constant rate from
// ignition to burnout. burnFraction therefore progresses linearly from 0 to 1.
// At or after burnout, only dry mass remains and mass stays constant.
double massAt(double time, const SimulationParameters& parameters) {
    // A zero-duration burn is treated as instantaneous burnout. This branch
    // also avoids dividing by zero when calculating burnFraction.
    if (parameters.burnTime == 0.0 || time >= parameters.burnTime) {
        return parameters.dryMass();
    }
    const double burnFraction =
        std::clamp(time / parameters.burnTime, 0.0, 1.0);
    return parameters.initialMass -
           parameters.propellantMass * burnFraction;
}
// The exponential atmosphere captures the main trend that aerodynamic forces
// weaken as altitude increases and air becomes thinner. A trial RK4 state may
// fall slightly below zero near touchdown, so density is evaluated at sea-level
// conditions rather than extrapolated to an unphysical below-ground value.
double airDensityAt(double altitude) {
    const double nonnegativeAltitude = std::max(0.0, altitude);
    return SEA_LEVEL_DENSITY *
           std::exp(-nonnegativeAltitude / ATMOSPHERE_SCALE_HEIGHT);
}
// Equations of motion
// The state derivatives are:
//   d(altitude)/dt = velocity
//   d(velocity)/dt = net force / mass
//   d(horizontal position)/dt = horizontal velocity
//   d(horizontal velocity)/dt = horizontal force / mass
// The force sign convention is positive upward and positive in the selected
// crosswind direction. Drag always opposes vertical velocity. Writing the
// velocity factor as v*abs(v) gives speed squared while preserving direction:
// the leading minus sign makes drag downward during ascent and upward during
// descent.
Derivative equationsOfMotion(double time,
                             const State& state,
                             double finCorrectionForce,
                             const SimulationParameters& parameters) {
    // Evaluate mass at this exact RK4 stage time. This matters during motor
    // burn because each stage may see a slightly different remaining mass.
    const double mass = massAt(time, parameters);
    // The simplified thrust curve is a constant positive force until burnout,
    // followed by zero thrust during coast and descent.
    const double thrustForce =
        (time < parameters.burnTime) ? parameters.thrust : 0.0;
    // Weight is mass times gravitational acceleration and acts downward, so it
    // will be subtracted when the signed vertical forces are combined.
    const double gravityForce = mass * GRAVITY;
    // Standard drag equation: Fd = 0.5*rho*v^2*Cd*A. airDensityAt supplies rho,
    // and the signed v*abs(v) form makes the result oppose rocket motion.
    const double dragForce =
        -0.5 * airDensityAt(state.altitude) * state.velocity *
        std::abs(state.velocity) * parameters.dragCoefficient *
        parameters.crossSectionalArea;
    // Upward thrust, downward weight, and signed drag form the total vertical
    // force. Newton's second law below divides this force by current mass.
    const double netForce = thrustForce - gravityForce + dragForce;
    // Crosswind is modeled as a constant external horizontal force. The fin
    // correction force has the sign chosen by the PID controller, so adding
    // both produces the net lateral disturbance/correction force.
    const double horizontalForce =
        parameters.crosswindForce + finCorrectionForce;
    // Position derivatives are velocities. Velocity derivatives are the
    // corresponding accelerations a = F/m in the vertical and lateral axes.
    return {
        state.velocity,
        netForce / mass,
        state.horizontalVelocity,
        horizontalForce / mass
    };
}
// Build a trial state by advancing from state along one derivative estimate.
// RK4 uses this helper to evaluate forces at its midpoint and endpoint states.
State addScaledDerivative(const State& state,
                          const Derivative& derivative,
                          double scale) {
    return {
        state.altitude + derivative.altitudeRate * scale,
        state.velocity + derivative.velocityRate * scale,
        state.horizontalPosition +
            derivative.horizontalPositionRate * scale,
        state.horizontalVelocity +
            derivative.horizontalVelocityRate * scale
    };
}
// PID fin controller
// A PID controller is a feedback law that turns measured path error into a fin
// command. It combines three complementary views of horizontal drift:
//   P (proportional): responds to current displacement. More drift immediately
//      requests more correction toward the desired vertical flight line.
//   I (integral): sums error over time. It removes steady offset that can remain
//      when a constant crosswind balances proportional control alone.
//   D (derivative): responds to the rate at which error changes. It anticipates
//      lateral motion and adds damping, reducing overshoot and oscillation.
// Here the target horizontal position is zero, so a positive wind-driven drift
// creates a negative error and therefore a correcting fin command. In a real
// rocket, fin deflection creates aerodynamic force and torque that restore a
// near-vertical attitude. This simplified model applies the equivalent lateral
// correction force directly, which demonstrates the stability control loop
// without introducing full rotational rigid-body dynamics.
ControlOutput updatePID(double horizontalPosition,
                        PIDState& pidState,
                        double dt,
                        const SimulationParameters& parameters) {
    // Navigation error is desired minus measured horizontal position.
    const double targetPosition = 0.0;
    const double error = targetPosition - horizontalPosition;
    // Numerical integration of error supplies the controller's memory of
    // persistent drift. Clamping is a simple anti-windup safeguard.
    pidState.integralError = std::clamp(
        pidState.integralError + error * dt,
        -MAX_INTEGRAL_ERROR,
        MAX_INTEGRAL_ERROR);
    // A finite difference estimates how quickly the error changed since the
    // previous millisecond. Save the current error for the next update.
    const double errorDerivative = (error - pidState.previousError) / dt;
    pidState.previousError = error;
    // Each gain converts its associated error quantity into radians of desired
    // fin deflection. Keeping the terms separate also makes CSV tuning easier.
    const double proportionalTerm = parameters.pidKp * error;
    const double integralTerm =
        parameters.pidKi * pidState.integralError;
    const double derivativeTerm = parameters.pidKd * errorDerivative;
    // The ideal PID output is the sum of P, I, and D actions.
    const double unrestrictedAngle =
        proportionalTerm + integralTerm + derivativeTerm;
    // Real fins have limited travel, so clip the ideal command to the modeled
    // mechanical limits before converting it to force.
    const double finAngle = std::clamp(
        unrestrictedAngle,
        -MAX_FIN_DEFLECTION_RAD,
        MAX_FIN_DEFLECTION_RAD);
    // FIN_FORCE_PER_RADIAN maps the limited angle to the lateral force passed
    // into the equations of motion during the next RK4 step.
    return {
        error,
        proportionalTerm,
        integralTerm,
        derivativeTerm,
        finAngle,
        FIN_FORCE_PER_RADIAN * finAngle
    };
}
// Fourth-order Runge-Kutta (RK4) integrator
// RK4 advances the state with much less numerical error than a single Euler
// slope. It samples the equations of motion four times across one timestep,
// then forms a weighted average of those slopes. Because time and trial state
// are passed to every evaluation, changing mass, air density, drag, and the
// thrust-to-coast transition are represented within the step.
// The PID controller is sampled once per millisecond. Its commanded fin force
// is held constant during k1 through k4, which represents a digital controller
// holding its output until the next sample.
State rk4Step(double time,
              const State& state,
              double dt,
              double finCorrectionForce,
              const SimulationParameters& parameters) {
    // k1: slope at the beginning of the timestep, using the known state.
    const Derivative k1 =
        equationsOfMotion(time, state, finCorrectionForce, parameters);
    // k2: slope halfway through the timestep, using a midpoint state predicted
    // by advancing half a step along k1.
    const Derivative k2 =
        equationsOfMotion(time + 0.5 * dt,
                          addScaledDerivative(state, k1, 0.5 * dt),
                          finCorrectionForce,
                          parameters);
    // k3: a refined midpoint slope. It predicts another half-step state using
    // k2 rather than k1, improving the estimate of conditions at mid-step.
    const Derivative k3 =
        equationsOfMotion(time + 0.5 * dt,
                          addScaledDerivative(state, k2, 0.5 * dt),
                          finCorrectionForce,
                          parameters);
    // k4: slope at the end of the timestep, using a full-step state predicted
    // from k3. It captures how forces look at the proposed endpoint.
    const Derivative k4 =
        equationsOfMotion(time + dt,
                          addScaledDerivative(state, k3, dt),
                          finCorrectionForce,
                          parameters);
    // Combine the four slopes with classical RK4 weights 1, 2, 2, 1. The
    // result advances all four state variables by one complete timestep.
    return {
        state.altitude +
            (dt / 6.0) * (k1.altitudeRate + 2.0 * k2.altitudeRate +
                          2.0 * k3.altitudeRate + k4.altitudeRate),
        state.velocity +
            (dt / 6.0) * (k1.velocityRate + 2.0 * k2.velocityRate +
                          2.0 * k3.velocityRate + k4.velocityRate),
        state.horizontalPosition +
            (dt / 6.0) *
                (k1.horizontalPositionRate +
                 2.0 * k2.horizontalPositionRate +
                 2.0 * k3.horizontalPositionRate +
                 k4.horizontalPositionRate),
        state.horizontalVelocity +
            (dt / 6.0) *
                (k1.horizontalVelocityRate +
                 2.0 * k2.horizontalVelocityRate +
                 2.0 * k3.horizontalVelocityRate +
                 k4.horizontalVelocityRate)
    };
}
// CSV output
// Write one complete telemetry sample. The CSV columns, in order, are:
//   time                         elapsed simulation time (s)
//   altitude                     vertical position above ground (m)
//   velocity                     signed vertical velocity (m/s)
//   acceleration                 signed vertical acceleration (m/s^2)
//   mass                         remaining vehicle mass (kg)
//   horizontal_drift             lateral position from target line (m)
//   horizontal_velocity          signed lateral velocity (m/s)
//   crosswind_force              applied disturbance force (N)
//   fin_correction_angle_deg     commanded fin deflection (deg)
//   fin_correction_force         modeled lateral fin force (N)
//   pid_error                    current position error (m)
//   pid_p_term                   proportional angle contribution (rad)
//   pid_i_term                   integral angle contribution (rad)
//   pid_d_term                   derivative angle contribution (rad)
void writeSample(std::ofstream& csv,
                 double time,
                 const State& state,
                 const ControlOutput& control,
                 const SimulationParameters& parameters) {
    // Re-evaluate the equations at the logged state so acceleration in the CSV
    // corresponds to the same time, state, mass, and controller output.
    const double acceleration =
        equationsOfMotion(time, state, control.finForce, parameters)
            .velocityRate;
    const double mass = massAt(time, parameters);
    // Six decimal places preserve the millisecond time resolution and enough
    // precision for plotting or comparing controller behavior between runs.
    csv << std::fixed << std::setprecision(6)
        << time << ',' << state.altitude << ',' << state.velocity << ','
        << acceleration << ',' << mass << ','
        << state.horizontalPosition << ',' << state.horizontalVelocity << ','
        << parameters.crosswindForce << ','
        << control.finAngle * 180.0 / PI << ',' << control.finForce << ','
        << control.error << ',' << control.proportionalTerm << ','
        << control.integralTerm << ',' << control.derivativeTerm << '\n';
}
// Format the logged-row count with thousands separators for readable terminal
// output; this has no effect on the CSV data.
std::string formatCount(std::size_t count) {
    std::string result = std::to_string(count);
    for (std::size_t position = result.length();
         position > 3;
         position -= 3) {
        result.insert(position - 3, ",");
    }
    return result;
}
// Find the first run number whose CSV and PNG names are both unused. Reserving
// a shared suffix keeps each telemetry file paired with its plot and prevents
// a new simulation from overwriting earlier test results.
int nextAvailableRunNumber(const std::filesystem::path& outputDirectory) {
    int runNumber = 1;
    while (true) {
        const std::string suffix = "_run" + std::to_string(runNumber);
        const std::filesystem::path csvCandidate =
            outputDirectory / ("flight_data" + suffix + ".csv");
        const std::filesystem::path pngCandidate =
            outputDirectory / ("flight_results" + suffix + ".png");

        if (!std::filesystem::exists(csvCandidate) &&
            !std::filesystem::exists(pngCandidate)) {
            return runNumber;
        }
        ++runNumber;
    }
}                                                            
int main(int argc, char* argv[]) {
    // Gather and validate the rocket, atmosphere interaction, wind, and control
    // parameters before any simulation state is created.
    const SimulationParameters parameters = readSimulationParameters();
    // Keep the CSV, plotting script, and graph beside the executable even if
    // the simulator is launched from a different working directory. The run
    // suffix associates this run's CSV and PNG without overwriting older files.
    const std::filesystem::path programDirectory =
        std::filesystem::absolute(argc > 0 ? argv[0] : "rocket_sim")
            .parent_path();
    const int runNumber = nextAvailableRunNumber(programDirectory);
    const std::string runSuffix = "_run" + std::to_string(runNumber);
    const std::filesystem::path csvPath =
        programDirectory / ("flight_data" + runSuffix + ".csv");
    const std::filesystem::path pngPath =
        programDirectory / ("flight_results" + runSuffix + ".png");
    const std::filesystem::path plotScriptPath =
        programDirectory / "plot_results.py";
    // Open this run's CSV before integration. The header names correspond
    // exactly to the values documented in writeSample above and to the column
    // names used by plot_results.py.
    std::ofstream csv(csvPath);
    if (!csv) {
        std::cerr << "Error: could not create "
                  << csvPath.filename().string() << '\n';
        return 1;
    }
    csv << "time,altitude,velocity,acceleration,mass,"
           "horizontal_drift,horizontal_velocity,crosswind_force,"
           "fin_correction_angle_deg,fin_correction_force,pid_error,"
           "pid_p_term,pid_i_term,pid_d_term\n";
    // Main simulation loop
    // Begin at rest on the launch pad and on the intended vertical flight line.
    // Continue in fixed 0.001-second RK4 steps until vertical velocity changes
    // from upward to downward (apogee) and altitude returns to ground level.
    double time = 0.0;
    State state{0.0, 0.0, 0.0, 0.0};
    PIDState pidState{0.0, 0.0};
    ControlOutput control =
        updatePID(
            state.horizontalPosition, pidState, TIME_STEP, parameters);
    double maxAltitude = state.altitude;
    double maxVelocity = state.velocity;
    double maxHorizontalDrift = std::abs(state.horizontalPosition);
    std::size_t loggedTimesteps = 0;
    bool apogeeReached = false;
    // Log the initial t = 0 state before taking the first integration step.
    writeSample(csv, time, state, control, parameters);
    ++loggedTimesteps;
    while (true) {
        // Advance physics with the fin force commanded at the current sample,
        // then update the digital PID controller for the new state.
        const double previousVelocity = state.velocity;
        state = rk4Step(
            time, state, TIME_STEP, control.finForce, parameters);
        time += TIME_STEP;
        control = updatePID(
            state.horizontalPosition, pidState, TIME_STEP, parameters);
        // Maintain summary statistics without storing the full trajectory in
        // memory; the detailed history is streamed directly to the CSV.
        maxAltitude = std::max(maxAltitude, state.altitude);
        maxVelocity = std::max(maxVelocity, state.velocity);
        maxHorizontalDrift =
            std::max(maxHorizontalDrift,
                     std::abs(state.horizontalPosition));
        // A positive-to-nonpositive vertical-velocity crossing marks apogee.
        if (previousVelocity > 0.0 && state.velocity <= 0.0) {
            apogeeReached = true;
        }
        // Report touchdown at exactly ground level rather than recording the
        // small negative altitude produced by the final fixed-size time step.
        if (apogeeReached && state.altitude <= 0.0) {
            state.altitude = 0.0;
            writeSample(csv, time, state, control, parameters);
            ++loggedTimesteps;
            break;
        }
        writeSample(csv, time, state, control, parameters);
        ++loggedTimesteps;
    }
    csv.close();
    std::cout << "\nSimulation complete. "
              << formatCount(loggedTimesteps)
              << " timesteps logged.\n\n"
              << "Flight summary\n"
              << "--------------\n"
              << std::fixed << std::setprecision(2)
              << "Max altitude:     " << maxAltitude << " m\n"
              << "Max velocity:     " << maxVelocity << " m/s\n"
              << "Max horizontal drift: "
              << maxHorizontalDrift << " m\n"
              << "Total flight time: " << time << " s\n"
              << "\nCSV generated: "
              << csvPath.filename().string() << '\n';
    // Automatically regenerate the graph only after the CSV is closed, which
    // guarantees that Python can read the complete dataset. Prefer the local
    // virtual environment because it contains pandas and matplotlib; otherwise
    // fall back to the system "python" command.
    const std::filesystem::path virtualEnvironmentPython =
        programDirectory / ".venv" / "Scripts" / "python.exe";
    const std::string plotCommand =
        std::filesystem::exists(virtualEnvironmentPython)
            ? "\"\"" + virtualEnvironmentPython.string() + "\" \"" +
                  plotScriptPath.string() + "\" \"" + csvPath.string() +
                  "\" \"" + pngPath.string() + "\"\""
            : "python \"" + plotScriptPath.string() + "\" \"" +
                  csvPath.string() + "\" \"" + pngPath.string() + "\"";
    // system() launches plot_results.py and passes this run's unique CSV input
    // and PNG output paths. Quoting preserves paths that contain spaces. The
    // return code determines whether the success or failure message is shown.
    const int plotStatus = std::system(plotCommand.c_str());
    if (plotStatus == 0) {
        std::cout << "Graph generated: "
                  << pngPath.filename().string() << '\n';
    } else {
        std::cerr << "Graph regeneration failed (plot_results.py returned "
                  << plotStatus << ").\n";
    }
    return 0;
}
