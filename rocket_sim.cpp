#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

// ============================================================================
// Constants and model-rocket parameters
// ============================================================================

constexpr double TIME_STEP = 0.001;            // RK4 and logging step (s)
constexpr double GRAVITY = 9.80665;             // Gravitational acceleration (m/s^2)
constexpr double SEA_LEVEL_DENSITY = 1.225;     // Air density at sea level (kg/m^3)
constexpr double ATMOSPHERE_SCALE_HEIGHT = 8500.0; // Exponential scale height (m)

constexpr double ROCKET_DIAMETER = 0.08;        // Body diameter (m)
constexpr double PI = 3.14159265358979323846;
constexpr double DEFAULT_CROSS_SECTIONAL_AREA =
    PI * ROCKET_DIAMETER * ROCKET_DIAMETER / 4.0; // Frontal area (m^2)

// Crosswind and fin-control parameters. The crosswind is represented by a
// small, constant horizontal disturbance force in the positive direction.
constexpr double FIN_FORCE_PER_RADIAN = 4.0;      // N/rad
constexpr double MAX_FIN_DEFLECTION_DEG = 12.0;   // Mechanical fin limit
constexpr double MAX_FIN_DEFLECTION_RAD =
    MAX_FIN_DEFLECTION_DEG * PI / 180.0;

constexpr double MAX_INTEGRAL_ERROR = 20.0; // Anti-windup limit (m*s)

// These defaults are displayed in the prompts and can be accepted with Enter.
struct SimulationParameters {
    double initialMass = 0.75;        // Rocket plus propellant (kg)
    double propellantMass = 0.25;     // Propellant consumed during burn (kg)
    double thrust = 25.0;             // Constant burn thrust (N)
    double burnTime = 2.0;            // s
    double dragCoefficient = 0.75;
    double crossSectionalArea = DEFAULT_CROSS_SECTIONAL_AREA; // m^2
    double crosswindForce = 0.30;      // N; sign selects wind direction
    double pidKp = 0.040;              // rad/m
    double pidKi = 0.005;              // rad/(m*s)
    double pidKd = 0.080;              // rad*s/m

    double dryMass() const {
        return initialMass - propellantMass;
    }
};

struct State {
    double altitude;           // m
    double velocity;           // m/s; upward is positive
    double horizontalPosition; // m; intended flight line is position zero
    double horizontalVelocity; // m/s
};

struct Derivative {
    double altitudeRate;           // m/s
    double velocityRate;           // m/s^2
    double horizontalPositionRate; // m/s
    double horizontalVelocityRate; // m/s^2
};

struct PIDState {
    double integralError;
    double previousError;
};

struct ControlOutput {
    double error;             // Current horizontal position error (m)
    double proportionalTerm; // P contribution to commanded angle (rad)
    double integralTerm;     // I contribution to commanded angle (rad)
    double derivativeTerm;   // D contribution to commanded angle (rad)
    double finAngle;          // Limited fin command (rad)
    double finForce;          // Resulting horizontal correction force (N)
};

// ============================================================================
// Interactive input handling
// ============================================================================
//
// Every prompt shows the current default, and a blank line accepts it. The
// parser requires one finite number with no extra text. If a value violates its
// physical validation rule, an error is shown and the same question is asked
// again.
template <typename Validator>
double promptForValue(const std::string& label,
                      double defaultValue,
                      Validator isValid,
                      const std::string& validationMessage) {
    while (true) {
        std::cout << label << " [default " << std::setprecision(8)
                  << defaultValue << "]: ";

        std::string line;
        if (!std::getline(std::cin, line)) {
            std::cout << "\nInput ended; using the default value.\n";
            return defaultValue;
        }

        if (line.find_first_not_of(" \t\r\n") == std::string::npos) {
            return defaultValue;
        }

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
    SimulationParameters parameters;

    std::cout << "Rocket Flight Simulation Setup\n"
              << "Press Enter at any prompt to accept its default.\n\n";

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

// Return the rocket's mass at a specified time. Propellant is consumed at a
// constant rate during the burn, after which the dry mass remains constant.
double massAt(double time, const SimulationParameters& parameters) {
    if (parameters.burnTime == 0.0 || time >= parameters.burnTime) {
        return parameters.dryMass();
    }

    const double burnFraction =
        std::clamp(time / parameters.burnTime, 0.0, 1.0);
    return parameters.initialMass -
           parameters.propellantMass * burnFraction;
}

// Use a simplified exponential atmosphere. Altitudes below ground are treated
// as ground level so that the model remains well behaved at touchdown.
double airDensityAt(double altitude) {
    const double nonnegativeAltitude = std::max(0.0, altitude);
    return SEA_LEVEL_DENSITY *
           std::exp(-nonnegativeAltitude / ATMOSPHERE_SCALE_HEIGHT);
}

// ============================================================================
// Equations of motion
// ============================================================================
//
// The state derivatives are:
//   d(altitude)/dt = velocity
//   d(velocity)/dt = net force / mass
//   d(horizontal position)/dt = horizontal velocity
//   d(horizontal velocity)/dt = horizontal force / mass
//
// Drag always points opposite the velocity. Writing it as
// -0.5*rho*v*abs(v)*Cd*A supplies the correct sign for both ascent and descent.
// The vertical equations are unchanged; the horizontal equation simply adds
// the fixed crosswind disturbance and the PID-commanded fin correction.
Derivative equationsOfMotion(double time,
                             const State& state,
                             double finCorrectionForce,
                             const SimulationParameters& parameters) {
    const double mass = massAt(time, parameters);
    const double thrustForce =
        (time < parameters.burnTime) ? parameters.thrust : 0.0;
    const double gravityForce = mass * GRAVITY;
    const double dragForce =
        -0.5 * airDensityAt(state.altitude) * state.velocity *
        std::abs(state.velocity) * parameters.dragCoefficient *
        parameters.crossSectionalArea;

    const double netForce = thrustForce - gravityForce + dragForce;
    const double horizontalForce =
        parameters.crosswindForce + finCorrectionForce;

    return {
        state.velocity,
        netForce / mass,
        state.horizontalVelocity,
        horizontalForce / mass
    };
}

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

// ============================================================================
// PID fin controller
// ============================================================================
//
// A PID controller combines three views of the horizontal drift error:
//   P (proportional): reacts immediately to the current displacement.
//   I (integral): accumulates persistent error, removing steady crosswind bias.
//   D (derivative): reacts to how quickly error changes and damps overshoot.
//
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
    const double targetPosition = 0.0;
    const double error = targetPosition - horizontalPosition;

    pidState.integralError = std::clamp(
        pidState.integralError + error * dt,
        -MAX_INTEGRAL_ERROR,
        MAX_INTEGRAL_ERROR);
    const double errorDerivative = (error - pidState.previousError) / dt;
    pidState.previousError = error;

    const double proportionalTerm = parameters.pidKp * error;
    const double integralTerm =
        parameters.pidKi * pidState.integralError;
    const double derivativeTerm = parameters.pidKd * errorDerivative;
    const double unrestrictedAngle =
        proportionalTerm + integralTerm + derivativeTerm;
    const double finAngle = std::clamp(
        unrestrictedAngle,
        -MAX_FIN_DEFLECTION_RAD,
        MAX_FIN_DEFLECTION_RAD);

    return {
        error,
        proportionalTerm,
        integralTerm,
        derivativeTerm,
        finAngle,
        FIN_FORCE_PER_RADIAN * finAngle
    };
}

// ============================================================================
// Fourth-order Runge-Kutta (RK4) integrator
// ============================================================================
//
// RK4 evaluates the time-dependent equations of motion four times per step,
// including the changing mass and the end of the thrust phase. The controller
// is sampled once per timestep, so its fin force is held constant through all
// four RK4 evaluations for that step.
State rk4Step(double time,
              const State& state,
              double dt,
              double finCorrectionForce,
              const SimulationParameters& parameters) {
    const Derivative k1 =
        equationsOfMotion(time, state, finCorrectionForce, parameters);
    const Derivative k2 =
        equationsOfMotion(time + 0.5 * dt,
                          addScaledDerivative(state, k1, 0.5 * dt),
                          finCorrectionForce,
                          parameters);
    const Derivative k3 =
        equationsOfMotion(time + 0.5 * dt,
                          addScaledDerivative(state, k2, 0.5 * dt),
                          finCorrectionForce,
                          parameters);
    const Derivative k4 =
        equationsOfMotion(time + dt,
                          addScaledDerivative(state, k3, dt),
                          finCorrectionForce,
                          parameters);

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

// ============================================================================
// CSV output
// ============================================================================

void writeSample(std::ofstream& csv,
                 double time,
                 const State& state,
                 const ControlOutput& control,
                 const SimulationParameters& parameters) {
    const double acceleration =
        equationsOfMotion(time, state, control.finForce, parameters)
            .velocityRate;
    const double mass = massAt(time, parameters);

    csv << std::fixed << std::setprecision(6)
        << time << ',' << state.altitude << ',' << state.velocity << ','
        << acceleration << ',' << mass << ','
        << state.horizontalPosition << ',' << state.horizontalVelocity << ','
        << parameters.crosswindForce << ','
        << control.finAngle * 180.0 / PI << ',' << control.finForce << ','
        << control.error << ',' << control.proportionalTerm << ','
        << control.integralTerm << ',' << control.derivativeTerm << '\n';
}

// Format the logged-row count with thousands separators for readable output.
std::string formatCount(std::size_t count) {
    std::string result = std::to_string(count);
    for (std::size_t position = result.length();
         position > 3;
         position -= 3) {
        result.insert(position - 3, ",");
    }
    return result;
}

// Find the first run number whose CSV and PNG names are both unused. This
// preserves every earlier result instead of overwriting it.
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
    const SimulationParameters parameters = readSimulationParameters();

    // Keep the CSV, plotting script, and graph beside the executable even if
    // the simulator is launched from a different working directory.
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

    // Open the output file and write the requested CSV column headings.
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

    // ========================================================================
    // Main simulation loop
    // ========================================================================
    //
    // Begin at rest on the launch pad. Continue in fixed 0.001-second RK4
    // steps until the rocket has passed apogee and returns to ground level.
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

    writeSample(csv, time, state, control, parameters);
    ++loggedTimesteps;

    while (true) {
        const double previousVelocity = state.velocity;
        state = rk4Step(
            time, state, TIME_STEP, control.finForce, parameters);
        time += TIME_STEP;
        control = updatePID(
            state.horizontalPosition, pidState, TIME_STEP, parameters);

        maxAltitude = std::max(maxAltitude, state.altitude);
        maxVelocity = std::max(maxVelocity, state.velocity);
        maxHorizontalDrift =
            std::max(maxHorizontalDrift,
                     std::abs(state.horizontalPosition));

        if (previousVelocity > 0.0 && state.velocity <= 0.0) {
            apogeeReached = true;
        }

        // Report touchdown at exactly ground level rather than displaying the
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

    // Automatically regenerate the graph after the completed CSV is closed.
    const std::filesystem::path virtualEnvironmentPython =
        programDirectory / ".venv" / "Scripts" / "python.exe";
    const std::string plotCommand =
        std::filesystem::exists(virtualEnvironmentPython)
            ? "\"\"" + virtualEnvironmentPython.string() + "\" \"" +
                  plotScriptPath.string() + "\" \"" + csvPath.string() +
                  "\" \"" + pngPath.string() + "\"\""
            : "python \"" + plotScriptPath.string() + "\" \"" +
                  csvPath.string() + "\" \"" + pngPath.string() + "\"";
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
