# Hodgkin-Huxley Model Simulation

This repository contains MATLAB scripts for simulating the Hodgkin-Huxley model of action potential generation in neurons. The model describes how action potentials in neurons are initiated and propagated. It is a set of nonlinear differential equations that approximates the electrical characteristics of excitable cells such as neurons and cardiac cells.

## Files

1. `q1p1.m`: This script simulates the full Hodgkin-Huxley model. It includes the ability to simulate without threshold & rebound behavior observance and to test Threshold & Rebound behavior.

2. `q1p1_reduced.m`: This script simulates a reduced version of the Hodgkin-Huxley model. It includes the ability to simulate without threshold & rebound behavior observance and to test Threshold & Rebound behavior.

3. `hh_diff_eq.m`: This function file contains the differential equations for the full Hodgkin-Huxley model.

4. `hhdiff_reduced.m`: This function file contains the differential equations for the reduced Hodgkin-Huxley model.

5. `q2.m`: This script simulates the behavior of two reciprocally inhibiting neurons under constant current injection and sinusoidal current injection. It also includes the ability to simulate the effect of different frequencies of sinusoidal input on the firing rate of the neurons.

## Usage

To run the simulation, open either `q1p1.m` or `q1p1_reduced.m` in MATLAB and run the script. The script will generate a series of plots showing the evolution of the membrane potential, the conductance of the sodium and potassium channels, the evolution of the gating variables, and the current evolution.

To run the reciprocal_inhibition simulation, open q2.m in MATLAB and run the script. The script will generate a series of plots showing the neuronal firing in response to the input current, the change in membrane potential over time, and the change in the threshold potential over time.

## Parameters

The scripts use the following parameters:

- Initial Membrane Voltage: -65 mV
- Sodium gate (m): 0.05
- Potassium gate (n): 0.32
- Leak gate (h): 0.59
- Sodium conductance (gNa): 120 mS/cm^2
- Potassium conductance (gK): 36 mS/cm^2
- Leak conductance (gL): 0.3 mS/cm^2
- Membrane capacitance (mC): 1 uF/cm^2
- Nernst potentials: nL=-61 mV, nNa=55 mV, nK=-77 mV

For the reciprocal inhibition, The scripts use the following parameters:

- Resting Membrane Potential: 0 mV
- Membrane Resistance: 10 MOhm
- Membrane Capacitance: 1 nF
- Threshold Potential: 5 mV
- Spike Potential: 70 mV
- Threshold Time Constant: 50 ms
- Synaptic Reversal Potential: -15 mV
- Synaptic Time Constant: 15 ms
- Peak Synaptic Conductance: 0.1 nS

## Requirements

MATLAB is required to run these scripts. The scripts have been tested on MATLAB R2023a.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License

[MIT](https://choosealicense.com/licenses/mit/)
