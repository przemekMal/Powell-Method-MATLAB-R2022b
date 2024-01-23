# Powell Method with Linear Constraints in MATLAB

This repository contains an implementation of Powell's optimization algorithm with linear constraints, written in MATLAB.

## Usage

To run the algorithm, use the file `powell_method_with_penalty_function.m`. Ensure you have MATLAB installed.

```matlab
% Example usage:
% Define the objective function
f = @(x)  x(1)^2+x(2)^2;

% Set initial conditions
x0 = [5.6, -2.6536];

% Set linear constraints (if any)
g = @(x) x(1)-x(2)-2;
G = {g};

% Run the optimization
[minimum, xes, iter] = powell_method_with_penalty_function(f, x0, [], 10, 10e-3, 12, G);

Test Files

Test files with names starting with 'test' are provided to help you validate the functionality of the algorithm. Feel free to explore and adapt these test files for your specific use case.
File Structure

    powell_method_with_penalty_function.m: Main script implementing Powell's method with linear constraints.
    test.m, testRastringera.m, ...: Example test scripts.

Contributors

    przemekMal
