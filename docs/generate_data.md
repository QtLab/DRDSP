# Generating Data

The method requires data samples from the attractors of the original model. These data sets can either be loaded in from a file or generated at run-time.

To generate a data set, the `DataGenerator` class will perform a simulation of the original model for multiple parameter values.

```cpp
auto parameters = ParameterList( 1.0, 2.0, 11 );   // Set of 11 parameter values in the range [1,2] (uniformly spaced)

DataGenerator<ExampleFamily> dataGenerator;
dataGenerator.initial.setRandom( 10 );   // Set the initial condition for the simulations
dataGenerator.tStart = 1000;             // Start recording data at this time (allow for transients to decay)
dataGenerator.tInterval = 100;           // Time interval over which to record data
dataGenerator.dtMax = 0.001;             // Maximum time-step to use in the simulation
dataGenerator.print = 1000;              // Number of data samples to take over this interval (uniformly spaced)

// Generate data for the given parameter values using 4 threads.
DataSystem data = dataGenerator.GenerateDataSystem( parameters, 4 ); 
```

To use a custom solver, use `DataGenerator<ExampleFamily,ExampleSolver>`.
