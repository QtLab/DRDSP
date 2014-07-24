# Obtain a Reduced Family

The `ReducedDataSystem` class takes the original family, the data sets, and an orthonormal matrix and computes the data required to find the reduced family:

```cpp
ReducedDataSystem reducedData;

// Compute data for the given family, data and projection using 4 threads.
reducedData.ComputeData( exampleFamily, data, projSecant.W, 4 );
```

The `RBFFamilyProducer` class takes this `ReducedDataSystem` and produces the reduced family:

```cpp
RBFFamilyProducer<RBF<RBFType>> producer( 30 );   // the number of rbfs to use
auto reducedFamily = producer.BruteForce( reducedData,
                                          data.parameters,
                                          1000,   // the number of iterations to perform
                                          4 );    // the number of threads
```

The `RBFType` is a radial basis function type, which can be one of the following:

* `ThinPlateSpline` -- r^2 log(r)
* `PolyharmonicSpline<3>` -- r^3
* `Multiquadratic` -- sqrt( 1 + r^2 )
* `InverseQuadratic` -- 1 / ( 1 + r^2 )
* `InverseMultiquadratic` -- 1 / sqrt( 1 + r^2 )
* `Gaussian` -- exp(-r^2)

The resulting `reducedFamily` is a family of `RBFModel`s parameterized by the original parameter space. This can be used in any of the methods that expect a family, such as `DataGenerator` and `BifurcationDiagramGenerator`.

The actual type of the reduced family will be `PMapFamily<RBFFamily<RBF<RBFType>>,AffineXd>`.
