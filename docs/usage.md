# Usage

In order to use the library with your example system, you need to implement the dynamics of the system.

### Implement your model

A model is essentially a vector field.

A base class `Model` is provided in `DRDSP/dynamics/model.h`. Inherit from this and pass it the dimension of your state space in the constructor.

```cpp
struct ExampleModel : Model<> {
	double a = 1.0,     // parameters and other data used to specify the model
	       b = 1.0,
	       c = 1.0;

	ExampleModel() : Model<>(10) {}   // set the dimension of the state space in the constructor
	
	// the vector field
	VectorXd operator()( const VectorXd& state ) const;
	
	// the partial derivatives of the vector field
	MatrixXd Partials( const VectorXd& state ) const;
};
```

### Implement a family

A family is a function object that takes a parameter and returns a model. This is where you choose which parameter(s) you're interested in investigating.

A base class `Family` is provided in `DRDSP/dynamics/model.h`. Inherit from this and pass it the dimension of your state space and parameter space in the constructor.

```cpp
struct ExampleFamily : Family<ExampleModel> {

	ExampleFamily() :
		Family<ExampleModel>(10,1)    // set the dimension of the state space and parameter space
	{}
	
	ExampleModel operator()( const VectorXd& parameter ) const {
		ExampleModel model;
		model.c = parameter[0];      // in this family we are investigating the c parameter.
		return model;
	}
};
```
