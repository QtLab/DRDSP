# Models with Non-Euclidean State Spaces

If the state space of your model is something other than R^n, you must embed the state space into R^n. Refer to the pendulum example to see this in action.

### Implement an embedding

```cpp
struct ExampleEmbedding : Embedding {

	ExampleEmbedding() :
		Embedding(10,12)       // the dimension of the state space and embedding space
	{}
	
	// Evaluate the embedding.
	VectorXd operator()( const VectorXd& state ) const;
	
	// The partial derivatives of the embedding.
	MatrixXd Derivative( const VectorXd& state ) const;
	
	// The adjoint of the derivative. If the coordinate bases are orthonormal this will be the matrix transpose of the derivative.
	MatrixXd DerivativeAdjoint( const VectorXd& state ) const;
	
	// The second partial derivatives of the ith component of the embedding.
	MatrixXd Derivative2( const VectorXd& state, uint32_t i ) const;
};
```

To use an embedding with the underlying model/family you can use these as your model/family: `ModelEmbedded<ExampleModel,ExampleEmbedding>` and `FamilyEmbedded<ExampleFamily,ExampleEmbedding>`.

### Implement a wrap function

A wrap function is applied after every time-step to modify the state vector. This allows you to keep angular variables in the range [-pi,pi) (for example).

```cpp
struct ExampleWrap {
	void operator()( VectorXd& x ) const {
		Wrap(x[3],-M_PI,M_PI);
	}
};
```

In order to use a custom wrap function with the provided RK4 integrator, you need to make a custom version of the solver:

```cpp
typedef RKDynamicalSystem<ExampleModel,ExampleWrap> ExampleSolver;
```
