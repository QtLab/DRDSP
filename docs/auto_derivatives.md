# Automatically Computing the Partial Derivatives

DRDSP features a technique called automatic differentiation that can compute the partial derivatives of your vector field for you. In order to use this, there are two things you need to do. First include the header

```cpp
#include <DRDSP/auto_diff.h>
```

and use it in the model definition like so:

```cpp
MatrixXd Partials( const VectorXd& state ) const {
	return AutoDerivative( *this, state );
}
```

For examples that have sparse derivatives, such as a discretized PDE, you can use the sparse version:

```cpp
SparseMatrix<double> Partials( const VectorXd& state ) const {
	return AutoDerivativeSparse( *this, state );
}
```

Secondly, you need to make sure the implemention of the vector field has a templated scalar type:

```cpp
template<typename Derived>
Matrix<typename Derived::Scalar,-1,1> operator()( const MatrixBase<Derived>& x ) const {
	// evaluate the vector field
}
```

DRDSP will compute the partial derivatives by evaluating the vector field with dual numbers. For more information on automatic differentiation and dual numbers see [here](https://cwzx.wordpress.com/2014/05/31/automatically-computing-the-derivatives-of-a-vector-field-using-dual-numbers/).
