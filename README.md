# Covariant_Derivative

A Mathematica package for symbolic general relativity calculations, designed to perform covariant derivative manipulations and compute various geometric tensors automatically.

## Overview

This package provides an elegant and efficient implementation of general relativity tensor calculations using Mathematica's built-in symbolic operations. The design philosophy emphasizes:

- **Symbolic elegance**: Avoids explicit index loops (like `Table`) in favor of Mathematica's native tensor operations
- **Performance optimization**: Leverages vectorized operations and built-in functions for faster computation
- **Mathematical clarity**: Uses standard general relativity conventions and notation

## Design Philosophy

The package is designed around Mathematica's powerful symbolic computation capabilities, making extensive use of:

- **Tensor operations**: `Transpose`, `TensorProduct`, `TensorContract` for index manipulations
- **Symbolic differentiation**: `D[]` for computing derivatives of metric tensors
- **Matrix operations**: `Inverse[]`, `Tr[]`, and dot products for efficient linear algebra
- **Higher-order functions**: `Through[]`, `Map[]` for functional programming patterns

This approach avoids the computational overhead of explicit index summations and makes the code more readable and maintainable.

## Tensor Index Notation

Following the xAct package convention, tensor indices are denoted as:
- **Upper indices**: Prefixed with minus sign (e.g., `-a`, `-b`)  
- **Lower indices**: Use the symbol directly (e.g., `c`, `d`)

For example, the tensor $G^{ab}{}_{cd}{}^{ef}{}_{hi}$ is represented as:
```mathematica
G[-a, -b, c, d, -e, -f, h, i]
```

## Mathematical Foundations

### Christoffel Symbols
The Christoffel symbols are computed using:
$$\Gamma^\mu_{\nu\rho} = \frac{1}{2}g^{\mu\sigma}\left(\partial_\nu g_{\sigma\rho} + \partial_\rho g_{\nu\sigma} - \partial_\sigma g_{\nu\rho}\right)$$

### Riemann Curvature Tensor
Using Wald's convention:
$$R^\rho{}_{\sigma\mu\nu} = \partial_\mu\Gamma^\rho_{\nu\sigma} - \partial_\nu\Gamma^\rho_{\mu\sigma} + \Gamma^\rho_{\mu\lambda}\Gamma^\lambda_{\nu\sigma} - \Gamma^\rho_{\nu\lambda}\Gamma^\lambda_{\mu\sigma}$$

### Ricci Tensor and Scalar
- Ricci tensor: $R_{\mu\nu} = R^\rho{}_{\mu\rho\nu}$
- Ricci scalar: $R = g^{\mu\nu}R_{\mu\nu}$

### Einstein Tensor
$$G_{\mu\nu} = R_{\mu\nu} - \frac{1}{2}Rg_{\mu\nu}$$

### Weyl Conformal Tensor
In $n$ dimensions:
$$C_{\mu\nu\rho\sigma} = R_{\mu\nu\rho\sigma} - \frac{2}{n-2}\left(g_{\mu[\rho}R_{\sigma]\nu} - g_{\nu[\rho}R_{\sigma]\mu}\right) + \frac{2R}{(n-1)(n-2)}g_{\mu[\rho}g_{\sigma]\nu}$$

## Functions

### Core Tensor Functions
- **`ChristoffelSymbol[g, xx]`**: Computes Christoffel symbols from metric tensor
- **`RiemannTensor[g, xx]`**: Computes Riemann curvature tensor
- **`RicciTensor[g, xx]`**: Computes Ricci tensor
- **`RicciScalar[g, xx]`**: Computes Ricci scalar
- **`EinsteinTensor[g, xx]`**: Computes Einstein tensor
- **`WeylTensor[g, xx]`**: Computes Weyl conformal tensor

### Covariant Derivative
- **`CovariantDerivative[a, ten, g, xx]`**: Computes covariant derivative of tensor `ten` with respect to coordinate `a`

The covariant derivative automatically handles:
- Metric compatibility: $\nabla_\mu g_{\nu\rho} = 0$
- Index raising and lowering using the metric tensor
- Proper tensor contractions for repeated indices

## Usage Example

```mathematica
(* Load the package *)
Get["gr_package.wl"]

(* Define coordinates and metric *)
xx = {t, r, θ, φ};
g = DiagonalMatrix[{-1, 1, r^2, r^2 Sin[θ]^2}]; (* Schwarzschild-like *)

(* Compute geometric quantities *)
christoffel = ChristoffelSymbol[g, xx];
riemann = RiemannTensor[g, xx];
ricci = RicciTensor[g, xx];
einstein = EinsteinTensor[g, xx];
```

## Installation and Testing

See `gr_test.ipynb` for detailed examples and test cases. To use with Jupyter notebooks:

1. Install Wolfram Engine
2. Set up Wolfram Language kernel in Jupyter
3. Load the package and run the test examples

## Performance Notes

The package is optimized for symbolic computation and scales well with the complexity of the metric tensor. The use of Mathematica's built-in tensor operations provides significant performance advantages over index-based implementations, especially for:

- High-dimensional spacetimes
- Metrics with complex symbolic entries  
- Repeated calculations with different coordinate systems

## Author

**Peng Liu** - Physicist working on holographic duality and black hole physics

## References

- Wald, R. M. *General Relativity* (University of Chicago Press, 1984)
- xAct tensor computer algebra system notation conventions

