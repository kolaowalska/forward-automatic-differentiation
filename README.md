This code implements automatic differentiation (FAD) using jet expansion to compute the value of a funtion as well as its first and second derivatives with respect to multiple variables. The use of operator overloading allows the propagation of derivatives seamlessly during arithmetic operations and function evaluations. This approach is optimal for calculating derivatives of complex functions in a numerically stable and computationally efficient way by avoiding symbolic differentiation's complexity and numerical differentiation's inacurracies, such as finite differences.

The `Jet` class represents a multivariate Taylor expansion of a function, storing the function's value, gradients, and Hessians, as well as provides operator overloading for basic arithmetic operations and mathematical functions. Each `Jet` object holds the function value (`f`), first-order partial derivatives (`dx`, `dy`), and second-order partial derivatives (`dxx`, `dxy`, `dyy`). The class also includes operator overloads and function implementations to propagate derivatives automatically for:
* Arithmetic operations `+`,  `-`, `*`, `/` between `Jet` objects and constants
* Basic mathematical functions, such as `sin`, `cos`, `exp`, applied to `Jet` objects

The integer `M` in the `source` function symbolizes the number of evaluations and then pais of inputs $(x_0, y_0)$, which are later used for creating input pairs.



