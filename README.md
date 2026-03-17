# Root-Finding-Project
Numerical root-finding toolbox in Python — implements Bisection, Newton, Secant, False Position, Fixed-Point, Muller, Steffensen, and Horner's methods via a single solve() interface. No external libraries.

# RootFindingProblem

A Python numerical root-finding toolbox implementing classical iterative methods for solving `f(x) = 0`. Built using only Python's standard library and `cmath` — no external solvers.

---

## Project Description

`RootFindingProblem` is a clean, object-oriented implementation of seven classical root-finding algorithms. It exposes a single public interface (`solve()`) that dispatches to the appropriate internal method based on a string argument, making it easy to switch between methods and compare results.

The project includes a worked example solving the polynomial:

```
600x⁴ - 550x³ + 200x² - 20x - 1 = 0   on [0.1, 1]
```

---

## Implemented Methods

| Method | Key | Description |
|---|---|---|
| Bisection | `bisection` | Bracket method; repeatedly halves the interval |
| Fixed-Point Iteration | `fixed_point` | Iterates `x = g(x)` until convergence |
| Newton's Method | `newton` | Uses function and derivative; quadratic convergence |
| Secant Method | `secant` | Derivative-free Newton variant using two prior points |
| False Position | `false_position` | Bracket method using linear interpolation |
| Horner's Method | `horner` | Efficient polynomial evaluation |
| Muller's Method | `muller` | Fits a parabola; can find complex roots |
| Steffensen's Method | `steffensen` | Accelerates fixed-point via Aitken's Δ² formula |

---

## Algorithm Details

### Bisection Method
Requires a bracket `[a, b]` where `f(a)` and `f(b)` have opposite signs. At each step, the midpoint `c = (a + b) / 2` is evaluated. The sub-interval containing the sign change is retained. Guaranteed to converge, but slowly (linear).

### Fixed-Point Iteration
Reformulates `f(x) = 0` as `x = g(x)`. Starting from an initial guess `x₀`, iterates `xₙ₊₁ = g(xₙ)`. Converges when `|g'(x)| < 1` near the root. A common reformulation is `g(x) = x - λ·f(x)` for a small λ.

### Newton's Method
Uses the tangent line at each iterate: `xₙ₊₁ = xₙ - f(xₙ) / f'(xₙ)`. Requires the derivative `df`. Exhibits quadratic convergence near a simple root but may fail if `f'(x) ≈ 0`.

### Secant Method
Approximates the derivative using two previous iterates:
```
xₙ₊₁ = xₙ - f(xₙ) · (xₙ - xₙ₋₁) / (f(xₙ) - f(xₙ₋₁))
```
No derivative needed. Superlinear convergence (order ≈ 1.618).

### False Position (Regula Falsi)
Like bisection, but uses linear interpolation instead of the midpoint:
```
c = (a·f(b) - b·f(a)) / (f(b) - f(a))
```
Maintains a bracket, so it is always convergent, but may be slow if one endpoint is stagnant.

### Horner's Method
Evaluates a polynomial `P(x) = a₀ + a₁x + a₂x² + ... + aₙxⁿ` using nested multiplication:
```
P(x) = (...((aₙ·x + aₙ₋₁)·x + aₙ₋₂)·x + ... + a₀)
```
Reduces multiplications from O(n²) to O(n). Used internally for all polynomial evaluations.

### Muller's Method
Fits a quadratic (parabola) through three points and solves for its root using the quadratic formula. Selects the root closer to `x₂` to avoid large jumps. Naturally handles complex roots via `cmath`. Order of convergence ≈ 1.839.

### Steffensen's Method
Accelerates fixed-point iteration using Aitken's Δ² formula:
```
xₙ₊₁ = xₙ - (g(xₙ) - xₙ)² / (g(g(xₙ)) - 2g(xₙ) + xₙ)
```
Achieves quadratic convergence from a linearly converging fixed-point scheme.

---

## Error Handling

All methods raise meaningful exceptions:

| Condition | Exception |
|---|---|
| Invalid bracket (same sign) in bisection / false position | `ValueError` |
| Missing derivative `df` in Newton's method | `ValueError` |
| Missing fixed-point function `g` in fixed-point / Steffensen | `ValueError` |
| Derivative or denominator too close to zero | `ZeroDivisionError` |
| No convergence after `max_iter` iterations | `RuntimeError` |

---

## File Structure

```
root-finding-project/
├── root_find.py      # RootFindingProblem class (all algorithms)
├── examples.py       # Interactive demo for all methods
└── README.md         # This file
```

---

## How to Run

**Requirements:** Python 3.x (no external packages needed)

```bash
python examples.py
```

An interactive menu will appear. Select a letter to run a method:

```
EXERCISE 24: POLYNOMIAL ROOT FINDING
  600x^4 - 550x^3 + 200x^2 - 20x - 1 = 0  on [0.1, 1]

  Available Methods:
    a) Bisection Method
    b) Newton's Method
    c) Secant Method
    d) False Position Method
    e) Muller's Method
    f) Fixed-Point Iteration
    g) Horner's Method (polynomial evaluation)
    h) Steffensen's Method
    i) Run all methods (comparison)
    q) Quit
```

Each method prompts for parameters (interval, initial guess, tolerance). Press Enter to accept defaults.

---

## Quick Code Example

```python
from root_find import RootFindingProblem

# Define the function and its derivative
def f(x):  return 600*x**4 - 550*x**3 + 200*x**2 - 20*x - 1
def df(x): return 2400*x**3 - 1650*x**2 + 400*x - 20

# Bisection
p = RootFindingProblem(f=f)
root = p.solve("bisection", a=0.1, b=1.0, tol=1e-4)
print(f"Bisection root: {root:.6f}")

# Newton's Method
p2 = RootFindingProblem(f=f, df=df)
root = p2.solve("newton", x0=0.5, tol=1e-4)
print(f"Newton root:    {root:.6f}")

# Secant Method
root = p.solve("secant", x0=0.1, x1=1.0, tol=1e-4)
print(f"Secant root:    {root:.6f}")

# Horner's Method (polynomial evaluation)
coeffs = [-1, -20, 200, -550, 600]   # [a0, a1, a2, a3, a4]
p3 = RootFindingProblem()
val = p3.solve("horner", coeffs=coeffs, x=0.5)
print(f"P(0.5) = {val:.6f}")

# Muller's Method (can find complex roots)
root = p.solve("muller", x0=0.1, x1=0.5, x2=1.0, tol=1e-4)
print(f"Muller root:    {root}")

# Fixed-Point Iteration  (g(x) = x - λ·f(x))
def g(x): return x - 0.001 * f(x)
p4 = RootFindingProblem(f=f, g=g)
root = p4.solve("fixed_point", x0=0.5, tol=1e-4)
print(f"Fixed-point:    {root:.6f}")

# Steffensen's Method
root = p4.solve("steffensen", x0=0.5, tol=1e-4)
print(f"Steffensen:     {root:.6f}")
```

---

## Constructor Reference

```python
RootFindingProblem(f=None, df=None, g=None)
```

| Parameter | Type | Required for |
|---|---|---|
| `f` | `Callable` | All methods except `horner` |
| `df` | `Callable` | `newton` only |
| `g` | `Callable` | `fixed_point`, `steffensen` |

```python
p.solve(method, **kwargs)
```

| Method string | Required kwargs |
|---|---|
| `"bisection"` | `a`, `b`, `tol` |
| `"false_position"` | `a`, `b`, `tol` |
| `"newton"` | `x0`, `tol` |
| `"secant"` | `x0`, `x1`, `tol` |
| `"muller"` | `x0`, `x1`, `x2`, `tol` |
| `"fixed_point"` | `x0`, `tol`, `max_iter` |
| `"steffensen"` | `x0`, `tol`, `max_iter` |
| `"horner"` | `coeffs`, `x` |
