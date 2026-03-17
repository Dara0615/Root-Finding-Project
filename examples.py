"""
Exercise 24: Polynomial Root Finding

Problem: Find a solution in [0.1, 1] accurate to within 10^-4 for:
    600x^4 - 550x^3 + 200x^2 - 20x - 1 = 0

Methods available:
a) Bisection method
b) Newton's method
c) Secant method
d) Method of False Position
e) Muller's method
f) Fixed-Point iteration
g) Horner's method (for polynomial evaluation)
h) Steffensen's method
i) Run all methods (comparison)
"""

from root_find import RootFindingProblem
import cmath


def print_header(title):
    """Print formatted section header."""
    print("\n" + "=" * 70)
    print(f"  {title}")
    print("=" * 70)


def print_result(method, result, error=None):
    """Print formatted result."""
    if error:
        print(f"  {method:25s}: ERROR - {error}")
    else:
        if isinstance(result, complex):
            print(f"  {method:25s}: {result}")
        else:
            print(f"  {method:25s}: {result:.10f}")


def define_polynomial():
    """
    Define the polynomial P(x) = 600x^4 - 550x^3 + 200x^2 - 20x - 1
    and its derivative using Horner's method.
    """
    # Coefficients for P(x) = -1 - 20x + 200x^2 - 550x^3 + 600x^4
    # Ordered as [a0, a1, a2, a3, a4] where P(x) = a0 + a1*x + a2*x^2 + ...
    coeffs = [-1, -20, 200, -550, 600]
    
    # Coefficients for P'(x) = -20 + 400x - 1650x^2 + 2400x^3
    deriv_coeffs = [-20, 400, -1650, 2400]
    
    return coeffs, deriv_coeffs


def evaluate_polynomial(coeffs, x):
    """Evaluate polynomial using Horner's method."""
    solver = RootFindingProblem()
    return solver.solve("horner", coeffs=coeffs, x=x)


#Bisection
def run_bisection():
  
    print_header("Bisection Method")
    print("  Polynomial: 600x^4 - 550x^3 + 200x^2 - 20x - 1 = 0")
    
    coeffs, _ = define_polynomial()
    
    def f(x):
        return evaluate_polynomial(coeffs, x)
    
    # Get user input
    try:
        a = float(input("  Enter left endpoint a [0.1]: ") or "0.1")
        b = float(input("  Enter right endpoint b [1.0]: ") or "1.0")
        tol = float(input("  Enter tolerance [1e-4]: ") or "1e-4")
    except ValueError:
        print("  Invalid input. Using defaults.")
        a, b, tol = 0.1, 1.0, 1e-4
    
    print(f"\n  Using interval [{a}, {b}] with tolerance {tol}")
    
    # Verify sign change
    f_a, f_b = f(a), f(b)
    print(f"  f({a}) = {f_a:.6f}")
    print(f"  f({b}) = {f_b:.6f}")
    
    if f_a * f_b > 0:
        print(f"  Warning: f(a) and f(b) have same sign. Method may fail.")
        return
    
    solver = RootFindingProblem(f=f)
    
    try:
        root = solver.solve("bisection", a=a, b=b, tol=tol)
        f_root = f(root)
        print_result("Bisection", root)
        print(f"  Verification: f({root:.6f}) = {f_root:.2e}")
    except Exception as e:
        print_result("Bisection", None, error=str(e))

#Newton
def run_newton():
   
    print_header("Newton's Method")
    print("  Polynomial: 600x^4 - 550x^3 + 200x^2 - 20x - 1 = 0")
    
    coeffs, deriv_coeffs = define_polynomial()
    
    def f(x):
        return evaluate_polynomial(coeffs, x)
    
    def df(x):
        return evaluate_polynomial(deriv_coeffs, x)
    
    # Get user input
    try:
        x0 = float(input("  Enter initial guess x0 [0.5]: ") or "0.5")
        tol = float(input("  Enter tolerance [1e-4]: ") or "1e-4")
    except ValueError:
        print("  Invalid input. Using defaults.")
        x0, tol = 0.5, 1e-4
    
    print(f"\n  Using initial guess x0 = {x0} with tolerance {tol}")
    print(f"  f({x0}) = {f(x0):.6f}")
    print(f"  f'({x0}) = {df(x0):.6f}")
    
    solver = RootFindingProblem(f=f, df=df)
    
    try:
        root = solver.solve("newton", x0=x0, tol=tol)
        f_root = f(root)
        print_result("Newton", root)
        print(f"  Verification: f({root:.6f}) = {f_root:.2e}")
    except Exception as e:
        print_result("Newton", None, error=str(e))

#Secant
def run_secant():
    
    print_header("Secant Method")
    print("  Polynomial: 600x^4 - 550x^3 + 200x^2 - 20x - 1 = 0")
    
    coeffs, _ = define_polynomial()
    
    def f(x):
        return evaluate_polynomial(coeffs, x)
    
    # Get user input
    try:
        x0 = float(input("  Enter first initial guess x0 [0.1]: ") or "0.1")
        x1 = float(input("  Enter second initial guess x1 [1.0]: ") or "1.0")
        tol = float(input("  Enter tolerance [1e-4]: ") or "1e-4")
    except ValueError:
        print("  Invalid input. Using defaults.")
        x0, x1, tol = 0.1, 1.0, 1e-4
    
    print(f"\n  Using initial guesses x0 = {x0}, x1 = {x1} with tolerance {tol}")
    print(f"  f({x0}) = {f(x0):.6f}")
    print(f"  f({x1}) = {f(x1):.6f}")
    
    solver = RootFindingProblem(f=f)
    
    try:
        root = solver.solve("secant", x0=x0, x1=x1, tol=tol)
        f_root = f(root)
        print_result("Secant", root)
        print(f"  Verification: f({root:.6f}) = {f_root:.2e}")
    except Exception as e:
        print_result("Secant", None, error=str(e))

#False_position
def run_false_position():
  
    print_header("False Position Method")
    print("  Polynomial: 600x^4 - 550x^3 + 200x^2 - 20x - 1 = 0")
    
    coeffs, _ = define_polynomial()
    
    def f(x):
        return evaluate_polynomial(coeffs, x)
    
    # Get user input
    try:
        a = float(input("  Enter left endpoint a [0.1]: ") or "0.1")
        b = float(input("  Enter right endpoint b [1.0]: ") or "1.0")
        tol = float(input("  Enter tolerance [1e-4]: ") or "1e-4")
    except ValueError:
        print("  Invalid input. Using defaults.")
        a, b, tol = 0.1, 1.0, 1e-4
    
    print(f"\n  Using interval [{a}, {b}] with tolerance {tol}")
    
    # Verify sign change
    f_a, f_b = f(a), f(b)
    print(f"  f({a}) = {f_a:.6f}")
    print(f"  f({b}) = {f_b:.6f}")
    
    if f_a * f_b > 0:
        print(f"  Warning: f(a) and f(b) have same sign. Method may fail.")
        return
    
    solver = RootFindingProblem(f=f)
    
    try:
        root = solver.solve("false_position", a=a, b=b, tol=tol)
        f_root = f(root)
        print_result("False Position", root)
        print(f"  Verification: f({root:.6f}) = {f_root:.2e}")
    except Exception as e:
        print_result("False Position", None, error=str(e))

#Muller
def run_muller():
   
    print_header("Muller's Method")
    print("  Polynomial: 600x^4 - 550x^3 + 200x^2 - 20x - 1 = 0")
    
    coeffs, _ = define_polynomial()
    
    def f(x):
        return evaluate_polynomial(coeffs, x)
    
    # Get user input
    try:
        x0 = float(input("  Enter first point x0 [0.1]: ") or "0.1")
        x1 = float(input("  Enter second point x1 [0.5]: ") or "0.5")
        x2 = float(input("  Enter third point x2 [1.0]: ") or "1.0")
        tol = float(input("  Enter tolerance [1e-4]: ") or "1e-4")
    except ValueError:
        print("  Invalid input. Using defaults.")
        x0, x1, x2, tol = 0.1, 0.5, 1.0, 1e-4
    
    print(f"\n  Using points x0 = {x0}, x1 = {x1}, x2 = {x2} with tolerance {tol}")
    print(f"  f({x0}) = {f(x0):.6f}")
    print(f"  f({x1}) = {f(x1):.6f}")
    print(f"  f({x2}) = {f(x2):.6f}")
    
    solver = RootFindingProblem(f=f)
    
    try:
        root = solver.solve("muller", x0=x0, x1=x1, x2=x2, tol=tol)
        if isinstance(root, complex) and abs(root.imag) < 1e-10:
            root = root.real
        
        f_root = f(root)
        print_result("Muller", root)
        print(f"  Verification: f({root:.6f}) = {f_root:.2e}")
    except Exception as e:
        print_result("Muller", None, error=str(e))

#Fixed_point
def run_fixed_point():
   
    print_header("Fixed-Point Iteration")
    print("  Polynomial: 600x^4 - 550x^3 + 200x^2 - 20x - 1 = 0")
    print("  Reformulation: x = x - lambda * P(x)")
    
    coeffs, _ = define_polynomial()
    
    def f(x):
        return evaluate_polynomial(coeffs, x)
    
    # Get user input
    try:
        x0 = float(input("  Enter initial guess x0 [0.5]: ") or "0.5")
        lam = float(input("  Enter relaxation parameter lambda [0.001]: ") or "0.001")
        tol = float(input("  Enter tolerance [1e-4]: ") or "1e-4")
        max_iter = int(input("  Enter max iterations [1000]: ") or "1000")
    except ValueError:
        print("  Invalid input. Using defaults.")
        x0, lam, tol, max_iter = 0.5, 0.001, 1e-4, 1000
    
    # Define g(x) = x - lambda * f(x)
    def g(x):
        return x - lam * f(x)
    
    print(f"\n  Using x0 = {x0}, lambda = {lam}, tolerance {tol}")
    print(f"  g(x) = x - {lam} * P(x)")
    print(f"  g({x0}) = {g(x0):.6f}")
    
    solver = RootFindingProblem(f=f, g=g)
    
    try:
        root = solver.solve("fixed_point", x0=x0, tol=tol, max_iter=max_iter)
        f_root = f(root)
        print_result("Fixed-Point", root)
        print(f"  Verification: f({root:.6f}) = {f_root:.2e}")
    except Exception as e:
        print_result("Fixed-Point", None, error=str(e))

#Stefensen
def run_steffensen():
  
    print_header("Steffensen's Method")
    print("  Polynomial: 600x^4 - 550x^3 + 200x^2 - 20x - 1 = 0")
    print("  Reformulation: x = x - lambda * P(x)")
    
    coeffs, _ = define_polynomial()
    
    def f(x):
        return evaluate_polynomial(coeffs, x)
    
    # Get user input
    try:
        x0 = float(input("  Enter initial guess x0 [0.5]: ") or "0.5")
        lam = float(input("  Enter relaxation parameter lambda [0.001]: ") or "0.001")
        tol = float(input("  Enter tolerance [1e-4]: ") or "1e-4")
        max_iter = int(input("  Enter max iterations [100]: ") or "100")
    except ValueError:
        print("  Invalid input. Using defaults.")
        x0, lam, tol, max_iter = 0.5, 0.001, 1e-4, 100
    
    # Define g(x) = x - lambda * f(x)
    def g(x):
        return x - lam * f(x)
    
    print(f"\n  Using x0 = {x0}, lambda = {lam}, tolerance {tol}")
    print(f"  g(x) = x - {lam} * P(x)")
    print(f"  g({x0}) = {g(x0):.6f}")
    
    solver = RootFindingProblem(f=f, g=g)
    
    try:
        root = solver.solve("steffensen", x0=x0, tol=tol, max_iter=max_iter)
        f_root = f(root)
        print_result("Steffensen", root)
        print(f"  Verification: f({root:.6f}) = {f_root:.2e}")
    except Exception as e:
        print_result("Steffensen", None, error=str(e))

#Horner
def run_horner():
   
    print_header("Horner's Method")
    print("  Polynomial: P(x) = 600x^4 - 550x^3 + 200x^2 - 20x - 1")
    print("  Coefficients: [-1, -20, 200, -550, 600]")
    
    coeffs, deriv_coeffs = define_polynomial()
    
    # Get user input
    try:
        x = float(input("  Enter x value to evaluate [0.5]: ") or "0.5")
    except ValueError:
        print("  Invalid input. Using default.")
        x = 0.5
    
    solver = RootFindingProblem()
    
    # Evaluate P(x)
    result = solver.solve("horner", coeffs=coeffs, x=x)
    print(f"\n  P({x}) = {result:.6f}")
    
    # Evaluate P'(x)
    deriv_result = solver.solve("horner", coeffs=deriv_coeffs, x=x)
    print(f"  P'({x}) = {deriv_result:.6f}")
    
    # Verify by direct computation
    direct = 600*x**4 - 550*x**3 + 200*x**2 - 20*x - 1
    print(f"\n  Direct computation: {direct:.6f}")
    print(f"  Difference: {abs(result - direct):.2e}")

#All Methods
def run_all_comparison():
   
    print_header("Comparison of All Methods")
    print("  Problem: 600x^4 - 550x^3 + 200x^2 - 20x - 1 = 0")
    print("  Using default parameters")
    
    coeffs, deriv_coeffs = define_polynomial()
    
    def f(x):
        return evaluate_polynomial(coeffs, x)
    
    def df(x):
        return evaluate_polynomial(deriv_coeffs, x)
    
    def g(x):
        return x - 0.001 * f(x)
    
    solver_full = RootFindingProblem(f=f, df=df, g=g)
    
    methods = [
        ("Bisection", "bisection", {"a": 0.1, "b": 1.0, "tol": 1e-4}),
        ("False Position", "false_position", {"a": 0.1, "b": 1.0, "tol": 1e-4}),
        ("Newton", "newton", {"x0": 0.5, "tol": 1e-4}),
        ("Secant", "secant", {"x0": 0.1, "x1": 1.0, "tol": 1e-4}),
        ("Muller", "muller", {"x0": 0.1, "x1": 0.5, "x2": 1.0, "tol": 1e-4}),
        ("Fixed-Point", "fixed_point", {"x0": 0.5, "tol": 1e-4, "max_iter": 1000}),
        ("Steffensen", "steffensen", {"x0": 0.5, "tol": 1e-4, "max_iter": 100}),
    ]
    
    print(f"\n  {'Method':<20} {'Root':<15} {'f(Root)':<12} {'Status'}")
    print("  " + "-" * 60)
    
    results = []
    for name, method, params in methods:
        try:
            root = solver_full.solve(method, **params)
            if isinstance(root, complex) and abs(root.imag) < 1e-10:
                root = root.real
            
            f_val = f(root)
            status = "OK" if abs(f_val) < 1e-3 else "CHECK"
            print(f"  {name:<20} {root:<15.6f} {f_val:<12.2e} {status}")
            results.append((name, root))
        except Exception as e:
            print(f"  {name:<20} {'ERROR':<15} {'N/A':<12} {str(e)[:15]}")
    
    # Calculate statistics
    real_roots = [r for _, r in results if isinstance(r, (int, float))]
    if len(real_roots) > 1:
        avg = sum(real_roots) / len(real_roots)
        print(f"\n  Average of successful methods: {avg:.6f}")


def display_menu():
    """Display the main menu."""
    print("\n" + "=" * 70)
    print("EXERCISE 24: POLYNOMIAL ROOT FINDING")
    print("  600x^4 - 550x^3 + 200x^2 - 20x - 1 = 0  on [0.1, 1]")
    print("=" * 70)
    print("\n  Available Methods:")
    print("    a) Bisection Method")
    print("    b) Newton's Method")
    print("    c) Secant Method")
    print("    d) False Position Method")
    print("    e) Muller's Method")
    print("    f) Fixed-Point Iteration")
    print("    g) Horner's Method (polynomial evaluation)")
    print("    h) Steffensen's Method")
    print("    i) Run all methods (comparison)")
    print("    q) Quit")
    print("-" * 70)


def main():
    """Main interactive loop."""
    menu_options = {
        'a': run_bisection,
        'b': run_newton,
        'c': run_secant,
        'd': run_false_position,
        'e': run_muller,
        'f': run_fixed_point,
        'g': run_horner,
        'h': run_steffensen,
        'i': run_all_comparison,
    }
    
    while True:
        display_menu()
        choice = input("  Select a method (a-i, q to quit): ").strip().lower()
        
        if choice == 'q':
            print("\n  Goodbye!")
            break
        elif choice in menu_options:
            print()  # Empty line for spacing
            try:
                menu_options[choice]()
            except Exception as e:
                print(f"\n  An error occurred: {e}")
            input("\n  Press Enter to continue...")
        else:
            print("\n  Invalid choice. Please select a letter from a to i, or q to quit.")
            input("  Press Enter to continue...")


if __name__ == "__main__":
    main()