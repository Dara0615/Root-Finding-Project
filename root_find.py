import cmath
import math

class RootFindingProblem:
    def __init__(self, f = None, df = None, g = None):
        self.f = f    # f : Callable - Function whose root is sought to (f(x) = 0)
        self.df = df  # df : Callable - Derivative of f(x) - required in Newton's method
        self.g = g    # g : Callable - Function used in fixed-point g(x) where x = g(x) - required in fixed-point iteration method
        
    def solve(self, method, **kwargs):
        methods = {
            'bisection' : self.__bisection,
            'fixed_point' : self.__fixed_point,
            'newton' : self.__newton,
            'secant' : self.__secant,
            'false_position' : self.__false_position,
            'horner' : self.__horner,
            'muller' : self.__muller,
            'steffensen' : self.__steffensen
        }
        
        if method not in methods:
            raise ValueError(
                f"Unknown method: '{method}'."
                f"Available methods: {list(methods.keys())}"
            )
        
        if self.f is None and method not in ['fixed_point', 'horner']:
            raise ValueError ("function f(x) must be provided for this method")
        
        # Validate g is provided for fixed-point methods
        if self.g is None and method in ['fixed_point', 'steffensen']:
            raise ValueError(
                f"method '{method}' requires fixed-point function g(x). "
                f"Initialize with RootFindingProblem(g=g)"
            )
        
        return methods[method](**kwargs)
    
    # Bisection Method
    
    def __bisection(self, a, b, tol = 1e-6, max_iter = 1000):
        fa, fb = self.f(a), self.f(b)
        
        if fa * fb > 0:
            raise ValueError(
                f"Invalid interval: f({a}) = {fa:.6e} and f({b}) = {fb:.6e}"
                f"must have opposite signs for bisection method"
            )
        
        for i in range (max_iter):
            c = (a + b) / 2.0
            fc = self.f(c)
            
            if abs (fc) < tol or (b - a) / 2 < tol:
                return c
            
            if fa * fc < 0:
                b = c
                fb = fc
            else:
                a = c
                fa = fc
                
        raise RuntimeError(
            f"Bisection method failed to converge after {max_iter} iterations. "
            f"Last interval: [{a}, {b}]"
        )
    
    # Fixed-point method
    
    def __fixed_point(self, x0, tol = 1e-6, max_iter = 1000):
        # Validation now handled in solve(), but keep as safety check
        if self.g is None:
            raise ValueError(
                "Fixed-point function g(x) not provided. "
                "Initialize with RootFindingProblem(f=f, g=g)"
            )
        
        x = x0
        
        for i in range (max_iter):
            x_new = self.g(x)
            
            if abs (x_new - x) < tol:
                return x_new
            
            x = x_new
            
        raise RuntimeError(
            f"Fixed-point iteration failed to converge after {max_iter} iterations. "
            f"Last value: x = {x}"
        )
        
        
    # Newton method
    
    def __newton(self, x0, tol = 1e-6, max_iter= 1000):
        if self.df is None:
            raise ValueError(
                "Derivative df(x) not provided for Newton's method. "
                "Initialize with RootFindingProblem (f=f, df=df)"
            )
        
        x = x0
        
        for i in range (max_iter):
            fx = self.f(x)
            
            if abs(fx) < tol:
                return x
            
            dfx = self.df(x)
            
            if abs(dfx) < 1e-15:
                raise ZeroDivisionError(
                    f"Derivative too close to zero at x = {x}. "
                    f"f'(x) = {dfx}"
                )
                
            x_new = x - fx / dfx
            
            if abs (x_new - x) < tol:
                return x_new
            
            x = x_new
            
        raise RuntimeError(
            f"Newton's method failed to converge after {max_iter} iterations. "
            f"Last value: x = {x}"
        )
        
    # Secant method
    
    def __secant(self, x0, x1, tol = 1e-6, max_iter = 1000):
        for i in  range (max_iter):
            fx0, fx1 = self.f(x0), self.f(x1)
            
            if abs(fx1) < tol:
                return x1
            
            denom = fx1 - fx0
            
            if abs(denom) < 1e-15:
                raise ZeroDivisionError(
                    f"Derivative by zero in secant method: f({x1}) = f({x0}) = {fx1}"
                )
            
            x_new = x1 - fx1 * (x1 - x0) / denom
            
            if abs (x_new - x1) < tol:
                return x_new
            
            x0, x1 = x1, x_new
            
        raise RuntimeError(
            f"Secant method failed to converge after {max_iter} iterations. "
            f"Last value: x0 = {x0}, x1 = {x1}"
        )
        
        
    # False_position
    
    def __false_position (self, a, b, tol = 1e-6, max_iter = 1000):
        fa, fb = self.f(a), self.f(b)
        
        if fa * fb > 0:
            raise ValueError(
                f"Invalid interval: f({a}) = {fa:.6e} and f({b}) = {fb:.6e} "
                f"must have opposite signs for false position method"
            )
            
        for i in range (max_iter):
            # Linear interpolation to find c where line crosses x-axis
            
            c = (a * fb - b * fa) / (fb - fa)
            fc = self.f(c)
            
            if abs(fc) < tol or abs(b - a) < tol:
                return c

            if fa * fc < 0:
                b = c
                fb = fc
            else:
                a = c
                fa = fc
            
        raise RuntimeError(
            f"False position method failed to converge after {max_iter} iterations. "
            f"Last interval: [{a}, {b}]"
        )
        
    # Horner method
    
    def __horner(self, coeffs, x):
        result = 0
        
        for coeff in reversed (coeffs):
            result = result *  x + coeff
        return result
    
    # Muller method
    
    def __muller (self, x0, x1, x2, tol = 1e-6, max_iter = 1000):
        
        x0, x1, x2 = complex(x0), complex(x1), complex(x2)
        
        for i in range (max_iter):
            f0, f1, f2 = self.f(x0), self.f(x1), self.f(x2)
            
            h0 = x1 - x0
            h1 = x2 - x1
            d0 = (f1 - f0) / h0 if h0 != 0 else 0
            d1 = (f2 - f1) / h1 if h1 != 0 else 0
            
            a = (d1 - d0) / (h1 + h0) if (h1 + h0) != 0 else 0
            b = a * h1 + d1
            c = f2
            
            discriminant = cmath.sqrt(b**2 - 4*a*c)
            
            if abs (b + discriminant) > abs (b - discriminant):
                denom = b + discriminant
            else:
                denom = b - discriminant
                
            if abs (denom) < 1e-15:
                raise ZeroDivisionError(
                    f"Division by zero in Muller's method at iteration {i}"
                )
            
            dx = -2 * c / denom
            x_new = x2 + dx
            
            if abs(dx) < tol:
                return x_new
            
            x0, x1, x2 = x1, x2, x_new
            
        raise RuntimeError(
            f"Muller's method failed to converge after {max_iter} iterations. "
            f"Last value: x = {x2}"
        )
    
    # Steffensen
    
    def __steffensen(self, x0, tol=1e-6, max_iter = 1000):
        # Validation now handled in solve(), but keep as safety check
        if self.g is None:
            raise ValueError(
                "Fixed-point function g(x) not provided for Steffensen's method. "
                "Initialize with RootFindingProblem(f=f, g=g)"
            )
            
        x = x0
        
        for i in range(max_iter):
            # Compute g(x) and g(g(x))
            gx = self.g(x)
            ggx = self.g(gx)
            
            # Aitken's Δ² formula
            denom = ggx - 2 * gx + x
            
            if abs(denom) < 1e-15:
                raise ZeroDivisionError(
                    f"Division by zero in Steffensen's method at x = {x}"
                )

            x_new = x - (gx - x)**2 / denom
            
            if abs(x_new - x) < tol:
                return x_new
            
            x = x_new
        
        raise RuntimeError(
            f"Steffensen's method failed to converge after {max_iter} iterations. "
            f"Last value: x = {x}"
        )