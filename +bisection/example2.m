% Find roots by varying the function coefficients.

[A, B] = meshgrid(linspace(1,2,6), linspace(4,12,10));
f = @(x) A.*x.^0.2 + B.*x.^0.87 - 15;
xstar = bisection(f,0,5)