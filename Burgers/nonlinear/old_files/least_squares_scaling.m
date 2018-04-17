function coefficients = least_squares_scaling(RHS,nonlin0,nonlin1,nonlin2,nonlin3,t,N)

resid = @(coeffs) sum((RHS-nonlin0-(coeffs(1)*(1/(N/2-coeffs(2)))*t*nonlin1...
                           +coeffs(1)*(1/(N/2-coeffs(2)))*(1/(N/2-coeffs(3)))*t^2/2*nonlin2...
                           -coeffs(1)*(1/(N/2-coeffs(2)))*(1/(N/2-coeffs(3)))*(1/(N/2-coeffs(4)))*t^3/6*nonlin3)).^2);
                       
coeffs = fminsearch(resid,[N/2;0;0;0]);

coefficients = [coeffs(1)*(1/(N/2-coeffs(2)));
                coeffs(1)*(1/(N/2-coeffs(2)))*(1/(N/2-coeffs(3)))
                coeffs(1)*(1/(N/2-coeffs(2)))*(1/(N/2-coeffs(3)))*(1/(N/2-coeffs(4)))];
           