function coeffs = renormalized_ROM_coeffs(epsilon,N)

coeffs = zeros(4,1);
coeffs(2) = -1.212317206140178*epsilon^-3.691035548147793*N^-5.735647633869574;
coeffs(4) = 0.350302790566717*epsilon^-7.388110167590722*N^-11.471931044430516;