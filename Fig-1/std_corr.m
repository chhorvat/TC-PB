function corrfac = std_corr(N)

corrfac = ((N/2).^(1/2)) .* gamma((N-1)/2) .* (1./gamma(N/2)); 
