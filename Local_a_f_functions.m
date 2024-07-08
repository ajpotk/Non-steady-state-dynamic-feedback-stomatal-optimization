function [a_f_func, da_f_dpi_L_func] = Local_a_f_functions(a_f_max, pi_L_star, beta)

a_f_func = @(pi_L) a_f_max * exp(-((pi_L/pi_L_star) .^ beta));
da_f_dpi_L_func = @(pi_L) -beta * (pi_L/pi_L_star).^beta .* a_f_func(pi_L)./pi_L;

end

