function [y] = camel6(xx)
    x1 = (xx(:, 1) - 0.5) * 6;  
    x2 = (xx(:, 2) - 0.5) * 4; 
    term1 = (4 - 2.1 * x1.^2 + (x1.^4) / 3) .* x1.^2;
    term2 = x1 .* x2;
    term3 = (-4 + 4 * x2.^2) .* x2.^2;
    y = term1 + term2 + term3;
end
