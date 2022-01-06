function fit = base_functions(x, R, o, o_local, idx)

% persistent fhd

% 1. Unimodal Functions
if     (idx ==  1) fhd = str2func('f1');   % Sphere Function
elseif (idx ==  2) fhd = str2func('f2');
elseif (idx ==  3) fhd = str2func('f3');
elseif (idx ==  4) fhd = str2func('f4');
% 2. Multi-Modal Functions
%    2.1. Strong global structure:
elseif (idx ==  5) fhd = str2func('f5');
elseif (idx ==  6) fhd = str2func('f6');
elseif (idx ==  7) fhd = str2func('f7');
elseif (idx ==  8) fhd = str2func('f8');
%    2.2. Moderate or weak global structure:
elseif (idx ==  9) fhd = str2func('f9');
elseif (idx == 10) fhd = str2func('f10');
elseif (idx == 11) fhd = str2func('f11');
elseif (idx == 12) fhd = str2func('f12');
    %3. Overlapping Functions
elseif (idx == 13) fhd = str2func('f13');
elseif (idx == 14) fhd = str2func('f14');
elseif (idx == 15) fhd = str2func('f15');
    % 4. Composition Functions
elseif (idx == 16) fhd = str2func('f16');
end

if idx ~= 16
    fit = feval(fhd, R, x);
else
    fit = feval(fhd, R, o, o_local, x);
end
end

%% ------------------------- BASE FUNCTIONS ------------------------

%--------------------------------------------------------------------------
% f1 Sphere Function 
%   The simplest base function
%--------------------------------------------------------------------------
function fit = f1(R, x)
    % Done
    fit = sum(x.*x, 1);
end


%--------------------------------------------------------------------------
% f2 Elliptic Function
%   In comparison to f1, the feature of variable scaling is introduced;
%--------------------------------------------------------------------------
function fit = f2(R, x)
    % Done
    [D, ps] = size(x);
    condition = 1e+4;
    coefficients = condition .^ linspace(0, 1, D); 
    fit = coefficients * x.^2; 
end

%--------------------------------------------------------------------------
% f3 Rotated Elliptic Function
%   In comparison to f2, the feature of non-separable variable is introduced;
%--------------------------------------------------------------------------
function fit = f3(R, x)
    % Done

    [D ps] = size(x);
    x = R * x;
    condition = 1e+4;
    coefficients = condition .^ linspace(0, 1, D); 
    fit = coefficients * x.^2; 
end

%--------------------------------------------------------------------------
% f4 Step Elliptic Function
%   In comparison to f3, the feature of plateaus is introduced;
%--------------------------------------------------------------------------
function fit = f4(R, x)
    % Done
    [D, ps] = size(x);
    condition = 1e+4;
    coefficients = condition .^ linspace(0, 1, D);
    
    z_hat = R * x;
    z = zeros(D, ps);
    z(abs(z_hat) > 0.5) = floor(z_hat(abs(z_hat) > 0.5) + 0.5);
    z(abs(z_hat) <= 0.5) = 0.1 * floor(10 * z_hat(abs(z_hat) <= 0.5) + 0.5);
    
    fit = max(abs(z(1, :))/1e6, coefficients(2 : end) * z(2 : end, :).^2);
end

%------------------------------------------------------------------------------
% f5 Rastrigin's Function, Original
%   In comparison to f2, the feature of multi-modal is introduced;
%------------------------------------------------------------------------------
function fit = f5(R, x)
    % Done
    [D, ps] = size(x);
    x = 0.0512 * x;
    fit = 10*(D - sum(cos(2*pi*x), 1)) + sum(x.^2, 1);

    fit = 250 * fit;
end

%------------------------------------------------------------------------------
% f6 Buche-Rastrigin's Function
%   In comparison to f6, the feature of asymmetric is introduced;
%------------------------------------------------------------------------------
function fit = f6(R, x)
    % Done
    [D, ps] = size(x);
    
    x = 0.0512 * x ;
    
    condition = sqrt(10);
    coefficients = condition .^ linspace(0, 1, D);
    coefficients = repmat(coefficients', 1, ps);

    b_1 = zeros(D, ps);
    b_1(1:2:end, :) = 1;
    b_2 = zeros(D, ps);
    b_2(x > 0) = 1;
    b_2(x <= 0) = 0;
    
    b = b_1 & b_2;
    
    coefficients(b) = 10 * coefficients(b);

    x = coefficients .* x;
    fit = 10*(D - sum(cos(2*pi*x), 1)) + sum(x.^2, 1);
    
    fit = 125 * fit;
end

%------------------------------------------------------------------------------
% f7 Griewank's Function
%------------------------------------------------------------------------------
function fit = f7(R, x)
    % Done
    [D ps] = size(x);
    
    x = R * x;
    x = T_diag(x, 100);
    i = repmat([1:D]', 1, ps);
    fit_1 = sum(x.^2/4000, 1);
    fit_2 = prod(cos(x./ sqrt(i) ), 1);
    
    fit = 1e4 * (fit_1 - fit_2 + 1);
end

%------------------------------------------------------------------------------
% f8 Rotated Rastrigin's Function
%   In comparison to f6, the feature of non-separable variable is introduced;
%------------------------------------------------------------------------------
function fit = f8(R, x)
    % Done
    [D, ps] = size(x);
    x = 0.0512 * x;
    x = R * x;
    x = T_diag(T_asy(T_irreg(x), 0.2), 10);
    fit = 10*(D - sum(cos(2*pi*x), 1)) + sum(x.^2, 1);
    
    fit = 250 * fit;
end

%------------------------------------------------------------------------------
% f9 Weierstrass's Function
%------------------------------------------------------------------------------
function fit = f9(R, x)
    % Done
    [D ps] = size(x);
    
    a = 0.5;
    b = 3;
    k_max = 20;
    
    x = R * (5e-3 * x);
    
    k = (0 : k_max)';
    
    temp_1 = 0;
    for i = 1 : D
        temp_1 = temp_1 + sum(a.^k .* cos(2 * pi * b.^k .* repmat((x(i, :) + 0.5), k_max+1, 1)));
    end
    temp_2 = D * sum(a.^k .* cos(pi * b.^k));
    
    fit = 1e4 * (temp_1 - temp_2);
end


%------------------------------------------------------------------------------
% f12 Ackley's Function
%------------------------------------------------------------------------------
function fit = f10(R, x)
    % Done
    [D ps] = size(x);
    a = 3e4;
    b = 0.5;
    
    x = R * 0.32 * x;
    
    fit = sum(x.^2,1);
    fit = a - a.*exp(-b.*sqrt(fit./D))-exp(sum(cos(2.*pi.*x),1)./D)+exp(1);
end


%------------------------------------------------------------------------------
% f9 Lunacek Bi-Rastrigin's Function
%------------------------------------------------------------------------------
function fit = f11(R, x)
    % Done
    [D ps] = size(x);
    
    x = T_diag(R * 0.0512 * x, 100);

    mu0 = 2.5;
    s = 1 - (1/(2*(sqrt(D+20))-8.2));
    
    d = 1;
    mu1 = -sqrt((mu0^2-d)/s);
    
    fit = min(sum((x - mu0).^2, 1), d*D + s*sum((x - mu1).^2, 1)) + ...
      10 * (D - sum(cos(2*pi*(x - mu0)), 1)) + sum(x.^2, 1);

    fit = 250 * fit;
end


%------------------------------------------------------------------------------
% f11 Schwefel's Function
%------------------------------------------------------------------------------
function fit = f12(R, x)
    % Done
    [D, ps] = size(x);
    
    z = R * 10 * x;         % TODO
    
    z = z + 420.9687462275036;
    
    fit = zeros(D, ps);
    
    fit(abs(z) <= 500) = z(abs(z) <= 500).* sin(sqrt(abs(z(abs(z) <= 500))));
    fit(z > 500) = (500 - mod(z(z > 500), 500)) .* sin(sqrt(abs(500 - mod(z(z > 500), 500)))) ...
        - ((z(z > 500)-500).^2)/(10000 *D);
    fit(z < -500) = (mod(z(z < -500), 500) - 500) .* sin(sqrt(abs(mod(z(z < -500), 500) - 500))) ...
        - ((z(z < -500) + 500).^2)/(10000 *D);
    
    fit = 30 * (418.982887 * D - sum(fit, 1));
end

%------------------------------------------------------------------------------
% f13 Rosenbrock's Function
%------------------------------------------------------------------------------
function fit = f13(R, x)
    % Done
    [D ps] = size(x);
    x = 0.05 * x + 1;      
    
    fit = sum(100.*(x(1:D-1,:).^2-x(2:D, :)).^2+(x(1:D-1, :)-1).^2, 1);
end

%------------------------------------------------------------------------------
% f14 The Expanded Rosenbrock's plus Grewangk's Function
%------------------------------------------------------------------------------
function fit = f14(R, x)
    % Done
    [D ps] = size(x);
    x = 0.05 * x + 1;
    fit = 0;
    R = 1;
    for i = 1 : D-1
        temp_fit = f13(R, x(i:i+1, :));
        fit = fit + f7(R, temp_fit);
    end
    
    temp_fit = f13(R, x([D, 1], :));
    fit = fit + f7(R, temp_fit);
    
    fit = 1e-4 * fit;
end

%------------------------------------------------------------------------------
% f15 The Expanded Schaffer's F6 Function
%------------------------------------------------------------------------------
function fit = f15(R, x)
    % Done
    [D ps] = size(x);
    fit = 0;
    x = x + 1;
    for i = 1 : D-1
        fit = fit + schaffer_F6(x(i, :), x(i+1, :));
    end
    fit = fit + schaffer_F6(x(1, :), x(D, :));
    
    fit = 4e4 * fit;
end

function g = schaffer_F6(x, y)

    g = 0.5 + ((sin(sqrt(x.^2 + y.^2))).^2 - 0.5) ./ (1+0.001*(x.^2 + y.^2)).^2;

end

%------------------------------------------------------------------------------
% f16 The Composition Function
%------------------------------------------------------------------------------
function fit = f16(R, o_global, o_local, x)
    % Done
    [D ps] = size(x);
    w = zeros(3, ps);
    sigma = [10, 20, 30];
    for i = 1 : 3
        if i == 1
            o = o_global;
        else
            o = o_local(:, i-1);
        end
        w(i, :) = (1 ./ sqrt((sum((x - o).^2, 1)))) .* exp(-(sum((x - o).^2, 1)) ./ ...
            (2 * D * sigma(i)^2));
    end
    
    if ~isempty(find(w == Inf, 1))
        w(w == Inf) = 1;
        w(w ~= 1) = 0;
    else
        w = w ./ sum(w, 1);
    end
    
%     fit = f13(R, x);
    fit = w(1, :).* f8(R, x) * 10 + w(2, :).* (f10(R, x + o_global - o_local(:, 1)).*1e-3 + 100)...
        +  w(3, :).* (2*f12(R, x + o_global - o_local(:, 2)) + 200);
end

%------------------------------------------------------------------------------
% This transformation function is used to break the symmetry of symmetric 
% functions.
%------------------------------------------------------------------------------
function g = T_asy(f, beta)
    [D popsize] = size(f);
    g = f;
    temp = repmat(beta * linspace(0, 1, D)', 1, popsize); 
    ind = f > 0;
    g(ind) = f(ind).^ (1 + temp(ind) .* sqrt(f(ind)));  
end


%------------------------------------------------------------------------------
% This transformation is used to create the ill-conditioning effect.
%------------------------------------------------------------------------------
function g = T_diag(f, alpha)
    [D popsize] = size(f);
    scales = repmat(sqrt(alpha) .^ linspace(0, 1, D)', 1, popsize); 
    g = scales .* f;
end


%------------------------------------------------------------------------------
% This transformation is used to create smooth local irregularities.
%------------------------------------------------------------------------------
function g = T_irreg(f)
   a = 0.1;
   g = f; 
   idx = (f > 0);
   g(idx) = log(f(idx))/a;
   g(idx) = exp(g(idx) + 0.49*(sin(g(idx)) + sin(0.79*g(idx)))).^a;
   idx = (f < 0);
   g(idx) = log(-f(idx))/a;
   g(idx) = -exp(g(idx) + 0.49*(sin(0.55*g(idx)) + sin(0.31*g(idx)))).^a;
end

