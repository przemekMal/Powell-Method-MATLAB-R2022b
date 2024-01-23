% Usage:
%   [minimum, xes, iter] = powell_method(f, x0, d, N, epsilon, gold_step)
%
% Parameters:
%   - f: Objective function to minimize
%   - x0: Initial guess for the optimum
%   - d: search direction matrix
%   - epsilon: possible error
%   - gold_step: one-dimensional search interval using the Golden Section method
%
function [minimum, xes, iter] = powell_method(f, x0, d, N, epsilon, gold_step);

    if nargin < 1
        error("No function!!!");
    end
    if nargin < 2
        error("No entry point"); 
    end
    if (nargin < 3) || isempty(d)
        n = length(x0);
        d = eye(n); % def search direction matrix
    end
    if nargin < 4
        N = 20; % def num iteration
    end
    if nargin < 5
        epsilon = 1e-03; % def epsilon
    end
    if nargin < 6
        gold_step = 10; % def gold_step
    else
        gold_step = gold_step * (1+sqrt(5))/2;
    end

    p0 = x0;
    n = length(p0);
    p = zeros(n+1, n);
    p(1,:) = p0;
    minimum = ones(1, n) * (1000000000);
    
    i = 0;

    %Do wydruku
    xes = [];
    iter = 1;
    %Koniec wydruku
    while i< N
        i = i + 1;
        %Do wydruku
        xes = vertcat(xes, x0);
        %Koniec wydruku
        for r = 1:n
            g = @(myalpha)  f(p(r, :) + myalpha .* d(r, :));     
            a = p(r,:) - gold_step * (sqrt(5)-1) / 2*i;       b = p(r,:) + gold_step * (sqrt(5)-1) / 2*i;
            alpha_min = zloty_podzial(g, a, b, epsilon / 1000);
            p(r + 1,:) = p(r, :) + alpha_min .* d(r, :);
        end
         % termination condition check
        if (norm(f(p(n+1,:)) - f(x0)) <= epsilon) || (n == 1)
            xes = vertcat(xes, p(n+1,:));
            iter = i;
            if f(minimum) > f(p(n+1,:))
                minimum = p(n + 1, :);
            end     
            return 
        end
        
        %elimination of linear dependence of basis vectors
        idxMax = 2; delta = f(p(1,:)) - f(p(2,:));
        for m = idxMax:n+1
            %if abs(f(p(m-1,:)) - f(p(m,:))) > delta
            if (f(p(m-1,:)) - f(p(m,:))) > delta
                delta = f(p(m-1,:)) - f(p(m,:));
                idxMax = m;
            end
        end
        idxMax = idxMax - 1;
        f_1 = f(x0);
        f_2 = f(p(n+1,:));
        f_3 = f(2 * p(n+1,:) - x0);

        if (f_3 >= f_1)  || ((f_1 - 2*f_2 + f_3)*(f_1 - f_2 - delta)^2 >= (0.5 * delta * (f_1 - f_3)^2))
            % Updating search directions
            if f(2* p(n+1,:) - x0) < f(p(n+1,:))
                x0 = 2* p(n+1,:) - x0;
            else
                x0 = p(n+1,:);
            end
            p(1,:) = x0;
        else
            temp_d = p(n+1,:) - x0;
            for u = idxMax:n-1
                d(u,:) = d(u+1,:);
            end
            d(n,:) = temp_d;
            %Determining the minimum of the function 𝑓 along the direction 𝒅𝑛(𝑖+1)
            g = @(myalpha)  f(p(n+1, :) + myalpha .* temp_d);
            a = p(n+1,:) - gold_step * (sqrt(5)-1) / 2*i;  b = p(n+1,:) + gold_step * (sqrt(5)-1) / 2*i;
            alpha_min = zloty_podzial(g, a, b, epsilon);
            x0 = p(n+1, :) + alpha_min .* d(n, :);
            p(1,:) = x0;
        end
        if f(minimum) > f(x0)
            minimum = x0;
        end
    end
    iter = N;
end
% Author: przemekMal
% Date: 23.01.2024
%
% GitHub repository: https://github.com/przemekMal
