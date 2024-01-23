function minimum_x = zloty_podzial(f, a, b, epsilon)
    if nargin < 1
        error('Brak funkcji'); % Wyświetl komunikat
    end
    if nargin < 2
        a = -1000000;
    end
    if nargin < 3
        b = 1000000;
    end
    if nargin < 4
        epsilon = 1e-6;
    end

    r = (sqrt(5)-1) / 2;
    
    delta1 = b - (b - a) * r;
    delta2 = a + (b - a) * r;
    
    iteration = 1;
    while (b - a >= epsilon)
        if iteration > 100
            minimum_x = (a + b) / 2;
            return
        end
        if f(delta1) < f(delta2)
            b = delta2;
            delta2 = delta1;
            delta1 = b - (b - a) * r;
        else
            a = delta1;
            delta1 = delta2;
            delta2 = a + (b - a) * r;
        end
        iteration = iteration + 1;
    end
    minimum_x = (a + b) / 2; % Minimum znajduje się w środku przedziału
end