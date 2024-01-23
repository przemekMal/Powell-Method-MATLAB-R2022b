%kod wykorszytany do funkcji Himmelbalu's
format short
clc
clear all

%Himmelblau's function     f(x,y)=(x^2+y-11)^2+(x+y^2-7)^2
%f(0,0) < epsilon
f = @(x) (x(1)^2+x(2)-11)^2+(x(1)+x(2)^2-7)^2;
x0 = [0, 0];

g = @(x) x(1)-x(2)-2;
G= {g};
%Nalezy usunac G jako argument poniżej, jeśli funkcja bez ograniczeń            ⇓⇓⇓
[minimum, xes, iter] = algorytm_powella_z_funkcja_kary(f, x0, [], 100, 10e-3, 6, G);
xes
minimum
f(minimum)

% Wykres funkcji trójwymiarowej
tytul = "Wykres siatki punktów w trójwymiarowej przestrzeni po liczbie: " + num2str(iter) + " iteracji.";

x = -6:0.2:6;
y =  -6:0.2:6;
[X, Y] = meshgrid(x, y);
Z = (X.^2 + Y - 11).^2 + (X + Y.^2 - 7).^2;
G = Y - X +2;

figure;
subplot(1, 2, 1);
contour(X, Y, Z, 50); % Dodanie wykresu konturowego
hold on;
scatter(minimum(:,1), minimum(:,2), 150, 'filled', 'MarkerFaceColor', 'g');

%Rsyowanie funkcji g(x)
contour(X, Y, G, [0 0], 'LineColor', 'blue', 'LineWidth', 1.5); % Dodanie wykresu ograniczenia g(x) = 0

plot(xes(:,1), xes(:,2), 'r.-', 'LineWidth', 1.5); % Dodanie ścieżki optymalizacji
xlabel('Oś X');
ylabel('Oś Y');
title('Wykres konturowy funkcji');
%Legenda z ogrniczeniem
legend('Funkcja', 'Minimum',  'Ograniczenie', 'Ścieżka optymalizacji');


subplot(1, 2, 2);
surf(X, Y, Z);
[odX, doX, odY, doY, odZ, doZ] = deal(-6, 6, -6, 6, -1, 2000);
axis([odX, doX, odY, doY, odZ, doZ]);
hold on;
for i = 1:size(xes, 1)
    scatter3(xes(i, 1), xes(i, 2), f([xes(i,1), xes(i,2)]), 'filled', 'MarkerFaceColor', 'r');
end

for i = 2:size(xes, 1)
    if (((xes(i-1, 1) < doX) && (xes(i-1, 1) > odX)) && ((xes(i-1, 2) < doY) && (xes(i, 2) > odY))) && (((xes(i, 1) < doX) && (xes(i, 1) > odX)) && ((xes(i, 2) < doY) && (xes(i, 2) > odY)))
        x_line = [xes(i-1, 1), xes(i, 1) ]; % Przykładowe punkty do narysowania linii
        y_line = [xes(i-1, 2), xes(i, 2) ];
        z_line = [f([xes(i-1,1), xes(i-1,2)]), f([xes(i,1), xes(i,2)])];
        line(x_line, y_line, z_line, 'Color', 'red', 'LineStyle', '-.'); % Rysuje czerwoną linię kropkowaną między punktami
    end
end
scatter3(minimum(:,1),minimum(:,2), f([minimum(:,1), minimum(:,2)]), 150, 'filled', 'MarkerFaceColor', 'g');
colorbar; % Dodanie paska kolorów dla wartości F_grid
xlabel('Oś X');
ylabel('Oś Y');
zlabel('Oś Z');
title(tytul);