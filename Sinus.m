clc; clear all; close all;

T = 0.02; dt = 0.0001; N = T / dt; 

Um = 537.4823; % Амплитуда напряжения
f = 50; % Частота (Гц)
omega = 2 * pi * f; % Угловая частота (рад/с)
phase1 = -2*pi/3;  phase2 = 0;  phase3 = 2*pi/3; % Начальная фаза (рад)

% Инициализация массивов
Ua = zeros(1, N+1); Ub = zeros(1, N+1); Uc = zeros(1, N+1);
t = zeros(1, N+1);

% Начальные условия
Ua(1) = Um * sin(phase1); Ub(1) = Um * sin(phase2); Uc(1) = Um * sin(phase3);
t(1) = 0;

fi(1)  = 0; fi_dt = 360/N;
Ud = zeros(1, N+1); Uq = zeros(1, N+1);
Ud(1) = ((Ub(1)-Uc(1))/sqrt(3))*sin(deg2rad(fi))+Ua(1)*cos(deg2rad(fi));
Uq(1) = ((Ub(1)-Uc(1))/sqrt(3))*cos(deg2rad(fi))+Ua(1)*sin(deg2rad(fi));

for i = 1:N
    fi(i+1) = fi(i) + fi_dt;
    Ud(i+1) = ((Ub(i)-Uc(i))/sqrt(3))*sin(deg2rad(fi(i+1)))+Ua(i)*cos(deg2rad(fi(i+1)));
    Uq(i+1) = ((Ub(i)-Uc(i))/sqrt(3))*cos(deg2rad(fi(i+1)))+Ua(i)*sin(deg2rad(fi(i+1)));

    t(i+1) = t(i) + dt;
    Ua(i+1) = Um * sin(omega * t(i+1) + phase1);
    Ub(i+1) = Um * sin(omega * t(i+1) + phase2);
    Uc(i+1) = Um * sin(omega * t(i+1) + phase3);
end

figure
plot(t, Ua, t, Ub, t, Uc);
xlabel('Время (с)'); ylabel('Напряжение (В)');
grid on;

figure
plot(t, Ud, t, Uq);
xlabel('Время (с)'); ylabel('Напряжение (В)');
legend('Ud', 'Uq');
grid on;