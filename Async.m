clear all;
clc;
close all;

% Параметры модели
dLs = 0.0016; % индуктивность рессеивания статора
dLr = 0.0025; % индуктивность рессеивания ротора
Lm = 0.086; % индуктивность намагничивания двигателя
Ls = Lm + dLs; % индуктивность на статоре
Rs = 0.258; % сопротивление статора
Rr = 0.145; % Сопротивление ротора
ZP = 1;
Lr = Lm + dLr;
KR = Lm/Lr; % коэффицент ротора
Tr = Lr/Rr;
J=1; % Момент инерции
kc = 0.1; % Коэффициент момента сопротивления
% Входные сигналы (примерные значения, заменить на реальные)
Usd = 0.9; % Напряжение на статоре по d-оси
Usq = 0.9; % Напряжение на статоре по q-оси
% Временные параметры
T = 200; % Время моделирования
dt = 0.0001; % Шаг интегрирования
N = T / dt; % Количество шагов

Isd(1) = 0; % Ток статора по d-оси
Isq(1) = 0; % Ток статора по q-оси
PsiR(1) = 0.1; % Поток ротора по q-оси
omega_e(1) = 0; % Электрическая угловая скорость
Mem(1) = 0;
omega_R(1) = 0; % угловая скорость ротора




Um = 220; % Амплитуда напряжения
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


    omega_R(i+1) = omega_R(i)+dt*((Mem(i) - kc*omega_R(i))/J);   
    Isd(i+1) =  Isd(i) + dt*(1/(Ls-Lm*KR) * (Ud(i+1) - Isd(i)*Rs  - KR*(((1/Tr) * (Lm*Isd(i) - PsiR(i))))) + omega_e(i)*Isq(i));
    Isq(i+1) = Isq(i) + dt*(1/(Ls-Lm*KR) * (Uq(i+1) -  Isq(i)*Rs - KR*omega_e(i)*PsiR(i))- omega_e(i)*Isd(i));
    PsiR(i+1) = PsiR(i) + dt*((1/Tr) * (Lm*Isd(i) - PsiR(i)));
    omega_e(i+1) = ZP*omega_R(i) + (Lm/(PsiR(i)*Tr))*Isq(i);
    Mem(i+1) = ((3*ZP*KR)/2)*PsiR(i)*Isq(i);
end

% Построение графиков
% figure;
% subplot(3,1,1);
% plot(Isd);
% title('Ток статора по d-оси');
% xlabel('Время (с)');
% ylabel('I_{sd} (A)');
% 
% subplot(3,1,2);
% plot(Isq);
% title('Ток статора по q-оси');
% xlabel('Время (с)');
% ylabel('I_{sq} (A)');
% 
% subplot(3,1,3);
% plot(PsiR);
% title('Поток ротора по q-оси');
% xlabel('Время (с)');
% ylabel('Ψ_{sq} (Wb)');

figure;
plot(omega_e);
title('omega_e');
grid on;

figure;
plot(omega_R);
title('omega_R');
grid on;

