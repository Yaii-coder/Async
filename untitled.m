clear all; clc; close all;

Lm = 0.086; % индуктивность намагничивания двигателя
Ls = 0.6223; % индуктивность статора
Lr = 3.118; % индуктивность ротора

Rs = 0.318; % сопротивление статора
Rr = 0.141; % Сопротивление ротора
ZP = 1; % Число пар полюсов двигателя 

Tr = Lr / Rr; 
KR = Lm / Lr; % коэффицент ротора

J = 1; % Момент инерции
k_c = 0.1; % Коэффициент момента сопротивления

T = 200; dt = 0.01; N = T / dt; 

Isd(1) = 0; Isq(1) = 0; Mem(1) = 0; omega_e(1) = 0; omega_R(1) = 0;
PsiR(1) = 0.1; 
I = 11.2;

Um = 380; % Амплитуда напряжения
f = 50; % Частота (Гц)
omega = 2 * pi * f; % Угловая частота (рад/с)
phase1 = -2*pi/3;  phase2 = 0;  phase3 = 2*pi/3; % Начальная фаза (рад)

% Инициализация массивов
Ua = zeros(1, N+1); Ub = zeros(1, N+1); Uc = zeros(1, N+1);
t = zeros(1, N+1);

% Начальные условия
t(1) = 0;

fi(1)  = 0; fi_dt = 360/N;
Usd = zeros(1, N+1); Usq = zeros(1, N+1);
% Usd(1) = ((Ub(1)-Uc(1))/sqrt(3))*sin(deg2rad(fi))+Ua(1)*cos(deg2rad(fi));
% Usq(1) = ((Ub(1)-Uc(1))/sqrt(3))*cos(deg2rad(fi))+Ua(1)*sin(deg2rad(fi));

for i = 1:N
    PsiR(i+1) = PsiR(i) + dt * ((1/Tr) * (Lm * I - PsiR(i)));
    omega_e(i+1) = ZP * omega_R(i) + (Lm / (PsiR(i) * Tr)) * I;
    Mem(i+1) = ((3 * ZP * KR) / 2) * PsiR(i) * I;

    omega_R(i+1) = omega_R(i) + dt * ((Mem(i) - k_c  * omega_R(i)) / J);

    fi(i+1) = fi(i) + fi_dt;
    % Usd(i+1) = ((Ub(i)-Uc(i))/sqrt(3))*sin(deg2rad(fi(i+1)))+Ua(i)*cos(deg2rad(fi(i+1)));
    % Usq(i+1) = ((Ub(i)-Uc(i))/sqrt(3))*cos(deg2rad(fi(i+1)))+Ua(i)*sin(deg2rad(fi(i+1)));

    t(i+1) = t(i) + dt;
    Ua(i+1) = Um * sin(omega * t(i+1) + phase1);
    Ub(i+1) = Um * sin(omega * t(i+1) + phase2);
    Uc(i+1) = Um * sin(omega * t(i+1) + phase3);
end

subplot(3,1,3);
plot(PsiR);
title('Поток ротора');
xlabel('Время (с)');
ylabel('ΨR (Wb)');

figure;
subplot(3,1,1);
plot(Mem);
title('Электромагнитный момент');
xlabel('Время (с)');
ylabel('M (Nm)');

subplot(3,1,2);
plot(omega_e);
title('omega e');
xlabel('Время (с)');
ylabel(' ');

subplot(3,1,3);
plot(omega_R);
title('omega R');
xlabel('Время (с)');
ylabel(' ');