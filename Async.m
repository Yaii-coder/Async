clear all; clc; close all;

Lm = 0.086; % индуктивность намагничивания двигателя
dls = 0.0016; % Индуктивность рассеяния статора 
dlr = 0.0025; % Индуктивность рассеяния ротора
Ls = Lm + dls; % индуктивность статора
Lr = Lm + dlr; % индуктивность ротора

Rs = 0.258; % сопротивление статора
Rr = 0.145; % Сопротивление ротора
ZP = 1;

Tr = Lr / Rr; % Число пар полюсов двигателя 
KR = Lm / Lr; % коэффицент ротора

% Usd = 0.9; % Напряжение на статоре по d-оси
% Usq = 0.9; % Напряжение на статоре по q-оси

J = 1; % Момент инерции
k_c = 0.1; % Коэффициент момента сопротивления

% T = 20; dt = 0.01; N = T / dt; 
T = 20; dt = 0.0001; N = T / dt; 

Isd(1) = 0; Isq(1) = 0; Mem(1) = 0; omega_e(1) = 0; omega_R(1) = 0;
PsiR(1) = 0.1; 



Um = 200; % Амплитуда напряжения
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
Usd = zeros(1, N+1); Usq = zeros(1, N+1);
Usd(1) = ((Ub(1)-Uc(1))/sqrt(3))*sin(deg2rad(fi))+Ua(1)*cos(deg2rad(fi));
Usq(1) = ((Ub(1)-Uc(1))/sqrt(3))*cos(deg2rad(fi))+Ua(1)*sin(deg2rad(fi));

for i = 1:N
    Isd(i+1) =  Isd(i) + dt * ((1 / (Ls - Lm * KR)) * (Usd(i) - Isd(i) * Rs - KR * (((1/Tr) * (Lm * Isd(i) - PsiR(i))))) + omega_e(i) * Isq(i));
    Isq(i+1) = Isq(i) + dt * (1/(Ls - Lm * KR) * (Usq(i) -  Isq(i) * Rs - KR * omega_e(i) * PsiR(i)) - omega_e(i) * Isd(i));
    PsiR(i+1) = PsiR(i) + dt * ((1/Tr) * (Lm * Isd(i) - PsiR(i)));
    omega_e(i+1) = ZP * omega_R(i) + (Lm / (PsiR(i) * Tr)) * Isq(i);
    Mem(i+1) = ((3 * ZP * KR) / 2) * PsiR(i) * Isq(i);

    omega_R(i+1) = omega_R(i) + dt * ((Mem(i) - k_c  * omega_R(i)) / J);

    fi(i+1) = fi(i) + fi_dt;
    Usd(i+1) = ((Ub(i)-Uc(i))/sqrt(3))*sin(deg2rad(fi(i+1)))+Ua(i)*cos(deg2rad(fi(i+1)));
    Usq(i+1) = ((Ub(i)-Uc(i))/sqrt(3))*cos(deg2rad(fi(i+1)))+Ua(i)*sin(deg2rad(fi(i+1)));

    t(i+1) = t(i) + dt;
    Ua(i+1) = Um * sin(omega * t(i+1) + phase1);
    Ub(i+1) = Um * sin(omega * t(i+1) + phase2);
    Uc(i+1) = Um * sin(omega * t(i+1) + phase3);
end

figure;
subplot(3,1,1);
plot(Isd);
title('Ток статора по d-оси');
xlabel('Время (с)');
ylabel('I_{sd} (A)');

subplot(3,1,2);
plot(Isq);
title('Ток статора по q-оси');
xlabel('Время (с)');
ylabel('I_{sq} (A)');

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