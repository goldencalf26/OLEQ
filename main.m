% Optimal Linear Estimator of Quaternion (OLEQ)
% Автор: Рогожин Игорь, студент: Самарский университет, Кафедра космических
% исследований
% Дата: 07.03.2024

clc; clear all;
format short

% Нормализированные векторы измерений (случайно сгенерированные) в орбитальной системе координат(ОСК) 
% Магнитометр | Акслерометр | GPS
Dr1 = [0.8190    0.5670    0.0875]'; Dr2 = [0.2449    0.4809    0.8419]'; Dr3 = [0.7004    0.1144    0.7045]';
Dr = [Dr1';Dr2';Dr3']';

% Углы ориентации(Эйлера)
psi = 5 * pi/180;
alph = 10 * pi/180;
phi = 5 * pi/180;
% Переход в связанную систему координат ССК
C_true = calculate_A_SSK(psi, alph, phi); % Истинная матрица перехода из ОСК в ССК

% Векторы измерений в ССК
Db1 = C_true*Dr1;
Db2 = C_true*Dr2;
Db3 = C_true*Dr3;
Db = [Db1';Db2';Db3']';


% Рассчитываем матрицы W для каждого вектора измерений
W_matrix1 = calculate_W_matrix(1, Db, Dr);
W_matrix2 = calculate_W_matrix(2, Db, Dr);
W_matrix3 = calculate_W_matrix(3, Db, Dr);


% СКО шума измерений для каждого вектора
sigma1 = 0.001; sigma2 = 0.001; sigma3 = 0.001;
totalStdDev = sqrt(1 / ((1 / sigma1^2) + (1 / sigma2^2) + (1 / sigma3^2)));
% Коэффициенты веса
a1 = totalStdDev^2 / sigma1^2;
a2 = totalStdDev^2 / sigma2^2;
a3 = totalStdDev^2 / sigma3^2;

% АЛГОРИТМ OLEQ
W = a1*W_matrix1 + a2*W_matrix2 + a3*W_matrix3;
q_prev = [1; 0; 0; 0];
q = [0; 0; 0; 1];
times = 0;
I = eye(4);
G = 0.5 * (W + I);

while (norm(q - q_prev) > 1e-8)
    q_prev = q;
    q = G * q_prev;
    q = q / norm(q);
    times = times + 1;
end

q = q / norm(q);

P1 = [q(1), q(2), -q(3), -q(4);
        -q(4), q(3), q(2), -q(1);
         q(3), q(4), q(1), q(2)];

P2 = [q(4), q(3), q(2), q(1);
         q(1), -q(2), q(3), -q(4);
        -q(2), -q(1), q(4), q(3)];

P3 = [-q(3), q(4), -q(1), q(2);
         q(2), q(1), q(4), q(3);
         q(1), -q(2), -q(3), q(4)];

% определений расчётой матрицы перехода в ССК
C = [P1*q P2*q P3*q];

% найдём ошибку определения матрицы поврота, в качестве оценки ошибки возьмём пространственный угол Эйлера
matrixError = C_true*C';% матрица ошибки
orientationError = 0.5*(trace(matrixError) - 1);
orientationError = rad2deg(acos(abs(orientationError))); % ошибка
disp(['Ошибка(пространственный угол Эйлера), градусы: ' num2str(orientationError)])
