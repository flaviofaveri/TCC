%% Modelagem conv. Buck em Espaços de Estado
%format long;
clear all;
close all;
clc;

%% Parâmetros de entrada
Vi = 36.95;
Vo = 20;
L = 437e-6;
C = 4.1e-6;
R = 11.4;

D = Vo/Vi;
D_ = 1-D;

%% Matrizes da análise de circuitos
% Matriz dos acumuladores de energia
K = [L 0; 0 C];

% Matrizes da etapa de chave fechada
A1 = [0 -1; 1 -1/R];
B1 = [1; 0];
C1 = [0 1];     % Matriz C para a saída em tensão
E1 = 0;

% Matrizes da etapa de chave aberta
A2 = [0 -1; 1 -1/R];
B2 = [0; 0];
C2 = [0 1];     % Matriz C para a saída em tensão
E2 = 0;

% Matrizes médias
A = A1*D + A2*D_;
B = B1*D + B2*D_;
C = C1*D + C2*D_;
E = E1*D + E2*D_;

%% Cálculos para linearização
U = Vi;
X = -inv(A)*B*U;
Y = (-C*inv(A)*B + E)*U;

Ap = inv(K)*A;
Bp = inv(K)*((A1-A2)*X + (B1-B2)*U);
Cp = C;
Ep = ((C1-C2)*X + (E1-E2)*U);

[num,den] = ss2tf(Ap,Bp,Cp,Ep);
Gps = tf(num,den);
% % step(Gps)

Tf = 0.001;
t = 0:Tf/999:Tf;
opt = stepDataOptions('InputOffset',0.49,'StepAmplitude',0.1);
sys = ss(Ap,Bp,Cp,Ep);
y = step(Gps,t,opt);
stepinfo(y,t,1,'SettlingTimeThreshold',0.05,'RiseTimeLimits',[0.1 0.9])
plot(t,y);
grid on;
title('Resposta à perturbação em torno de D');
xlabel('t [s]');
ylabel('Tensão de saída [V]');