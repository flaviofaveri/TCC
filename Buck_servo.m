%% Projeto Controlador p/ conv. Buck com metodologia de Servosistemas
format long;
clear all;
close all;
clc;

%% Par�metros da planta original
delta1 = 4.2;       % Tens�o pico a pico medida no sobressinal
delta2 = 20;        % Tens�o pico a pico medida de 0 a regime
tp = 143.5e-6;      % Tempo de instante do pico

%% C�lculos para determina��o da FT da planta a ser controlada
Mp_p = delta1/delta2;                         % Sobressinal
zeta_p = -log(Mp_p)/sqrt(pi^2 + log(Mp_p)^2); % Fator de amortecimento
wn_p = pi/(tp*sqrt(1 - zeta_p^2));            % Frequ�ncia natual
Vi = 36.95;

% Lembrando que G(s) = Vi*wn^2 / 1 + 2*zeta*wn + wn^2
A_p = Vi*wn_p^2;                    % Numerador da FT
B_p = 2*zeta_p*wn_p;                % Fator central da FT
Gps = tf(A_p,[1 B_p A_p]);          % FT da planta original

%% Pr�-requisitos do projeto:
%  Ts5% igual a metade do valor obtido em malha aberta
%  Erro nulo em regime permanente para resposta ao degrau
%  Sobressinal menor que 10%
%  Estabilidade

Mp = 0.05;      % Novo sobressinal
ts5 = 100e-6;   % Novo tempo de acomoda��o 5%

% C�lculo de zeta e wn a partir do pr�-requisitos supra citados
zeta = log(Mp)/sqrt(pi^2 + log(Mp)^2);
wn = 3/(zeta*ts5);

% Polos Alocados
p1 = -zeta*wn +(wn*sqrt(1-zeta^2))*i;
p2 = -zeta*wn -(wn*sqrt(1-zeta^2))*i;
p3 = -abs(p1)*3;

% Matriz de polos para determina��o da matriz de ganho K de retroa��o
% de estados
pK = [p1 p2 p3];

%% Discretiza��o do integrador do Observador de Estados
% T = ts5/10;

%% Projeto da matriz de ganho K de retroa��o de estados
% Par�metros de entrada
Vi = 36.95;
Vo = 20;
L = 437e-6;
C = 4.1e-6;
R = 11.4;

D = Vo/Vi;
D_ = 1-D;

%% Matrizes da an�lise de circuitos
% Matriz dos acumuladores de energia
K = [L 0; 0 C];

% Matrizes da etapa de chave fechada
A1 = [0 -1; 1 -1/R];
B1 = [1; 0];
C1 = [0 1];     % Matriz C para a sa�da em tens�o
E1 = 0;

% Matrizes da etapa de chave aberta
A2 = [0 -1; 1 -1/R];
B2 = [0; 0];
C2 = [0 1];     % Matriz C para a sa�da em tens�o
E2 = 0;

% Matrizes m�dias
A = A1*D + A2*D_;
B = B1*D + B2*D_;
C = C1*D + C2*D_;
E = E1*D + E2*D_;

%% C�lculos para lineariza��o
U = Vi;
X = -inv(A)*B*U;
Y = (-C*inv(A)*B + E)*U;

Ap = inv(K)*A;
Bp = inv(K)*((A1-A2)*X + (B1-B2)*U);
Cp = C;
Ep = ((C1-C2)*X + (E1-E2)*U);

% [num,den] = ss2tf(Ap,Bp,Cp,Ep);
% Gps = tf(num,den);

sys = ss(Ap,Bp,Cp,Ep);
Gps = tf(sys);

% Matrizes auxilires para transforma��es
H=[0;0];
I=0;

% Matrizes A e B para c�lculo da matriz de ganho K
Ach = [Ap H; -Cp I];
Bch = [Bp; I];

% Matriz de ganho K de retroa��o de estados
Kch = place(Ach,Bch,pK);
Teste = eig(Ach-Bch*Kch);

%% Cria��o de sistema para teste no matlab
% AA = (Ach-Bch*Kch);
% BB = [0;0;1];
% CC = [C 0];
% DD = 0;
% 
% % Gera��o da resposta controlada ao degrau
% t = 0:0.000001:0.001;
% opt = stepDataOptions('InputOffset',0,'StepAmplitude',20);
% sys = ss(AA,BB,CC,DD);
% y = step(sys,t,opt);
% stepinfo(y,t,1,'SettlingTimeThreshold',0.05,'RiseTimeLimits',[0.1 0.9])
% plot(t,y);
% grid on;
% title('Uni-Step Response');
% xlabel('t [s]');
% ylabel('Output y');

%% Verifica��o da controlabilidade e observabilidade do sistema 
sys = ss(Ap,Bp,Cp,Ep);
Co = ctrb(sys);
Ob = obsv(sys);
k1 = rank(Co);
k2 = rank(Ob);
