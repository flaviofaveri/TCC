%% Projeto Controlador com Observador de Ordem Mínima p/ Conv. Buck
format long;
clear all;
close all;
clc;

%% Parâmetros da planta original
delta1 = 4.2;       % Tensão pico a pico medida no sobressinal
delta2 = 20;        % Tensão pico a pico medida de 0 a regime
tp = 143.5e-6;      % Tempo de instante do pico

%% Cálculos para determinação da FT da planta a ser controlada
Mp_p = delta1/delta2;                         % Sobressinal
zeta_p = -log(Mp_p)/sqrt(pi^2 + log(Mp_p)^2); % Fator de amortecimento
wn_p = pi/(tp*sqrt(1 - zeta_p^2));            % Frequência natual
Vi = 36.95;

% Lembrando que G(s) = Vi*wn^2 / 1 + 2*zeta*wn + wn^2
A_p = Vi*wn_p^2;                    % Numerador da FT
B_p = 2*zeta_p*wn_p;                % Fator central da FT
Gps = tf(A_p,[1 B_p A_p]);          % FT da planta original

%% Pré-requisitos do projeto:
%  Ts5% igual a metade do valor obtido em malha aberta
%  Erro nulo em regime permanente para resposta ao degrau
%  Sobressinal menor que 10%
%  Estabilidade

Mp = 0.05;      % Novo sobressinal
ts5 = 100e-6;   % Novo tempo de acomodação 5%

% Cálculo de zeta e wn a partir do pré-requisitos supra citados
zeta = log(Mp)/sqrt(pi^2 + log(Mp)^2);
wn = 3/(zeta*ts5);

% Polos Alocados
p1 = -zeta*wn +(wn*sqrt(1-zeta^2))*i;
p2 = -zeta*wn -(wn*sqrt(1-zeta^2))*i;
p3 = -abs(p1)*3;

% Matriz de polos para determinação da matriz de ganho K de retroação
% de estados
pK = [p1 p2 p3];

% Matriz de polos para determinação da matriz de ganho Ke do observador
% de estados de ordem mínima
% Escolhida a parte real de p1. Deve ser pelo menos 2x maior que o polo
% de K
pKe = [real(p1)*8];

%% Discretização do integrador do Observador de Estados
T = ts5/50;

%% Projeto da matriz de ganho K de retroação de estados
% Parâmetros de entrada
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

% [num,den] = ss2tf(Ap,Bp,Cp,Ep);
% Gps = tf(num,den);

sys = ss(Ap,Bp,Cp,Ep);
Gps = tf(sys);

% Matrizes auxilires para transformações
H=[0;0];
I=0;

% Matrizes A e B para cálculo da matriz de ganho K
Ach = [Ap H; -Cp I];
Bch = [Bp; I];

% Matriz de ganho K de retroação de estados
Kch = place(Ach,Bch,pK);
Teste = eig(Ach-Bch*Kch);

%% Projeto da matriz de ganho Ke do observador de ordem mínima
% Separação da matriz A
Aaa = Ap(2,2);
Aab = Ap(2,1);
Aba = Ap(1,2);
Abb = Ap(1,1);

% Separação da matriz B
Ba = Bp(2,1);
Bb = Bp(1,1);

% Matriz de ganho Ke do observador de ordem mínima
Ke = (place(Abb',Aab',pKe))';

% Verifica se os autovalores de A-KeC são os polos escolhidos
Teste2 = eig(Abb-Ke*Aab);

% Cálculo das matrizes dos blocos do observador de ordem mínima
Achp = Abb - Ke*Aab
Bchp = Achp*Ke + Aba - Ke*Aaa
Fchp = Bb - Ke*Ba
Dchp = [1;Ke]
Cchp = [0;eye(1)]
