G = [ 4/3   -1/3    0       0       1;      %[9]
      -1/3   4/3    -1      0       0;      %[0]
       0    -1      5/2     -1      -1/2;   %[v3]3
       0     0      -1      2       -1;     %[v4]4.5
       1     0      -1/2    -1      5/2];   %[v5]6
   
G = [   
    1.2357      -0.7786     0       -0.4571
    -0.7786     1.25        -0.4571 -0.0143
    0           -0.4571     0.8238  -0.3667
    -0.4571     -0.0143     -0.3667 0.838
];
% 
% f = [2 4];
% p = [1 3];
% Vp = [0 10]';
% Gff = G(f, f);
% Gfp = G(f, p);
% B = Gfp*Vp;
% Vf = Gff\B;
% val = -3.7076, -4.4392


   
f = [3 4 5];
p = [1 2];
Vp = [9 0]';
Gff = G(f, f);
Gfp = G(f, p);
B = Gfp*Vp;
Vf = Gff\B;
D = Gff.*eye(3);

[bicgSeq, ~, itSeq] = bicg(Gff, [0;0;0], [0;0;9], D, 1000, 0.001);

%EBE
%e = 1...n_el 
n_el= 6;
%E_e = [n_e x n_e]; Matriz de coeficiente do elemento e
n_e = 2; %numero de nós locais
n_nod = 5; %Número de nós globais

%L_e = matriz booleana [n_e x n_nod]: Mapeia um nó local para um nó global
L1 = [1 0 0 0 0; 0 1 0 0 0]; 
L2 = [0 1 0 0 0; 0 0 1 0 0];
L3 = [0 0 1 0 0; 0 0 0 1 0];
L4 = [0 0 0 1 0; 0 0 0 0 1];
L5 = [0 0 1 0 0; 0 0 0 0 1];
L6 = [0 0 0 0 1; 1 0 0 0 0];
L = [L1; L2; L3; L4; L5; L6];

N_e = n_e * n_el;
%[E_e] =  diagonal de blocos
%b = [N_e x 1]

b1 = [0; 0]; %[ne, nnod]
b2 = [0; 0];
b3 = [0; 0];
b4 = [0; 0];
b5 = [0; 0];
b6 = [0; 0];
b = [b1; b2; b3; b4; b5; b6];

M = [1 -1;-1 1];

E1 = M/3;
E2 = M;
E3 = M;
E4 = M;
E5 = M/2;
E6 = M;
E = blkdiag(E1, E2, E3, E4, E5, E6);

A = L'*E*L;
bb = L'*b;
