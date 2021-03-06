clear;
clc;

max_it = 1000;
tol = 0.001;

cont = [1 3; 10 0];

n = [
    0.5    2.3  5     3;
    2.6    0.8  1.3     3.2;
];

t = [
    1   2;
    2   3;
    4   4;
];

qtdT = size(t, 2);
qtdN = size(n, 2);

[A, adj] = matEl(n, t, qtdN, qtdT);

G = matG(A, t, qtdN, qtdT);

%obtém as ondições de x
[x, pts] = defContorno(n, qtdN);
r_v = [10 -1 0 -1]
f = find(r_v == -1);
p = find(r_v ~= -1);
Vp = r_v(p)';
Gff = G(f, f);
Gfp = G(f, p);
b = -Gfp*Vp;
Vf = Gff\b;
r_v(r_v == -1) = Vf;
m_v = reshape(r_v, r_qtdNos)';

%Direto
sol = G\b;

%inicialização do rhs
b = cell(qtdT,1);
for e = 1:qtdT
    b{e} = zeros(3,1);
end

%atribuição das condições de contorno
for c = cont
    no = c(1);
    val = c(2);
    els = adj{no};
    qtd = size(els,2);
    for i = els
        e = i(1);
        ni = i(2);
        b{e} = b{e} - A{e}(:,ni)*val;
        b{e}(ni) = val/qtd;          
        A{e}(:,ni) = 0;
        A{e}(ni,:) = 0;
        A{e}(ni, ni) = 1/qtd;
    end
end

  flag = 0;                                 % initialization
  iter = 0;
  
  x = zeros(qtdN,1);
  
  rhs = zeros(qtdN,1);
  for e = 1:qtdT
      map = t(:,e)';
      rhs(map) = rhs(map) + b{e};
  end

  bnrm2 = norm( rhs );
  if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

  r = rhs;

  error = norm( r ) / bnrm2;
  if ( error < tol ) return, end

  for iter = 1:max_it                       % begin iteration

      z = zeros(qtdN,1);
      for e = 1:qtdT
          map = t(:,e)';
          M = diag(A{e});
          m = power(M,-1);          
          z(map) = z(map) + m.*r(map);
      end     
     rho = (r'*z);

     if ( iter > 1 ),                       % direction vector
        beta = rho / rho_1;
        p = z + beta*p;
     else
        p = z;
     end

      q = zeros(qtdN,1);
      for e = 1:qtdT
          map = t(:,e)';
          q(map) = q(map) + A{e}*p(map);
      end
     alpha = rho / (p'*q );
     x = x + alpha * p;                    % update approximation vector

     r = r - alpha*q;                      % compute residual
     error = norm( r ) / bnrm2;            % check convergence
     if ( error <= tol ), break, end 

     rho_1 = rho;

  end

  if ( error > tol ) flag = 1; end         % no convergence