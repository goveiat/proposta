function [x, error, iter, flag] = ebe_bicg(A, x, b, M, max_it, tol)


for u = 1:length(U)
    r(u) = 0;
    D(u) = 0;
end

for c = 1:length(colors)
    for e = 1:length(els(c))
        r(e) = r(e) + b(e) -A(e)*v(e);
        D(e) = D(e) + A(e)*eye(2);
    end
end

for u = 1:length(U)
    D(u) = inv(D(u));
    r_til(u) = r(u);
    d(u) = dot(D(u),r(u));
    d_til(u) = dot(D(u), r_til(u));
    q(u) = 0;
    q_til = 0;
end


del = 0;
for u = 1:length(U)
    del = del + r(u)*d(u);
end

while del > prec
    for c = 1:length(colors)
        for e = 1:length(els(c))
            q(e) = q(e) + A(e)*d(e);
            q_til(e) = q_til(e) + A(e)'*d_til(e);
        end
    end
    
    alpha = 0;
    for u = 1:length(U)
        alpha = alpha + dot(d_til(u),q(u));
    end 
    
    for u = 1:length(U)
        x(u) = x(u) + dot((del/alpha), d(u));
        r(u) = r(u) - dot((del/alpha), q(u));
        r_til(u) = r-til(u) - dot((del/alpha), q_til(u));
    end   
    
    del_til = del;
    del = 0;
    for u = 1:length(U)
        del = del + r(u)*D(u)*r_til(u);
    end
    
    alpha = del/del_til;
    
    for u = 1:length(U)
        d(u) = dot(D(u), r(u)) + dot(alpha, d(u));
        d_til(u) = dot(D(u), r_til(u)) + dot(alpha, d_til(u));
        q(u) = 0;
        q_til = 0;
    end     
end


%  -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
%  [x, error, iter, flag] = bicg(A, x, b, M, max_it, tol)
%
% bicg.m solves the linear system Ax=b using the 
% BiConjugate Gradient Method with preconditioning.
%
% input   A        REAL matrix
%         M        REAL preconditioner matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         max_it   INTEGER maximum number of iterations
%         tol      REAL error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm
%         iter     INTEGER number of iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it
%                          -1 = breakdown

  iter = 0;                              % initialization
  flag = 0;

  bnrm2 = norm( b );
  if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

  r = b - A*x;
  error = norm( r ) / bnrm2;
  if ( error < tol ) return, end

  r_tld  = r;

  for iter = 1:max_it                    % begin iteration

     z = M \ r;
     z_tld = M' \ r_tld;
     rho   = ( z'*r_tld );
     if ( rho == 0.0 ),
        break
     end

     if ( iter > 1 ),                    % direction vectors
        beta = rho / rho_1;
        p   = z  + beta*p;
        p_tld = z_tld + beta*p_tld;
     else
        p = z;
        p_tld = z_tld;
     end

     q = A*p;                            % compute residual pair
     q_tld = A'*p_tld;
     alpha = rho / (p_tld'*q );

     x = x + alpha*p;                    % update approximation
     r = r - alpha*q;
     r_tld = r_tld - alpha*q_tld;

     error = norm( r ) / bnrm2;          % check convergence
     if ( error <= tol ), break, end

     rho_1 = rho;

  end

  if ( error <= tol ),                   % converged
     flag =  0;
  elseif ( rho == 0.0 ),                 % breakdown
     flag = -1;
  else
     flag = 1;                           % no convergence
  end

% END bicg.m
