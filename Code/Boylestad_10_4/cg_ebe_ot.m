function [x, error, iter, flag, rhs] = cg_ebe_ot(A, b, t, qtdN, qtdT, max_it, tol, mapa)

    flag = 0;                                 % initialization
    iter = 0;
  
  x = zeros(qtdN,1);
  rhs = zeros(qtdN,1);
  
  nW = 2;
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool('mesh',nW);
    end

  for e = 1:qtdT
      map = t([1 2 3],e)';
      rhs(map) = rhs(map) + b{e};
  end
    
  bnrm2 = norm( rhs );
  if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

  r = rhs;

  error = norm( r ) / bnrm2;
  if ( error < tol ) return, end

  for iter = 1:max_it                       % begin iteration

      z = zeros(qtdN,1);
      w = zeros(1,nW);
      w(1:nW-1) = idivide(int32(qtdT),nW);
      w(nW) = w(1) + rem(qtdT,nW);
      
 
      
      for e = 1:qtdT
          L(e,1:3,1:3) = A{e};
      end  
      %tic;
      for e = 1:qtdT
          R(:,e) = r(mapa(:,e));
      end      
      parfor e = 1:qtdT
          M = diag(A{e});
          m = power(M,-1);
          kk(:,e) = m.*R(:,e);
      end  
      for e = 1:qtdT
          z(mapa(:,e)) = z(mapa(:,e)) + kk(:,e);
      end      
      %toc;

%       spmd
%         codistL = codistributor('1d',1,3*w,[3*qtdT 3]);  
%         codistMAP = codistributor('1d',2,w,[3*qtdT 3]); 
%         cLd = codistributed(L,codistL);
%         cRd = codistributed(R,codistR);          
%           map = t([1 2 3],e)';
%           M = diag(A{e});
%           m = power(M,-1);    
%           h = r(map);
%           z(map) = z(map) + m.*h;
%       end
      
     rho = (r'*z);

     if ( iter > 1 ),                       % direction vector
        beta = rho / rho_1;
        p = z + beta*p;
     else
        p = z;
     end

      q = zeros(qtdN,1);
      for e = 1:qtdT
          map = t([1 2 3],e)';
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
  
end

