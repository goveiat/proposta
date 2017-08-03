function [x, error, iter, flag, tempo] = cg_ebe_par(A, MAP, rhs, qtdCor, qtdN, max_it, tol, model)

    flag = 0;                                 % initialization
    iter = 0;
  
    x = zeros(qtdN,1);
    qtdCores = length(qtdCor);
    tempo = zeros(qtdCores,1);
    
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool('mesh',2);
    end
      

    bnrm2 = norm( rhs );
    if  ( bnrm2 == 0.0 )
        bnrm2 = 1.0; 
    end

    r = rhs;

    error = norm( r ) / bnrm2;
    if ( error < tol ) 
        return
    end

    %figure
    for iter = 1:max_it                       % begin iteration
        

        z = zeros(qtdN,1);
        for c = 1:qtdCores          
            dof = MAP{1,c}(:);
            L = A(c,:);
            s = qtdCor(c);
            R = reshape(r(dof),3,[]);
            Z = zeros(3,s);
            
            startTime = tic;
            parfor j = 1:s
                m = power(diag(L{j}),-1); 
                h = R(:,j);
                Z(:,j) = m.*h;
            end
            toc(startTime)

            z(dof) = z(dof) + reshape(Z,[],1);  
        end


        rho = (r'*z);

        if ( iter > 1 )                       % direction vector
            beta = rho / rho_1;
            p = z + beta*p;
        else
            p = z;
        end
        
        q = zeros(qtdN,1);
        for c = 1:qtdCores            
            dof = MAP{1,c}(:);
            L = A(c,:);
            s = qtdCor(c);
            P = reshape(p(dof),3,[]);
            Q = zeros(3,s);
            parfor j = 1:s
                Q(:,j) = L{j}*P(:,j);
            end
            q(dof) = q(dof) + reshape(Q,[],1);  
        end         

        alpha = rho / (p'*q );
        x = x + alpha * p;
        r = r - alpha*q; 


        error = norm( r ) / bnrm2;

        if ( error <= tol )
            break
        end 

        rho_1 = rho;

         

%         pdeplot(model,'XYData',x','ColorMap','jet');
%         drawnow
    end

    if ( error > tol ) 
        flag = 1; 
    end        
  
end

