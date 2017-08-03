function [x, error, iter, flag, tempo] = cg_ebe_spmd(A, MAP, rhs, qtdCor, qtdN, max_it, tol, model)

    flag = 0;                                 % initialization
    iter = 0;
  
    x = zeros(qtdN,1);
    qtdCores = length(qtdCor);
    tempo = zeros(qtdCores,1);
    
    nW = 2;
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool('mesh',nW);
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
     
%     figure
    for iter = 1:max_it 
        

        z = zeros(qtdN,1);
        for c = 1:qtdCores            
            dof = MAP{1,c}(:);
            L = A(c,:);
            L = cell2mat(L');
            s = qtdCor(c);
            R = reshape(r(dof),3,[]);
            w = zeros(1,nW);
            w(1:nW-1) = idivide(int32(s),nW);
            w(nW) = w(1) + rem(s,nW);  
            
            
            startTime = tic;
            spmd
                codistL = codistributor('1d',1,3*w,[3*s 3]);                
                codistR = codistributor('1d',2,w,[3 s]); 
                cLd = codistributed(L,codistL);
                cRd = codistributed(R,codistR);
                Ld = getLocalPart(cLd);
                Rd = getLocalPart(cRd);
                s = size(Ld,1);
                Z = zeros(1,s);
                for i=1:3:s
                    j = 1+(i-1)/3;
                    k = [i i+1 i+2];
                    m = power(diag(Ld(k,:)),-1); 
                    Z(:,k) = m.*Rd(:,j);
                end                  
            end
            tempo(c) = toc(startTime);

            Z = [Z{:}];
            z(dof) = z(dof) + Z';  
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
            L = cell2mat(L');
            s = qtdCor(c);
            P = reshape(p(dof),3,[]);
            w(1:nW-1) = idivide(int32(s),nW);
            w(nW) = w(1) + rem(s,nW);  
            spmd
                codistL = codistributor('1d',1,3*w,[3*s 3]);                
                codistP = codistributor('1d',2,w,[3 s]); 
                cLd = codistributed(L,codistL);
                cPd = codistributed(P,codistP);
                Ld = getLocalPart(cLd);
                Pd = getLocalPart(cPd);
                s = size(Ld,1);
                Q = zeros(1,s);
                
                for i=1:3:s
                    j = 1+(i-1)/3;
                    k = [i i+1 i+2];
                    Q(:,k) = Ld(k,:)*Pd(:,j);
                end                  
            end  
            Q = [Q{:}];
            q(dof) = q(dof) + Q';  
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

