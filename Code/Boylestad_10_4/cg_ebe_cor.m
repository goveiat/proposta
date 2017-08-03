function [x, error, iter, flag, tempo] = cg_ebe_cor(A, b, cores, qtdN, max_it, tol, model)

    flag = 0;                                 % initialization
    iter = 0;
  
    x = zeros(qtdN,1);
    rhs = zeros(qtdN,1);
    qtdCores = length(cores);
    tempo = zeros(qtdCores,1);
    
    for c = 1:qtdCores
      t = cores{c};
      for e = t
            map = e([1 2 3]);
            rhs(map) = rhs(map) + b{e(5)};
      end       
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
%     fig = figure;
    for iter = 1:max_it                       % begin iteration
        
        z = zeros(qtdN,1);
        for c = 1:qtdCores
            t = cores{c};
            
            startTime = tic;
            for e = t
                map = e([1 2 3]);
                M = diag(A{e(5)});
                m = power(M,-1);   
                h = r(map);
                z(map) =  m.*h;
            end 
            tempo(c) = toc(startTime);
            
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
            t = cores{c};
            for e = t
                map = e([1 2 3]);
                q(map) = q(map) + A{e(5)}*p(map);
            end      
        end   

        alpha = rho / (p'*q );
        r = r - alpha*q; 
        
        x = x + alpha * p;


        error = norm( r ) / bnrm2;

        if ( error <= tol )
            break
        end 

        rho_1 = rho;
%         pdeplot(model,'XYData',x','ColorMap','jet');
%         drawnow
        % Capture the plot as an image 
%       frame = getframe(fig); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
%       % Write to the GIF File 
%       if iter == 1 
%           imwrite(imind,cm,'cg_ebe_cor2','gif', 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,'cg_ebe_cor2','gif','WriteMode','append'); 
%       end        
    end

    if ( error > tol ) 
        flag = 1; 
    end        
  
end

