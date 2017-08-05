clc
clear

cont = [1 3; 0 10];

n = [
    0.8    1.4  2.1     1.2;
    1.8    1.4  2.1     2.7;
];

t = [
    1   2;
    2   3;
    4   4;
];

qtdT = size(t, 2);
qtdN = size(n, 2);

Q = {
    [
        1 0 0 0;
        0 1 0 0;
        0 0 0 1;
    ]
    [
        0 1 0 0;
        0 0 1 0;
        0 0 0 1;
    ]    
};

%obtenção das matrizes elementares
[A, adj] = matEl(n, t, qtdN, qtdT);


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

%montagem do rhs global
B = zeros(qtdN,1);
for e = 1:qtdT
    ng = t(:,e)';
    B(ng,1) = B(ng,1) +b{e};
end

adjEl = cell(qtdT,1);
for no = 1:qtdN
    es = adj{no}(1,:);
    for e = es
        for ad = es
            if ad ~= e
                if isempty(find(adjEl{e} == ad,1))
                    adjEl{e}(1,end+1) = ad;
                end
            end            
        end
    end
end


r = B;
x = zeros(qtdN,1);

%Obtenção dos precondicionadores
for e = 1:qtdT   
    m{e} = tril(A{e});   
end

%Acumulação dos adjacentes
for e = 1:qtdT   
    sum = zeros(3);
    for ad = adjEl{e}
        sum = sum + m{ad};
    end    
    m_{e} = m{e} + sum;   
end

%sol. do sistema de precond
for e = 1:qtdT         
    h{e} = (m_{e} \ r{e});    
end

for e = 1:qtdT       
    sum = zeros(3, 1);
    for ad = adjEl{e}
        sum = sum + h{ad};
    end    
    h_{e} = h{e} + sum;   
    p{e} = h_{e};
end


lambda_0 = 0;
for e = 1:qtdT
    lambda_0 = lambda_0 + r{e}'*h_{e} ; 
end 


k=1;
while k <100
    
    %2
    r_r = 0;
    for e = 1:qtdT
        sum = zeros(3,1);
        for ad = adjEl{e}
            sum = sum + r{ad};
        end
        v = r{e} + sum;
        r_r = r_r + r{e}'*v;
    end
    
    p_Ap = 0;    
    for e = 1:qtdT
        p_Ap = p_Ap + p{e}'*A{e}*p{e};
    end  
    
    alpha = r_r/p_Ap;
    
    %3
    for e = 1:qtdT
        x{e} = x{e} + alpha*p{e};
        r{e} = r{e} - alpha*A{e}*p{e};
    end  
    
    for e = 1:qtdT         
        h{e} = (m_{e} \ r{e});    
    end  
    
    for e = 1:qtdT       
        sum = zeros(3, 1);
        for ad = adjEl{e}
            sum = sum + h{ad};
        end    
        h_{e} = h{e} + sum;   
    end    
  
    
    lambda_1 = 0;
    for e = 1:qtdT
        lambda_1 = lambda_1 + r{e}'*h{e} ; 
    end   
    
    beta = lambda_1/lambda_0;
    
    if k>1 && beta < 1
        x
    else
        for e = 1:qtdT
            p{e} = h{e}+beta*p{e};
        end 
    end
        
    lambda_0 = lambda_1;
    
    k = k+1;
end
