clc
clear

G = [ 4/3   -1/3    0       0       1;      
      -1/3   4/3    -1      0       0;      
       0    -1      5/2     -1      -1/2;   
       0     0      -1      2       -1;     
       1     0      -1/2    -1      5/2];

cont = [1 2; 9 0];
   
numEl = 6;

adj = {
    [2 6; 2 1; 1 2] 
    [1 3 4; 1 2 2; 2 1 1]
    [2 4 5 6; 1 1 2 2; 2 1 2 1]
    [2 3 5; 1 1 2; 2 1 1]
    [3 4 6; 2 1 2 ; 2 2 1]
    [1 3 5; 2 1 1; 1 2 2];
};

Q1 = [
     1 0 0 0 0;
     0 1 0 0 0;
     0 0 1 0 0;
     0 0 1 0 0;
     0 0 0 0 1;
     0 0 0 1 0;
];

Q2 = [
     0 1 0 0 0;
     0 0 1 0 0;
     0 0 0 1 0;
     0 0 0 0 1;
     0 0 0 1 0;
     1 0 0 0 0;
];


zeros_col = {
    [0; 0]
    [0; 0]
    [0; 0]
    [0; 0]
    [0; 0]
    [0; 0]
};

M = [1 -1;-1 1];

A = {
    [0.25 0; 0 0.25]
    [0.75 0; 0 1]
    M/2
    M
    M
    [1 0; 0 0.75]
};


m = {
    cell(2)
    cell(2)
    cell(2)
    cell(2)
    cell(2)
    cell(2)    
};

h = {
    cell(2, 1)
    cell(2, 1)
    cell(2, 1)
    cell(2, 1)
    cell(2, 1)
    cell(2, 1)    
};

b = {
    [2.25; 0]
    [0; 0]
    [0; 0]
    [0; 0]
    [0; 0]
    [9; 6.75]
}

% for c = cont   
%     e = c(1);
%     u = c(2);
%     el1 = find(Q1(:,e) == 1);
%     el2 = find(Q2(:,e) == 1);
%     els = [el1' el2'];
%     k = 0;
%     temp = els;
%     for i = els
%         for j = adj{i}
%             if any(temp == j(1))
%                 temp(temp == j(1)) = [];
%                 k = k+1;
%                 J{k} = j;                
%             end
%         end        
%     end
%     
%     den = 0;
%     for i = 1:k
%         e = J{i}(1);
%         no = J{i}(3);
%         den = den + A{e}(no,no);
%     end
%     
%     b{e} = b{e} - A{e}(:,no)*u;
%     for i = 1:k
%         e = J{i}(1);
%         no = J{i}(3);
%         v = A{e}(no,no)/den;        
%         b{e}(no) = v*u;        
%         A{e}(no,:) = 0;
%         A{e}(:,no) = 0;
%         A{e}(no,no) = v;
%     end    
% end

x = zeros_col;
r = b;

for e = 1:numEl   
    m{e} = tril(A{e});   
end



for e = 1:numEl   
    sum = zeros(2);
    for ad = adj{e}(1,:)
        sum = sum + tril(A{ad});
    end    
    m_{e} = m{e} + sum;   
end

for e = 1:numEl         
    h{e} = (m_{e} \ r{e});    
end

for e = 1:numEl       
    sum = zeros(2, 1);
    for ad = adj{e}(1,:)
        sum = sum + h{ad};
    end    
    h_{e} = h{e} + sum;   
    p{e} = h_{e};
end


lambda_0 = 0;
for e = 1:numEl
    lambda_0 = lambda_0 + r{e}'*h_{e} ; 
end 

k=1;
while k <100
 
    %2
    r_h = 0;
    for e = 1:numEl
        v = [0; 0];
        for ad = adj{e}(1,:)
            v = v + h{ad};
        end    
        r_h = r_h + r{e}'*(h{e} + v);
    end
    
    p_Ap = 0;    
    for e = 1:numEl
        p_Ap = p_Ap + p{e}'*A{e}*p{e};
    end  
    
    alpha = r_h/p_Ap;
    
    %3
    for e = 1:numEl
        x{e} = x{e} + alpha*p{e};
        r{e} = r{e} - alpha*A{e}*p{e};
    end      

    
    %4
    for e = 1:numEl         
        h{e} = (m_{e} \ r{e});    
    end    
    
    %5
    sum = zeros(2, 1);
    for e = 1:numEl       
        sum = zeros(2, 1);
        for ad = adj{e}(1,:)
            sum = sum + h{ad};
        end    
        h_{e} = h{e} + sum;    
    end    
    
    lambda_1 = 0;
    for e = 1:numEl
        lambda_1 = lambda_1 + r{e}'*h_{e} ; 
    end   
    
    beta = lambda_1/lambda_0;
    
    if 0
    else
        for e = 1:numEl
            p{e} = h_{e}+beta*p{e};
        end               
    end
        
%     lambda_0 = lambda_1;
    
    k = k+1;
end
