function dydt = cc_oscillator2(t,y,p)

Vcc = p(1);
RL = p(2);
L = p(3);
C1 = p(4);
C2 = p(5);
Ree = p(6);
Vee = p(7);
beta = p(8);
Ron = p(9);
Vth = p(10);
Vb = p(11);
lambda = p(12);
n = p(13);

Ib = zeros(n,1);
Ic = zeros(n,1);

Vc = zeros(n,1);
Ve = zeros(n,1);
IL = zeros(n,1);

% [Ib,Ic] = oscillator_param(t,y,beta,Ron,Vth,Vb,n);

for j = 1:n %one direction
    
    if (Vb-y((3*j)-1)) <= Vth
        Ib(j) = 0;
    else
        Ib(j) = (Vb-y((3*j)-1)-Vth)/Ron;
    end
    Ic(j) = beta*Ib(j);
    
    k1 = mod((3*j),(3*n)) + 1;

    Vc(j) = ((y(3*j)-Ic(j))/C1) + ((((Vee+Vb-y((3*j)-1))/Ree)+y(3*j)+Ib(j))/C2) + (lambda/C1)*...
    ((y(k1)-y((3*j)-2))/Ree);
    Ve(j) = (((Vee+Vb-y((3*j)-1))/Ree)+y(3*j)+Ib(j))/C2;
    IL(j) = (Vcc-y((3*j)-2)+Vb-(y(3*j)*RL))/L;
end

% for j = 1:n %bidirectional
%     
%     if (Vb-y((3*j)-1)) <= Vth
%         Ib(j) = 0;
%     else
%         Ib(j) = (Vb-y((3*j)-1)-Vth)/Ron;
%     end
%     Ic(j) = beta*Ib(j);
%     
%     k1 = mod((3*j),(3*n)) + 1;
%     k2 = mod((3*(j+(n-2))),(3*n)) + 1;
% 
%     Vc(j) = ((y(3*j)-Ic(j))/C1) + ((((Vee+Vb-y((3*j)-1))/Ree)+y(3*j)+Ib(j))/C2) + (lambda/C1)*...
%     ((y(k1)-2*y((3*j)-2)+y(k2))/Ree);
%     Ve(j) = (((Vee+Vb-y((3*j)-1))/Ree)+y(3*j)+Ib(j))/C2;
%     IL(j) = (Vcc-y((3*j)-2)+Vb-(y(3*j)*RL))/L;
% end

dydt = [Vc(1); Ve(1); IL(1)];
for k = 2:n
    dydt = [dydt; Vc(k); Ve(k); IL(k)];
end


end