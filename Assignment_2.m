%studentID: 4372182;
clear all; clc;
A = [4 7.5; 0 -2];
B = [0;1];
%% Question 3
K = [5.6 7];
eA =@(h) [exp(4*h) 7.5/6*(exp(4*h)-exp(-2*h)); 0 exp(-2*h)];

F_cl =@(h,hl) eA(hl)- (eA(hl) - eA(hl-h))/A*B*K;

h = 0.2;
delta = 6;
Q3_save_eig = zeros(1,delta);
for i = 1:delta
    Q3_save_eig(i) = max(abs(eig(F_cl(h,i*h))));
end
%max 3 consegutive package losses

figure();
plot([h,(delta)*h],[1,1],'r--','LineWidth',1.5);
 hold on;
 plot((1:delta)*h,Q3_save_eig,'LineWidth',2);
 hold off;
axis([h,delta*h,0,4]);
xlabel('h_l');
ylabel('max absolute eigenvalue');
title('eigenvalues of the cl system for various h_l');
legend('max abs eigenvalue','stable bound')

%% 3.2
p = sdpvar(1);
P0 = sdpvar(2);
P1 = sdpvar(2);
A0 = F_cl(h,h);
A1 = eA(h);
eps = 1e-8;
Constraints = [p<=1;p>=0;P0>=eps*eye(2);P1>=eps*eye(2)];
Constraints = [Constraints; P0 - A0'*(p*P1 + (1-p)*P0)*A0 >= eps*eye(2)];
Constraints = [Constraints; P1 - A1'*(p*P1 + (1-p)*P0)*A1 >= eps*eye(2)];
Objective = -p;
optimize(Constraints,Objective);
%upper bound found 0.2
p_star = value(p);
%% 3.3
x0 = 1000*rand(2,100); %100 random start vectors
x0_1 = x0;
x0_2 = x0;
x0_3 = x0;
N = 100;
Q3_save_eig_0 = zeros(1,N);
Q3_save_eig_1 = zeros(1,N);
Q3_save_eig_2 = zeros(1,N);
Q3_save_eig_3 = zeros(1,N);
for i = 1:N %N time steps
    if rand<= p_star
        x0 = eA(h) * x0;
    else 
        x0 = F_cl(h,h) * x0;
    end
    if rand<= 2*p_star
        x0_1 = eA(h) * x0_1;
    else 
        x0_1 = F_cl(h,h) * x0_1;
    end
    if rand<= 3* p_star
        x0_2 = eA(h) * x0_2;
    else 
        x0_2 = F_cl(h,h) * x0_2;
    end
    if rand<= 4*p_star
        x0_3 = eA(h) * x0_3;
    else 
        x0_3 = F_cl(h,h) * x0_3;
    end
    Q3_save_eig_0(i) = max(diag(x0' * x0)); %every ||x_k||^2
    Q3_save_eig_1(i) = max(diag(x0_1' * x0_1));
    Q3_save_eig_2(i) = max(diag(x0_2' * x0_2));
    Q3_save_eig_3(i) = max(diag(x0_3' * x0_3));
end

figure();
 plot(0:(N-1),Q3_save_eig_0,'LineWidth',2);
 hold on;
 plot(0:(N-1),Q3_save_eig_1,'LineWidth',2);
 plot(0:(N-1),Q3_save_eig_2,'LineWidth',2);
 plot(0:(N-1),Q3_save_eig_3,'LineWidth',2);
 hold off;
%axis([0,N,0,4]);
xlabel('time steps k');
set(gca, 'YScale', 'log')
ylabel('norm');
title('Max norm of 100 random x_0 start vectors over time');
legend('p^*','2p^*','3p^*','4p^*')