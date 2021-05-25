%studentID: 4372182;
clear all; clc;
A = [4 7.5; 0 -2];
B = [0;1];
%% Question 1
%contious time, constant sampling time h, no delays tau=0
%K_bar places poles on {-2,-3}
Q1_K_bar = [5.6 7];

%closed loop eigenvalues check
Q1_cl_eig = eig(A-B*Q1_K_bar);

eA =@(h) [exp(4*h) 7.5/6*(exp(4*h)-exp(-2*h)); 0 exp(-2*h)];
N = 50;
dt = 100;
Q1_save_eig = zeros(1,N);
for i = 1:N
    h = (i-1)/dt;
    Q1_cl_A = (eA(h) - (eA(h) - eye(2))/A*B*Q1_K_bar);
    Q1_save_eig(i) = max(abs(eig(Q1_cl_A)));
end

figure();
plot([0,(N-1)/dt],[1,1],'r--','LineWidth',1.5);
 hold on;
 plot((0:N-1)/dt,Q1_save_eig,'LineWidth',2);
 hold off;
axis([0,N/dt,0,4]);
xlabel('h');
ylabel('max absolute eigenvalue');
title('eigenvalues of the cl system for various h');
legend('max abs eigenvalue','stable bound')

%% Question 2
%h  tau; extended state xe_k = [x_k u_k-1]
Q2_F =@(h,tau) [eA(h) (eA(h)-eA(h-tau))/A*B; zeros(1,3)];
Q2_G =@(h,tau) [(eA(h-tau)-eye(2))/A*B; 1];
Q2_cl_A =@(h,tau) Q2_F(h,tau) - Q2_G(h,tau)*[Q1_K_bar 0];

N = 50; %values of h
M = 10; %values of tau
dt = 100;
Q2_save_eig = zeros(M,N);
for i = 1:N
    h = (i-1)/dt;
    for j = 1:M
        tau = h/j;
        Q2_save_eig(j,i) = max(abs(eig(Q2_cl_A(h,tau))));
    end
end

figure();
X = 0:1/dt:(N-1)/dt;
plot([0,(N-1)/dt],[1,1],'r--','DisplayName',"stable bound",'LineWidth',1.5)
hold on
for j = 1:M
   plot(X,Q2_save_eig(j,:),'DisplayName',"tau = h/"+num2str(j),'LineWidth',2)
end
hold off
legend;
xlabel('h');
ylabel('max absolute eigenvalue');
axis([0,(N-1)/dt,0,4]);
title('eigenvalues of the cl system for various h and delays tau');

%% improved gain Q2
%choice of h=0.3
h = 0.3; %without delays, controller is stable
Q2_K_bar_new = [5.6 7 0.3]; 
Q2_save_eig_2 = zeros(1,10);
range = zeros(1,10);
for j = 1:10
    tau = h/j; 
    range(j) = tau;
    Q2_save_eig_2(j) =  max(abs(eig(Q2_F(h,tau) - Q2_G(h,tau) * Q2_K_bar_new)));
end

figure;
plot(range,ones(1,10),'--r','DisplayName',"stable bound",'LineWidth',1.5 );
hold on
plot(range,Q2_save_eig(:,31),'DisplayName',"original gain",'LineWidth',2 );
plot(range,Q2_save_eig_2,'DisplayName',"new gain",'LineWidth',2 );
hold off
xlabel('tau');
ylabel('max absolute eigenvalue');
title('improving controller for h=0.3 w.r.t. delays');
legend;
axis([h/10,h,0,4]);

%% Question 3
%extended state xe_k = [x_k u_k-1 u_k-2]
%Q3_F =@(h,tau) [eA(2*h) eA(h)*(eA(2*h-tau)-eye(2))/A*B+(eA(h)-eA(2*h-tau))/A*B];
%Q3_F =@(h,tau) [eA(2*h) (eA(3*h-tau)-eA(2*h-tau))/A*B (eA(2*h)-eA(3*h-tau))/A*B
%                0 0 0 0; 0 0 1 0];
Q3_F =@(h,tau) [eA(h) (eA(2*h-tau)-eye(2))/A*B (eA(h)-eA(2*h-tau))/A*B
                0 0 0 0; 0 0 1 0];
Q3_G =@(h,tau) [0;0 ; 1; 0];
Q3_cl_A =@(h,tau) Q3_F(h,tau) - Q3_G(h,tau)*[Q1_K_bar 0 0];

N = 75; %values of h
M = 10; %values of tau
dt = 500;
Q3_save_eig = zeros(M,N);
for i = 1:N
    h = (i-1)/dt;
    for j = 1:M
        tau = h+h/j;
        Q3_save_eig(j,i) = max(abs(eig(Q3_cl_A(h,tau))));
    end
end

figure();
X = 1/dt:1/dt:N/dt;
plot([1/dt,N/dt],[1,1],'r--','DisplayName',"stable bound",'LineWidth',1.5)
hold on
for j = 1:M
   plot(X,Q3_save_eig(j,:),'DisplayName',"tau = h+h/"+num2str(j),'LineWidth',2)
end
hold off
legend;
axis([0,(N-1)/dt,0.7,2]);
xlabel('h');
ylabel('max absolute eigenvalue');
title('eigenvalues of the cl system for various h and delays tau');

%% improve gain Q3
%choice of h = 0.1
h = 0.1;
Q3_K_bar = zeros(10,4);
for i = 1:10
    tau = h+h/i;
    Q3_K_bar(i,:) = place(Q3_F(h,tau),Q3_G(h,tau),[-0.2,0.2,0.1,-0.1]);
end

N = 10;
Q3_K_bar_old = [5.6 7 0 0]; 
Q3_K_bar_new = [5.6 7 0.3 0.2]; 
Q3_save_eig_old = zeros(1,N);
Q3_save_eig_2 = zeros(1,N);
range = zeros(1,N);
for j = 1:N
    tau = h;
    if j<N, tau = h+h/j; end 
    range(j) = tau;
    Q3_save_eig_old(j) =  max(abs(eig(Q3_F(h,tau) - Q3_G(h,tau) * Q3_K_bar_old)));
    Q3_save_eig_2(j) =  max(abs(eig(Q3_F(h,tau) - Q3_G(h,tau) * Q3_K_bar_new)));
end

figure;
plot(range,ones(1,N),'--r','DisplayName',"stable bound",'LineWidth',1.5 );
hold on
plot(range,Q3_save_eig_old,'DisplayName',"original gain",'LineWidth',2 );
plot(range,Q3_save_eig_2,'DisplayName',"new gain",'LineWidth',2 );
hold off
xlabel('tau');
ylabel('max absolute eigenvalue');
title('improving controller for h=0.1 w.r.t. delays');
legend;
axis([h,2*h,0.7,2]);

%% Question 4
%unfortunatly I did not manage to obtain a working optimization problem

%{
eps = 1e-5 ;
h = sdpvar(1);
N = sdpvar(1,4); %N = KP^-1
M = sdpvar(4); %M = P^-1
Objective = -h;
Constraints = [M>=eps*eye(4); h>=eps];
Constraints = [Constraints; [M ([Q2_F(h,0.2*h) zeros(3,1); 0 0 1 0]*M-[Q2_G(h,0.2*h) ;zeros(1)]*N)';[Q2_F(h,0.2*h) zeros(3,1); 0 0 1 0]*M-[Q2_G(h,0.2*h) ;zeros(1)]*N M]>=eps*eye(8)];
Constraints = [Constraints; [M ([Q2_F(h,0.5*h) zeros(3,1); 0 0 1 0]*M-[Q2_G(h,0.5*h) ;zeros(1)]*N)';[Q2_F(h,0.5*h) zeros(3,1); 0 0 1 0]*M-[Q2_G(h,0.5*h) ;zeros(1)]*N M]>=eps*eye(8)];
Constraints = [Constraints; [M (Q3_F(h,h)*M-Q3_G(h,h)*N)';Q3_F(h,h)*M-Q3_G(h,h)*N M]>=eps*eye(8)];
Constraints = [Constraints; [M (Q3_F(h,1.5*h)*M-Q3_G(h,1.5*h)*N)';Q3_F(h,1.5*h)*M-Q3_G(h,1.5*h)*N M]>=eps*eye(8)];        
optimize(Constraints,Objective);

P = inv(value(M));
K = value(N)*P;
h_ans = value(h);
%}
