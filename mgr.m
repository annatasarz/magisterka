beep off; clc; clear all; close all;
model=0;
iter=1;

N_in=1
% N_in_m=1;
T_d=0.1
alpha_e=log(2)/T_d

alpha=normrnd(0.0398, 0.0168)
beta=alpha/10
gamma=0.0093
%%
K=20;
ro= 0.004;

%todo
lambda=0.01; 
beta_g=1;
delta=0;
%%%%todo: add stop times to therapies (idea: step from 0 to 1 and vice versa)

t_m=normrnd(500, 100);

start_time=0;
stop_time=20000;
step=0.001;
chemo = timeseries(0,(start_time:step:stop_time)');
radio=chemo;

chemo_t=100;
radio_t=normrnd(200, 50);
radio_t=1000000;
delta=chemo;
for i=1:iter
    i
simOut=sim('mgr_sim','StartTime', num2str(start_time),'StopTime',num2str(stop_time),'FixedStep',num2str(step));
NT=simOut.get('NT');
N0=simOut.get('N0');
N1=simOut.get('N1');
end;

plot(NT.Time, NT.Data);
hold on 
plot(N1.Time, N1.Data);
hold on 
plot(N0.Time, N0.Data);

% plot(radio_t/10, 1,'*');

% plot(T_m);
