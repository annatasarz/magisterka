beep off; clc; clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MASTER THESIS %%                                                             %
% Comparison of selected approaches to model radio- and chemotherapy effects    % 
%                                                                               %
% author: A. Tasarz                                                             %
% supervisor: prof. J. Œmieja                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initial parameters %%

model=0;                                        % growth model choice, 0 - Gompertz, 1 - exponential
radio_on=1;                                     % radiotherapy (0 - not applied)
chemo_on=1;                                     % chemotherapy (0 - not applied)

iter=1000;                                      % number of generated patients

N_in_diam=6;                                    % initial tumour diameter
K_diam=30;                                      % maximum tumour diameter
N_death=13;                                     % death condition diameter
gamma=0.0093;                                   % sensitivity to chemotherapy parameter
ro_u=7.00e-05;                                  % growth parameter mean
ro_s=7.23e-03;                                  % growth parameter std
alpha_u=0.0398;                                 % 1st radiosensitivity parameter mean
alpha_s=0.0168;                                 % 1st radiosensitivity parameter std
met_scale=1.5;                                  % parameter increasing growth rate of the metastasis

step=1;                                         % step - 1 day
duration=5;                                     % finish time - years   

radio_start=1;                                  % radiotherapy start day
radio_duration=7*7;                             % it lasts for 7 weeks
radio_stop=radio_start+radio_duration;          % radiotherapy finish time
radio_dose=63;                                  % dose which is applied during treatment in Gy

surface=1.7;                                    % average body surface   
cis_dose=0.1*surface;                           % dose of cisplatin    
vin_dose=0.05*surface*0.4;                      % dose of vinblastine

vin_days=[0 7 14 21 28];                        % days at which vinblastine is administered
cis_days=[0 28];                                % days at which cisplatin is administered
cis_hl=10;                                      % cisplatine half-life


%% Calculation and initialization of necessary parameters %%

start_time=0;                                   % starting time for the simulation 
stop_time=duration*365/step;                    % finish time for the simulation 

N_in=(N_in_diam/2)^3*pi*4/3;                    % calculating volume of initial tumour
K=(K_diam/2)^3*pi*4/3;                          % calculating volume of maximum possible tumour 
death_cond=(N_death/2)^3*pi*4/3;

N_in=N_in*5.8e8;                                % calculating no of cells of initial tumour
K=K*5.8e8;                                      % calculating no of cells of maximum possible tumour 
death_cond=death_cond*5.8e8;                    % calculating no of cells of death condition

radio_amp=radio_dose/radio_duration;            % calculating radiotherapy control amplitude
cis_det=step/cis_hl;                            % calculating cisplatine deterioration rate

check_start=radio_stop;                         % setting time of checking whether patient was cured
check_stop=check_start+step;


scale=10000;                                    % necessary scaling
mu = [alpha_u*sqrt(scale) ro_u*sqrt(scale)];    % means of alpha and rho
sigma = [alpha_s^2*scale 0.87                   % covariance matrix
        0.87 ro_s^2*scale];
R = mvnrnd(mu,sigma,iter);                      % generating parameters
R=R./sqrt(scale);                               % rescaling parameters
alpha_v=R(:,1);
ro_v=R(:,2);
corr(alpha_v, ro_v)

k_m=[];
k=0;


%% Data generation and simulation %%
for i=1:iter
    t_m   = normrnd(20, 10);                    % the time at which metastasis appears
%     t_m=10000000000;                          % ommiting metastasis
                               
    ro=ro_v(i);                                 % growth rate for the gompertz model
    ro_m=met_scale*ro;                          % growth rate for the metastasis (gompertz)

    ro_e=3*ro;                                  % growth rate for the exponential model
    ro_m_e=met_scale*ro_e;                      % growth rate for the metastasis (exponential)
    
    alpha=alpha_v(i);                           % 1st raiosensitivity parameter
    beta  = alpha/10;                           % 2nd raiosensitivity parameter
  
    if(ro>=0)                                   % eliminating patients with tumour shrinking
        k=k+1;
        simOut=sim('mgr_sim','StartTime', num2str(start_time),'StopTime',num2str(stop_time),'FixedStep',num2str(step));
        NT=simOut.get('NT');
        N0=simOut.get('N0');
        N1=simOut.get('N1');
      
        k_m=[k_m max(NT.Time)];
    end;
end;

%% Plot tumor size %%
% figure
% plot(NT.Time./365, NT.Data);
% hold on
% plot(N1.Time./365, N1.Data);
% xlabel('Time [years]') 
% ylabel('Tumor size') 
% xlim([0 stop_time]./(365))
% ylim([0 K]);
% hold on

%% Kaplan-Meier curve %%
figure
patients=(1:k)./k;
k_m=flip(sort(k_m));
plot(k_m./(365), patients, 'k');
xlim([0 stop_time]./(365));
xlabel('Time [years]') 
ylabel('Overall survival') 
hold on
  
%% notification %%
% [y,Fs]=audioread('paint it black.mp3'); 
% PO=audioplayer(y,Fs);
% play(PO);

