%% Trial of Apodization

Fs = details_para.Fs;
ppm = details_para.ppm_referenced;
[N, N_vox ] = size(ppm');
n_t = ((1:N)')/Fs;

 FID_TD = real(data_value.FID_TD);
 FID_TD = FID_TD(:,23);
 
 %FID_TD_Ref
%  FID_TD_ref = real(data_value.FID_TD_ref)
%  FID_TD = FID_TD_ref
%  plot(ppm,FID_TD)
%% Apodization

% Exponential Apodization
a1 = pi*analysis_pipeline.process(4).para(end,1);
w1 = exp(-n_t*a1);
sm_signal_exp = w1.*FID_TD;
p = plot(ppm,sm_signal_exp,ppm,FID_TD)
p(1).LineStyle = '--'
p(2).LineStyle = ':' 
ax = gca;
ax.XDir = 'reverse';
legend('exp apodized','FID TD')

% Gaussian Apodization
a2 = (pi*analysis_pipeline.process(4).para(end,2))/(2*sqrt(log(2)));
w2 = exp(-(n_t.^2)*a2.^2);
sm_signal_gauss = w2.*FID_TD;

p = plot(ppm,sm_signal_gauss,ppm,FID_TD)
ax = gca;
ax.XDir = 'reverse';
p(1).LineStyle = '--'
p(2).LineStyle = ':' 
% hold on
% plot(ppm,FID_TD)
legend('Gaussian apodized','FID TD')

% Mixed (Exponential + Gaussian) Apodization
w3 = exp(-(n_t*a1 + (n_t.^2)*(a2.^2)));
sm_signal_gauss_exp = w3.*FID_TD;
p = plot(ppm,sm_signal_gauss_exp,ppm,FID_TD)
ax = gca;
ax.XDir = 'reverse';
p(1).LineStyle = '--'
p(2).LineStyle = ':' 
legend('Mix apodized','FID_TD')

% Effects of Negative and positive exponential
neg_ex = FID_TD*exp(-0.75);
pos_ex = FID_TD*exp(0.75);
p = plot(ppm,neg_ex,'r',ppm,pos_ex,'b',ppm,FID_TD,'k')
ax = gca;
p(1).LineStyle = '--'
p(2).LineStyle = ':'
p(3).LineStyle = '-.'
ax.XDir = 'reverse';
legend('Negative Apodized', 'Positive Apodized', 'Original TD Signal')


%Plots of all Apodization outputs
plot(ppm,neg_ex)
hold on
plot(ppm,pos_ex)
hold on
plot(ppm,FID_TD)
hold on
plot(ppm,sm_signal_exp)
hold on
plot(ppm,sm_signal_gauss)
hold on
plot(ppm,sm_signal_gauss_exp)
ax = gca;
ax.XDir = 'reverse';
legend('Negative Exponential: -0.1','Positive Exponential: 0.1','Original','Exponential: 0.1','Gaussian: 1','Mix:0.1 1')

FID_FD_Apodized = fftshift(fft(sm_signal_gauss_exp));
FID_FD_Apodized_real = real(FID_FD_Apodized);
FID_FD_Apodized_imag = imag(FID_FD_Apodized);
plot(FID_FD_Apodized_real)

% Phase Correction (Zero Order and First Order)
phi_0 = analysis_pipeline.process(5).para(23);
phi_1 = analysis_pipeline.process(5).para(46);
%[N,~] = size(sig);
N = 2048;
fres = details_para.fres;
phase_mat = repmat(phi_0,[N,1]) + fres'*(2*pi*-phi_1*10^-3);  %get the phase correction values
FD_FID = FID_FD_Apodized_real
af_ap_pc = (exp(1i*phase_mat).*FD_FID); 

%Window limiting by ppm (User defined limits)
prompt = 'Input starting ppm index';
start_ppm = input(prompt);

prompt = 'Input ending ppm index';
stop_ppm = input(prompt);

start_ppm_idx = find(ppm>start_ppm);
start_ppm_idx = start_ppm_idx(1);
stop_ppm_idx = find(ppm<stop_ppm);
stop_ppm_idx = stop_ppm_idx(end);

ppm = ppm(start_ppm_idx:stop_ppm_idx);
real_FID_FD = real(af_ap_pc);

X = ppm;
real_FID_FD = real_FID_FD(start_ppm_idx:stop_ppm_idx)';

[N, N_vox ] = size(real_FID_FD');
n_t = ((1:N)')/Fs;
plot(ppm,real_FID_FD)
ax = gca;
ax.XDir = 'reverse';
legend('Real part of FID FD')

%% Baseline correction

BL = ALS(real_FID_FD', 10^3, 0.01);
BLC = real_FID_FD' - BL;
% if min(BLC)<0
%     BLC = BLC + abs(min(BLC));
% end
figure, plot(ppm, BLC)
ax = gca;
ax.XDir = 'reverse';
hold on
plot(ppm,real_FID_FD)
plot(ppm,BL)
legend('BLC','Original','Baseline')

figure, plot(ppm,BLC), title('Baseline corrected signal')
ax = gca;
ax.XDir = 'reverse';
legend('Baseline Corrected Signal')


