% test of SPA in the ULA case

clear all
% close all
% clc


K = 2;          % source number
M = 22;         % array length
B = 6;          % beamspace dimension
N = 100;        % snapshot number
sigma = 1;      % noise variance
rmes = 0;
degrad = pi/180;
% true frequency and power
theta = [-10;5];
% f = [.1; .11; .5];
f = 1/2*sin(theta*degrad);
p = [1; 1];
% noise variance mode, 1: equal; 2: different
mode = 1;
rmse_spa = 0;
rmse_cam = 0;
rmse_epuma = 0;
rmse_sparrow = 0;
rmse_music = 0;
rmse_bsw = 0;
rmse_gl = 0;
rmse_bsadmm = 0;
snr_index = 1;
mse_spa = zeros(1,11);
mse_cam = zeros(1,11);
mse_epuma = zeros(1,11);
mse_sparrow = zeros(1,11);
mse_music = zeros(1,11);
mse_bsw = zeros(1,11);
mse_gl = zeros(1,11);
mse_bsadmm = zeros(1,11);
alpha = 0;
A = exp(1i*2*pi*kron((0:M-1)',f'));
S = sqrt(diag(p))*exp(1i*2*pi*rand(K,N));
% load signal0714
load sig1020;
S(2,:) = alpha*S(1,:)+(1-alpha)*S(2,:);


iter = 100;
SNR = linspace(-15,10,11);
Beam_w = sqrt(1/M)*exp(1i*2*pi/M*(-1:B-2)'*(-(M-1)/2:(M-1)/2));

for iS = 1:11
    snr = SNR(iS);
    parfor iter_num=1:iter
        warning off;
        Yt        = A*S;
        Y = awgn(Yt,snr,'measured');
        
%         [freq_spa, pow_spa] = SPA(Y, Beam_w,K);
        [freq_music,pow_music] = BS_MUSIC(Y,Beam_w,K);
        [freq_GL,pow_music] = GL_ANM(Y,Beam_w,K);        
        [freq_cam, pow_cam] = Beam_ANM(Y, Beam_w,K);
        [freq_bsadmm]=BS_admm(Y,Beam_w, K,K);

        % MSE
        rmse_cam = rmse_cam+sum((theta - sort(freq_cam(1:K))).^2) ;
        rmse_music = rmse_music+sum((theta - sort(freq_music(1:K))).^2);
%         rmse_spa = rmse_spa+sum((theta' - sort(freq_spa(1:K))).^2) ;        
        rmse_bsadmm = rmse_bsadmm +sum((theta - sort(freq_bsadmm(1:K))).^2) ;
        rmse_gl = rmse_gl +sum((theta - sort(freq_GL(1:K))).^2) ;
    end
    mse_cam(snr_index)=sqrt(rmse_cam/iter)/ K;
%     mse_spa(snr_index)=sqrt(rmse_spa/iter)/ K;
    mse_bsadmm(snr_index)=sqrt(rmse_bsadmm/iter)/ K;
    mse_music(snr_index)=sqrt(rmse_music/iter)/ K;
    mse_gl(snr_index)=sqrt(rmse_gl/iter)/ K;
    snr_index = snr_index + 1;
    rmse_spa = 0;
    rmse_cam = 0;
    rmse_epuma = 0;
    rmse_sparrow = 0;
    rmse_bsadmm = 0;
    rmse_music = 0;
    rmse_bsw = 0;
    rmse_gl = 0;
end
snr = SNR;
mz = 4;
lw = 0.5;
load crb;
figure
semilogy(snr,mse_music, '-p', 'markersize', mz, 'linewidth', 2); hold on;
% semilogy(snr,mse_spa, '-o', 'markersize', mz, 'linewidth', 2); 
semilogy(snr,mse_gl, '-d', 'markersize', mz, 'linewidth', 2);
semilogy(snr,mse_cam, '->', 'markersize', mz, 'linewidth', 2);
semilogy(snr,mse_bsadmm, '-s', 'markersize', mz, 'linewidth', 2)
semilogy(snr,crb, '-k', 'markersize', mz, 'linewidth', 2)
grid on;
legend('MUSIC','GL-ANM','BS-ANM','BS-ADMM','CRB')
xlabel('SNR (dB)')
ylabel('RMSE (\circ)')