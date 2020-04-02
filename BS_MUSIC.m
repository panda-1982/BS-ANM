function [freq, pow] = BS_MUSIC(Y, beam_w, K)
%
% Input
% Y: observation
% beam_w: beamspace weight
% K: the number of sources
%
% Output:
% freq: frequency estimate
% Written by Pan Jie, 2019   E-mail:panjie@yzu.edu.cn
degrad = pi/180;


[M, N] = size(Y);

if N > 1
    Rhat = beam_w*Y * Y'*beam_w' / N;
end


[U S V] = svd(Rhat);
En = U(:,K+1:end);
GG=En*En';
index = 1;
power = zeros(1,1251);
for ag = -20:0.02:10
    f = 1/2*sin(ag*degrad);
    A_t = beam_w*exp(1i*2*pi*kron((0:M-1)',f'));
    power(index) =log(1/abs(A_t'*En*En'*A_t));
    index = index +1;
end
power = power/max(power(:));
fi = Peaksearch_m(power, [-20:0.02:10]', K);
fi=sort(fi);
% postprocessing and parameter estimation
freq = fi;
pow = 0;
end

