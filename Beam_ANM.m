function [freq, pow] = Beam_ANM(Y, beam_w, K)

% [freq, pow] = Beam_ANM(Y, beam_w, Omega)
%
% Beam_ANM implements the Beamspace Atomic Norm Minimization approach (BS-ANM) for linear array
% signal processing.
%
% Input
% Y: observation
% beam_w: beamspace weight
% K: the number of sources
%
% Output:
% freq: frequency estimate
% Written by Pan Jie, 2019   E-mail:panjie@yzu.edu.cn


[M, N] = size(Y);


if N > 1
    Rhat = Y * Y' / N;
end

cvx_quiet true
cvx_precision default

[freq, pow] = BSANM(Rhat, beam_w, N,K);




end



function [freq, pow] = BSANM(Rhat,beam_w, N,K)
% SPA when Rhat is nonsingular
    M = size(Rhat, 1);
    B = size(beam_w,1);
    Rhat = beam_w*Rhat*beam_w';
    Rhat =Rhat/sum(real(diag(Rhat)));
    delta = svd(Rhat);
    delta = sqrt(N/3600)*B*delta(end);
    eta = (delta);
    degrad = pi/180;
    
    cvx_solver sdpt3
    cvx_begin sdp
    variable x(B,B) hermitian,
    variable u(M) complex,
    variable R(B,B) hermitian
    
    [x R; R beam_w*toeplitz(u)*beam_w'] >= 0,
    norm(R-Rhat,'fro')<=eta
    minimize trace(x) + real(trace(beam_w*toeplitz(u)*beam_w'));
    cvx_end
    sig = 0;
    R_cov = beam_w*toeplitz(u)*beam_w';
    [U SS V] = svd(R_cov);
    En = U(:,K+1:end);
    GG=beam_w'*En*En'*beam_w;
    MM=M+1;
    a = zeros(2*MM-1,1);
    for j=-(MM-1):(MM+1)
        a(j+MM) = sum( diag(GG,j) );
    end
    
    ra=roots([a]);
    rb=ra(abs(ra)<1);
    freq_bsanm=-asin(angle(rb)/pi)/degrad;
    rc_i = freq_bsanm<10&freq_bsanm>-20;
    rc = rb(rc_i);
    % pick the n roots that are closest to the unit circle
    [~,I]=sort(abs(abs(rc)-1));
    if size(I,1)>=K
        w=angle(rc(I(1:K)));
    else
        w=angle([rc ;zeros(K-size(I,2),1)]);
    end
    freq=sort(-asin(w/pi)/degrad);
    pow = 0;
end





