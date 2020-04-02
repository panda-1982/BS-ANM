function [freq] = BS_admm(Z,beam_w, tau,K)
% BS_admm implements the Beamspace Atomic Norm Minimization approach
% (BS-ANM) via ADMM
%
% Input
% Y: observation
% beam_w: beamspace weight
% K: the number of sources
%
% Output:
% freq: frequency estimate
% Written by Pan Jie, 2019   E-mail:panjie@yzu.edu.cn
%size of mesurements
B = size(beam_w,1);
% K = 2;
Z = beam_w*(Z*Z')*beam_w';
[U S V] = svd(Z);
Z = U(:,1:tau);
Z = Z*pinv(Z);
Z_size = size(Z);
n = Z_size(1,1);
L = Z_size(1,2);
N = size(beam_w,2);
degrad = pi/180;
%parameters
maxIter = 500; % maximum number of ADMM steps
rho = 0.75; % penalty parameter in augmented Lagrangian
tol_abs = 1e-3; %absolute tolerance 1e-3
converged = 0;

%initialization
YOld = zeros(n+L,n+L);
Lambda = zeros(n+L,n+L);
e1 = toeplitz_adjoint(eye(n),beam_w);
[T_r, T_i]= D_adjoint(beam_w);

for count = 1:1:maxIter  
    
    %update the variables W,X, and u
    W = YOld(n+1:n+L,n+1:n+L)...
        + (1/rho)*(Lambda(n+1:n+L,n+1:n+L) - (tau/2)*eye(L));
    
    T = [T_r -T_i*1i];
    u_v = (2/rho)*pinv([T_r -T_i*1i])*(toeplitz_adjoint(Lambda(1:n,1:n),beam_w)... 
         + rho/2*toeplitz_adjoint(YOld(1:n,1:n),beam_w) - tau/2*e1);
    u = u_v(1:N)+1i*u_v(N+1:2*N);
    X = Z;
    u_old = [real(u);imag(u)];
    %temp is the matrix that should be psd.
    temp = [beam_w*toeplitz((u))*beam_w', X; X', W];
    %projection of Q onto the semidefinite cone
    Q = temp - Lambda/rho;
    [V,E] = eig(Q);   
    e = real(diag(E));
    idx = ((e)>0);
    Y = V(:,idx)*diag(e(idx))*V(:,idx)';
    Y = (Y+Y')/2;   
    sum(diag(temp(1:B,1:B)));
    %stop criteria.
    pri_res = temp - Y;
    res = norm(Y(1:B,1:B)-YOld(1:B,1:B),'fro')/norm(Y(1:B,1:B),'fro');
    
    converged = (res<tol_abs);
    if converged, break; end
   
    Lambda = Lambda + rho*(Y - temp);
    YOld = Y;    
end

if count == maxIter
%     fprintf('ADMM fails!');
end
Tu = Y(1:B,1:B);
R_cov = Tu;
M = size(beam_w,2);
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



function T = toeplitz_adjoint(A,beam_w)
N = size(beam_w,2);
T = zeros(N,1);
T(1) = sum(diag(beam_w'*A*beam_w));
for n = 1:(N-1)
    T(n+1) = sum(diag(beam_w'*A*beam_w,n));
end

function [T_r, T_i] = D_adjoint(beam_w)
N = size(beam_w,2);

T_r = zeros(N,N);
T_i = zeros(N,N);
for m = 0:(N-1)
    for n = 0:(N-1)
        eta1 =diag(ones(N-m,1), m);
        eta2 =(diag(ones(N-abs(-n),1),-n)+diag(ones(N-abs(-n),1),n));
        if n ==0
            eta2 = eta2/2;
        end
        T1 = eta1*(beam_w'*beam_w).';
        T2 = eta2*(beam_w'*beam_w).';
        T_r(m+1,n+1) = sum(diag(T1*T2));
    end
end
for m = 0:(N-1)
    for n = 0:(N-1)
        eta1 =diag(ones(N-m,1), m);
        eta2 = -diag(ones(N-abs(-n),1),-n)+diag(ones(N-abs(-n),1),n);
        if n ==0
            eta2 = eta2/2;
        end
        T1 = eta1*(beam_w'*beam_w).';
        T2 = eta2*(beam_w'*beam_w).';
        T_i(m+1,n+1) = sum(diag(T1*T2));
    end
end


