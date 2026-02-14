function [thetaHD,dxf,dphif,dvarphif,dP,dQ,Finfo] = identifier(thetaH,u,t,x, ...
    xf,phif,varphif,P,Q, k,l,GammaTheta)
persistent P_s Q_s te

[phi,p] = SIDBasis(x);

g=[0; 1];
varphi=g*u; % define

%--------------------------------------------------------------------------
n=length(x);
if isempty(P_s)
    P_s = zeros(p,p);
    Q_s = zeros(p,n);
    te = 0;
end

dxf=(x-xf)/k;
dphif=(phi-phif)/k;
dvarphif=(varphi-varphif)/k;

% normalization
phifbar = phif;
b = ((x-xf)/k-varphif);

% 辅助矩阵
dP=-l*P + phifbar*phifbar';
dQ=-l*Q + phifbar*b';

if rank(P)==p && te==0
    te=t 
    out=F_info(P)
end


% 足够信息矩阵选择性更新
if F_info(P_s) <= F_info(P)
    P_s = P;
    Q_s = Q;
    % fprintf('t=%.4f, F_info(P_s)=%.4f \n',t,F_info(P_s))
end
Finfo = F_info(P_s);

thetaHD=-GammaTheta*(P_s*thetaH-Q_s);
thetaHDtmp=reshape(thetaHD,p*n,1);
thetaHD = proj_rectangle(reshape(thetaH,p*n,1), thetaHDtmp, -5, 5, 1);

% meP = real(min(eig(P_s)));
end

function min_lambda = F_info(X)
    % X = (X+X')/2; %强制对称化？
    lambda = eig(X);
    lambda(lambda < 0) = 0;
    min_lambda = min(lambda);
end