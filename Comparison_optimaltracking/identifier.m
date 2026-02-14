function [thetaHD,dxf,dphif,dvarphif,dP,dQ,meP] = identifier(thetaH,u,t,x, ...
    xf,phif,varphif,P,Q, k,l,GammaTheta)

[phi,p] = SIDBasis(x);

g=[0; 1];
varphi=g*u; % define

%--------------------------------------------------------------------------
n=length(x);
% if t==0 %no use
%     xf=x;
%     phif=phi;
%     varphif=varphi;
% end
dxf=(x-xf)/k;
dphif=(phi-phif)/k;
dvarphif=(varphi-varphif)/k;

% normalization
phifbar = phif;
b = ((x-xf)/k-varphif);
% phifbar = phif/(1+phif'*phif);
% b = ((x2-x2f)/k-varphif)/(1+phif'*phif);


% 辅助矩阵
dP=-l*P + phifbar*phifbar';
dQ=-l*Q + phifbar*b';

thetaHD=-GammaTheta*(P*thetaH-Q);
thetaHD=reshape(thetaHD,p*n,1);

meP = min(eig(P));