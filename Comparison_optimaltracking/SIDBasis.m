function [sigth,p]=SIDBasis(x)

% sigth=[tanh(x)];
% sigth=[x(1);
%     x(2)
%     x(1)^3;
%     x(2)*abs(x(2));
%     tanh(10*x(2));
%     % tanh(50*x(2));
%     abs(x(2))*tanh(10*x(2));
%     exp(-0.5*abs(x(2)))*tanh(20*x(2));
%     exp(-2.0*abs(x(2)))*tanh(20*x(2))];

sigth=[x(1);
    x(2);
    tanh(50*x(2));
    exp(-abs(x(2)));
    exp(-abs(x(2)))*tanh(50*x(2))];

p=size(sigth,1);