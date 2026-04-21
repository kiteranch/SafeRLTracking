function [sigth,p]=SIDBasis(x)

sigth=[x(1);
    x(2);
    tanh(50*x(2));
    exp(-abs(x(2)));
    exp(-abs(x(2)))*tanh(50*x(2))];

p=size(sigth,1);
end