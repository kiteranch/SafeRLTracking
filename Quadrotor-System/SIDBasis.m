function [sigth,p]=SIDBasis(x)
phi = x(1);
theta = x(2);
psi = x(3);
dphi = x(4);
dtheta = x(5);
dpsi = x(6);

sigth= [phi*cos(dtheta) theta*sin(phi) psi*sin(phi) dphi*cos(theta) ...
        dtheta*cos(dpsi) dpsi*sin(theta)]';

p=size(sigth,1);
end