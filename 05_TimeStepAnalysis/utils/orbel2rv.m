function [r,v] = orbel2rv(a,e,i,W,w,nu,MU_BODY)

%Calculates pos/vel vectors from classical orbital elements
%input: a, e, i(rad), W(rad), w(rad), nu(rad), MU_BODY

% Ref: Fundamentals of Astrodynamcis and Applications
%      D. A. Vallado, pg 151 (algorithm 6)

p = a*(1-e^2);

r = [p*cos(nu)/(1+e*cos(nu))
     p*sin(nu)/(1+e*cos(nu))
    		0             ];

v = [-sqrt(MU_BODY/p)*sin(nu)
      sqrt(MU_BODY/p)*(e + cos(nu))
                0                   ];

cW = cos(W);  sW = sin(W);
cw = cos(w);  sw = sin(w);
ci = cos(i);  si = sin(i);

T = [cW*cw-sW*sw*ci   -cW*sw-sW*cw*ci   sW*si
     sW*cw+cW*sw*ci   -sW*sw+cW*cw*ci  -cW*si
     sw*si		cw*si            ci   ];

r = T*r;
v = T*v;
