% convert from mean anomaly to true anomaly
function f = M2f(M,e)
E = keplers_eqn(M,e,M);
f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
end