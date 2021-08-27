function x = keplers_eqn(Eguess,e,M)
tol = 1e-12;
error = 1;
Eold = Eguess;
while error > tol
    Enew = Eold - kep_eqn(Eold,e,M)/dkep_eqn(Eold,e);
    error = abs(Enew-Eold);
    Eold = Enew;
end
x = Enew;
end
function f = kep_eqn(E,e,M)
f = E - e*sin(E) - M;
end
function df = dkep_eqn(E,e)
df = 1 - e*cos(E);
end