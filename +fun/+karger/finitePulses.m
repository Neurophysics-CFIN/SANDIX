function S = finitePulses(f,D,R,delta,Delta,q)

persistent options
if isempty(options)
    options = odeset();
    options.AbsTol = 1e-12;
    options.RelTol = 1e-12;
end

q2 = q^2;
factor = (q/delta)^2;
q2_fun = @(t) (factor*t.^2).*(t<=delta) + q2*(t>delta & t<=Delta) + (factor*(delta+Delta-t).^2).*(t>Delta);
dSdt = @(t,S) (R - q2_fun(t)*D) * S;
[~,S] = ode45(dSdt,[0 delta Delta Delta+delta],f',options);
S = S(4,:);

end