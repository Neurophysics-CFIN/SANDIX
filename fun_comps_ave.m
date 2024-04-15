function [S,J,f] = fun_comps_ave(pars,y,info,f)
% compartment keys:
% 0: dot
% 1: stick
% 2: ball
% 3: ellipsoid
% 4: Gaussian sphere
% 5: karger model with sticks and ball

comps = info.comps;
b(:,1) = info.b;
Delta(:,1) = info.Delta;
delta(:,1) = info.delta;
sqw(:,1) = sqrt(info.weights);
% clearvars info0
persistent info0 R2D

% get R2D function, which transforms radius to diffusivity
if any(comps==4)
    D0 = info.D0_sphere;
    if ~isequal(info,info0)
        info0 = info;
%         eqn = chebfun(@(x)besselj(3/2,x)-x.*besselj(5/2,x),[1 100]);
%         am = roots(eqn)';
        am = [2.08157597781806	5.94036999057271	9.20584014293666	12.4044450219020	15.5792364103872	18.7426455847748	21.8996964794928	25.0528252809929	28.2033610039524	31.3520917265645	34.4995149213670	37.6459603230864	40.7916552312719	43.9367614714198	47.0813974121542	50.2256516491831	53.3695918204908	56.5132704621986	59.6567290035279	62.8000005565198	65.9431119046552	69.0860849466452	72.2289377620155	75.3716854092873	78.5143405319308	81.6569138240367	84.7994143922025	87.9418500396598	91.0842274914688	94.2265525745684	97.3688303629010];
        B = 2/D0./delta.^2./(Delta-delta/3)*2.*delta./am.^4./(am.^2-2);
        C = 2/D0./delta.^2./(Delta-delta/3)./am.^4./(am.^2-2)./am.^2/D0;
        D1 = -D0*(Delta-delta).*am.^2;
        D2 = -D0*delta.*am.^2;
        D3 = -D0*Delta.*am.^2;
        D4 = -D0*(Delta+delta).*am.^2;
        R2D = @(r) r^4*sum( B-r^2*C.*(2+exp(D1/r^2)-2*exp(D2/r^2)-2*exp(D3/r^2)+exp(D4/r^2)) ,2);
    end
end

Ncomps = length(comps);
Sn = ones(length(b),Ncomps);
parsCount = 1;
for n = 1:Ncomps
    if comps(n)==0 % dot
        
    elseif comps(n)==1 % stick
        Dpar = pars(parsCount);
        Sn(:,n) = sqrt(pi/4./(Dpar*b+1e-16)).*erf(sqrt(Dpar*b+1e-16));
        parsCount = parsCount + 1;
        
    elseif comps(n)==2 % ball
        D = pars(parsCount);
        Sn(:,n) = exp(-D*b);
        parsCount = parsCount + 1;
        
    elseif comps(n)==3 % ellipsoid
        Dpar = pars(parsCount);
        Dper = pars(parsCount+1);
        Sn(:,n) = exp(-Dper*b).*sqrt(pi/4./((Dpar-Dper)*b+1e-16)).*erf(sqrt((Dpar-Dper)*b+1e-16));
        parsCount = parsCount + 2;
        
    elseif comps(n)==4 % sphere in Gaussian phase approximation
        R = pars(parsCount);
        Sn(:,n) = exp(-b.*R2D(R));
        parsCount = parsCount+1;
        
    elseif comps(n)==5 % karger model with isotropic sticks and ball
        Sn(:,n) = fun.karger.powderaverage_finitePulses(pars(parsCount:parsCount+3),info);
        parsCount = parsCount + 4;
        
    end
end

if nargin~=4
    f = lsqnonneg(sqw.*Sn,sqw.*y);
end
S = Sn*f;
J = [];


if nargout==2

    assert(all(comps~=4) && all(comps~=5), "Analytic gradient not implemented for Gaussian sphere or Kärger model")
    
    J = zeros(length(b),length(pars));
    parsCount = 1;
    for n = 1:Ncomps
        if comps(n)==0 % dot
            
        elseif comps(n)==1 % stick
            Dpar = pars(parsCount);
            J(:,parsCount) = f(n)/2/Dpar*(exp(-Dpar*b) - Sn(:,n));
            parsCount = parsCount + 1;
            
        elseif comps(n)==2 % ball
            J(:,parsCount) = -f(n)*b.*Sn(:,n);
            parsCount = parsCount + 1;
            
        elseif comps(n)==3 % ellipsoid
            Dpar = pars(parsCount);
            Dper = pars(parsCount+1);
            J(:,parsCount) = f(n)/2/(Dpar-Dper).*(exp(-Dpar*b) - Sn(:,n));
            J(:,parsCount+1) = -f(n)*b.*Sn(:,n) - J(:,parsCount);
            parsCount = parsCount + 2;
            
        end
    end
    
end


