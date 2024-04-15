function [pars,sse,f] = fit_comps_ave(data,pars,model_info,A_con,b_con,display,options)
% constraint x <- A_con*x + b_con
if ~exist("A_con","var")
    A_con = [];
end

if ~exist("display","var") || isempty(display)
    display = true;
end

% weights
if ~isfield(model_info,"weights")
    model_info.weights = ones(size(data,1),1);
end
sqw(:,1) = sqrt(model_info.weights);

% specify options as normally to lsqcurvefit
if ~exist("options","var") || isempty(options)
    if all(model_info.comps~=4) && all(model_info.comps~=5)
        options = optimoptions(@lsqcurvefit,'Display','off','SpecifyObjectiveGradient',true);
    else
        options = optimoptions(@lsqcurvefit,'Display','off','SpecifyObjectiveGradient',false);
    end
end

% compartment codes:
% 0: dot
% 1: stick
% 2: ball
% 3: ellipsoid
% 4: Gaussian sphere
% 5: karger model with sticks and ball 
lb = [];
ub = [];
Dlb = 0;
Dub = 3;
comps = model_info.comps;
for n = 1:length(comps)
    if comps(n)==0
    elseif any(comps(n)==[1 2])
        lb = [lb Dlb];
        ub = [ub Dub];
    elseif comps(n)==3
        lb = [lb Dlb Dlb];
        ub = [ub Dub Dub];
    elseif comps(n)==4
        lb = [lb -20];
        ub = [ub 20];
    elseif comps(n)==5
        lb = [lb 0.01 Dlb Dlb 0];
        ub = [ub 0.99 Dub Dub 1];
    end
end
lb = lb';
ub = ub';

dims = size(pars);
dims(end+1) = 1;
N = prod(dims(2:end));
sse = zeros(dims(2:end));
f = zeros([length(comps) dims(2:end)]);
if display
    WaitMessage = parfor_wait(N, 'Waitbar', true);
else
    WaitMessage = [];
end
if isempty(A_con) % without linear parameter constraints
    parfor i = 1:N
        y = data(:,i);
        [pars(:,i),sse(i)] = lsqcurvefit(@(x,~)S_fun(x,y,model_info,sqw), pars(:,i), [],sqw.*y, lb,ub, options);
        [~,~,f(:,i)] = fun_comps_ave(pars(:,i), y, model_info);
        if display
            WaitMessage.Send;
        end
    end
else
    if isempty(b_con)
        b_con = zeros(size(A_con,1),1);
    end
    lb = A_con\(lb-b_con);
    ub(~isfinite(ub)) = 1e6;
    ub = A_con\(ub-b_con);
    ub(ub>=1e6) = inf;
    parfor i = 1:N
        y = data(:,i);
        [pars(:,i),sse(i)] = lsqcurvefit(@(x,~)con_S_fun(x,y,model_info,sqw,A_con,b_con), pars(:,i), [],sqw.*y, lb,ub, options);
        [~,~,f(:,i)] = fun_comps_ave(A_con*pars(:,i)+b_con, y, model_info);
        if display
            WaitMessage.Send;
        end
    end
end


% the purpose of these helper functions are only to apply weights to the fit
function [S,J] = S_fun(x, y, model_info, sqw)
if nargout==2
    [S,J] = fun_comps_ave(x, y, model_info);
    S = sqw.*S;
    J = sqw.*J;
else
    S = fun_comps_ave(x, y, model_info);
    S = sqw.*S;
end

function [S,J] = con_S_fun(x, y, model_info, sqw, A_con,b_con)
if nargout==2
    [S,J] = fun_comps_ave(A_con*x+b_con, y, model_info);
    J = J*A_con;
    S = sqw.*S;
    J = sqw.*J;
else
    S = fun_comps_ave(A_con*x+b_con, y, model_info);
    S = sqw.*S;
end

