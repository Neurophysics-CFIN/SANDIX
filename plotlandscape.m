function pars_best = plotlandscape(y,in)
nbins = in.nbins;

switch in.model
    case 'eSANDIX'
        in.comps = [1 4 5]; % eSANDIX
        if ~isfield(in,'D0_sphere')
            in.D0_sphere = 2.0;
        end
        S_fun = @(x) fun_comps_ave(x,y,in);
        labels = ["$D_{stick}$","$R$","$f_n$", "$D_n$", "$D_e$", "$r_n$"];
        if all(isfield(in,{'ub','lb'}))
            bins = arrayfun(@(a,b)(linspace(a,b,nbins)),in.lb,in.ub,'UniformOutput',false);
        else
            bins = {linspace(0.1,0.9,nbins)
                linspace(5,15,nbins)
                linspace(0.1,0.9,nbins)
                linspace(0,3,nbins)
                linspace(0,3,nbins)
                linspace(0.05,0.5,nbins)};
        end
        pars_best = landscape(y, S_fun, bins, [], labels);
    case 'SANDIX'
        in.comps = [4 5]; % SANDIX
        if ~isfield(in,'D0_sphere')
            in.D0_sphere = 2.0;
        end
        S_fun = @(x) fun_comps_ave(x,y,in);
        labels = ["$Dstick$","$R$","$f_n$", "$D_n$", "$D_e$", "$r_n$"];
        if all(isfield(in,{'ub','lb'}))
            bins = arrayfun(@(a,b)(linspace(a,b,nbins)),in.lb,in.ub,'UniformOutput',false);
        else
            bins = {
                linspace(0,3,nbins)
                linspace(5,15,nbins)
                linspace(0.1,0.9,nbins)
                linspace(0,3,nbins)
                linspace(0,3,nbins)
                linspace(0.05,0.5,nbins)};
        end
        pars_best = landscape(y, S_fun, bins, [], labels);
    case 'SMEX'
        in.comps = 5; % SMEX
        S_fun = @(x) fun_comps_ave(x,y,in);
        labels = ["$f_n$", "$D_n$", "$D_e$", "$r_n$"];
        if all(isfield(in,{'ub','lb'}))
            bins = arrayfun(@(a,b)(linspace(a,b,nbins)),in.lb,in.ub,'UniformOutput',false);
        else
            bins = {
                linspace(0.1,0.9,nbins)
                linspace(0,3,nbins)
                linspace(0,3,nbins)
                linspace(0.02,0.5,nbins)};
        end
        pars_best = landscape(y, S_fun, bins, [], labels);
    case 'SMEXdot'
        in.comps = [0 5]; %
        S_fun = @(x) fun_comps_ave(x,y,in);
        labels = ["$f_n$", "$D_n$", "$D_e$", "$r_n$"];
        if all(isfield(in,{'ub','lb'}))
            bins = arrayfun(@(a,b)(linspace(a,b,nbins)),in.lb,in.ub,'UniformOutput',false);
        else
            bins = {
                linspace(0.1,0.9,nbins)
                linspace(0,3,nbins)
                linspace(0,3,nbins)
                linspace(0.05,0.5,nbins)};
        end
        pars_best = landscape(y, S_fun, bins, [], labels);
    otherwise
        error("Invalid model specified")
end
