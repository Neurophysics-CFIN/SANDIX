function out = fitwrapper(data,in)

% fits signal using various GM models with exchange as in Olesen et al,
% https://doi.org/10.1016/j.neuroimage.2022.118976. Input "in" is a
% structure with as a minimun fiels corresponding to b, delta and Delta,
% i.e., in.b, in.delta and in.Delta, and in.model specifying the model.
% Out contains the fitted modelparameters with descriptive names, as well
% as the predicted signal in out.signal. Starting guesses
% must be specified in in.fn, etc, see below for naming conventions for the fields:
%
% fstick: stick fraction (impermeable neurites)
% Dstick stick diffusivity
% fsphere: soma fraction
% R  soma radius
% fn neurite fraction
% Dn neurite diffusivity
% De extra-cellular diffusivity
% rn neurite exchange rate.

if ~isfield(in,'weights')
    in.weights = ones(size(in.b));
end
if ~isfield(in,'use_ODF')
    in.use_ODF = false; % solve rate equations(ODF) explicitely (true) or use neural network lookup table (false)
end
fit_options = optimoptions(@lsqcurvefit,'Display','final-detailed');
switch in.model
    case 'eSANDIX'
        in.comps = [1 4 5]; % eSANDIX
        if ~isfield(in,'D0_sphere')
            in.D0_sphere = 2.0;
        end
        in.pars0 =  [in.Dstick,in.R,in.fn,in.Dn,in.De,in.rn]';
        [out.pars,out.sse,out.f] = fit_comps_ave(data, in.pars0, in, [],[], [], fit_options); % S = f(1) * S_KM + f(2) * Sdot -->
        out.fstick = out.f(1)/(out.f(1) + out.f(2) + out.f(3));
        out.Dstick = out.pars(1);
        out.fsphere = out.f(2)/(out.f(1) + out.f(2) + out.f(3));
        out.R = out.pars(2);
        out.fn = out.pars(3)*out.f(3)/(out.f(1) + out.f(2) + out.f(3));
        out.Dn = out.pars(4);
        out.De = out.pars(5);
        out.rn = out.pars(6);
        out.signal = fun_comps_ave(out.pars,[],in,out.f);
    case 'SANDIX'
        in.comps = [4 5]; % SANDIX
        if ~isfield(in,'D0_sphere')
            in.D0_sphere = 2.0;
        end
         in.pars0 =  [in.R,in.fn,in.Dn,in.De,in.rn]';
        [out.pars,out.sse,out.f] = fit_comps_ave(data, in.pars0, in, [],[], [], fit_options); % S = f(1) * S_KM + f(2) * Sdot -->
        out.fsphere = out.f(1)/(out.f(1) + out.f(2));
        out.R = out.pars(1);
        out.fn = out.pars(2)*out.f(2)/(out.f(1) + out.f(2));
        out.Dn = out.pars(3);
        out.De = out.pars(4);
        out.rn = out.pars(5);
        out.signal = fun_comps_ave(out.pars,[],in,out.f);
    case 'SMEX'
        in.comps = 5; % SMEX
        in.pars0 = [in.fn,in.Dn,in.De,in.rn]';
        [out.pars,out.sse,out.f] = fit_comps_ave(data, in.pars0, in, [],[], [], fit_options); % S = f(1) * S_KM + f(2) * Sdot -->
        out.fn = out.pars(1);
        out.Dn = out.pars(2);
        out.De = out.pars(3);
        out.rn = out.pars(4);
        out.signal = fun_comps_ave(out.pars,[],in,out.f);
    case 'SMEXdot'
        in.comps = [0 5]; %
        in.pars0 = [in.fn,in.Dn,in.De,in.rn]';
        [out.pars,out.sse,out.f] = fit_comps_ave(data, in.pars0, in, [],[], [], fit_options); % S = f(1) * S_KM + f(2) * Sdot -->
        out.fn = out.pars(1)*out.f(2)/(out.f(1) + out.f(2));
        out.fdot = out.f(1)/(out.f(1) + out.f(2));
        out.Dn = out.pars(2);
        out.De = out.pars(3);
        out.rn = out.pars(4);
        out.signal = fun_comps_ave(out.pars,[],in,out.f);
    otherwise
        error("Invalid model specified")
end

