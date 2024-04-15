close all
clearvars
load("example_data.mat","data","mask","b","Delta","delta")
%% fit single voxel
y = data(:,587); % some voxel in cortex

% set up  structure with model name and starting guesses etc.
info.model = 'eSANDIX' %Possibilities SMEX, SMEXdot, SANDIX, eSANDIX
info.Dn=1;
info.De = 2;
info.rn = .3;
info.fn =0.7;
info.D0_sphere = 2;
info.R = 10;
info.Dstick = 1.1;

% dMRI protocol
info.b = b;
info.Delta = Delta;
info.delta = delta;
info.weights = ones(size(b));
info.use_ODF = false; % solve rate equations explicitely (true) or use neural network lookup table (false)

out = fitwrapper(y,info)