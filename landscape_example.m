clear all
close all

load 4Sune.mat

info.b = b/1e3;
info.Delta = Delta;
info.delta = delta;
info.weights = ones(size(b));
info.use_ODF = false; %
info.model = 'eSANDIX';

%You can supply lower and upper bounds for the region to explore. See
%plotlandscape for the order of the parameters for the different models. If
% none are given, default values are used, see plotlandscape.m
% info.lb = [0  5   0.05 0  0 .1 ];
% info.ub = [1  16  0.95 1  1 .5 ];
%Number of grids points along each parameter dimension
info.nbins = 20;

tic
pars_best = plotlandscape(y,info); 
toc
