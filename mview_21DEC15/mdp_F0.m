function data = mdp_F0(data, varargin)
%MDP_F0  - MVIEW dataproc that computes running F0 estimate using PRAAT
%
% appends F0; assumes audio is DATA(1)
% dependencies:  ComputeF0, PraatF0, ep.praat, copy of Praat.app (or praatcon.exe)

% mkt 02/08

[F0,sr] = ComputeF0(data);
data(end+1) = struct('NAME','F0','SRATE',sr,'SIGNAL',F0);
