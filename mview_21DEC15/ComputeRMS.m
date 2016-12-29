function [RMS,sr] = ComputeRMS(s, wl);
%COMPUTERMS  - find smoothed windowed RMS signal
%
%	usage:  RMS = ComputeRMS(s, wl)
%
% This procedure computes the Root-Mean-Square of signal S using overlapping
% Hamming windows of length WL msecs
%
% S may be an AUDIO object, a MAVIS-compatible array of structs (first element 
% assumed to be audio), a {S,SRATE} cell object, or just a vector of samples
% (in this case WL is interpreted as a sample-length window)
%
% WL is optional and defaults to 20
%
% e.g. compute RMS using a 25 msec window
% rms = ComputeRMS({s,sr}, 25);
%
% see also ComputeZC

% mkt 04/05

if nargin < 1,
	eval('help ComputeRMS');
	return;
end;
if nargin<2 | isempty(wl), wl = 20; end;

sr = [];
if isa(s,'audio'),
	sr = s.SRATE;
	s = double(s);
elseif isstruct(s),
	sr = s(1).SRATE;
	s = s(1).SIGNAL;
elseif iscell(s),
	sr = s{2};
	s = s{1};
end;
s = s(:);
if ~isempty(sr), wl = round(wl*sr/1000); end;
if wl > floor(length(s)/3), wl = floor(length(s)/3); end;

spad = [s(1:wl);s;s(length(s)-wl+1:end)];
RMS = sqrt(abs(filtfilt(2*hamming(wl)/wl,1,spad.^2)));
RMS([1:wl length(spad)-wl+1:end]) = [];
