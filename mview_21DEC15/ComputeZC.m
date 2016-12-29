function [zc,sr] = ComputeZC(s, wl, dontFilter)
%COMPUTEZC  - find windowed zero-crossing signal
%
%	usage:  zc = ComputeZC(s, wl, dontFilter)
%
% This procedure computes zero crossings for signal S using overlapping rectangular
% windows of length WL msecs
%
% S may be an AUDIO object, a MAVIS-compatible array of structs (first element 
% assumed to be audio), a {S,SRATE} cell object, or just a vector of samples
% (in this case WL is interpreted as a sample-length window)
%
% WL is optional and defaults to 20
%
% output is smoothed with a WIN length moving average filter unless optional
% nonzero DONTFILTER is specified
%
% see also ComputeRMS

% mkt 04/05

if nargin < 1,
	eval('help ComputeZC');
	return;
end;
if nargin<2 || isempty(wl), wl = 20; end;
if nargin<3 || isempty(dontFilter), dontFilter = 0; end

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
wl2 = ceil(wl/2);
s = [zeros(wl2,1);s;zeros(wl2,1)];
zc = filter(rectwin(wl),1,[0;abs(diff(s>=0))]);
zc = zc(wl2*2+1:end);

if dontFilter, return; end
zc = filtfilt(2*hamming(wl)/wl,1,zc);
