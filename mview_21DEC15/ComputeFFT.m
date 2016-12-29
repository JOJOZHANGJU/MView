function [p,f] = ComputeFFT(s, offset, mu, win, frame)
%COMPUTEFFT  - compute spectrum
%
%	usage:  [p,f] = ComputeFFT(s, offset, mu, win, frame)
%
% given signal S this procedure computes the magnitude power spectrum centered on
% the specified OFFSET (msecs) and returns the result P in dB
%
% optionally returns corresponding frequency bins F (Hz)
%
% S may be an AUDIO object, a MAVIS-compatible array of structs (first element 
% assumed to be audio), or a {S,SRATE} cell object
%
% These arguments are optional:
%	MU specifies pre-emphasis factor (default .98)
% 	WIN specifies analysis window length (default 30 msecs)
%	FRAME specifies the number of FFT analysis points (def 1024)
%
% see also COMPUTEAOS, COMPUTELPC

% mkt 01/08

% parse args

if nargin < 2, 
	eval('help ComputeFFT');
	return;
end;
nargs = nargin;
if isa(s,'audio'),
	sr = s.SRATE;
	s = double(s);
elseif isstruct(s),
	sr = s.SRATE;
	s = s.SIGNAL;
elseif iscell(s),
	sr = s{2};
	s = s{1};
else,
	error('argument error (signal)');
end;
if nargin<3 || isempty(mu), mu = .98; end;
if nargin<4 || isempty(win), win = 30; end;
if nargin<5 || isempty(frame), frame = 1024; end;

% remove DC bias
s = s - mean(s);

% pre-emphasize
if mu > 0,
	s = filter([1 -mu], 1, s);				% s[n] = s[n] - mu*s[n-1]
	s = [s(2:end);s(end)];					% fix filter delay
end;

% get analysis frame
wins = round(win*sr/1000);
wins = wins + 1-mod(wins,2);				% make sure it's odd
offs = floor(offset*sr/1000)-1;
ht = [-1 1]*floor(wins/2)+offs;
if ht(1) < 1,
	ht = [1 wins];
elseif ht(2) > length(s),
	ht = length(s) - [wins-1 0];
end;
if diff(ht)+1 > length(s),
	ht = [1 length(s)];
	wins = length(s);
end;
s = s(ht(1):ht(2));
s = s(:)';

% compute (magnitude) power spectrum
p = abs(fft(hamming(wins)' .* s, frame));
p = p(2:frame/2+1);						% drop upper reflection
p = 20*log10(abs(p/20+eps)');
f = linspace(0,sr/2,size(p,1)+1)';
f(1) = [];
