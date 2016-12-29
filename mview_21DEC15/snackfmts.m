function [fmt,bw] = snackfmts(s, sr, varargin)
%SNACKFMTS  - interface to the KTH SNACK library for formant tracking
%
%	usage:  [fmt,bw] = snackfmts(s, sr, ...)
%
% returns the formants and associated bandwidths of signal S with sampling rate SR in Hz
% each [nPoints x nFmts]
%
% supported 'NAME',VALUE options (defaults in {brackets}):
%   OVERLAP  - time between successive frames in msecs {1}
%   WINDUR   - duration of analysis window in msecs {49}
%   WINTYPE  - window type: 0=rectangular; {1=Hamming}; 2=cos**4; 3=Hanning
%   PREEMP   - preemphasis factor {0.98}
%   NFORM    - number of formants to track (1-5) {4}
%   LPCORD   - LPC order {12}
%   THRESH   - suppress formants below this percent of max RMS-ZC ([] disables) {.1}
%   RMS      - pre-computed RMS
%   ZC       - pre-computed ZC
%
% provides an interface the SNACK library of Kare Sjolander (http://www.speech.kth.se/snack/)

% mkt 03/10 revision of de-TCLed version

if nargin < 1,
	eval('help snackfmts');
	return;
end;

if exist('snack') ~= 3,
	error('missing SNACK MEX interface for this architecture');
end;

% defaults
overlap = 1;
windur = 49;
wintype = 1;
preemp = .98;
nform = 4;
lpcord = 12;
thresh = .1;
dsfreq = 10000;		% downsampling frequency
rms = [];
zc = [];

% parse args
for ai = 2 : 2 : length(varargin),
	switch upper(varargin{ai-1}),
		case 'OVERLAP',overlap = varargin{ai};
		case 'WINDUR', windur = varargin{ai};
		case 'WINTYPE',wintype = varargin{ai};		
		case 'PREEMP', preemp = varargin{ai};
		case 'NFORM',  nform = varargin{ai};
		case 'LPCORD', lpcord = varargin{ai};
		case 'THRESH', thresh = varargin{ai};
		case 'RMS', rms = varargin{ai};
		case 'ZC', zc = varargin{ai};
		otherwise, error('unrecognized SNACKFMTS parameter (%s)', varargin{ai-1});
	end;
end;

% map to 16 bit range
s = 2^15 * (s/max(abs(s)));

% compute formants
[fmt,bw] = snack(s,sr,overlap/1000,windur/1000,wintype,preemp,nform,lpcord,dsfreq);
fmt = fmt';	bw = bw';	% [nPoints x nFmts]

% suppress below threshold
if ~isempty(thresh),
	if isempty(rms), rms = ComputeRMS({s,sr}); end
	if isempty(zc), zc = ComputeZC({s,sr}); end
	q = rms - zc;
	q = interp1(q,linspace(1,length(s),length(fmt)));
	k = find(q < max(q)*thresh);
	fmt(k,:) = NaN;
	bw(k,:) = NaN;
end;
