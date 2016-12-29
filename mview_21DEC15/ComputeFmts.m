function [fmts,bw,sr] = ComputeFmts(s, alg, order, isF, outSR, varargin);
%COMPUTEFMTS  - compute formant tracks
%
%	usage:  [fmts,bw,sr] = ComputeFmts(s, alg, order, isF, outSR, ...)
%
% This procedure computes continuous formant tracks through signal S
% using one of two ALGorithms:  'SNACK' (default) | 'PRAAT'
%
% S may be an AUDIO object, a MAVIS-compatible array of structs 
% (first element assumed to be audio), a string interpreted as 
% a MS WAV filename, or a {S,SRATE} cell object
%
% optional ORDER defaults to sampling rate in kHz + 4
%
% specify nonzero ISF for female talkers
%
% values resampled to OUTSR Hz rate if specified
%
% returns fmts [nSamps x nFmts] and optional sampling rate SR
% for the SNACK method formant bandwidths BW may also be returned
%
% additional supported 'NAME',VALUE parameters passed to SNACKMEX:
%	NFORM  - number of formants; default is 5
%	PREEMP - pre-emphasis; default is .98
%	RMS    - precomputed RMS matching S
%	ZC     - precompute ZC matching S
%
% see also TRACKFMTS, SNACKMEX, PRAATFMTS

% mkt 10/08

% snack parameters
nform = 5;		% default number of formants
pe = .98;		% default pre-emphasis
rms = [];
zc = [];

% parse args
if nargin < 1,
	eval('help ComputeFmts');
	return;
end;
fName = '';
if isa(s,'audio'),
	sr = s.SRATE;
	s = double(s);
elseif isstruct(s),
	sr = s.SRATE;
	s = s.SIGNAL;
elseif ischar(s),
	fName = s;
	[p,f,e] = fileparts(fName);
	if isempty(e), fName = fullfile(p,[f,'.wav']); end;
	[s,sr] = wavread(fName);
elseif iscell(s),
	sr = s{2};
	s = s{1};
else,
	error('argument error (signal)');
end;
if nargin<2 || isempty(alg), alg = 'SNACK'; end;
if nargin<3, order = []; end
if nargin<4 || isempty(isF), isF = 0; end;
if nargin<5, outSR = []; end;
for k = 2 : 2 : length(varargin),
	switch upper(varargin{k-1}),
		case 'NFORM', nform = varargin{k};
		case 'PREEMP', pe = varargin{k};
		case 'RMS', rms = varargin{k};
		case 'ZC', zc = varargin{k};
		otherwise, error('unrecognized SNACKMEX parameter (%s)', varargin{k-1});
	end;
end;

% switch by algorithm
switch upper(alg),
	case 'SNACK',
		if isempty(order), order = round(sr/1000) + 4; end;
		[fmts,bw] = snackfmts(s,sr,'preemp',pe,'lpcord',order,'nform',nform,'rms',rms,'zc',zc);
	case 'PRAAT',
		doDel = 0;
		if isempty(fName),
			if max(abs(s)) > 1, s = s / max(abs(s)); end;
			fName = 'TempPraatData.wav';
			if verLessThan('matlab','8.3.0'),
				wavwrite(s,sr,fName);
			else,
				audiowrite(fName,s,sr);
			end
			doDel = 1;
		end;
		fmts = PraatFmts(fName, isF);
		bw = [];
		if doDel, delete(fName); end;
	otherwise,
		error('unrecognized algorithm (%s)', alg);
end;

% set low amplitude regions to NaN
if isempty(rms), 
	rms = ComputeRMS({s,sr}); 
	rms = rms - min(rms); 
end;
% if isempty(zc) zc = ComputeZC({s,sr}); end;
% q = (rms/max(rms)) - (zc/max(abs(zc)));
q = (rms/max(rms));
q = interp1(q,linspace(1,length(s),size(fmts,1)));
fmts(find(q<.05),:) = NaN;
if size(fmts,1) == 1,
	outSR = sr;
	return;
end;

% resample if necessary
if isempty(outSR),
	dur = 1000*length(s)/sr;
	sr = sr * size(fmts,1) / length(s);
end;
for fi = 1 : size(fmts,2),
	f = fmts(:,fi);
	idx = isnan(f);
	f(idx) = 0;
	k = linspace(1,length(f),round(length(s)/sr*outSR))';
	f = interp1(f,k);
	idx = interp1(double(idx),k);
	f(find(idx>0)) = NaN;
	if fi == 1,
		ff = zeros(size(f,1),size(fmts,2));
	end;
	ff(:,fi) = f;
end;
fmts = ff;
sr = outSR;

