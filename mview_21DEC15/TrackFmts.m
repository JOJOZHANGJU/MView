function [fmts,bw,offs,bnds] = TrackFmts(varargin)
%TRACKFMTS  - compute and optionally plot formant tracks
%
%	usage:  [fmts,bw,offs,bnds] = TrackFmts(s, ...)
%
% Given speech data S this procedure computes LPC formants over voiced regions 
% (determined by RMS thresholding).  Starting from the local RMS peak formants
% are then tracked bidirectionally using a minimum distance principle.
%
% S can be AUDIO object, a MAVIS compatible array-of-structs (of which the first
% specifies audio data), or a vector of samples followed by a scalar sampling rate;
% in this case use TrackFmts(s, sr, ...)
%
% optionally returns formant tracks FMTS, associated bandwidths BW [nFrames x NFMTS],
% msec OFFSets of each formant frame [nFrames], 
% and voiced region analysis BNDS (msecs) [nRegions x onset,offset]
%
% The following optional 'NAME',VALUE argument pairs are supported 
% (defaults shown as {VALUE}):
%	'BWCLIP' - discard formant peaks with bandwidth >= BWCLIP cutoff (Hz):  {600}
%	'F0'	 - compute F0 track (first col of returned FMTS): {[](disabled)}|'M'(ale)|'F'(emale)
%	'F0LIM'  - F0 cutoff frequency (Hz; values above this ignored):  {200 (M), 300 (F)}
%	'HICUT'	 - spectrogram display cutoff:  {Nyquist}
%	'NFMTS'  - number of formants to return (Hz):  {3}
%	'OLAP'   - analysis window overlap (msecs):  {2}
%	'OFFSET' - offset within analysis region for start of tracking (percentage):  {RMS peak}
%	'ORDER'  - number of LPC coefficients (default is Vallabha & Tuller k<.15 heuristic)
% 	'PLOT'	 - plot results:  {'T' if nargout < 1; 'F' if nargout>0}
%	'PREEMP' - pre-emphasis (0<=mu<=1):  {.98}
%	'RMSWIN' - RMS smoothing window (msecs):  {20}
%	'RMSTHR' - RMS threshold (percent of max):  {.3}
%	'SPECRES'- displayed spectrogram freq resolution factor (1 broad .. 4 narrow):  {2}
%	'WSIZE'  - analysis window size (msecs):  {25}
%
% These parameters may also be specified by defining a variable called TFParams
% in the base workspace; e.g. 
%	>> TFParams = {'ORDER',3,'BWCLIP',400} 
%
% Because the algorithm relies on propagation of the formant structure found at local 
% RMS peaks, manipulation of the RMSTHR parameter to adjust the number of analysis regions 
% and thus the number of starting formant templates can improve problematic cases.
%
% see also TRACKFMTSBINS

% mkt 11/04

%	defaults

nFmts = 3;				% number of returned formants
rmsWin = 20;			% RMS smoothing window (msecs)
rmsThr = .2;			% RMS threshold (percent of max)
mu = .98;				% pre-emphasis
wSize = 25;				% analysis frame size (msecs)
offset = [];			% starting point for analysis
oLap = 2;				% overlap (msecs)
order = [];				% LPC order
bwClip = 600;			% formant peak cutoff (Hz)
hiCut = [];				% high end cutoff (Hz)
doPlot = (nargout<1);	% plot
specRes = 2;			% spectrogram multiplier
doF0 = '';				% compute F0
f0Lim = [];				% F0 limit

mSize = 10;				% plotting marker size
grayTweak = 4;			% spectrogram plotting factor

%mSize = 2; grayTweak = 8;

%	parse args

if nargin < 1,
	eval('help TrackFmts'); 
	return;
end;

if isa(varargin{1}, 'audio'),
	s = varargin{1};
	sr = s.SRATE;
	s = double(s);
	a = 2;
elseif isstruct(varargin{1}),
	s = varargin{1};
	sr = s(1).SRATE;
	s = s(1).SIGNAL;
	a = 2;
elseif nargin<2,
	error('need sampling rate');
else,
	s = varargin{1};
	sr = varargin{2};
	a = 3;
end;

params = [evalin('base','TFParams','{}') , varargin(a:end)];
if mod(length(params),2),
	error('TrackFmts argument error:  parameters must be specified in ''NAME'',VALUE format');
end;

for ai = 1 : 2 : length(params),
	switch upper(params{ai}),
		case 'BWCLIP', bwClip = params{ai+1};
		case 'F0', doF0 = upper(params{ai+1}(1));
		case 'F0LIM', f0Lim = params{ai+1};
		case 'HICUT', hiCut = params{ai+1};
		case 'NFMTS', nFmts = params{ai+1};
		case 'OFFSET', offset = params{ai+1};
		case 'OLAP', oLap = params{ai+1};
		case 'ORDER', order = params{ai+1};
		case 'PLOT', doPlot = strcmpi(params{ai+1}(1),'T');
		case 'PREEMP', mu = params{ai+1};
		case 'RMSWIN', rmsWin = params{ai+1};
		case 'RMSTHR', rmsThr = params{ai+1};
		case 'SPECRES', specRes = params{ai+1};
		case 'WSIZE', wSize = params{ai+1};
		otherwise, error(sprintf('argument error:  %s', params{ai}));
	end;
end;
if isempty(doF0),
	doF0 = 0;
else,
	sex = strcmpi(doF0(1),'F');
	doF0 = 1;
	if isempty(f0Lim),
		if sex, f0Lim=300; else, f0Lim=200; end;
	end;
end;
if nFmts > 5, nFmts = 5; end;
if isempty(hiCut), hiCut = floor(sr/2); end;

%	compute signal RMS

win = round(rmsWin*sr/1000);
nSamps = length(s);
if win > nSamps/3, win = floor(nSamps/3) - 1; end;
b = 2 * hamming(win)/win;
rms = sqrt(abs(filtfilt(b,1,s.^2)));
rms = rms - min(rms);
rmsThr = max(rms)*rmsThr;				% threshold for formant analysis

%	find regions exceeding RMS threshold for analysis

k = find(diff(rms > rmsThr));			% above-threshold regions to analyze
if rms(1) > rmsThr, k = [1 ; k]; end;
if rms(end) > rmsThr, k = [k ; nSamps]; end;
bnds = reshape(k, 2, length(k)/2)';			% [nRegions x onset,offset] (speech sample units)

% 	pre-emphasize

urS = s;									% for plotting
urBnds = bnds;
dur = 1000*(nSamps-1)/sr;
if mu > 0,
	s = filter([1 -mu], 1, s);				% s[n] = s[n] - mu*s[n-1]
	s = [s(2:end);s(end)];					% fix filter delay
end;

% 	compute formants

s = s - mean(s);							% remove DC bias
wSizeS = floor(wSize*sr/1000);				% analysis frame length (samps)
win = hamming(wSizeS);						% window
shift = oLap*sr/1000;						% overlap (fractional samples)
nFrames = round(nSamps/shift);				% number of analysis frames
fmts = repmat(NaN,nFrames,5);				% [nFrames x F1..F5]
bw = fmts;
k = find(diff(bnds,1,2) < shift*2);
bnds(k,:) = [];								% kill excessively short regions
urBnds(k,:) = [];

% 	loop over voiced sections

s = [zeros(wSizeS,1) ; s ; zeros(wSizeS,1)];	% pad
si = ceil(wSizeS/2);							% fractional sample index
bnds = bnds + si - 1;							% padding offset
bi = 1;											% bounds index
sf = [];										% analysis region starting frame
for fi = 1 : nFrames,							% frame index

% process in-range frame
	if si >= bnds(bi,1) & si <= bnds(bi,2),			% in range
		if isempty(sf), sf = fi; end;				% remember starting frame

% compute formants this frame
		k = round(si);
		sx = s(k:k+wSizeS-1) .* win;
		R = flipud(fftfilt(conj(sx),flipud(sx)));	% unbiased autocorrelation estimate

% if order unspecified, compute optimal LPC order using a modification of the Vallabha & Tuller heuristic:
% 	Pmin = round(sr/1000)
%	Pmax = round(Pmin*1.5)
%	Popt such that reflection coefficient Kopt+1,+2,+3 all < .15
%
% cf. Vallabha, G. & Tuller, B. (2002) "Systematic errors in the formant analysis of
%	steady-state vowels," Speech Comm., 38, p141-160

		if isempty(order),
			Pmin = round(sr/1000);
			Pmax = round(Pmin*1.5);
			a = levinson(R,Pmax);
			k = poly2rc(a);				% reflection coefficients
			k = k(:)';
			q = (abs(k) < .15);			% find first sequence of 3 within tolerance
			q = [q , 0 , 0] + [0 , q , 0] + [0 , 0 , q];
			q = find(q>2);
			if isempty(q),
				Popt = Pmax;						% maxed out
			else,
				Popt = q(1) - 2;
				if Popt < Pmin, Popt = Pmin; end;	% minned out
				a = levinson(R, Popt);				% recompute
			end;
		else,
			a = levinson(R, order);					% LPC
		end;
		r = flipud(sort(roots(a)));					% roots
		f = angle(r(find(imag(r) > 0)))' * sr/(2*pi);	% formants
		b = -log(abs(r(find(imag(r) > 0))))' * sr/pi;	% bandwidths
		[f,idx] = sort(f);								% sort formants
		b = b(idx);
		k = find(b >= bwClip);						% clip formants with excessive BW
		f(k) = [];
		b(k) = [];
		nf = min(5,length(f));						% number of valid formants
		fmts(fi,1:nf) = f(1:nf);					% record formants this frame
		bw(fi,1:nf) = b(1:nf);

% collate formants for completed region
		if (si+shift>bnds(bi,2) & fi>sf+1) | (fi == nFrames & sf<nFrames-1),
			if isempty(offset),
				[v,rmsMax] = max(rms(urBnds(bi,1):urBnds(bi,2)));	% RMS max offset is starting point for formant tracking
			else,
				rmsMax = diff(urBnds(bi,:))*offset;					% explicit offset
			end;
			rmsMax = rmsMax + urBnds(bi,1) - 1;
			rmsMX(bi) = rmsMax;								% for plotting
			rmsMax = floor(((rmsMax-1)/sr)*1000/oLap)+1;	% convert to frame units
			if rmsMax <= sf, rmsMax = sf+1; elseif rmsMax >= fi, rmsMax = fi-1; end;

% starting from RMS peak, track formants bidirectionally through successive frames using 
% minimum distance between corresponding formants for collation
% assumption:  formants at RMS peak are valid

			fmts(sf:rmsMax-1,:) = flipud(track(flipud(fmts(sf:rmsMax-1,:))));	% first half
			fmts(rmsMax+1:fi,:) = track(fmts(rmsMax+1:fi,:));					% second half
			sf = [];							% reset starting frame
			bi = bi + 1;							% next bounds pair
			if bi > size(bnds,1), break; end;		% no more analysis regions
		end;	% (completed region)		
	end;	% (in-range frame)
	si = si + shift;							% shift window
end;
fmts = fmts(:,1:nFmts);							% clip to desired number of output formants
fmts(find(fmts>hiCut)) = NaN;

%	compute F0 over voiced intervals

if doF0,
	wSizeS2 = floor(20*sr/1000);
	wSizeS = wSizeS2 * 2;
	F0 = ones(nFrames,1) * NaN;
	s = [urS ; zeros(wSizeS,1)];				% pad original signal
	si = 1;										% fractional sample index
	bi = 1;										% bounds index
	for fi = 1 : nFrames,
	
% process in-range frame
		if si >= bnds(bi,1) & si <= bnds(bi,2),			% in range
			k = round(si);
			sx = s(k:k+wSizeS-1);
			F0(fi) = ComputeF0(sx, sr, sex);
			if si+shift > bnds(bi,2),
				bi = bi + 1;							% next bounds pair
				if bi > size(bnds,1), break; end;		% no more analysis regions
			end;	% (completed region)
		end;	% (in-range frame)
		si = si + shift;							% shift window
	end;	% (voiced sections)
	F0(find(F0 > f0Lim)) = NaN;						% kill values above limit
end;

%	plot

if doPlot,
	set(figure,'name',inputname(1));
	if doF0, nPan = 4; else, nPan = 3; end;
	subplot(nPan,1,1);
	ylim = max(abs(urS)) * 1.1;
	plot([urS,rms*4-ylim*.9]); 
	set(gca,'xtick',[],'ytick',[],'xlim',[1 nSamps],'ylim',[-1 1]*ylim); 
	title(inputname(1),'interpreter','none');
	h1 = line([1;1]*urBnds(:)',repmat([-1;1]*ylim,1,size(urBnds,1)*2),'color','r');
	h2 = line([1;1]*rmsMX,repmat([-1;1]*ylim,1,length(rmsMX)),'color','g');
	h3 = line(get(gca,'xlim'),[1 1]*rmsThr*4-ylim*.9,'color','m','linestyle',':');
	legend([h1(1),h2(1),h3(1)],'analysis bnds','RMS peak','RMS thresh',1);
	p1 = get(gca,'position');
	subplot(nPan,1,3);
	pos = get(gca,'position');
	pos(4) = p1(2) - pos(2);
	set(gca,'position',pos);	
	[b,f] = ComputeSpecgram(urS, sr, specRes, mu~=0);
	f(find(f > hiCut)) = [];
	b = b(1:length(f),:);
	ns = size(b,2);
	set(gca, 'xlim', [1 ns], 'ylim', [0 hiCut]);
	imagesc(1:ns, f, b);
	set(gca, 'ydir', 'normal', 'xtick',[], 'box','on');
	ylabel('Hz');
	colormap(flipud(gray(256).^grayTweak));
	axes('position', pos, 'ytick',[], 'xtick',[], 'color','none');
	line([1:nFrames],fmts,'marker','.','markerSize',mSize,'linestyle','none');
	set(gca,'xlim',[1 nFrames],'ylim',[0,hiCut]);
	axes('position', pos, 'ytick',[], 'xlim', [0 dur],'color','none');
	if doF0, 
		subplot(nPan,1,4);
		line([1:nFrames],F0, 'marker','.','markerSize',mSize/2,'linestyle','none','color',[0 .5 .5]);
		set(gca,'xlim',[1 nFrames],'ylim',[0 f0Lim],'xtick',[]);
		ylabel('F0  (Hz)');
		box on;
		axes('position', get(gca,'position'), 'ytick',[], 'xlim', [0 dur],'color','none');
	end;
	xlabel('msecs');
end;

if doF0, fmts = [F0 , fmts]; end;
if nargout < 1, clear fmts; end;
if nargout > 2,
	offs = linspace(0, dur, nFrames)';
	bnds = 1000*(bnds-1)/sr;
end;


%=============================================================================
% COMPUTESPECGRAM  - compute spectrogram
%
%	returns [nFFTs x nFrames] image, freq axis

function [b,f] = ComputeSpecgram(s, sRate, specRes, doPremp);

if doPremp, s = diff(s); end;		% 1st difference unless preemphasis set to 0
frame = 512;						% analysis frame
wSize = 6*specRes;					% analysis window (msecs)
shift = 1;							% window shift (msecs)
nSamps = length(s);
wSize = floor(wSize*sRate/1000);	% window size (samples)
wSize = wSize + mod(wSize,2);		% make sure it's even
shift = shift * sRate/1000;			% overlap (fractional samples)
nFrames = round(nSamps/shift);
w = hanning(wSize);
b = zeros(frame, nFrames);
sx = wSize/2 + 1;					% fractional sample index
s = [zeros(wSize/2,1) ; s ; zeros(wSize,1)];
for fi = 1 : nFrames,
	si = round(sx);
	pf = abs(fft(w .* s(si:si+wSize-1),frame*2));
	b(:,fi) = pf(2:frame+1);		% drop DC, upper reflection
	sx = sx + shift;
end;
b = filter(ones(3,1)/3,1,abs(b),[],2);	% clean holes
f = linspace(0,sRate/2,frame);


%=============================================================================
% TRACK  - track formants through successive frames using minimum distance
%			assumes formants at starting point (RMS peak) are valid

function fmts = track(allFmts)

[nFrames,nFmts] = size(allFmts);
fmts = repmat(NaN, size(allFmts));
fmts(1,:) = allFmts(1,:);
lastF = allFmts(1,:);
lastD = zeros(1,nFmts);
for fi = 2 : nFrames,

% find distance between all formants in this frame to good formants in previous one
	f0 = ones(nFmts,1) * lastF;
	f = allFmts(fi,:)' * ones(1,nFmts);
	[d,idx] = min(abs(f-f0));					% distances and indices of closest formants from this frame to last
	fmts(fi,:) = f(sub2ind(size(f),idx,[1:nFmts]));
if all(isnan(d)), break; end;
	d(find(isnan(d))) = 1e6;					% force exclusion for missing formants	

% check for gaps:  exclude candidates whose between frame distance jumped more than 100Hz from the last good frame
	k = find(abs(d-lastD) > 100);
% fprintf('%d   %d %d %d %d %d   %.1f %.1f %.1f %.1f %.1f   %d %d %d %d %d\n', fi, round(fmts(fi,:)), d, (abs(d-lastD) > 20*lastD));
	if isempty(k),
		lastD = d;
	else,
		fmts(fi,k) = NaN;
		notK = [1:nFmts];
		notK(k) = [];
		lastD(notK) = d(notK);
		lastF(notK) = fmts(fi,notK);
	end;
end;

