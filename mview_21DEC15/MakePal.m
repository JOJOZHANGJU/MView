function pal = MakePal(fl, varargin)
%MAKEPAL  - construct palate trace from specified files
%
%	usage:	pal = MakePal(fl, ...)
%
% This procedure scans MAVIS-compatible array-of-struct data stored in mat  
% files specified by cellstr array FL.  Data are assumed to be oriented such
% that X:post>ant, Y:left>right, Z:inf>sup.  The highest inferior/superior (Z) 
% point at each posterior/anterior (X) offset sampled in RESolution steps 
% (default 1mm) is used to obtain a series of points through which an 
% interpolating spline is fit.
%
% returns PAL [nPoints x X,Y,Z]
%
% The following optional 'NAME',VALUE argument pairs are supported 
% (defaults shown as {VALUE}):
%
%	NAMES	- cellstr array of sensor names to scan (empty for all):  {'T*'} 
%	NOISY   - display diagnostics; plot distribution with palate trace:  {'T'}
%	RES     - sampling resolution (mm):  {.5}
%	TOL     - smoothing tolerance (mm; [] disables):  {5}
%
% see also RESAMPLEPAL, PALSELECTOR, EXTRACTPAL, SURVEYDATA

% mkt 07/15

% parse args
if nargin<1, eval('help MakePal'); return; end;
noisy = 1;
tNames = {'T*'};	% trajectory names to match
res = .5;
tol = 5;
for ai = 2 : 2 : length(varargin),
	switch upper(varargin{ai-1}),
		case 'NAMES', tNames = varargin{ai};
		case 'NOISY', noisy = (upper(varargin{ai}(1)) == 'T');
		case 'RES', res = varargin{ai};
		case 'TOL', tol = varargin{ai};
		otherwise, error(sprintf('unrecognized MAKEPAL option (%s)', varargin{ai-1}));
	end;
end;
if ischar(tNames), tNames = {tNames}; end;

% expand wildcards
q = LoadMAT(fl{1});
sNames = {q.NAME};		% all sensor names
if isempty(tNames), 
	tNames = sNames; 
	k = cellfun('size',{q.SIGNAL},2);
	tNames(k<3) = [];	% kill non-mvt trajectories
end
nn = {};
for ni = 1 : length(tNames),
	nn = [nn , ExpWildCards(tNames{ni},sNames)];
end;
tNames = nn;

% glom data
data = [];
bad = 0;
for fi = 1 : length(fl),
	q = LoadMAT(fl{fi});
	for ti = 1 : length(tNames),
		si = find(strcmpi(tNames{ti}, {q.NAME}));
		if isempty(si), continue; end;
		s = q(si).SIGNAL(:,1:3);			% ignore non-spatial components
		k = find(~isreal(s(:,1)) | isnan(s(:,1)));
		if ~isempty(k),
			bad = bad + length(k);
			s(k,:) = [];					% strip imaginary and missing data
		end;
		data = [data ; s];
	end;
	fprintf('.');
end;
fprintf('\n');

% bin
dataXZ = data(:,[1 3]);
N = floor(max(range(dataXZ))/res);			% number of bins
xz = linspace(min(dataXZ(:)),max(dataXZ(:)),N);
Dxz = interp1(xz,1:N,dataXZ,'nearest');
Dxz = accumarray(Dxz,1,[N N])';				% 2D histogram, x:z

% make palate
x = NaN(N,1); z = x;
for k = 1 : N, 
	q = find(Dxz(:,k),1,'last'); 
	if ~isempty(q), 
		x(k) = xz(k); 
		z(k) = xz(q); 
	end; 
end;
k = find(isnan(x));
x(k) = []; z(k) = [];
y = ones(length(x),1)*nanmean(data(:,2));

% smooth palate
if ~isempty(tol), [~,z] = spaps(x,z,tol); z = z'; end;

% resample palate
pal = ResamplePal([x , y , z]);

% report/plot
if noisy,
	fprintf('\n%d files, %d samples, %d bad/missing; sensors: ', length(fl), size(data,1), bad); 
	fprintf('%s,',tNames{1:end-1});
	fprintf('%s\n        ',tNames{end});
	xyz = {'X','Y','Z'};
	fprintf('%7s',xyz{:});
	fprintf('  (mm)\n   Range:'); fprintf('%7.1f',range(data));
	fprintf('\n     Min:'); fprintf('%7.1f',min(data));
	fprintf('\n     Max:'); fprintf('%7.1f',max(data));
	fprintf('\n    Mean:'); fprintf('%7.1f',mean(data)); 
	fprintf('\n     Std:'); fprintf('%7.1f',std(data)); 
	fprintf('\n');

	figure;
	set(surf(xz,xz,Dxz),'linestyle','none');
	colormap(1-sqrt(gray));
	xlabel('X');
	ylabel('Z');
	zlabel('COUNT');
	view(2);
	grid off;
	hold on;
	line(pal(:,1),pal(:,3),zeros(size(pal,1),1),'color','r','linewidth',2);
end;

if nargout < 1, clear pal; end;
