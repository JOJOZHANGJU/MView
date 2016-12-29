function fileList = gfl(mask,showExt,mostRecentOnly)
%GFL  - return a list of files as a cell array
%
%	usage:  fileList = gfl(mask, showExt, mostRecentOnly)
%
% Returns a cell string array of filenames matching MASK (default = '*')
%
% The extension of matched files is not included by default; specify 
% non-zero SHOWEXT to include it
%
% By default all files matching MASK are returned.  To return only the
% most recent in cases where repeated tokens were retained specify
% non-zero MOSTRECENTONLY
%
% Examples:
%	return all files matching "M01*.mat" w/o the ".mat" extension
% >> fl = gfl('M01*.mat');
%
%	return all files matching "M01*.emd" including the ".emd" extension
% >> fl = gfl('M01*.emd',1);
%
%	return most recent files matching "M01*.daq"
% >> fl = gfl('M01*.daq',0,2)

% mkt 05/01

if nargin<1 | isempty(mask), mask = '*'; end;
if nargin<2 | isempty(showExt), showExt = 0; end; 
if nargin<3 | isempty(mostRecentOnly), mostRecentOnly = 0; end; 

% build list
fileList = [];
files = dir(mask);
fileList = {files.name}';

% strip extensions
[p,f,ext] = fileparts(mask);
fl = fileList;
for k = 1 : length(fl),
	[p,f] = fileparts(fl{k});
	fl{k} = f;
end;

% delete all but most recent
%	old format:  foo.emd, foo_1.emd, foo_2.emd, etc.
%	new format:  foo_01.daq, foo_02.daq, etc.
% ASSUMPTION:  all reps end in digit

if mostRecentOnly,
	fl2 = {};
	if showExt, fl1 = fileList; else, fl1 = fl; end;
	fi = 1;
	while fi <= length(fl),
		k = fl{fi}(end);
		if k>='0' && k<='9',	% number -- assume first of possible new format series
			n = findstr(fl{fi},'_');
			k = strmatch(fl{fi}(1:n(end)),fl);
		else,				% not a number -- assume first of possible old format series
			k = strmatch(fl{fi},fl);
		end;
		fl2{end+1,1} = fl1{k(end)};
		fi = k(end) + 1;
	end;
	fileList = fl2;

elseif ~showExt,
	fileList = fl;
end;

