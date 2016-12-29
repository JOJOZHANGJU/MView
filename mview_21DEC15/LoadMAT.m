function var = LoadMAT(fileName, varName)
%LOADMAT  - load variable from MAT file
%
%	usage:  var = LoadMAT(fileName, varName)
%
% FILENAME holds the path to the MAT file
% if FILENAME is empty file selected using dialog box
% if FILENAME contains the '*' character it is interpreted as a mask
%
% if file holds a single variable it is returned immediately, else
% if varName is empty returned variable selected interactively, else
% returns variable matching VARNAME
%
% see also PEEKMAT

% mkt 12/98

if nargout < 1,
	eval('help LoadMAT')
	return;
end;
if nargin < 1, fileName = []; end;
if nargin < 2, varName = []; end;

%	check for mask

if findstr(fileName, '*'),
	mask = fileName;
	fileName = [];
else,
	mask = '*.mat';
end;

%	if fileName is [] obtain interactively

if isempty(fileName)
	[file, pathName] = uigetfile(mask, 'Select MAT file');
	if file == 0, return; end;
	fileName = [pathName file];
end;

%	select variable

vars = who('-file', fileName);
if isempty(vars),					% try prepending cwd
	fileName = [pwd fileName];
	vars = who('-file', fileName);
end;
if isempty(vars),
	error(sprintf('Can''t find %s', fileName));
end;
if length(vars) > 1,
	if isempty(varName),
		i = GetVar(vars);			% get selected variable
		if isempty(i), return; end;		% cancel
	else,
		i = strmatch(varName, vars, 'exact');
		if isempty(i),
			error(sprintf('Can''t find %s in %s', varName, fileName));
		end;
	end;
else,
	i = 1;
end;
name = vars{i};
var = load(fileName, name);				% load specified variable
if isstruct(var),
	eval(['var=var.',name,';']);
end;
	

%=============================================================================
% GETVAR  - return index of selected variable within vars list

function index = GetVar(vars)

figPos = get(0, 'ScreenSize');

width = 180;
height = 175;
figPos = [figPos(1)+(figPos(3)-width)/2, figPos(2)+(figPos(4)-height)/2, width, height];

cfg = dialog('Name', 'Select Variable', ...
	'Position', figPos, ...
	'ButtonDownFcn', '', ...
	'UserData', 0);

% variable list
list = uicontrol(cfg, ...
	'Position',[(width-150)/2 height-115 150 106], ...
	'String', vars, ...	
	'ListBoxTop', 1, ...
	'Style','listbox', ...
	'Value', 1, ...
	'Callback', 'if strcmp(get(gcbf,''SelectionType''),''open''),set(gcbf,''UserData'',1);uiresume;end');
	
% OK, cancel buttons
uicontrol(cfg, ...		% buttons
	'Position',[width/2-70 15 60 25], ...
	'String','OK', ...
	'Callback','set(gcbf,''UserData'',1);uiresume');
uicontrol(cfg, ...
	'Position',[width/2+10 15 60 25], ...
	'String','Cancel', ...
	'Callback','uiresume');

% wait for input
if ismac,				% focus X
	[s,r] = unix('osascript -e ''tell application "MATLAB" to activate''');
end;
uiwait(cfg);

% process response
if get(cfg, 'UserData'),
	index = get(list, 'Value');
else,
	index = [];
end;
delete(cfg);
