function mOut = lp_template(lpState, action, varargin)
%LP_TEMPLATE  - MELBA/MVIEW user labelling procedure
%
% This is an example MELBA/MVIEW labelling procedure.  It mimics default labelling
% behavior, but draws labels with dashed lines and in a user-specified color (either
% single letter or RGB triplet).  Cursor movement is tracked while the mouse button 
% is down (see the MOVE handler), and the label is created on release (see UP).
%
% Clone & tailor this file to create new user labelling procedures.
%
% User labelling procedures must supply four handlers:
%	CONFIG	- process any configuration necessary
%	DOWN	- mouseDown handler
%	EXPORT	- label export handler
%	PLOT	- label plotting handler
%
% User procedures may store internal state information in a variable which is
% passed as the first argument supplied to each handler; the CONFIG handler is
% called when the user procedure is first selected or subsequently reconfigured
% and can be used to initialize and modify this.  The DOWN handler is called upon 
% a modified mouse click or when the user chooses menu option 'Make Label.'  A
% label can be created immediately, or on the associated mouseUp event (see below
% for an example of button-down tracking).  The EXPORT handler is used to write
% label values to a tab-delimited spreadsheet file.  The PLOT handler is called 
% whenever label plotting objects need to be created.
%
% To create a label, user labelling procedures must at some point call
%
%	feval(GetCaller, 'MAKELBL', label);	% MELBA/MVIEW generic version
%	melba('MAKELBL', label);			% MELBA-specific
%
% where LABEL is a struct with at least the following fields defined
%   NAME     - string specifying displayed label name
%   OFFSET   - label offset from time 0 in msecs
%	VALUE	 - data written to spreadsheet on export
%   HOOK     - any additional information
%
% all of which may be empty.  If the OFFSET is empty MELBA/MVIEW sets it to the
% current cursor offset.  If the HOOK field is a string it is displayed (as 
% the "Note" field) and may be edited, else it is stored in the form passed.
%
% Access to the internal state of MELBA/MVIEW itself is via the following mechanism:
%
%	state = get('gcbf', 'userData');
%
% This returns a struct variable with the following fields of interest:
%
%	state.HEAD		- start of currently displayed selection (msecs)
%	state.TAIL		- end of displayed selection
%	state.CURSOR	- current cursor location (msecs)
%	state.LPSTATE	- user procedure internal state
%
% (MELBA)
%	state.SIGNAL	- the full data vector
%	state.SRATE		- its sampling rate in Hertz
%
% (MVIEW)
%	state.DATA		- concurrently sampled dataset
%
% Note that arbitrarily storing anything back into the userData field of 
% the current callback figure (gcbf) will lobotomize it....

% mkt 06/04

%	1st argument is internal state

idString = 'LP_TEMPLATE';

%	branch by action (2nd argument)

switch upper(action),
		
%-----------------------------------------------------------------------------
% CONFIG:  handle configuration
%
% 	returns MOUT = new internal state
%
% for a minimal handler just set mOut = 1 and return
% mOut = [] flags cancelled

	case 'CONFIG',
		mOut = DoConfig(lpState, idString);
		return;
		
		
%-----------------------------------------------------------------------------
% DOWN:  handle mouseDown
%
%	arg(1)	- cursor loc (msecs)
%	arg(2)	- immediate (nonzero if label set by menu command)
%
%	return MOUT = new internal state
%				= [] no update
%
% here we elect to track cursor movement while the mouse button is down, and set
% the label on its release; alternatively a label could be created here immediately

	case 'DOWN',
		mOut = [];						% no state update needed
		if nargin>3 & varargin{2},		% menu label creation -- jump to mouseUp handler
			lp_template([], 'UP', 1);
			return;
		end;

% init motion handlers (unnecessary on immediate label creation)		
		set(gcbf,  ...					
			'WindowButtonMotionFcn', sprintf('%s MOVECUR;',lpState.CALLER), ...			% default motion handler
			'WindowButtonUpFcn', 'lp_template([],''UP'');', ...	% intercept mouseUp
			'pointer', 'crosshair');							% set crosshair cursor


%-----------------------------------------------------------------------------
% EXPORT:  label export handler
%
%	arg(1)	- labels
%	arg(2)	- signal name
%
% here we write each label's name, offset, and color (hook); value is ignored
% default label export behavior may be called using
%	feval(lpState.CALLER, 'LEXPORT');

	case 'EXPORT',
		labels = varargin{1};
		name = varargin{2};
		if isempty(labels), return; end;
		[fileName, pathName] = uiputfile([name '.lab'], 'Save labels as');
		if fileName == 0, return; end;		% cancelled
		fileName = [pathName, fileName];

% open the file
		[fid, msg] = fopen(fileName, 'wt');
		if fid == -1							
			error(sprintf('error attempting to open %s', fileName));
		end;

% write headers, data
		fprintf(fid, 'LABEL\tOFFSET\tCOLOR\n');
		for n = 1 : length(labels),
			fprintf(fid, '%s\t%.1f\t%s\n', labels(n).NAME, labels(n).OFFSET, labels(n).HOOK);
		end;

% clean up
		fclose(fid);
		fprintf('Labels written to %s\n', fileName);
		
		
%-----------------------------------------------------------------------------
% PLOT:  plot label into the current axes
%
%	arg(1)	- label to plot
%	arg(2)	- yLim of current axis
%
%	return MOUT = updated label
%
% MELBA/MVIEW uses a vector of plotted object handles attached to each label for 
% figure updating while editing.  Although any object may be plotted as part of a 
% user label, each handle should be tagged 'LABEL' and stored in the label.HANDS
% vector, which is created when the label is initialized.
%
% In this example labels with an empty HOOK field are plotted as normal labels, 
% otherwise the HOOK field is assumed to hold the color of the label

	case 'PLOT',
		label = varargin{1};
		ylim = varargin{2};
		if isempty(lpState), lpState = DefCfg(idString); end;
		if isempty(label.HOOK),		% use default handler
			mOut = feval(lpState.CALLER, 'LPLOT', label, ylim);
			return;
		end;

% plot the current offset value
		switch lpState.CALLER,
			case 'melba', ylim = ylim * .98;
			otherwise,
		end;
		x = [1 1] * label.OFFSET;
		if length(label.HOOK) > 1,
			c = str2num(label.HOOK);
			if length(c) ~= 3 | any(c)<0 | any(c)>1, c = 'y'; end;
		else,
			c = label.HOOK;
		end;
		label.HANDS = line(x, ylim, ...
				'tag', 'LABEL', 'eraseMode', 'xor', ...
				'buttonDownFcn', sprintf('%s(''LMOVE'',''DOWN'');',lpState.CALLER), ...
				'lineStyle','--', 'color', c);

% plot the label name if any
		if ~isempty(label.NAME),
			label.HANDS(end+1) = text(x(1), ylim(2), [' ', label.NAME], ...
				'tag', 'LABEL', ...
				'eraseMode', 'xor', ...
				'verticalAlignment', 'top', ...
				'fontname', 'geneva', ...
				'fontsize', 9);
		end;
		mOut = label;
		

%-----------------------------------------------------------------------------
% UP:  handle mouseUp
%
%	note that we have called this handler via callback, and so LPSTATE must 
% 	be recovered from figure USERDATA

	case 'UP',
		set(gcbf, ...						% clear motion handlers
			'WindowButtonMotionFcn', '', ...
			'WindowButtonUpFcn', '', ...
			'pointer', 'arrow');
		
		state = get(gcbf, 'userData');		% retrieve state
		lpState = state.LPSTATE;			% retrieve user proc state
		label = struct('NAME', '', ...		% name supplied by dialog if necessary
						'OFFSET', [], ...	% use current cursor offset
						'VALUE', 0, ...		% extent
						'HOOK', lpState.PLOTCOLOR);			% color
		feval(lpState.CALLER,'MAKELBL', label);		% create the label
		
		
%-----------------------------------------------------------------------------
% error

	otherwise,
		error(['LP_TEMPLATE:  unrecognized action (', varargin{2}, ')']);
	
end;


%=============================================================================
% DEFCFG  - set default configuration
%
%	returns default values for lpState

function lpState = DefCfg(idString)

lpState = struct('SOURCE', idString, ...
					'PLOTCOLOR', '1 .2 .8', ...
					'CALLER', GetCaller);
				

%=============================================================================
% DOCONFIG  - config handler
%
%   returns non-empty lpState on OK, [] on cancel
%
% in this example we configure the plotting color

function lpState = DoConfig(lpState, idString)

figPos = get(0, 'ScreenSize');
width = 300; height = 170;
figPos = [figPos(1)+(figPos(3)-width)/2, figPos(2)+(figPos(4)-height)/2, width, height];

% initialize if necessary
if isempty(lpState) | ~strcmp(lpState.SOURCE, idString),
	lpState = DefCfg(idString);
end;

cfg = dialog('Name', idString, ...
	'Tag', upper(lpState.CALLER), ...
	'menubar', 'none', ...
	'Position', figPos, ...
	'KeyPressFcn', 'set(gcbf,''UserData'',1);uiresume', ...
	'UserData', 0);

% about
blurb = ['This is an example MELBA/MVIEW labelling procedure.  It draws labels ' ...
		'in the specified color (either single letter or RGB triplet).'];

uicontrol(cfg, ...
	'Style', 'frame', ...
	'Position', [10 height-70 width-20 60]);
uicontrol(cfg, ...
	'Style', 'text', ...
	'HorizontalAlignment', 'left', ...
	'String', blurb, ...
	'Position', [13 height-67 width-26 54]);

% editable text field and label
uicontrol(cfg, ...
	'Style','text', ...
	'HorizontalAlignment', 'right', ...
	'String','Label Plotting Color:', ...
	'Units', 'characters', ...
	'Position', [1 4 21 1.5]);
colorField = uicontrol(cfg, ...
	'Style', 'edit', ...
	'HorizontalAlignment', 'left', ...
	'String', lpState.PLOTCOLOR, ...
	'Units', 'characters', ...
	'Position', [23 4.2 20 1.5], ...
	'Callback', 'set(gcbf,''UserData'',1);uiresume');

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
uiwait(cfg);
if get(cfg, 'UserData'),
	c = get(colorField, 'string');
	if length(c) == 1 & ~isempty(findstr(c,'ymcrgbwk')),
		lpState.PLOTCOLOR = c;
	else,
		n = str2num(c);
		if ~isempty(n) & length(n)==3 & max(n)<=1 & min(n)>=0,
			lpState.PLOTCOLOR = c;
		end;
	end;
else,
	lpState = [];
end;
delete(cfg);


%=============================================================================
% GETCALLER  - get calling procedure name

function caller = GetCaller

stack = dbstack;
[p, caller] = fileparts(stack(end).name);
