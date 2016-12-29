function fmts = PraatFmts(fName, isF)
%PRAATFMTS  -- compute formants using Praat
%	
%	usage:  fmts = PraatFmts(fName, isF)
%
% specify nonzero ISF for female voice
%
% N.B. This procedure calls PRAAT.  For this to work properly the environment must be prepared as follows:
%
% - Windoze
% download a copy of praatcon.exe from http://www.fon.hum.uva.nl/praat/download_win.html
% and place it in the same location as this file (PraatF0)
%
% - OSX
% (assumes that Praat.app has been installed into the default /Applications folder)
% open a terminal window and cd to the directory location of this file (PraatF0)
% at the terminal command line enter the following to create a symbolic link to Praat:
%  ln -s /Applications/Praat.app/Contents/MacOS/Praat Praat
%
% - Linux
% create a symbolic link to your location of the Praat executable in the same location as this file (PraatF0)
%
% for all platforms the scripts eff.praat and efm.praat must also be present in the same location 
% as this file (PraatF0)

% mkt 08/08

if nargin < 1,
	eval('help PraatFmts');
	return;
end;
if nargin < 2 || isempty(isF), isF = 0; end;

p = fileparts(which('PraatFmts'));
fp = fileparts(fName);
if isempty(fp), fName = fullfile(pwd,fName); end;

if isF, script = 'eff.praat'; else script = 'efm.praat'; end;
if ispc,
	[s,r] = dos(sprintf('"%s" "%s" "%s"', fullfile(p,'praatcon.exe'), fullfile(p,script), fName));
else,
	[s,r] = unix(sprintf('"%s" "%s" "%s"', fullfile(p,'Praat'), fullfile(p,script), fName));
end;

if s,
	error('error attempting to compute formants for %s', fName);
end;

q = regexp(r,'([0-9X.]+)\s','tokens');
v = zeros(1,length(q));
for k = 1 : length(q),
	v(k) = str2num(cell2mat(q{k}));
end;
v = reshape(v,[4 length(v)/4])';
step = min(diff(v(:,1))) * 1000;		% time step in ms
if isempty(step),
	fmts = NaN(1,4);
	return;
end;

% convert [time offset, value] pairs to equispaced vector
if verLessThan('matlab','8.3.0'),
	[s,sr] = wavread(fName);
else,
	[s,sr] = audioread(fName);
end
dur = 1000*length(s)/sr;
fmts = NaN * zeros(round(dur/step),3);
t = floor(1000*v(:,1)/step) + 1;
for k = 1 : length(t),
	fmts(t(k),:) = v(k,2:4);
end;
