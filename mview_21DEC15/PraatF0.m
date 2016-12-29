function f0 = PraatF0(fName)
%PRAATF0  -- compute F0 contour using Praat
%	
%	usage:  f0 = PraatF0(fName)
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
% for all platforms the script ep.praat must also be present in the same location as this file (PraatF0)

% mkt 08/08

if nargin < 1,
	eval('help PraatF0');
	return;
end;

p = fileparts(which('PraatF0'));
fp = fileparts(fName);
if isempty(fp), fName = fullfile(pwd,fName); end;

if ispc,
	[s,r] = dos(sprintf('"%s" "%s" "%s"', fullfile(p,'praatcon.exe'), fullfile(p,'ep.praat'), fName));
else,
	[s,r] = unix(sprintf('"%s" "%s" "%s"', fullfile(p,'Praat'), fullfile(p,'ep.praat'), fName));
end;

if s,
	error('error attempting to compute pitch for %s', fName);
end;

q = regexp(r,'([0-9X.]+)\s','tokens');
if length(q) < 4, 
	error('error attempting to compute pitch for %s', fName);
end;
v = zeros(1,length(q));
for k = 1 : length(q),
	v(k) = str2num(cell2mat(q{k}));
end;
v = reshape(v,[2 length(v)/2])';
step = min(diff(v(:,1))) * 1000;		% time step in ms

% convert [time offset, value] pairs to equispaced vector
[s,sr] = wavread(fName);
dur = 1000*length(s)/sr;
f0 = NaN * zeros(round(dur/step),1);
t = floor(1000*v(:,1)/step) + 1;
for k = 1 : length(t),
	f0(t(k)) = v(k,2);
end;

