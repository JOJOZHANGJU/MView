function mat2wav(mask, scale, idx)
%MAT2WAV  - save audio track from array-of-structs MAT variable to WAV format file
%
%	usage:  mat2wav(mask, scale, idx)
%
% Use this procedure to save the audio track from an MVIEW-compatible array-of-structs
% variable stored in a MAT file to a separate MS WAVE format file with the same name.
%
% All files matching string MASK are converted.
%
% All files are scaled to +/- SCALE range.
%
% Optional IDX (into MAT variable) defaults to 1
%
% Example:  convert files matching JW11* to a +/-5 range
% >> mat2wav('JW11*',5)

% mkt 07/06

if nargin < 2,
	eval('help mat2wav');
	return;
end;
if nargin<3 | isempty(idx), idx = 1; end;

fileList = [];
files = dir(mask);
fileList = {files.name}';
[p,f,ext] = fileparts(mask);
fl = fileList;
for k = 1 : length(fl),
	[p,f] = fileparts(fl{k});
	fl{k} = f;
end;

for k = 1 : length(fl),
	try,
		data = load(fl{k},fl{k});
		data = eval(['data.',fl{k}]);
	catch,
		error(sprintf('unable to load array-of-structs data from %s.mat', fl{k}));
	end;
	s = data(idx).SIGNAL ./ scale;
	try,
		wavwrite(s,data(idx).SRATE,16,[fl{k},'.wav']);
	catch,
		error(sprintf('unable to write %s.wav', fl{k}));
	end;
	fprintf('wrote %s.wav\n', fl{k});
end;
