function WritePraatTextGrid(fName, tiers, dur)
%WRITEPRAATTEXTGRID  - write label data to PRAAT short TextGrid format file
%
%	usage:  WriteShortTextGrid(fName, tiers, dur)
%
% FNAME is the filename to write to (optional extension defaults to ".TextGrid")
%
% TIERS is an array of structs, one per desired PRAAT tier, with fields
%   NAME    - tier name
%   OFFSETS - label offsets (secs); either [nLabels x HEAD,TAIL] (interval) or [nLabels x 1] (point)
%	LABELS  - label strings {nLabels x 1}
%
% DURation is the labeled signal data duration in seconds
%
% see also READPRAATTIER

% mkt 11/14

% parse params
if nargin < 3, eval('help WritePraatTextGrid'); return; end
nTiers = length(tiers);

% open file
[p,f,e] = fileparts(fName);
if isempty(e), fName = fullfile(p,[f,'.TextGrid']); end
try,
	fid = fopen(fName, 'wt');
catch,
	error('error attempting to open %s', fName);
end;

% write header
fprintf(fid,'File type = "ooTextFile"\n');
fprintf(fid,'Object class = "TextGrid"\n');
fprintf(fid,'\n');
fprintf(fid,'0\n');
fprintf(fid,'%f\n', dur);
fprintf(fid,'<exists>\n');
fprintf(fid,'%d\n', nTiers);

% write tiers
for ti = 1 : nTiers,
	t = tiers(ti);
	nLabels = size(t.OFFSETS,1);
	if size(t.OFFSETS,2) == 1,		% point tier
		fprintf(fid,'"TextTier"\n');
		fprintf(fid,'"%s"\n', t.NAME);
		fprintf(fid,'0\n');
		fprintf(fid,'%f\n', dur);
		fprintf(fid,'%d\n', nLabels);
		for k = 1 : nLabels,
			fprintf(fid,'%f\n', t.OFFSETS(k));
			fprintf(fid,'"%s"\n', t.LABELS{k});
		end
	else,							% interval tier
		fprintf(fid,'"IntervalTier"\n');
		fprintf(fid,'"%s"\n', t.NAME);
		fprintf(fid,'0\n');
		fprintf(fid,'%f\n', dur);
		fprintf(fid,'%d\n', nLabels);
		for k = 1 : nLabels,
			fprintf(fid,'%f\n%f\n', t.OFFSETS(k,:));
			fprintf(fid,'"%s"\n', t.LABELS{k});
		end
	end	
end

% clean up
fclose(fid);
fprintf('wrote %s\n', fName);

