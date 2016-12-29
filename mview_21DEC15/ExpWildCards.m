function ms = ExpWildCards(ws,s,matchCase)
%EXPWILDCARDS  - expand wildcards in string
%
%	usage:  ms = ExpWildCards(ws, s, matchCase)
%
% Use this procedure to test expanded wildcard string WS against cellstr S using regexp
% returns matched instances MS
%
% by default matching is regardless of case; specify nonzero MATCHCASE for case sensitivity
%
% "*" is mapped to "\w*" (0 or more instances of arbitrary alphanumeric character)
% "!" is mapped to "\w" (single instance of arbitrary alphanumeric character
% any other REGEXP recognized sequences matched as is

% mkt 07/15

if nargin < 2, eval('help ExpWildCards'); return; end
if nargin < 3 || isempty(matchCase), matchCase = 0; end

ws = regexprep(ws,'*','\\w+');			% expand "*"
ws = regexprep(ws,'!','\\w');			% expand "!"

ms = {};
for si = 1 : length(s),
	if matchCase,
		q = regexp(s{si},['^',ws,'$'],'match');
	else,
		q = regexpi(s{si},['^',ws,'$'],'match');
	end
	if ~isempty(q), ms{end+1} = q{1}; end
end
