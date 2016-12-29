function KillAll
%KILLALL  - deletes all open windows regardless of protection

% mkt 10/99

show = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles', 'on');
h = findobj(0);							
h = h(find(ishandle(h)));			% find valid handles
h = h(strmatch('figure', get(h,'Type')));
delete(h);							% delete figures
set(0,'ShowHiddenHandles', show);
