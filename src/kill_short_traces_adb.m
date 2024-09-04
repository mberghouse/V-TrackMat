function traces=kill_short_traces_adb(traces)

minL = traces.opt.minLength;                 
disMin = traces.opt.minDisplacement;  

dx = cellfun(@(x) abs(x(1,1)-x(end,1)), traces.coos, 'UniformOutput',0);
dy = cellfun(@(x) abs(x(1,2)-x(end,2)), traces.coos, 'UniformOutput',0);
dist = cellfun(@(x,y) sqrt(x.^2+y.^2), dx,dy, 'UniformOutput',0);

N = cellfun(@(x) size(x,1), traces.coos, 'UniformOutput',0);
N = cell2mat(N);

dist = cell2mat(dist);
ind = find(dist>=disMin & N>=minL);

traces.coos = traces.coos(ind);
traces.opt.validIndex = ind;
traces.opt.InvalidIndex = setxor(1:numel(traces.coos),ind);