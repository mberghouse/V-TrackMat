function newtraces=glue_all_traces_new_adb(traces)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that calls the glue_function_adb iteratively
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrTraces=length(traces.coos); % number of original trajectories
nrNTraces=0; % number of glued trajectories
skipInd=[];

listPairID = traces.opt.pairsID_o;

for jg = 1:size(listPairID,1)
    if ~ismember(jg,skipInd) % skipInd is list of indices that has been glued
        sprintf('... glueing %d out of %d pairs' ,jg,length(traces.opt.pairsID(:,1)))
        thereIsNext=1; 
        % continue running the while loop as long as its 1
        nrNTraces=nrNTraces+1;
        % identify the pair for jg
        pair=listPairID(jg,1:2);
        % assign first pair
        coos1 = traces.coos{pair(1)}; % coos of trace1
        
        while thereIsNext==1
            coos2 = traces.coos{pair(2)}; % coos of trace2
            % update trace1 with glued trace
            coos1 = glue_function_adb(coos1,coos2, 0);
            
            % find next index in the first column of the list 
            nextInd=find(pair(2)==listPairID(:,1));
            % make sure this index is unique
            [aa,~,~] = unique(skipInd);
            if (~isempty(nextInd) && (length(aa)==length(skipInd)))
                skipInd=[skipInd; nextInd];
                pair_temp=listPairID(nextInd, :);
                % only update the second index in pair
                pair(2)=pair_temp(2);
                thereIsNext=1;
            else
                thereIsNext=0;
            end
        end
        newtraces.coos{nrNTraces} = coos1;
    end
end


%%now add the ones that were not candidates to the new list of longer trajectories
ind=setxor([traces.opt.pairsID(:,1);traces.opt.pairsID(:,2)],1:nrTraces);
newtraces.coos=[newtraces.coos, traces.coos(ind)];
newtraces.opt=traces.opt;

