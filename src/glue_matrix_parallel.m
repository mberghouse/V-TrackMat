function [glue_matrix,traces] = glue_matrix_parallel(traces)

DistMax = traces.opt.dxMax;           %make this criterion 10x more strict since we are glueing anachronistic trajs
n = 5;                                %how far back and forward in frames we look at average velocity components to match
if traces.opt.minLength-1 < n
    warning('shortest trajectories are too short to look n=5 steps forward and back for matchig velocities...')
end

%%build velocity components
traces.velacc=[];
nrTraces=length(traces.coos);
glue_matrix = zeros(nrTraces,6);

for jn=1:nrTraces
    
    if(mod(jn,1000) == 0)
        sprintf('created %d out of %d total',jn,nrTraces)
    end
    
    traces.coos{jn}(:,4) = 0;
    traces.velacc{jn}(:,1) = smooth(gradient(traces.coos{jn}(1:2*n,1))./gradient(traces.coos{jn}(1:2*n,3)),0.18,'rloess'); %u
    traces.velacc{jn}(:,2) = smooth(gradient(traces.coos{jn}(1:2*n,2))./gradient(traces.coos{jn}(1:2*n,3)),0.18,'rloess'); %v
    
    glue_matrix(jn,1)=traces.coos{jn}(1,1);
    glue_matrix(jn,2)=traces.coos{jn}(1,2);
    glue_matrix(jn,3)=mean(traces.velacc{jn}(1:n,1));  
    glue_matrix(jn,4)=mean(traces.velacc{jn}(1:n,2));  
    glue_matrix(jn,5)= abs(traces.coos{jn}(2,3) - traces.coos{jn}(1,3));
    glue_matrix(jn,6)= 0.5*(std(traces.velacc{jn}(1:n,1)) + std(traces.velacc{jn}(1:n,2)));
    
end
'...created a matrix... to glue them'

end