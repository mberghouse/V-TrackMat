%%%find_candidate_traces_4d function for UNIL
function traces=find_candidate_traces_2D2D_parallel_FM_2_adb(traces)
%% NB1: traces.coos{} contains 4 columns, [x y time ID=1 for possible gluable trace candidate or ID=0 otherwise]

DistMax = traces.opt.dxMax;           %make this criterion 10x more strict since we are glueing anachronistic trajs
n = 1;                                %how far back and forward in frames we look at average velocity components to match
if traces.opt.minLength-1 < n
    warning('shortest trajectories are too short to look n=5 steps forward and back for matching velocities...')
end

%%build velocity components
traces.velacc=[];
nrTraces=length(traces.coos);
glue_matrix = zeros(nrTraces,6);

% for each trajectory
for j=1:nrTraces

    traces.coos{j}(:,4) = 0;
    % velocity in x and y
    traces.velacc{j}(:,1) = gradient(traces.coos{j}(:,1))./gradient(traces.coos{j}(:,3)); %u
    traces.velacc{j}(:,2) = gradient(traces.coos{j}(:,2))./gradient(traces.coos{j}(:,3)); %v
    traces.dt{j} = gradient(traces.coos{j}(:,3)); % dt
    
    glue_matrix(j,1)=traces.coos{j}(1,1); % first point in traj (x)
    glue_matrix(j,2)=traces.coos{j}(1,2); % first point in traj (y)
    glue_matrix(j,3)=mean(traces.velacc{j}(1:n,1)); % mean velocity x first n points
    glue_matrix(j,4)=mean(traces.velacc{j}(1:n,2)); % mean velocity y first n points
    glue_matrix(j,5)= abs(traces.coos{j}(2,3) - traces.coos{j}(1,3)); % first dt 
    glue_matrix(j,6)= 0.5*(std(traces.velacc{j}(1:n,1)) + std(traces.velacc{j}(1:n,2))); % std in velocity first n
    
end
'...created a matrix... to glue them'
traces.opt.pairsID=[];          %gets filled with base traj ID, candidate traj ID, 2D Eucledian distance tail-end
SkipInd = [];
tic
%DistMax = 8;
%
clc
for jt1= 1:nrTraces
%     sprintf('... analyzing %d out of %d candidates' ,jt1,nrTraces)
    
    Ind =1:nrTraces;
    % SkipInd = [jt1,SkipInd];
    % Ind(SkipInd) = [];
    % Ind = Ind_t(Ind_t~=jt1);
    %search for candidates if second traj is not its own pair and don't
    %worry about time since we are glueing anachronistic trajectories together
    
    
    %%============modified from 6D distance from Xu et al 2008 we use in 3D-PTV
    xie = traces.coos{jt1}(end,1);  %takes last position x in trajectory jt1
    yie = traces.coos{jt1}(end,2);  %takes last position y in trajectory jt1
    
    uie = mean(traces.velacc{jt1}(end-n:end,1));                   %takes average vel component of last n frames
    vie = mean(traces.velacc{jt1}(end-n:end,2));
    sigma_v_cand = 0.5*(std(traces.velacc{jt1}(end-n:end,1)) + std(traces.velacc{jt1}(end-n:end,2)));
    dte = traces.dt{jt1}(end);
    
    xjs = glue_matrix(Ind,1);                                      %takes first non-NaN position in trajectory jt2
    yjs = glue_matrix(Ind,2);
    ujs = glue_matrix(Ind,3);                       %takes average vel component of first n frames
    vjs = glue_matrix(Ind,4);
    
    
    DistMax2 = 2*n*sqrt(uie^2+vie^2) * dte; % you can use DistMax2 for major axis
    
    %     return
    %     direction_ie = atand(vie./uie);
    %     traj_pt = [xie'; yie'];
    %     [ind_temp,~] = arrayfun(@(p) directionaldist(p, [xjs';yjs'],DistMax,1,5,0.5,'none'), traj_pt, 'UniformOutput',0);
    %
    %     numel(ind_temp)
    %
    %     pause;
    %dtjs = glue_martix(Ind,5);
    %          dist = sqrt((abs(xie-xjs)).^2 + (abs(uie-ujs)).^2 + ...
    %             (abs(yie-yjs)).^2 + (abs(vie-vjs)).^2);
    %     dist_space = sqrt((abs(xie-xjs)).^2 + (abs(yie-yjs)).^2);
    %     dist_space(jt1) = 3e3;
    %     ind_temp = find(dist_space <= DistMax);
    %     ind_cand = Ind(ind_temp)
%     figure(1); clf
%     quiver(xie, yie, uie, vie, 0.01); hold on
%     [DistMax DistMax2]
    ind_cand = find_particle_inside(xie, yie, xjs, yjs, atan2(vie,uie), [DistMax DistMax2]);
    
    
%     [uie vie DistMax2 DistMax]
%         return
    %%
    if(~isempty(ind_cand))
        dist_vel = sqrt((abs(uie-ujs(ind_cand))).^2 + (abs(vie-vjs(ind_cand))).^2);
        [~,ind_dv] = min(dist_vel);
        real_ind = ind_cand(ind_dv);
        
        %%select and keep track if: a reasonable match is found that meets maximum dist criterion
        %if  dist<DistMax
        %if((dist-DistMax) <= 0)
        
        if(dist_vel(ind_dv) <= (sigma_v_cand + glue_matrix(real_ind,6)))
            %             disp('real ind')
            %             jt1
            %%
            jt2_ind = real_ind;
            dist_pairs = dist_vel(ind_dv);
            traces.opt.pairsID = [traces.opt.pairsID;jt1,jt2_ind,dist_pairs,length(traces.coos{jt1})];
            traces.coos{jt2_ind}(:,4)=1;
            traces.coos{jt1}(:,4)=1;
            %   disp('here3')
            %    pause
            %     SkipInd = jt2_ind;
            
            %                plot to check
            %                             figure(3)
            %                             plot(traces.coos{jt1}(:,1),traces.coos{jt1}(:,2),'.-b')
            %                             hold on
            %                             plot(traces.coos{jt2_ind}(:,1),traces.coos{jt2_ind}(:,2),'.-r')
            % %                             hold off
            %                             axis equal
            %                             pause
            %             return
        end
        
    end
end
% return
%%put extra pairs away (we don't want one trace to be glued to several other ones, so we just take the ones with the shortes distance)
[~,ia,ic]=unique(traces.opt.pairsID(:,1));
uniqueID=traces.opt.pairsID(ia,:);
if size(traces.opt.pairsID(:,1))~=size(uniqueID(:,1))
    warning('Check size of pairsID bc they dont match size of uniqueID!!!');
end

sprintf('... we found %d pairs', size(traces.opt.pairsID,1))
% pause
for j=1:size(uniqueID,1)
    ind=find(traces.opt.pairsID(:,1)==uniqueID(j,1));
    tmp=traces.opt.pairsID(ind,:);
    [~,imin]=min(tmp(:,3));
    uniqueID(j,:)=tmp(imin,:);
end

fprintf('time: %1.1f min \n',toc/60);

% just to be safe
traces.opt.pairsID=uniqueID;

traces.opt.pairsID_o=traces.opt.pairsID;