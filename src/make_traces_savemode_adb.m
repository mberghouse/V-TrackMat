%%%make_traces function for UNIL
function [distMax]=make_traces_savemode_adb(pos,time,exp_param,opt,folder_save,xp_name,yp_name,tp_name,ds,no, N)
%%
tic
pos = unique(pos, 'rows', 'stable');
%%options for tracking
ghost = opt.ghost;                      %yes = 1, no = 0 for activating memory function
n_loop = opt.n_loop;                    %time loop to search candidates in the n*dx distances when memory function is on

%%define some vectors and function controls
%timeMax = 500;
timeMax = length(time)-1;

x = [];
y = [];
t = [];
pos1 = pos;
%%experimental parameters
Q = exp_param.Q;                        %uL/min = mm^3/min
w = exp_param.w;                        %mm
h = exp_param.h;                        %mm
n = exp_param.n;                        %mm^3/mm^3
pix_size = exp_param.pix_size;          %um
dt = 1;                      %s, inverse of frame rate

q = Q/(w*h)/60;                         %mm/s
v_mean = 10;                           %mm/s
dx_mean = v_mean*dt;                    %mm
%dx_mean = dx_mean;        %pixels
max_dx = dx_mean*1.0;                    %pixels we assume max displacement is 10x the average displacement

%%find out what the inter particle distance is at a given frame in time
ipd = [];
for ii=1:timeMax
    ind_c = find(pos(:,3)== ii);       %current given time
    tracking = ipdm([pos(ind_c,1),pos(ind_c,2)],'Subset','NearestNeighbor','result','struct');
    ipd = [ipd;tracking.distance];
end
% ipd(ipd==0)=NaN;
max_dist = max(ipd);
mean_dist = mean(ipd)
median_dist = median(ipd);
fprintf('mean ipd: ',mean_dist)

%%by default choose maximum distance (distMax) as the smallest of either median distance between particles in any frame or the maximum particle displacement
distMax = min(median_dist,max_dx)*4.2;
% return
% distMax = distMax*3;
%
%matObj = matfile('traces_temporar.mat','Writable',true);
%save('traject.mat','-struct','traces');
np = 0;
clear tracking
%% begin tracking code timeMax
pos = pos1;
opt.ghost=0;
clear ghost
ighost = [];
for ii = 1:timeMax
    
    fprintf('Analyzing frame %d over %d...\n',ii,timeMax);
    ind_c = find(pos(:,3)== ii);       %current time
    ind_n = find(pos(:,3)== ii+1);     %next time
    
    % first find all the particles in the first frame
    if ii == 1
        x(1:length(ind_c),ii) = pos(ind_c,1);
        y(1:length(ind_c),ii) = pos(ind_c,2);
        t(1:length(ind_c),ii) = time(pos(ind_c,3));
        % ghosts are 0 when there is a particle
        ighost(1:length(ind_c),ii) = 0;
        
        %padding for positions at next time (column)
        x(:,ii+1)=NaN;
        y(:,ii+1)=NaN;
        t(:,ii+1)=NaN;
        
        % make ghost as 1
        ighost(:,ii+1) = 1;
        dx = zeros(size(x));
        dy = zeros(size(x));
    end
    
    % save dead trajectories 
    if(mod(ii,20)==0 && ii>6)
%         return
        % death index= sum of how many NaN in the last 6 points of trajectory
        death = sum(isnan(x(:,end-5:end)),2);
        % find index of those trajectories
        ind_death = find(death > 5);
        
        fprintf('saving data for %d trajectories...\n',length(ind_death));
        % save those trajectories
%         pause
        np = save_temp_tracks_adb(x(ind_death,:),y(ind_death,:),t(ind_death,:),folder_save,xp_name,yp_name,tp_name,ind_death,ds,np,no, N, ii);
        
        % remove those trajectories from the main arrays
        x(ind_death,:) = [];
        y(ind_death,:) = [];
        t(ind_death,:) = [];
        ighost(ind_death,:) =[];
        
    end
    
    if ii > 1
        %padding for positions at next time (column)
        x(:,ii+1)=NaN;
        y(:,ii+1)=NaN;
        t(:,ii+1)=NaN;
        % padding ghost 
        ighost(:,ii+1) = 1;
        
        % find displacement, if NaN make it 0
        dx = gradient(x(:,ii-1:ii));dx=dx(:,2);dx(isnan(dx))=0;
        dy = gradient(y(:,ii-1:ii));dy=dy(:,2);dy(isnan(dy))=0;
    end
    
    cipdm=0;
    
    % search candidate particle in the next frame
    while sum(isnan(x(:,ii+1)))>0
        % find particles in the ii+1 frame that is NaN
        ind = find(isnan(x(:,ii+1)));
        
        %%find interparticle distance with correction for displacement in prior step
        tracking = ipdm([x(ind,ii)+dx(ind), y(ind,ii)+dy(ind)],[pos(ind_n,1),pos(ind_n,2)],'Subset','NearestNeighbor','result','struct');
        
        %%kick out redundant points that are not the minimum distance and remove them from the tracking structure entirely
        uniq_p = unique(tracking.columnindex);
        for r = 1:length(uniq_p)
            ind_red = find(tracking.columnindex == uniq_p(r));
            if length(ind_red)>1
                [~,I]=min(tracking.distance(ind_red));
                tracking.columnindex(ind_red(ismember(ind_red,ind_red(I))==0))=[];
                tracking.rowindex(ind_red(ismember(ind_red,ind_red(I))==0))=[];
                tracking.distance(ind_red(ismember(ind_red,ind_red(I))==0))=[];
            end
        end
        
        %%find those points in next frame that meet distance criterion and save their x and y positions
        j = find(tracking.distance <= distMax);
        
        % break when there is no particle to be updated
        if isempty(j)
            break
        end
        
        x(ind(tracking.rowindex(j)),ii+1) = pos(ind_n(tracking.columnindex(j)),1);
        y(ind(tracking.rowindex(j)),ii+1) = pos(ind_n(tracking.columnindex(j)),2);
        t(ind(tracking.rowindex(j)),ii+1) = time(pos(ind_n(tracking.columnindex(j)),3));
        
        % make ghosts 0 when particle is real
        ighost(ind(tracking.rowindex(j)),ii+1) = 0;
        
        %%clear out already matched points in pos at t+1 to prevent redundant assignment
        pos(ind_n(tracking.columnindex(j)),1:2)=NaN; % adb: changed from 1:2 to 1:3
        cipdm=cipdm+1;
        
    end
    
    % fill the ghost locations with projected displacement
    ign = find(ighost(:,ii+1)==1);
    x(ign, ii+1) = x(ign,ii)+dx(ign)*.4;
    y(ign, ii+1) = y(ign,ii)+dy(ign)*.4;
    t(ign, ii+1) = t(ign,ii)+dt;

    
    if ii>=3
        % indices of ghost particles in the previous frame
        ignn = find(ighost(:, ii-1)==1);
        % indices of those ghosts that found a real particle in next 2 frames
        ig2r = find(prod(ighost(ignn,ii-1:ii+1),2)==0);
        % remove those particles from ghost list
        ighost(ignn(ig2r), ii-1) = 0;
        % find indices of ghosts that didnt find real particle
        irnn = setxor(ignn,ignn(ig2r));
        % make those NaN
        x(irnn,ii-1:ii+1)=NaN;
        y(irnn,ii-1:ii+1)=NaN;
        t(irnn,ii-1:ii+1)=NaN;
    end
    
    %%start new trajectories for points too far from previous frame in x and y matrices as new rows
    k = find(~isnan(pos(ind_n,1))==1);
    if size(k)>0
        sprintf('making %d new trajs',length(k));
        x(end+1:end+length(k),:)=NaN;       %padding for new trajs (rows)
        y(end+1:end+length(k),:)=NaN;
        t(end+1:end+length(k),:)=NaN;
        ighost(end+1:end+length(k),:) = 1; %cellfun(@(x) add_ghost(x, length(k), 1), ghost, 'UniformOutput',0);
        
        x(end+1-length(k):end,ii+1) = pos(ind_n(k),1);
        y(end+1-length(k):end,ii+1) = pos(ind_n(k),2);
        t(end+1-length(k):end,ii+1) = time(pos(ind_n(k),3));
        ighost(end+1-length(k):end,ii+1) = 0;
        clear k
    end
    
end
% save survival tracks
ind_surv = 1:size(x(:,1:end),1);%len(x)
%size(x(:,1:end),1)
np_end = save_temp_tracks_adb(x(:,1:end),y(:,1:end),t(:,1:end),folder_save,xp_name,yp_name,tp_name,ind_surv,ds,np,no, N);

sprintf('... fraction of particles not matched and converted to new trajectories is %f',...
    sum(~isnan(pos(pos(:,3)<timeMax,1)))/length(pos(pos(:,3)<timeMax,1)))
