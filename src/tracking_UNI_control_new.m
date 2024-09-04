script;clc;
clear all;warning off;
fh = findall(0,'type','figure');
for indY=1:length(fh)
    clf(fh(indY));
end
pause(0);




%% I. load pos file containing a list of particle coordinates for each frame recorded
% tic
%%% Redo 500 16x, 1000 16x, 1000 4x, 500 32x, 
to_be_tracked = 1;
if(to_be_tracked)
    pos = readmatrix('C:\Users\marcb\Desktop\OF_PTsim\homo_spots2\tmate_poslist_1000part_32xspeed.csv');
    pos = sort(pos(2:end,5:8),8);
    pos(:,3)=[];
    %pos_t = load('C:\Users\marcb\Desktop\OF_PTsim\PT Simulation\trackpy_poslist_2000part.csv');
    %pos = load('C:\Users\marcb\Desktop\lbm_sim_600particles\pos_lst_600particles.dat');
    
end

S=1;
dt = 1;
N = max(pos(:,3))+1;

folder = 'C:\Users\marcb\Desktop\filippo_PT_tempfolder\';
dir = [folder, 'dataSampled_dt\'];
if isfolder(dir)==0
    mkdir(dir);
end

folder_save = dir;
time = (1:N);
%time = linspace(0,N*dt,N);
%time = time - time(1) + diff(time(1:2));


'... step I complete: pos'
% II. make traces
clc
clear traces

cd 'C:\Users\marcb\Desktop\Filippo_new\'

%%experimental parameteres
exp_param.Q = 0.7;                                                          %uL/min = mm^3/min, flow rate
exp_param.w = 7;                                                            %mm, width of cell
exp_param.h = 11;                                                        %mm, height of cell
exp_param.n = 0.46;                                                          %mm^3/mm^3, porosity
exp_param.pix_size = 1;                                                 %um, pixel size (px_size/magnification)
exp_param.dt = 1;                                                       %s, inverse of frame rate
nx = 512;
ny = 512;
st = 1000;
ds = linspace(nx/st,nx,st);


%%tracking code options
opt.ghost = 1;                                                              %yes = 1, no = 0 for activating memory function
opt.n_loop = 5;                                                             %time loop to search candidates in the n*dx distances when memory function is on
opt.plotTrajs = 1;                                                          %yes = 1, no = 0 for the option of plotting trajectories
opt.step_save = 0;

% save temp trajectories parameters
save_temp_track_opt = 1;
no = ['Filippo_', sprintf('%d',S)];


%save temporary trajects
if(save_temp_track_opt==1)
    param_name = ['param_',sprintf('%s',no),'.dat'];                          %open writtable .bin files
    xp_name = ['xp_',sprintf('%s',no),'.bin'];
    fileIDx = fopen([folder_save,xp_name],'w');fclose(fileIDx);
    yp_name = ['yp_',sprintf('%s',no),'.bin'];
    fileIDy = fopen([folder_save,yp_name],'w');fclose(fileIDy);
    tp_name = ['tp_',sprintf('%s',no),'.bin'];
    fileIDt = fopen([folder_save,tp_name],'w');fclose(fileIDt);
    
    
    %create and save trajectories in binary files
    [distMax]=make_traces_savemode_adb(pos,time,exp_param,opt,folder_save,xp_name,yp_name,tp_name,ds,no, N);
else
    % file names to be read
    param_name = ['param_',sprintf('%s',no),'.dat'];
    xp_name = ['xp_',sprintf('%s',no),'.bin'];
    yp_name = ['yp_',sprintf('%s',no),'.bin'];
    tp_name = ['tp_',sprintf('%s',no),'.bin'];
end
distMax

%figure(1);[~,distMax]=make_traces_savemode(pos,time,exp_param,opt,folder_save,xp_name,yp_name,tp_name,ds,no);
%figure(1);[traces_o,distMax]=make_traces(pos,exp_param,opt); %traces.coos{j}(:, [x y time])

'... step II complete: tracks saved in binary files'
%return
% pause
%% II.b construct a struct file from binary saved trajectories (if the binary file already exists)
read_binary_coos = 1;
% no = 02;

param = load([folder_save,'param_save_',sprintf('%s',no),'.dat']);
if(read_binary_coos)
    [traces_o] = read_binary_tracks_bis_new(folder_save,no,param);
end

'... step II.b complete: struct file traces_o created...'

%% III. kill too short and practically immobile traces and plot those we will keep
clc
clear traces
traces=traces_o;                                                           %loads traces from previous section
traces.opt.minLength=14;                                                  %minimum number of frames to keep a trajectory (this number should be at least the size of n in find_candidate_traces)
traces.opt.minDisplacement=.001;                                            %particle must be displaced at least this distance during its existence (here ~ 3 particle diameters)
traces_s=kill_short_traces_adb(traces);                                    %traces.coos{j}(:, [x y time])
% 
%figure(1); hold on; cellfun(@(x) plot(x(:,1), x(:,2), 'r'), traces_o.coos(traces_s.opt.InvalidIndex), 'UniformOutput',0); axis equal xy
%figure(1); hold on; cellfun(@(x) plot(x(:,1), x(:,2)), traces_s.coos, 'UniformOutput',0); axis equal xy

'... step III complete: traces_s'

%% IV. find ANACHRONISTIC candidates for gluing in same field of view & selects among multiple by shortest distance with matching tail-head position and velocity continuity
%distMax = ds(1);
clc
clear traces
traces=traces_s;                                                          %load traces from previous section (it needs tr.opt struct)
dx = cellfun(@(x) diff(x(:,1)), traces.coos, 'UniformOutput',0);
dx = cell2mat(dx');
dy = cellfun(@(x) diff(x(:,2)), traces.coos, 'UniformOutput',0);
dy = cell2mat(dy');

traces.opt.dxMax=sqrt(mean(dx)^2+mean(dy)^2);

c_loop = 1;
n_glue = numel(traces.coos);
while c_loop<2
    traces_c_fm = find_candidate_traces_2D2D_parallel_FM_2_adb(traces);
    if size(traces_c_fm.opt.pairsID,1)<n_glue
        n_glue = size(traces_c_fm.opt.pairsID,1);
    else
        break;
    end
    glue_trajectories = 1;
    if(glue_trajectories)
        traces_new=glue_all_traces_new_adb(traces_c_fm);
    end
    traces = traces_new;
    c_loop = c_loop+1
end
fig = figure(9); 
fig.Position =  [30 30 700 700];
hold on; cellfun(@(x) plot(x(:,2), x(:,1),'LineWidth',1), traces_new.coos, 'UniformOutput',0); axis equal xy
xlim([600,1000])
ylim([600,1000])
% figure(1); %clf;
% indexPair = traces_new.opt.pairsID(:,1);
% indexPair = indexPair(find(indexPair<numel(traces_new.coos)));
% % toc
% hold on; cellfun(@(x) plot(x(:,1), x(:,2)), traces_new.coos(indexPair), 'UniformOutput',0); axis equal xy
% %axis([1400 1540 1180 1280])
% hold off;
'... step IV complete: traces_c'
%%
% N = cellfun(@(x) x(end,1)-x(1,1), traces_s.coos, 'UniformOutput',0);
% N=cell2mat(N)/nx;
% 
% fprintf('Number of initial trajectories: %d ...\n', numel(N))
% 
% % Nn = cellfun(@(x) x(end,1)-x(1,1), traces_new.coos, 'UniformOutput',0);
% % Nn=cell2mat(Nn)/nx;
% 
% [pdfs,Xs,~,~] = PDF_function3(N,ones(size(N)),25,'lin',5);
% % [pdfn,Xn,~,~] = PDF_function3(Nn,ones(size(Nn)),25,'log',5);
% 
% figure(1);
% clf
% loglog(Xs,pdfs, 'k'); hold on
% % loglog(Xn,pdfn, 'b'); hold on
% 
% %
% P =0.25;
% x0=cellfun(@(x) x(1,1), traces_new.coos, 'UniformOutput',0);
% x0 = cell2mat(x0);
% t0=cellfun(@(x) x(1,3), traces_new.coos, 'UniformOutput',0);
% t0 = cell2mat(t0);
% 
% ind = find(x0 <= P*nx & t0<= 0.1*max(t0));
% traces = traces_new.coos(ind);
% Nn = cellfun(@(x) x(end,1)-x(1,1), traces, 'UniformOutput',0);
% Nn=cell2mat(Nn)/nx;
% [pdfn,Xn,~,~] = PDF_function3(Nn,ones(size(Nn)),25,'lin',5);
% %
% x0=cellfun(@(x) x(1,1), traces_s.coos, 'UniformOutput',0);
% x0 = cell2mat(x0);
% t0=cellfun(@(x) x(1,3), traces_s.coos, 'UniformOutput',0);
% t0 = cell2mat(t0);
% ind = find(x0 <= P*nx & t0<= 0.1*max(t0));
% traces = traces_s.coos(ind);
% Ns = cellfun(@(x) x(end,1)-x(1,1), traces, 'UniformOutput',0);
% Ns=cell2mat(Ns)/nx;
% [pdfs,Xs,~,~] = PDF_function3(Ns,ones(size(Ns)),25,'lin',5);
% 
% 
% figure(1);
% clf
% plot(Xs,pdfs, 'r'); hold on
% plot(Xn,pdfn, 'b'); hold on
% %%
% figure(2); clf
% histogram(Ns); hold on
% alpha 0.5
% histogram(Nn); hold on
% alpha 0.5
% %% IV.b find ANACHRONISTIC candidates for gluing BTW 2 DIFFERENT field of view & selects among multiple by shortest distance with matching tail-head position and velocity continuity
% glue_consequent_FOV = 0;
% clear traces
% traces=traces_s;
% if(glue_consequent_FOV)
%     %loads traces from previous section
%     traces.opt.dxMax=distMax/1;                                               %max distance to glue traces together
%     traces_c=find_candidate_traces_4D_parallel(traces);                       %traces.coos{j}(:, [x y time 0or1(if candidate)]) & creates traces.opt.pairsID (IDparent, IDcandidate, Eucledian dist tail-head)
%     
%     traces_new=glue_all_traces_nw_adb(traces_c);                                             %traces.coos{j}(:, [x y time 0or2(for nothingORglued)])
% end
% '... step IV complete: traces_new'
% 
% %% V. plot Candidates
% 
% clear traces
% traces=traces_c_fm;                                                            %matched candidate spaghetti
% figure(1);clf;
% count_candidates=0;
% for j=1:length(traces.coos)
%     p1 = plot(traces.coos{j}(1,1),traces.coos{j}(1,2),'k*','markersize',4,'DisplayName','Spaghetti start points');     %plots start of trace in black
%     hold on
%     p2 = plot(traces.coos{j}(:,1),traces.coos{j}(:,2),'b.-','markersize',2,'DisplayName','Spaghetti references');      %plots all traces in blue
%     if traces.coos{j}(:,4)==1
%         count_candidates=count_candidates+1;
%         p3 = plot(traces.coos{j}(:,1),traces.coos{j}(:,2),'r.-','markersize',2,'DisplayName','Spaghetti candidate');   %plots candidates traces.coos{j}(:,4)==1
%     end
% end
% axis equal;axis tight;
% hold off;
% 
% legend([p1 p2 p3]);
% title(sprintf('%d candidate traces to be glued',count_candidates));
% xlabel('x');ylabel('y');
% clear j p1 p2 p3 count_candidates
% 
% '... step V complete: plot of traces_c'

% %% VI. glue all
% 
% clear traces
% traces=traces_c_fm;
% %traces=traces_c;
% traces_new=glue_all_traces_new_adb(traces);                                             %traces.coos{j}(:, [x y time 0or2(for nothingORglued)])
% %traces_new = glue_all_traces_new_overlapped(traces);
% '... step IV complete: traces_new'
% 
% %% VI.b plot the glued tracks
% traces_new = traces_c_fm;
% show_longest_trajs = 1;                                                         %yes = 1, no = 0 to plot longest trajs as 0.5 max length in frames and displacement
% figure(2);[long,displacement]= plot_glued_trajs(traces_new,show_longest_trajs);   %plots new longer glued trajs and their histogram length
% 
% [N,X] = hist(displacement,100);
% figure(5)
% plot(X,N,'DisplayName','length in displacement(px)');
% xlabel 'trajectory length (px)'; ylabel 'count';
% legend SHOW;
% 
% '... step VI.b complete: traces_new plotted'
% pause
% %% VII.  save the traces_new.coos1 file for each position
% to_save = 0;
% if(to_save)
%     dir_save = folder_save;
%     
%     coos = traces_new.coos;
%     opt = traces_new.opt;
%     %opt = traces_s.opt;
%     %coos = traces_s.coos;
%     save([dir_save,'track_glued12.mat'],'coos');
%     save([dir_save,'opt_track_glued12.mat'],'opt');
%     %save([dir,'displacement_2_full_set.dat'],'displacement','-ASCII');
%     
%     '... step VII complete: saved in .mat file...'
% end
%% VIII. stitch fields of view together ... TBD

% see gluing_FOVS_def.m code

% Structure: traces.coos{i} contains trajectory coordinates (x, y, frame #)

all_speeds = []; % Initialize array to store all speeds

% Loop through each trajectory in traces.coos
for i = 1:length(traces.coos)
    current_traj = traces.coos{i};
    
    % Calculate speeds for current trajectory
    speeds = zeros(size(current_traj, 1) - 1, 1);
    for j = 1:size(current_traj, 1) - 1
        % Calculate displacement in x and y
        displacement = current_traj(j+1, 1:2) - current_traj(j, 1:2);
        distance = norm(displacement);
        
        % Calculate time difference using frame numbers
        dt = (current_traj(j+1, 3) - current_traj(j, 3)) * 0.1; % Assuming 10 fps (0.1s per frame)
        
        % Calculate speed
        speed = distance / dt;
        speeds(j) = speed;
    end
    
    % Append speeds to the all_speeds array
    all_speeds = [all_speeds; speeds];
end
all_speeds = all_speeds(all_speeds ~= 0);
% Save the array of all speeds



save('all_speeds001.mat', 'all_speeds');

xc = readmatrix('C:\Users\marcb\Desktop\OF_PTsim\homo_trajectories2\xc_homo_1000part_32xspeed.csv');
yc = readmatrix('C:\Users\marcb\Desktop\OF_PTsim\homo_trajectories2\yc_homo_1000part_32xspeed.csv');
total_ipdm=0;
xc = xc(1:50,:);
yc = yc(1:50,:);
for i=1:size(xc,1)
   x=xc(i,:)';
   y=yc(i,:)';
   x = rmmissing(x);
   y = rmmissing(y);
   d=ipdm([x,y],'Subset','NearestNeighbor','result','struct');
   idist = d.distance;
   idist = rmmissing(idist);
   mean_ipdm=mean(idist);
   total_ipdm = total_ipdm+mean_ipdm;
end
total_ipdm/size(xc,1)

