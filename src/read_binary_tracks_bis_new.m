%%Read binary saved tracks and create a "struct" file
function [traces] = read_binary_tracks_bis_new(folder,no,param)

%clear all;clc;close all;

%folder_save = '/media/fmiele/USB DISK/lagrangian_tracks/';
%folder = '/media/fmiele/Elements/lagr_tra_0808/pos_2_tracks_ds_600/';
%folder = '/media/fmiele/Elements/lagr_tra_0808/mask_traject/DATA_pos2_fullset_lower_conc/';
%folder_save = folder;

np = param(2);
nt = param(1);
 
fileID = fopen([folder,'xp_',sprintf('%s',no),'.bin']);
xp = fread(fileID,[np,nt],'double','s');
fclose(fileID);
fileID = fopen([folder,'yp_',sprintf('%s',no),'.bin']);
yp = fread(fileID,[np,nt],'double','s');
fclose(fileID);
fileID = fopen([folder,'tp_',sprintf('%s',no),'.bin']);
tp = fread(fileID,[np,nt],'double','s');
fclose(fileID);

disp('file loaded')

%% save trajectories as objects, if desired plot constructed trajectories and PDF of their length
long = [];
displacement = [];
traces.coos = {};    %coordinates (x y time)

for ii = 1:size(xp,2)
    sprintf('... restructuring trajs as objects %d of %d',ii,size(xp,2))
    
    ind = find(~isnan(xp(:,ii)));
    long = [long;sum(~isnan(xp(:,ii)))];
    traces.coos{ii} = [xp(ind,ii) yp(ind,ii) tp(ind,ii)]; %% attenzione...
    displacement = [displacement; ipdm(traces.coos{ii}(1,1:2),traces.coos{ii}(end,1:2))];
    
    clear xp(ind,ii) yp(ind,ii) tp(ind,ii) 
    
end

disp('ended')
%clear xp yp tp
% %%
% clc
% 
% vv = find(sum(isnan(xp)) < 100);
% fprintf('tot. long tracks %d \n',size(vv));
% 
% for k=1:size(vv,2)
% %     figure(1)
% %     clf
% %     plot(tp(:,vv(k)),xp(:,vv(k)),'.-')
% %     %hold on
% %     drawnow
% %     figure(2) 
% %     clf
% %     plot(tp(:,vv(k)),yp(:,vv(k)),'.-')
% %     %hold on
% %     drawnow
%     figure(3)
%    % clf;
%     plot(xp(:,vv(k)),yp(:,vv(k)),'.-')
%     axis equal
%     hold on
%     drawnow
% end