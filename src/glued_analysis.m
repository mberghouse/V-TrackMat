%% Stitching two images together by maching the mask

nn = 600;
cm2 = hot(nn);

buildingScene = imageDatastore('/media/fmiele/Elements/lagr_tra_0808/mask_traject/','FileExtensions',{'.png'});
a = readimage(buildingScene, 1);
b = readimage(buildingScene, 2);

d = readimage(buildingScene, 4);

e = imread([dir,'mask_123.png']);
c = imread([dir,'mask_12.png']);

f = imread([dir,'mask_pos5.png']);

figure(11)
clf;
imagesc(d)
colormap(cm2)
colorbar
set(gca,'Clim',[-0.3/nn,0.3])
axis equal
hold on

figure(7)
clf;
imagesc(a)
colormap(cm2)
colorbar
set(gca,'Clim',[-0.3/nn,0.3])
axis equal
hold on

figure(8)
clf;
imagesc(b)
colormap(cm2)
colorbar
set(gca,'Clim',[-0.3/nn,0.3])
axis equal
hold on

 if(1)
  a_shift = imtranslate(a,[33,24]);  %shift pos1 to match pos2
 end
 
c = [a_shift,b(:,2:end)];  %create mask pos12

figure(9)
imagesc(c)
colormap(cm2)
colorbar
set(gca,'Clim',[-0.3/nn,0.3])
axis equal
hold on

figure(10)
imagesc(c)
colormap(cm2)
colorbar
set(gca,'Clim',[-0.3/nn,0.3])
axis equal
hold on

c_shift = imtranslate(a,[33,24]);   % shift pos 12 to match pos3
imtranslate(a,[33,24]);          %create mask pos123

figure(12)
clf;
imagesc(e)
colormap(cm2)
colorbar
set(gca,'Clim',[-0.3/nn,0.3])
axis equal
hold on


%g = [e,f(:,2:end)];
g_shift = imtranslate(e,[36,23]);
g = [g_shift,f(:,2:end)];

figure(14)
clf;
imagesc(g)
colormap(cm2)
colorbar
set(gca,'Clim',[-0.3/nn,0.3])
axis equal
hold on

%%
%pair=traces.opt.pairsID(1,1:2); 

% figure(50)
% clf;
% plot(coos1(:,1),coos1(:,2),'rx-','DisplayName','traj 1');
% hold on
% plot(coos2(:,1),coos2(:,2),'bx-','DisplayName','traj 2');
% hold on
% plot(intCoos(:,1),intCoos(:,2),'gx-','DisplayName','missing bit');

dir = '/media/fmiele/Elements/lagr_tra_0808/pos_2_tracks_ds_600/';
X1 = load([dir,'pos_2_full_set.mat']);
mask = imread([dir,'mask_pos2.png']);

%X2 = load([dir 'new_traces_pos2.mat']);
Dx_shift = 34;
%Dy_shift = 22;

aa = 0;
figure(7)
clf;
imagesc(mask)
colormap(cm2)
colorbar
set(gca,'Clim',[-0.3/nn,0.3])
axis equal
hold on

for ii=1:length(X1.coos)

 x1 = X1.coos{ii}(:,:); 
 
 if(length(x1) > 450)
  %figure(7)
  %plot(x1(:,1),x1(:,2),'.-b')
  %hold on
  %drawnow
  aa = aa +1;
 end
end 
disp('done')
%% gluing criteria

x_max = length_x(a) + dx_shift; %to add for each fov

%find possible candidates

%% check on traject
  
Xlim = [];
figure(20)
clf;
for jj=1:length(traces_s.coos)
    
figure(20)
plot(traces_s.coos{jj}(:,1),traces_s.coos{jj}(:,2),'.-')

%x_lim = traces_new.coos{jjj}(end,1);
%y_lim = traces_new.coos{jjj}(end,2);

%Xlim =[Xlim,x_lim];
hold on
axis equal
drawnow
%pause(0.3)
end

%% saving new pos each time
%coos1 = traces_123.coos;
%save([dir,'traces_123_new.mat'],'coos1');

%%

% images posA and posB are stuck to create pos AB. A precedes B. B must be shifted for the
% total x-length of A (Dx). A must be shifted to match the edges with B (dx,dy).

traces_1 = load([dir 'traces_23_new.mat']);
traces_2 = load([dir 'traces_4_new.mat']);

dx = 33; %shift for A
dy = 24;
Dx_shift = 34; %shift for B
Dx_shift = length(c);

for mm=1:length(traces_1.coos1)
    traces_1.coos1{mm}(:,1) = traces_1.coos1{mm}(:,1) + 33;
    traces_1.coos1{mm}(:,2) = traces_1.coos1{mm}(:,2) + 24;
end


for nn =1:length(traces_2.coos1)
     traces_2.coos1{nn}(:,1) = traces_2.coos1{nn}(:,1) +Dx_shift;
end

%% gluing 123->4 (pos 2, 3,4 to 5)

traces_1234.coos = {};
%traces_123.coos = [traces_1.coos1,traces_2.coos1];
traces_13 = load([dir,'traces_123_def.mat']);
traces_5 = load([dir,'traces_new_5.mat']);

h = imread([dir,'mask_1234.png']);
gg = imread([dir,'mask_123.png']);

Dx = length(gg);

figure(14)
clf;
nn = 600;
cm2 = hot(nn);
imagesc(h)
colormap(cm2)
colorbar
set(gca,'Clim',[-0.3/nn,0.3])
axis equal
hold on

%shifting traces_123 dx,dy
for mm = 1:length(traces_13.coos1)
    
    traces_13.coos1{mm}(:,1) = traces_13.coos1{mm}(:,1) + 36;
    traces_13.coos1{mm}(:,2) = traces_13.coos1{mm}(:,2) + 23;
    
end

%shifting traces_5
for nn=1:length(traces_5.coos5)
    traces_5.coos5{nn}(:,1) =  traces_5.coos5{nn}(:,1) + Dx;
end

% creating traces_1234
traces_1234.coos = [traces_13.coos1,traces_5.coos5];

%plotting
for ii=1:length(traces.coos)

  %x1 = traces_new.coos{ii}(:,:);
   x1 = traces.coos{ii}(:,:);
  
  figure(14)
  plot(x1(:,1),x1(:,2),'.--')
  %hold on
  pause
 
end 
disp('done')

%plotting the new

for ii=1:20:length(traces_new.coos)

  x1 = traces_new.coos{ii}(:,:);
   %x1 = traces_1234.coos{ii}(:,:);
  
  figure(14)
  plot(x1(:,1),x1(:,2),'.--')
  hold on
  pause(0.2)
 
end 
disp('done')


%imwrite(g,[dir,'mask_1234.png'],'png','BitDepth',8);
%%
Xlim = [];
for jjj=1:length(traces_1.coos)
    
%figure(20)
%plot(traces_new.coos{jj}(:,1),traces_new.coos{jj}(:,2),'.r')

x_lim = traces_1.coos{jjj}(end,1);
y_lim = traces_1.coos{jjj}(end,2);

Xlim =[Xlim,x_lim];
%hold on
%axis equal
%pause(0.3)
end
%% printing figure

figure(14)
axis([0,8200,0,2000])

set(gca,'TickLabelInterpreter','latex')
%set(gca,'fontsize',6)
set(gcf,'papersize',[40,10],'PaperPosition',[0,0,40,10]);
print('-dpdf','-r1200',[dir,'glued1234_bis.pdf']);

ff =  imread([dir,'glued1234.png']);

figure(21)
imagesc(ff)

%%

traces = traces_new;
long = [];
displacement = [];

clf;

for ii = 1:length(traces.coos)
    ii/length(traces.coos)
    
    long = [long;length(traces.coos{ii}(:,1))];    
    displacement = [displacement; ipdm(traces.coos{ii}(1,1:2),traces.coos{ii}(end,1:2))];
    
    %%plot of all final trajectories 
    subplot(1,2,1);hold on;
    plot(traces.coos{ii}(:,1),traces.coos{ii}(:,2),'.-');
    xlabel 'x [pix]'; ylabel 'y [pix]';
    title 'New longer glued trajectories';
    axis tight;axis equal;
    
end

%rand
%pause
%%PDF of trajectory length (frames)
%subplot(1,2,2);hold on;
%[N,X] = hist(long);plot(X,N,'DisplayName','length in frames'); % linear plot
%[N,X] = hist(displacement);plot(X,N,'DisplayName','length in displacement');
[N,X] = hist(long);
loglog(X,N,'DisplayName','length in frames'); % log-log plot

[N,X] = hist(displacement);
figure(5)
loglog(X,N,'DisplayName','length in displacement(px)');
xlabel 'trajectory length (px)'; ylabel 'count';
legend SHOW;


%% save the displacement pdf

%reshape traces_new
clear traces_bis
traces_bis.coos = {};

for ll = 1:length(traces_new.coos)
    traces_bis.coos{ll}(:,1) = traces_new.coos{ll}(:,1);
    traces_bis.coos{ll}(:,2) = traces_new.coos{ll}(:,2);
    traces_bis.coos{ll}(:,3) = traces_new.coos{ll}(:,3);
end

disp_3 = displacement;
disp_1 = load([dir,'pos_4_displacement_round1.dat']);
disp_2 = load([dir,'pos_4_displacement_round2.dat']);

[N2,X2] = hist(disp_2);
[N1,X1] = hist(disp_1);
[N3,X3] = hist(disp_3);

figure(5)
clf;
semilogy(X1,N1,'DisplayName','length in displacement(px) round1');
hold on
semilogy(X2,N2,'DisplayName','length in displacement(px) round2');
hold on
semilogy(X3,N3,'DisplayName','length in displacement(px) round3');
xlabel 'trajectory length (px)'; ylabel 'count';
legend SHOW;

%%

tracks = load('/media/fmiele/Elements/lagr_tra_0808/mask_traject/traces_2_full_set.mat');
figure(3)
%clf;

for jj= 1:tracks(end,4)
    
   ind = find(tracks(:,4) == jj);
   figure(1)
   plot(tracks(ind,1),tracks(ind,2))
   hold on
   axis equal
   
   pause(0.3)
    
end

%disp = load('/media/fmiele/Elements/exp_pt_lower_conc/');

figure(14)
clf;
for ii=1:50:length(tracks.coos)

  x1 = tracks.coos{ii}(:,:);
  %x1 = traces_1234.coos{ii}(:,:);
  
  figure(14)
  plot(x1(:,1),x1(:,2),'.--')
  hold on
  pause(0.2)
  axis equal
 
end 
disp('done')
%%

if(save_param)
  param_name = ['param_',sprintf('%d',no),'.dat'];
  xp_name = ['xp_',sprintf('%d',no),'.bin'];
  fileIDx = fopen([folder_save,xp_name],'w');fclose(fileIDx);
  yp_name = ['yp_',sprintf('%d',no),'.bin'];
  fileIDy = fopen([folder_save,yp_name],'w');fclose(fileIDy);
  tp_name = ['tp_',sprintf('%d',no),'.bin'];
  fileIDt = fopen([folder_save,tp_name],'w');fclose(fileIDt);
end

p0 = 1; 
step_save = 5e2;

%x(ii,ind)' y(ii,ind)' t(ind)

if(mod(p,step_save) == 0)
    for ii = 1:step_save
      if(save_param)
        param = [ii,length(x(ii,ind))];
        save([folder_save,param_name],'param','-ascii');
        fileIDx = fopen([folder_save,xp_name],'a');
        fwrite(fileIDx,x(ii,ind),'double','s');
        fclose(fileIDx);
        fileIDy = fopen([folder_save,yp_name],'a');
        fwrite(fileIDy,y(ii,ind),'double','s');
        fclose(fileIDy);
        fileIDt = fopen([folder_save,tp_name],'a');
        fwrite(fileIDt,t(ii,ind),'double','s');
        fclose(fileIDt);
      end
    end
    time_elps = toc;
    fprintf('\n * Particle %d tracked over %d - time elapsed: %6.1f s... \n',p,Np,time_elps)
    toc
end

%% Read binary saved tracks

no = 2;
folder_save = '/media/fmiele/Elements/lagr_tra_0808/mask_traject/DATA_pos2_fullset/';
param = load([folder_save,'param_save_',sprintf('%02d',no),'.dat']);

np = param(1);
nt = param(2);

fileID = fopen([folder_save,'xp_',sprintf('%d',no),'.bin']);
xp = fread(fileID,[np,nt],'double','s');
fclose(fileID);
fileID = fopen([folder_save,'yp_',sprintf('%d',no),'.bin']);
yp = fread(fileID,[np,nt],'double','s');
fclose(fileID);
 
figure(4)
clf;
for k=1:10:6000
 figure(4)
 plot(xp(:,k),yp(:,k),'.--')
 axis equal
 hold on
 pause(0.2)
end

tracks = load('/media/fmiele/Elements/exp_pt_lower_conc/traces_test.mat');

ds_cum_1 = [];
ds_cum_2 = [];

for ii=1:length(tracks.coos)

  x1 = tracks.coos{ii}(:,:);
  ds_1 = sqrt(gradient(x1(:,1)).^2 + gradient(x1(:,2)).^2);
  ds_cum = cumsum(ds_1);
  ds_cum_1 = [ds_cum_1,ds_cum'];
  
end 

disp('done')

for jj=1:size(xp,1)
    
  ds_2 = sqrt(gradient(xp(:,jj)).^2 + gradient(yp(:,jj)).^2);
  ds_cum = cumsum(ds_2);
  ds_cum_2 = [ds_cum_2,ds_cum'];
    
end

disp('done')

[N1,X1] = hist(ds_cum_1);
figure(5)
loglog(X1,N1,'DisplayName','length in displacement(px)');
xlabel 'trajectory length (px)'; ylabel 'count';
legend SHOW;

[N2,X2] = hist(ds_cum_2);
figure(6)
loglog(X2,N2,'DisplayName','length in displacement(px)');
xlabel 'trajectory length (px)'; ylabel 'count';
legend SHOW;


fileID = fopen([folder_save,'tp_',sprintf('%d',no),'.bin']);
tp = fread(fileID,[nt,np],'double','s');
fclose(fileID);

%%
%tracks = load('/media/fmiele/Elements/exp_pt_lower_conc/tracking_results.dat')

figure(1)
clf;
for jj=1:length(tracks)
    
    ind = find(tracks(:,4) == jj);
    figure(1)
    plot(tracks(ind,1),tracks(ind,2))
    axis equal
    hold on
    pause
    
end

%% time ordering after gluing tracks
nn = 600;
cm2 = hot(nn);

script;clc;
clear all;warning off;
fh = findall(0,'type','figure');
for indY=1:length(fh)
  clf(fh(indY));
end
pause(0);

dir = '/media/fmiele/Elements/lagr_tra_0808/mask_traject/';
e = imread([dir,'mask_1234.png']);
traces_o_full = load('/media/fmiele/Elements/lagr_tra_0808/mask_traject/glued_tracks_fullset.mat');

% figure(8)
% clf;
% imagesc(e)
% colormap(cm2)
% colorbar
% set(gca,'Clim',[-0.3/nn,0.3])
% axis equal
% hold on

%% test tracks

traces_o_full = traces_new;

for ii=1:length(traces_o_full.coos)
    
    if(mod(ii,1000) == 0)
       ii
    end

  x1 = traces_o_full.coos{ii}(:,1);
  y1 = traces_o_full.coos{ii}(:,2);
  t1 = traces_o_full.coos{ii}(:,3);
  
  dt = diff(t1);
  ind_cut1 = find(dt <= 0)
  
  if(~isempty(ind_cut1))
      ii
      disp('negative dt')
      t1(ind_cut1)
      
      t1(ind_cut1)
      dt(ind_cut1)
      
      
      figure(1)
      clf;
      plot(x1,y1)
      hold on
      plot(x1(ind_cut1),y1(ind_cut1),'og')
      axis equal
      
      figure(2)
      clf;
      plot(x1,t1)
      hold on
      plot(x1(ind_cut1),t1(ind_cut1),'og')
      pause
  end
  
  %dt = abs(gradient(t1));
  ddt = abs(gradient(dt));
  
  ind_cut2 = find(abs(ddt)./abs(dt) >= 0.8);
  
  if(~isempty(ind_cut2))
      
      ind_cut2
      disp('jump too big')
      ii
      t1(ind_cut2)
      pause

%     ind_cut = ind_cut2(1) + 2;
%     traces_o_full.coos{ii}(:,4) = 0;
%     traces_o_full.coos{ii}(ind_cut2,4) = 1;
%     t1(ind_cut) = t1(ind_cut-1) + dt(ind_cut +1);
%     t1(ind_cut+1:end) = t1(ind_cut) + cumsum(dt(ind_cut+1:end));
%   
%     traces_o_full.coos{ii}(:,3) = t1;
  end
  
  %figure(8)
  %plot(x1,y1,'.--')
  %hold on
  %axis equal
  %drawnow
end 
disp('done')


traces_new.opt{1}

for mm = 1:length(traces_new.coos)
    
     pair=traces.opt.pairsID(mm,1:2);
     
     x_p = traces_c.coos{pair(1)};
     
     x_c = traces_c.coos{pair(2)};
    
     x_g = traces_new.coos{mm}; 
     
     dt_p = diff(x_p(:,3));
     dt_g = diff(x_g(:,3));
     
     ind_cut_p = find(dt_p <= 0)
     ind_cut_g = find(dt_g <= 0)
     
     
     if(dt_p(ind_cut_p) == 0)
         disp('dt zero')
     else
         disp('dt negative')
     end
     
     figure(5)
     clf;
     plot(x_p(:,1),x_p(:,2),'.b')
     hold on
     plot(x_c(:,1),x_c(:,2),'.r')
     hold on
     plot(x_p(ind_cut_p,1),x_p(ind_cut_p,2),'og')
     axis equal
     
     sqrt((x_p(ind_cut_p-1,1) - x_p(ind_cut_p,1))^2 + (x_p(ind_cut_p-1,2) - x_p(ind_cut_p,2))^2)
     
     figure(6)
     clf;
     plot(x_g(:,1),x_g(:,2),'.g')
     axis equal
     
%      x_old = traces_temp.coos{pair(1)};
%      
%      figure(7)
%      clf;
%      plot(x_old(:,1),x_old(:,2),'.m')
%      axis equal
     
     pause 
end


%% b method

dir = '/media/fmiele/Elements/lagr_tra_0808/mask_traject/';
e = imread([dir,'mask_1234.png']);
%traces_o_full = load('/media/fmiele/Elements/lagr_tra_0808/mask_traject/glued_tracks_fullset.mat');


for ii=1:length(traces_o_full.coos)

  x1 = traces_o_full.coos{ii}(:,1);
  %y1 = traces_o_full.coos{ii}(:,2);
  t1 = traces_o_full.coos{ii}(:,3);
  dt = diff(t1);
  dt = [dt(1),dt'];
  aa = find(abs(dt) >= 2);
  bb = find(dt < 0 );
  
  if(~isempty(aa) && ~isempty(bb))
      
   ind_cut_temp = [aa,bb];
   ind_cut = unique(ind_cut_temp);
   
    if(length(ind_cut) == 1)
        ind_cut = [ind_cut,length(t1)+1];
    end
    
    if (ind_cut(1) == 1)
        ind_cut(1) = [];
    end
      
      for k =1:length(ind_cut)-1  
        traces_o_full.coos{ii}(:,4) = 0;
        traces_o_full.coos{ii}(ind_cut(k),4) = 1;
        
        dt_temp = abs(dt(ind_cut(k)-1));
        t1(ind_cut(k):ind_cut(k+1)-1) = t1(ind_cut(k):ind_cut(k+1)-1) - t1(ind_cut(k)) + t1(ind_cut(k)-1) + dt_temp;
        traces_o_full.coos{ii}(:,3) = t1;
      end

  figure(8)
  plot(x1,t1,'.--')
  axis equal
  drawnow
  
  [cc,dd] = find(diff(t1)<=0) 
  size(x1)
  pause
  end  
  
%   figure(8)
%   plot(x1,y1,'.--')
%   hold on
%   axis equal
%   drawnow
end 
disp('done')

%%

t1 = [linspace(1,20,10),linspace(4,50,15),linspace(70,80,7),linspace(3,7,10)];
dt = abs(diff(t1));
dt = [dt(1),dt];

ddt = abs(diff(dt));
ddt = [ddt(1),ddt];
  
ind_cut2 = find(abs(ddt)./abs(dt) >= 0.8);
ind_cut = ind_cut2 + 2;

[pks,loc] = findpeaks(dt)
[pks1,loc1] = findpeaks(ddt)
[pk2,lo2] = findpeaks(ddt./dt)

a = 1:42;
figure(4)
plot(a,dt)
figure(5)
plot(a,ddt)
figure(6)
plot(a,t,'-')
%% save struct file

coos = traces_o_full.coos;
%opt = traces_o_full.opt;
save([dir,'traces_full_length_t_ordered.mat'],'coos');
%save([dir,'opt_traces_full_set_t_ordered.mat'],'opt');
%save([dir,'displacement_2_full_set.dat'],'displacement','-ASCII');

disp('done!')

%% plotting

dir = '/media/fmiele/Elements/lagr_tra_0808/mask_traject/';
e = imread([dir,'mask_1234.png']);

traces_o_full = traces_2345;

nn = 600;
cm2 = hot(nn);

figure(8)
clf;
imagesc(e)
colormap(cm2)
colorbar
set(gca,'Clim',[-0.3/nn,0.3])
axis equal
hold on

Np = length(traces_o_full.coos);

count = 0;
for ii=1:Np

  t1 = traces_o_full.coos{ii}(:,3)
  pause
  %aa = find(t1<0)
  
  if(length(traces_o_full.coos{ii}) > 1)
    count = count +1;
    fprintf('plotted %d over %d...\n',count,Np)
    x1 = traces_o_full.coos{ii}(:,1);
    y1 = traces_o_full.coos{ii}(:,2);      
    %ind_g = find(traces_o_full.coos{ii}(:,4));
    figure(8)
    plot(x1,y1,'.--')
    hold on
    %plot(x1(ind_g),y1(ind_g),'*g')
    %hold on
    axis equal
    drawnow
  end
    
end

folder_print = dir;

figure(8)
set(gcf,'PaperSize',[8,6],'PaperPosition',[0,0,8,6]);
print('-dpdf','-r600',[folder_print,'tracks2loop.pdf']);





