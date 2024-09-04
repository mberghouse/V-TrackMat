%% loading data
clear all;clc;close all;

nn = 600;
cm2 = hot(nn);

dir1 = '/media/fmiele/Elements/lagr_tra_0808/glued_tracks/';

%dir2 = 'G:\lagr_tra_0808\pos_2_tracks_ds_600\';
X2 = load([dir1,'glued_tracks_pos2_t_ordered.mat']);
mask2 = imread([dir1,'mask_pos2.png']);

%dir3 = 'G:\lagr_tra_0808\pos_3_tracks_ds_600\';
X3 = load([dir1,'glued_tracks_pos3_t_ordered.mat']);
%mask3 = imread([dir3,'mask_pos3.png']);

%dir4 = 'G:\lagr_tra_0808\pos_4_tracks_ds_600\';
X4 = load([dir1,'glued_tracks_pos4_t_ordered.mat']);
%mask4 = imread([dir4,'mask_pos4.png']);

%dir5 = 'G:\lagr_tra_0808\pos_5_tracks_ds_600\';
X5 = load([dir1,'glued_tracks_pos5_t_ordered.mat']);
%mask5 = imread([dir5,'mask_pos5.png']);

% Shift and reshape

dxdy = [33 24
    33 24
    36 23];

dx5 = length(mask2)*3 - 3; 
dx4 = length(mask2)*2 - 2;
dx3 = length(mask2) - 1;

'... data loaded ...'

%% shift the tracks

% general parameters to glue: I found them by hand. 

alpha = 1.15; % tollerance on the overlap area (in px)
dx_shift = 15; % shift plus in px

%pos 5
for jj = 1:length(X5.coos)
    X5.coos{jj}(:,1) =  X5.coos{jj}(:,1) + dx5-11;
    
    tt = diff(X5.coos{jj}(:,3));
    tt =[tt(1),tt'];    
    
    bb = find(tt < 0);
  
     if(~isempty(bb))
         disp('hmmm.... dt negative')
        if(length(bb) == 1 && bb(1)~=length(X5.coos{jj}))
         dt = cumsum(tt(bb+1:end));
         dt = [dt(1),dt];
         X5.coos{jj}(bb:end,3) = X5.coos{jj}(bb-1,3) + dt';
        elseif (bb(1) == 1)
            X5.coos{jj} = [];
            continue
        else
          X5.coos{jj}(bb(1):end,:) = [];
          continue
        end
     end  
    
    cc = find(tt == 0); 
     if(~isempty(cc))
         disp('dt equal to zero')
     end
    
end

disp('pos5 ended');

%% pos4
for jj = 1:length(X4.coos)
    X4.coos{jj}(:,1) =  X4.coos{jj}(:,1) + dx4 + dxdy(3,1)-5.5;
    X4.coos{jj}(:,2) =  X4.coos{jj}(:,2) + dxdy(3,2);
     
    tt = diff(X4.coos{jj}(:,3));
    tt =[tt(1),tt'];
    
    bb = find(tt < 0);
      if(~isempty(bb))
         disp('hmmm.... dt negative')
        if(length(bb) == 1 && bb(1)~=length(X4.coos{jj}))
         dt = cumsum(tt(bb+1:end));
         dt = [dt(1),dt];
         X4.coos{jj}(bb:end,3) = X4.coos{jj}(bb-1,3) + dt';
         elseif (bb(1) == 1)
            X4.coos{jj} = [];
            continue
        else
          X4.coos{jj}(bb(1):end,:) = [];
          continue
        end
      end
      
    cc = find(tt == 0); 
     if(~isempty(cc))
         disp('dt equal to zero...')
     end
    
    %add crop over the traject that cross the overlap
    if (X4.coos{jj}(end,1) > dx5 + dx_shift )
        ind = find(X4.coos{jj}(:,1) >= dx5 + dx_shift,1);
        X4.coos{jj}(ind:end,:) = [];
    end
end

disp('pos4 ended');

%% pos3

for jj = 1:length(X3.coos)
    X3.coos{jj}(:,1) =  X3.coos{jj}(:,1) + dx3 + sum(dxdy(1:2,1))-2;
    X3.coos{jj}(:,2) =  X3.coos{jj}(:,2) + sum(dxdy(1:2,2));
    
    tt = diff(X3.coos{jj}(:,3));
    tt =[tt(1),tt'];
    
    bb = find(tt < 0);
   if(~isempty(bb))
        disp('hmmm.... dt negative')
        if(length(bb) == 1 && bb(1)~=length(X3.coos{jj}))
         dt = cumsum(tt(bb+1:end));
         dt = [dt(1),dt];
         X3.coos{jj}(bb:end,3) = X3.coos{jj}(bb-1,3) + dt';
         elseif (bb(1) == 1)
            X3.coos{jj} = [];
            continue
        else
          X3.coos{jj}(bb(1):end,:) = [];
          continue
        end
   end
    
    cc = find(tt == 0); 
     if(~isempty(cc))
         disp('dt equal to zero...')
     end
    
    %add crop over the traject that cross the overlap
    if (X3.coos{jj}(end,1) > dx4 + alpha*dxdy(3,1) + dx_shift)
        ind = find(X3.coos{jj}(:,1) >= dx4 + alpha*dxdy(3,1) + dx_shift,1);
        X3.coos{jj}(ind:end,:) = [];
    end
end

disp('pos3 ended');

%% pos2

for jj = 1:length(X2.coos)
    X2.coos{jj}(:,1) =  X2.coos{jj}(:,1) + sum(dxdy(:,1));
    X2.coos{jj}(:,2) =  X2.coos{jj}(:,2) + sum(dxdy(:,2));
     
    tt = diff(X2.coos{jj}(:,3));
    tt =[tt(1),tt'];
    
    bb = find(tt < 0);
 
     if(~isempty(bb))
         disp('hmmm.... dt negative')
        if(length(bb) == 1 && bb(1)~=length(X2.coos{jj}))
         dt = cumsum(tt(bb+1:end));
         dt = [dt(1),dt];
         X2.coos{jj}(bb:end,3) = X2.coos{jj}(bb-1,3) + dt';
        else
          X2.coos{jj}(bb(1):end,:) = [];
          continue
        end
     end
    
     cc = find(tt == 0);
     
     if(~isempty(cc))
         disp('dt equal to zero')
     end
    
    %add crop over the traject that cross the overlap
     if (X2.coos{jj}(end,1) > dx3 + alpha*sum(dxdy(1:2,1)))
         ind = find(X2.coos{jj}(:,1) >= dx3 + alpha*sum(dxdy(1:2,1)),1);
         X2.coos{jj}(ind:end,:) = [];
     end
end

disp('pos2 ended');
%% collect tracks in one single struct file

traces_23.coos = X2.coos;
j2 = length(X2.coos);
for k=1:length(X3.coos)
    traces_23.coos{k+j2} = X3.coos{k};
end

traces_234.coos = traces_23.coos;
j3 = length(traces_23.coos);
for kk =1:length(X4.coos)
    traces_234.coos{kk + j3} = X4.coos{kk};
end

traces_2345.coos = traces_234.coos;
j4 = length(traces_234.coos);
for kkk = 1:length(X5.coos)
    traces_2345.coos{j4 + kkk} = X5.coos{kkk};
end

traces_temp = traces_2345;
dd = [];
disp('here')
for m = 1:length(traces_temp.coos)
    if(numel(traces_temp.coos{m}) < 1)
     dd = [m,dd];
    end
end
traces_temp.coos(dd) = [];

disp('done')

%%

% load the full mask
c = imread([dir1,'mask_1234.png']);

% plot the trajectories with a different marker for different FoV

figure(20)
clf;
imagesc(c)
colormap(cm2)
colorbar
set(gca,'Clim',[-0.3/nn,0.3])
axis equal
hold on
for jjj=1:500:length(traces_2345.coos)    
 figure(20)
 plot(traces_2345.coos{jjj}(:,1),traces_2345.coos{jjj}(:,2),'.r')
 hold on
 
 if (jjj>j2)
  plot(traces_2345.coos{jjj}(1,1),traces_2345.coos{jjj}(1,2),'*g')
  hold on
 end
 
 if (jjj>j3)
  plot(traces_2345.coos{jjj}(1,1),traces_2345.coos{jjj}(1,2),'*k')
  hold on
 end
 
 if (jjj>j4)
  plot(traces_2345.coos{jjj}(1,1),traces_2345.coos{jjj}(1,2),'*m')
  hold on
 end
 
 axis equal
 hold on
 drawnow
end

coos = traces_2345.coos;
%opt = traces_o.opt;
save([dir1,'traces_2345_FOVs_t_ordered.mat'],'coos');
%save([dir,'opt_2_full_length_t_ordered.mat'],'opt');
%save([dir,'displacement_2_full_set.dat'],'displacement','-ASCII');

disp('done!')







