script;clc;
clear all;warning off;
fh = findall(0,'type','figure');
for k=1:length(fh)
    clf(fh(k));
end
pause(0);

%% THIS CODE compute the background image for an experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
level = 2^8 - 1;

S=1;
% folder = 'Z:\colloid_density\20_01_23_het_05uL_B\';
% prefix = ['20_01_23_het_05uL_B_s',sprintf('%d',S),'t'];
folder = 'Z:\Training Experiment\03042020\PT\pos1\test01_04032020_1um_DAPI_0.6uLpm_pos_1\';
prefix = 'dapi_0.6ulpm_002t';

IMAGES = double(imread([folder,prefix,sprintf('%04d',S),'.tif']) )./ level;
[Ny,Nx,~] = size(IMAGES);
% N = 1000;
N = 7199;
% return
back = zeros(Ny,Nx);

% 7 - background calculation loop:
for k=1:N
    im = double(imread([folder,prefix,sprintf('%04d',k),'.tif'])) ./ level;
    back = back + im ./ N;
    fprintf('image no. %d \n',k)
end

imwrite(back,[folder,'background',sprintf('%d',S),'.png'],'png','BitDepth',8);
disp('done.');

%% ------------------------------------------
% this script analyze all the images
% and save particles coordinate for each time
% so that those can be analyszed with the track_part function
% ------------------------------------------
%

script;clc;
clear all;warning off;
fh = findall(0,'type','figure');
for k=1:length(fh)
    clf(fh(k));
end
pause(0);
% cd 'C:\Users\dscheidw\Dropbox\Matlab_scripts\Windows_office\Particle_Tracking/';

% nn = 300;
% cm = hot(nn);
% cm = cm(end:-1:1,:);

% 1 - physical parameters
level = 2^8 - 1;         % [bit]
S=1;

% 2 - Particle tracking parameters
% N = 100;
N = 7199;
visu = 1;
% tr2=0.10; % for  motile
tr2=0.09; % for beads
lobj = 50;
% folder = 'S:\scheidwe\EXPORT\Pietro_chambers\190521_motile_diffusion-05\';
% prefix = '190521_motile_diffusion-05_t';

folder = 'Z:\Training Experiment\03042020\PT\pos1\test01_04032020_1um_DAPI_0.6uLpm_pos_1\';
prefix = 'dapi_0.6ulpm_002t';

folder_save = [folder,'data_raw',sprintf('%d',S),'\'];mkdir(folder_save);
back = double(imread([folder,'background',sprintf('%d',S),'.png'])) ./ level;
% return
% back = back(1200:1800,1200:1800);
pos_lst=[];
xold = [];
yold = [];
ANGLE=[];
count=0;
%%
for k=1%:1:N
    
    %%
    imo = double(imread([folder, prefix, sprintf('%04d.tif',k)])) ./ level;
    %     imo = imo(1200:1800,1200:1800);
    
    %     imo = imgaussfilt(imo,1);
    im = imo-back;
    %     im=rgb2gray(im);
    im(im < 0) = 0;
    
    [gradx,grady] = ((gradient(im)));
    grad = sqrt(gradx.^2 + grady.^2);
    %     grad = (imgaussfilt(grad,1));
    % grad = stdfilt(im);
    
    %     b = bpass(grad,1,5);
    
    %     b = b ./ (max(b(:)));
    b = im ./ (max(im(:)));
    bb=b;
    b(b < 0.05) = 0;
    %     b(b>=tr2)=1;
    %     b = watershed(b);
    %     b = b(200:2001,1:1801);  % crop for local analysis
    b=bwareaopen(b,5);
    %     pk = pkfnd(b,0,lobj);
    %     cnt = cntrd(b,pk,lobj);
    %
    %     if    isempty(cnt)==1
    %         %         k=k+1;
    %         continue
    %     end
    
    
    %     cc(k,1)=length(cnt);
    %
    %     if   cc(k,1) > 100
    %
    %         cc(k,1)=nan;
    %         continue
    %     end
    
    if(visu)
        
        
        %                 hold off
        %                 drawnow
        %                 saveas(gcf,[folder_save,'nonmotile_tracked',sprintf('%d',k),'.png'])
        %         SE= strel('disk',10);
        %         c=bwareaopen(im,lobj);
        %         c = imdilate(c,SE);
        %         c = imerode(c,SE);
        %         c = bwpropfilt(c,'Area',[0 300]);
        %          %         c=b;
        %         c = im2bw(b,0);
        
        %%
        if    max(b(:))==0
            %         k=k+1;
            continue
        end
        
        %         c(c >0) = 1;
        st = regionprops(b,'Centroid','Orientation','MajorAxisLength','Perimeter');
        MajorAxis = cat(1, st.MajorAxisLength)*2;
        centroids = cat(1, st.Centroid);
        angle = cat(1, st.Orientation);
        cnt=length(centroids(:,1));
        xCentre = centroids(:,1);
        yCentre = centroids(:,2);
        
        corr_angle=[];
        
        
        xold = [xold;centroids(:,1)];
        yold = [yold;centroids(:,2)];
        %         figure(1)
        %         clf
        %         imagesc(b)
        %         axis equal xy
        %         colormap(cm)
        %         hold on
        %         plot(cnt(:,1),cnt(:,2),'ro','markersize',4)
        %         plot(xold,yold,'b.','markersize',8)
        %         colorbar
        %         drawnow
        
        figure(2)
        clf
        %         imagesc((bb))
        
        %         imagesc(imadjust(imo,[0 0.12]))
        imshowpair(imo,b,'falsecolor')
        axis equal xy
        colormap(gray)
        hold on
        plot(centroids(:,1),centroids(:,2),'ro','markersize',4)
        plot(xold,yold,'b.','markersize',3);
        axis([1500 1750 1100 1300])
        %         axis equal xy
        
    end
    return
    fprintf('image no. %d particle found %d \n',k, (cnt))
    count(k) = cnt;
    pos_lst = [pos_lst ; centroids(:,1),centroids(:,2),ones(length(centroids(:,1)),1) .* k];
    
end
[Ny,Nx] = size(imo);
mask = zeros(Ny,Nx);
C=nanmean(count);

save([folder_save,'pos_lst.dat'],'pos_lst','-ASCII');
% save([folder_save,'ANGLE.dat'],'ANGLE','-ASCII');
disp('Done! ');
