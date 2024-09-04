%% Define magnification and pixel size
% cd('E:\Ankur\Filtration Project\Filippo data\');
% load('traces_data.mat');
M=6;
px_size = 0.0065 / M;

% Calculate the velocity cells

ux = cellfun(@(x) px_size*gradient(x(:,1))./gradient(x(:,3)), traces_new.coos, 'UniformOutput',0);
uy = cellfun(@(x) px_size*gradient(x(:,2))./gradient(x(:,3)), traces_new.coos, 'UniformOutput',0);



% Convert cells of coos and velocities to vector arrays

xcoos = cellfun(@(x) [x(:,1)], traces_o.coos, 'UniformOutput',0);
ycoos = cellfun(@(x) [x(:,2)], traces_o.coos, 'UniformOutput',0);

xp = cell2mat(xcoos');
yp = cell2mat(ycoos');
umag = cellfun(@(x,y) sqrt(x.^2+y.^2), ux, uy, 'UniformOutput',0);
Vel = cell2mat(umag');
ux = cell2mat(ux'); uy = cell2mat(uy');
Velnorm = Vel/max(Vel(:));

% define colorbar
col = hot(numel(Vel));
% col = flipud(col);

%% Discretize velocity field
nx = 2048;
ny = 2044;
delx = 4; dely = 4;
[X,Y] = meshgrid(1:delx:nx,1:dely:ny);
uxbin = NaN(size(X)); uybin = NaN(size(Y)); partc = NaN(size(X));
for i = 1:size(X,1)-1
    for j = 1:size(X,2)-1
        xv = [X(i,j) X(i,j+1) X(i+1,j+1), X(i+1,j)];
        yv = [Y(i,j) Y(i,j+1) Y(i+1,j+1), Y(i+1,j)];
        [in,~] = inpolygon(xp, yp, xv, yv);
        if isempty(find(in==1, 1))==0
            uxbin(i,j) = nanmean(ux(in==1));
            uybin(i,j) = nanmean(uy(in==1));
            partc(i,j) = nnz(in);
            %         else
            %             uxbin(i,j)=NaN;
            %             uybin(i,j)=NaN;
        end
    end
end
%%
delx1 = unique(diff(unique(X))) * px_size;
delx2 = unique(diff(unique(Y))) * px_size;
[du1dx2, du1dx1] = gradCentralDifference_adb(uxbin, delx1, delx2);
[du2dx2, du2dx1] = gradCentralDifference_adb(uybin, delx1, delx2);

%% Calculate stretching field based on eigen value of Cauchy-Green tensor

% Deformation gradient tensor
F = NaN(4, size(du1dx1,1), size(du1dx1,2));

F(1,:,:) = 1+du1dx1;
F(2,:,:) = du1dx2;
F(3,:,:) = du2dx1;
F(4,:,:) = 1+du2dx2;
stretchingField = NaN(size(du1dx1));
for i=1:size(du1dx1,1)
    for j=1:size(du1dx1,2)
        Cij = reshape(F(:,i,j),2,2)*reshape(F(:,i,j),2,2)';
        if isnan(sum(Cij(:))) == 0
            eigMax = max(eig(Cij));
            stretchingField(i,j) = sqrt(eigMax);
        end
        
    end
end

%%
% change the colorbar matrix based on velocity
data = sortrows([Velnorm xp/4 yp/4 col]);
xp = data(:,2);
yp = data(:,3);
col1 = data(:,4:6);

%% plot scatter plot colored with velocity
figure(4); 
scatter(xp, yp, 3, col); axis equal
title([dir])
% axis([0 2000 0 2000])
% nc = 600;
% cm2 = hot(nc);
% cm2 = [[190,190,170]./255;cm2(end:-1:1,:)];
% colormap(cm2)
% caxis([0 50]);
box on
