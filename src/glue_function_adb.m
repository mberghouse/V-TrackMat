function coos_new = glue_function_adb(coos1,coos2, option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that takes two trajectories and joins them based on:
% overlap 1 (there is overlapping region)
% overlap 0 (there is no overlapping region)
% Input:
% coos1: trajectory 1 coordinates
% coos2: trajectory 2 coordinates
% option: 0= equal dt, 1=equal ds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% distance from end of trace1 to beginning of trace2
dist_s = sqrt((coos1(end,1)-coos2(1,1))^2 + (coos1(end,2)-coos2(1,2))^2);

% mid point between the two ends
c1 = [coos1(end,1) coos1(end,2)];
c2 = [coos2(1,1) coos2(1,2)];
c(1) = (c1(1)+c2(1))/2;
c(2) = (c1(2)+c2(2))/2;
theta = linspace(0,2*pi,101);
xc = c(1)+dist_s/2*cos(theta);
yc = c(2)+dist_s/2*sin(theta);

% find points that overlaps the two trajectories
IN1 = inpolygon(coos1(:,1), coos1(:,2), xc,yc);
IN2 = inpolygon(coos2(:,1), coos2(:,2), xc,yc);

% visualize the trajectories with common points
% figure(1); plot(coos1(:,1),coos1(:,2), 'or', coos2(:,1), coos2(:,2), 'sb'); axis equal
% hold on
% plot(xc,yc,'k')
% plot(coos1(IN1,1),coos1(IN1,2), '.', coos2(IN2,1), coos2(IN2,2), '.'); axis equal
% return

if ~isempty([IN1;IN2])
    overlap = 1;
else
    overlap = 0;
end

% compute velocity to do time interpolation
if overlap == 1
    vel_x= gradient(coos2(:,1))./gradient(coos2(:,3));
    vel_y= gradient(coos2(:,2))./gradient(coos2(:,3));
    vel_abs = sqrt(vel_x.^2+vel_y.^2);
    % take the first finite velocity at the end of trajectory 1
    vel_abs = vel_abs(find(vel_abs>0,1, 'first'));
elseif overlap == 0
    vel_x= gradient(coos1(:,1))./gradient(coos1(:,3));
    vel_y= gradient(coos1(:,2))./gradient(coos1(:,3));
    vel_abs = sqrt(vel_x.^2+vel_y.^2);
    % take the last finite velocity at the end of trajectory 1
    vel_abs = vel_abs(find(vel_abs>0,1, 'last'));
end
dist_t = dist_s/(vel_abs);

if overlap==1
    % remove points in the overlapping region to later interpolated
    % remove last point in traj-1 and first point in traj-2
    ind1 = [find(IN1); numel(coos1(:,1))];
    ind2 = [1; find(IN2)];
    coos1(ind1,:) = [];
    coos2(ind2,:) = [];
    if isempty(coos1)
        keep=coos2;
    elseif isempty(coos2)
        keep=coos1;
    end
    % change overlap to 0 after removing points
    overlap = 0;
end

if ~isempty(coos2) && ~isempty(coos1)
    
    % move t=0 time to the beginning of trace1
    t1 = coos1(:,3); t1 = t1-t1(1);
    t2 = coos2(:,3); t2 = t2-t2(1)+t1(end)+dist_t;
    
    if overlap ==0
        if option ==0
            t = [t1;t2];
            x = [coos1(:,1); coos2(:,1)];
            y = [coos1(:,2); coos2(:,2)];
            
            dtt = diff(t);
            dt = median(dtt);
            
            t_int = t(1):dt:t(end);
            t_int = t_int(:);
            
            [~,ia,~] = unique(t);
            t = t(ia);
            x = x(ia);
            y = y(ia);
            
            if isnan(x)
                disp('Porc')
                pause
            end
            
            coos_new(:,1) = interp1(t,x, t_int);
            coos_new(:,2) = interp1(t,y, t_int);
            coos_new(:,3) = t_int;
        end
        
        if option ==1
            % move s=0 to the beginning of trace 1
            s1 = sqrt(diff(coos1(:,1)).^2+diff(coos1(:,2)).^2); s1 = [0; s1];
            s2 = sqrt(diff(coos2(:,1)).^2+diff(coos2(:,2)).^2); s2 = [0; s2];
            
            s1 = cumsum(s1);
            s2 = cumsum(s2); s2 = s2+s1(end)+dist_s;
            
            ds = median(diff([s1;s2]));
            s_int = s1(1):ds:s2(end);
            s_int = s_int(:);
            coos_new(:,1) = interp1([s1;s2],[coos1(:,1); coos2(:,1)], s_int, 'linear');
            coos_new(:,2) = interp1([s1;s2],[coos1(:,2); coos2(:,2)], s_int, 'linear');
            coos_new(:,3) = interp1([s1;s2],[t1; t2], s_int, 'linear');
        end
        if max(coos_new(:,1))>2048
            vis=1;
            disp('Shit!')

        else
            vis=0;
        end
        if vis==1
            figure(1); clf
            subplot(1,2,1)
            plot(t,x,'.r'); hold on
            plot(coos_new(:,3), coos_new(:,1), 'o');
            subplot(1,2,2)
            plot(t,y,'.r'); hold on
            plot(coos_new(:,3), coos_new(:,2), 'o');
            %         plot(coos1(:,3),s1, '.r'); hold on
            %         plot(coos2(:,3),s2, '.b'); hold on
            %         iind1 = find(diff(coos1(:,3))<0);
            %         iind2 = find(diff(coos2(:,3))<0);
            %         plot(s1(iind1),coos1(iind1,3), 'o'); hold on
            %         plot(s2(iind2),coos2(iind2,3), 's'); hold on
            
%             figure(2); clf
%             plot(coos1(:,1),coos1(:,2), 'or', coos2(:,1), coos2(:,2), 'sb'); axis equal;
%             hold on
%             figure(2); plot(coos_new(:,1), coos_new(:,2), '.g')
%             
%             figure(3);
%             plot(diff(t), '.'); hold on
%             plot(diff(coos_new(:,3)), '.');

            close all
        end
    end
else
    disp('Bad trajectory')
    vis=1;
    coos_new =keep;
end
