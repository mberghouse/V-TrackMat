function coos_new = glue_all_traces_adb(jg, traces)

% Find indiviudual pair
% pair=traces.opt.pairsID(jg,1:2);
% coos1 = traces.coos{pair(1)}; % coos of trace1
% coos2 = traces.coos{pair(2)}; % coos of trace2

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
figure; plot(coos1(:,1),coos1(:,2), 'or', coos2(:,1), coos2(:,2), 'sb'); axis equal
hold on
plot(xc,yc,'k')
plot(coos1(IN1,1),coos1(IN1,2), '.', coos2(IN2,1), coos2(IN2,2), '.'); axis equal

if ~isempty([IN1;IN2])
    overlap = 1;
else
    overlap = 0;
end
overlap
%
clear coos_new
% compute velocity to do time interpolation
if overlap == 1
    vel_x= gradient(coos2(:,1))./gradient(coos2(:,3));
    vel_y= gradient(coos2(:,2))./gradient(coos2(:,3));
    vel_abs = sqrt(vel_x(1).^2+vel_y(1).^2);
elseif overlap == 0
    vel_x= gradient(coos1(:,1))./gradient(coos1(:,3));
    vel_y= gradient(coos1(:,2))./gradient(coos1(:,3));
    vel_abs = sqrt(vel_x(end).^2+vel_y(end).^2);
end
dist_t = dist_s/(vel_abs);

if overlap==1
    % remove points in the overlapping region to later interpolated
   ind1 = [find(IN1); numel(coos1(:,1))];
   ind2 = [1; find(IN2)];
   coos1(ind1,:) = [];
   coos2(ind2,:) = [];
   overlap = 0; 
end
% move t=0 time to the beginning of trace1
t1 = coos1(:,3); t1 = t1-t1(1);
t2 = coos2(:,3); t2 = t2-t2(1)+t1(end)+dist_t;

% move s=0 to the beginning of trace 1
s1 = sqrt(diff(coos1(:,1)).^2+diff(coos1(:,2)).^2); s1 = [0; s1];
s2 = sqrt(diff(coos2(:,1)).^2+diff(coos2(:,2)).^2); s2 = [0; s2];

s1 = cumsum(s1);
s2 = cumsum(s2); s2 = s2+s1(end)+dist_s;

ds = mean(diff([s1;s2]));
if overlap ==0
    %dist_s = s2(1)-s1(end);
    n_int = ceil(dist_s/ds)-1;
    % interpolate with equally spaced points separated by ds
%     s_int = s1(end)+ds: ds: s1(end)+n_int *ds;
    s_int = s1(1):ds:s2(end);
    s_int = s_int(:);
    coos_new(:,1) = interp1([s1;s2],[coos1(:,1); coos2(:,1)], s_int, 'spline');
    coos_new(:,2) = interp1([s1;s2],[coos1(:,2); coos2(:,2)], s_int, 'spline');
    coos_new(:,3) = interp1([s1;s2],[t1; t2], s_int, 'spline');
%     coos_new(:,1) = [coos1(:,1);
%         interp1([s1;s2],[coos1(:,1); coos2(:,1)], s_int, 'pchip');
%         coos2(:,1)];
%     coos_new(:,2) = [coos1(:,2);
%         interp1([s1;s2],[coos1(:,2); coos2(:,2)], s_int, 'pchip');
%         coos2(:,2)];
%     coos_new(:,3) = [t1;
%         interp1([s1;s2],[t1; t2], s_int, 'pchip');
%         t2];
end
% figure; plot(t1, coos1(:,1), '+m', t2, coos2(:,1), 'xg'); 
% figure; plot(s1, coos1(:,1), '+m', s2, coos2(:,1), 'xm'); axis equal
figure(2); plot(coos1(:,1),coos1(:,2), 'or', coos2(:,1), coos2(:,2), 'sb'); axis equal;
hold on
figure(2); plot(coos_new(:,1), coos_new(:,2), 'xg')