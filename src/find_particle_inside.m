function p_ind = find_particle_inside(xie, yie, xjs, yjs, theta, distMax)
% searches points inside an ellipse with major axis alined with theta
% theta has to be in radian

a=distMax(2); % major axis
b=distMax(1); % minor axis

% if major axis is too large make it circle
if a/b>20
    a=b;
end

t = linspace(-pi,pi,100);
x = xie + a*cos(t)*cos(theta) - b*sin(t)*sin(theta);
y = yie + b*sin(t)*cos(theta) + a*cos(t)*sin(theta);


IN = inpolygon(xjs,yjs,x,y);

p_ind = find(IN==1);
% plot(x,y); hold on
% axis equal xy
% pause


