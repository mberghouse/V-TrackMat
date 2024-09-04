%%%plot_glued_trajs function for UNIL
function [long,displacement]=plot_glued_trajs(traces,show_longest_trajs)

long = [];
displacement = [];

% clf;

for ii = 1:length(traces.coos)
    
    if rem(ii,1000)
        ii/length(traces.coos)
    end
    long = [long;size(traces.coos{ii},1)];
%     size(traces.coos{ii},1)
    displacement = [displacement; ipdm(traces.coos{ii}(1,1:2),traces.coos{ii}(end,1:2))];
    
%     if (mod(ii,200) == 0)
%      % plot of all final trajectories 
%      subplot(1,2,1);hold on;
%      plot(traces.coos{ii}(:,1),traces.coos{ii}(:,2),'.-');
%      xlabel 'x [pix]'; ylabel 'y [pix]';
%      title 'New longer glued trajectories';
%      axis tight;axis equal;
%      drawnow
%     end
end

%%PDF of trajectory length (frames)
% subplot(1,2,2);hold on;
% [N,X] = hist(long);plot(X,N,'DisplayName','length in frames'); % linear plot
% [N,X] = hist(displacement);plot(X,N,'DisplayName','length in displacement');
% %[N,X] = hist(long);loglog(X,N,'DisplayName','length in frames'); % log-log plot
% %[N,X] = hist(displacement);loglog(X,N,'DisplayName','length in displacement(px)');
% xlabel 'trajectory length (px)'; ylabel 'count';
% legend SHOW;

%%for fun plot of final trajectories of at least 0.5x the maximum length (frames)
% if show_longest_trajs == 1
%     ind = find(long > max(long)/2 | displacement > max(displacement)/2);
%     figure;hold on;
%     
%     for ii = 1:length(ind)
%         plot(traces.coos{ind(ii)}(:,1),traces.coos{ind(ii)}(:,2),'.-');
%         xlabel 'x [pix]'; ylabel 'y [pix]';
%         title('New longer glued trajectories of length (frames or displacment) > 0.5 the maximum');
%         axis tight;axis equal;
%     end
% end
