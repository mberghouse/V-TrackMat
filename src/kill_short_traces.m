%%%kill_short_traces function for UNIL
function traces=kill_short_traces(traces)

lengthTraces=[];
ind=[];
indminD=[];
anyp4found=0;

clf   %clear current figure

subplot(1,2,1);
p3 = [];
p4 = [];

for jj=1:length(traces.coos)
jj
    if(numel(traces.coos{jj}) < 1)
        disp('ecchime...');
        continue
    end
    
    if(mod(jj,500) == 0)
        %disp('500 less...')
        jj
    end
    
   % p1 = plot(traces.coos{j}(1,1),traces.coos{j}(1,2),'k*','markersize',3,'DisplayName','Start of spaghetti');          %label start of spaghetti
   % hold on
   % p2 = plot(traces.coos{j}(:,1),traces.coos{j}(:,2),'b.-','markersize',3,'DisplayName','All spaghetti');              %plot all spaghetti
    lengthTraces(jj)=length(traces.coos{jj});
    
     if length(traces.coos{jj})<traces.opt.minLength
            p3 = plot(traces.coos{jj}(:,1),traces.coos{jj}(:,2),'r.-','markersize',3,'DisplayName', 'Short spaghetti');   %plot short spaghetti
            ind=[ind;jj];
     end
     dist=ipdm(traces.coos{jj}(1,1:2),traces.coos{jj}(end,1:2));
     if dist < traces.opt.minDisplacement
            p4 = plot(traces.coos{jj}(:,1),traces.coos{jj}(:,2),'m.-','markersize',3,'DisplayName', 'Immobile spaghetti'); %plot immobile spaghetti
            ind=[ind;jj];
            indminD=[indminD;length(traces.coos{jj})];
            anyp4found=1;
     end
     
end
axis equal;axis tight;
hold off

% if length(ind)>0 & anyp4found == 1
%     legend([p1 p2 p3 p4],'Location','southeast')
% elseif length(ind)>0 & anyp4found == 0
%     legend([p1 p2 p3],'Location','southeast')
% else
%     legend([p1 p2],'Location','southeast')
% end

title(sprintf('minLength = %d',round(traces.opt.minLength)));
xlabel('x');ylabel('y');

ind=unique(ind);

traces.coos(ind)=[]; %deletes short or spaghetti (okay to delete with parenthesis)






