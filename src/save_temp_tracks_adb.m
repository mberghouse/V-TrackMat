function np = save_temp_tracks_adb(x,y,t,folder_save,xp_name,yp_name,tp_name,ind_death,ds,np,no, Nframes, ii)
% Modified to keep the uninterpolated x y t (choose option = 0)
option = 0;

if option == 0
    short_tr = 0; void_tr = 0;
    if isempty(ind_death)==0
        for jj = 1:length(ind_death)
            %size(~isnan(x(jj,:)))
            %length(ind_death)
                      
            ind_a = find(~isnan(x(jj,:)));
            ind_g = ind_a;
%             [~,ind_g] = unique(x(jj,ind_a),'stable'); % Make sure we dont undersample
%             if ~isempty(ind_g)
%                 ind_g = ind_g + (ind_a(1)-1);
%             end
            n_short = length(ind_g);
            
            if numel(find(isnan(x(jj,:)))) == numel(x(jj,:))
                void_tr=void_tr+1;
                %             pause
                continue;
            end
            
            if(n_short<5)
                short_tr = short_tr +1;
                continue
            end
            
            xs(1:n_short) = x(jj,ind_g);
            xs(n_short+1:Nframes) = NaN;
            ys(1:n_short) = y(jj, ind_g);
            ys(n_short+1:Nframes) = NaN;
            ts(1:n_short) = t(jj, ind_g);
            ts(n_short+1:Nframes) = NaN;
            deltat = diff(ts(:));
            
            if max(deltat)>5
                clc
                deltat
                figure(1); hold on
                plot(xs, ys); axis equal xy
                i_fault = find(deltat>5)
                fprintf('Trajectory index: %d', jj)
                pause
            end
            
            %         ts
            %         pause;
            fileIDx = fopen([folder_save,xp_name],'a');
            fwrite(fileIDx,xs,'double','s');
            fclose(fileIDx);
            fileIDy = fopen([folder_save,yp_name],'a');
            fwrite(fileIDy,ys,'double','s');
            fclose(fileIDy);
            fileIDt = fopen([folder_save,tp_name],'a');
            fwrite(fileIDt,ts,'double','s');
            fclose(fileIDt);
            
            
        end
        
        [void_tr short_tr length(ind_death)]
        p = length(ind_death) - short_tr -void_tr;
        % pause
        np = np + p;
        param = [np,Nframes];
        save([folder_save,'param_save_',sprintf('%s',no),'.dat'],'param','-ascii', '-v7.3');
        disp('Executed save function');
    end
    
end
if option == 1
    %%
    short_tr = 0; void_tr = 0;
    for jj=1:length(ind_death)
        %         jj
        if numel(find(isnan(x(jj,:)))) == numel(x(jj,:))
            void_tr=void_tr+1;
            %             pause
            continue;
        end
        
        ind_a = find(~isnan(x(jj,:)));
        xy = [x(jj,ind_a); y(jj,ind_a)]';
        [~,ind_g] = unique(xy, 'rows', 'stable');
        ind_g = ind_g + (ind_a(1)-1);
        ds_m = sqrt(diff(x(jj,ind_g)).^2 + diff(y(jj,ind_g)).^2);
        %         pause
        ds_m = [0,ds_m];
        dss_m = cumsum(ds_m(:));
        %dss_m = unique(dss_m);
        
        if((numel(ds_m) < 2) || dss_m(end) < ds(30))
            short_tr = short_tr +1;
            
            
            continue
        end
        fprintf('Short trajectories= %d\n', numel(short_tr))
        aa = find(ds < dss_m(end));
        a = aa(end);
        
        ts(1:a) = interp1(dss_m,t(jj,ind_g),ds(1:a));
        if(min(diff(ts(1:a))) < 0)
            disp('porco dio:');
            disp('====');
            ds(1:a)
            % pause
        end
        xs(1:a) = interp1(t(jj,ind_g),x(jj,ind_g),ts(1:a));
        ys(1:a) = interp1(t(jj,ind_g),y(jj,ind_g),ts(1:a));

        dt = diff(ts);
        cc = find(dt == 0);
        
        if(~isempty(cc))
            jj
            disp('dagane')
            pause
        end
        
        xs(a+1:length(ds)) = NaN;
        ys(a+1:length(ds)) = NaN;
        ts(a+1:length(ds)) = NaN;
        
        clc
        fileIDx = fopen([folder_save,xp_name],'a');
        fwrite(fileIDx,xs,'double','s');
        fclose(fileIDx);
        fileIDy = fopen([folder_save,yp_name],'a');
        fwrite(fileIDy,ys,'double','s');
        fclose(fileIDy);
        fileIDt = fopen([folder_save,tp_name],'a');
        fwrite(fileIDt,ts,'double','s');
        fclose(fileIDt);
        %         fprintf('Saved trajectory %d\n', jj);
    end
    
    [void_tr short_tr length(ind_death)]
    p = length(ind_death) - short_tr -void_tr;
    % pause
    np = np + p;
    param = [np,length(ds)];
    save([folder_save,'param_save_',sprintf('%s',no),'.dat'],'param','-ascii');
    disp('Executed save function');
end
end

