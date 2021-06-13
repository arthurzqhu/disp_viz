clear
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
load clouds.mat

%%
close all

do_test = false;
% Put all data into a single data structure
campaigns={'vocalspdi','masepdi','postpdi'};

% fractional height below which is considered cloud base
CB_thres = 0.05; 

load dp_pdi.mat

for c= 1:length(campaigns)
    camp=campaigns{c};
    clear fb
    
    % load flight start time, end time, cloud base/top height
    try
        fb = load([camp,'_flight_basics.mat']);
        fbvar = [camp,'_flight_basics'];
    catch
    end
    
    clear gen
    
    var='s_lwc_pdi'; %LWC variable name
    Nvar='s_ntot_pdi'; %Number concentration variable name
    nvar='s_conc_pdi';
    dvar='s_disp_pdi'; %Relative dispersion variable name
    
    switch camp
        case 'vocalspdi'
            bindata=load('dp_pdi.mat'); 
            binlims=bindata.pdi_bins(:,1);
            binlims(end+1)=bindata.pdi_bins(end,2);
            dlogD = log10(binlims(2)/binlims(1));
            filedir='vocals/';
            files=dir([filedir,'*.mat']);
            vocals_aer_files = dir('VOCALS_pcasp/PCASP_08*');
            for ifile=1:length(files)
                gen(ifile)=load([filedir,files(ifile).name]); %Read all of the files
                gen(ifile).s_disp_pdi(gen(ifile).s_disp_pdi==0)=NaN;
                gen(ifile).s_disp_pdi(gen(ifile).s_ntot_pdi<25)=NaN;
            end
        case 'masepdi'
            bindata=load('dp_pdi.mat');
            binlims=bindata.pdi_bins(:,1);
            binlims(end+1)=bindata.pdi_bins(end,2);
            dlogD = log10(binlims(2)/binlims(1));
            filedir='mase/';
            files=dir([filedir,'*.mat']);
            
            % calculate rhoa that is missing in MASE
            mm_v = .018;
            mm_d = .029;
            for ifile=1:length(files)
                gen(ifile)=load([filedir,files(ifile).name]);
                if ifile==1
                    gen=rmfield(gen,'s_wz'); %This is an extra field not found in the other days
                end
                %Remove suspicious RelDisp Values
                y=17.6*gen(ifile).s_ntot_pdi.^(-.8);
                gen(ifile).s_disp_pdi(gen(ifile).s_disp_pdi>y)=NaN;
                gen(ifile).s_disp_pdi(gen(ifile).s_disp_pdi==0)=NaN;
                gen(ifile).s_disp_pdi(gen(ifile).s_ntot_pdi<25)=NaN;
            end
            
            for ifile=1:length(files)
                Ps = gen(ifile).s_ps;
                T = gen(ifile).s_ta;
                rv_m = gen(ifile).s_mr/1000;
                rv_n = rv_m*mm_d/mm_v;
                fv = rv_n./(rv_n+1);
                fd = 1-fv;
                R0 = 8.314;
                R_tot = R0./(fv*mm_v + fd*mm_d);
                gen(ifile).s_rhoa = Ps*100./(R_tot.*(T+273.15));
            end
            
        case 'postpdi'
            bindata=load('dp_pdi.mat');
            binlims=bindata.pdi_bins(:,1);
            binlims(end+1)=bindata.pdi_bins(end,2);
            dlogD = log10(binlims(2)/binlims(1));
            filedir='post/';
            files=dir([filedir,'*.mat']);
            for ifile=1:length(files)
                tempdata=load([filedir,files(ifile).name]);
                try 
                    tempdata=rmfield(tempdata,'s_mr2'); %This is an extra field not found in the other campaigns
                catch
                end
                gen(ifile)=tempdata;
                gen(ifile).s_disp_pdi(gen(ifile).s_disp_pdi==0)=NaN;
                gen(ifile).s_disp_pdi(gen(ifile).s_ntot_pdi<25)=NaN;
            end
    end
    
    % filter statistically insignificant data
    
    binmean=(binlims(1:end-1)+binlims(2:end))/2; %Bin center diameter
    %%
    for ifile=1:length(files)
        %Calculate RH
        gen(ifile).s_rh = gen(ifile).s_mr/1000./wsat_tp(gen(ifile).s_ta,gen(ifile).s_ps)*100; 
        %Calculate normalized droplet concentration as the droplet concentration / max droplet concentration
        gen(ifile).s_normN = gen(ifile).(Nvar)/max(gen(ifile).(Nvar)); 
        ql = gen(ifile).s_lwc_pdi./gen(ifile).s_rhoa;
        qv = gen(ifile).s_mr/1000./(1+gen(ifile).s_mr/1000)*1000;
        gen(ifile).s_qt = ql + qv;
        
        %calculate the normalized altitude in cloud
        z = gen(ifile).s_ap;
        normAC = zeros(size(z));
        s_t = gen(ifile).s_t;
        ti = fb.(fbvar)(ifile).ti;
        tf = fb.(fbvar)(ifile).tf;
        
        nleg=length(fb.(fbvar)(ifile).ti);
        for jleg = 1:nleg
            z_CB = fb.(fbvar)(ifile).z_CB(jleg);
            z_CT = fb.(fbvar)(ifile).z_CT(jleg);
            
            if jleg == 1
                be4_1st_cld = s_t<ti(jleg);
                normAC(s_t<ti(jleg)) = (z(be4_1st_cld)-z_CB)/(z_CT-z_CB);
            end
            
            incloud_idx = s_t>=ti(jleg)&s_t<=tf(jleg);
            normAC(incloud_idx) = (z(incloud_idx)-z_CB)/(z_CT-z_CB);
            % set regions before the first cloud and after the last
            % cloud to nan
%             normAC(t<ti(1) | t>tf(end))=nan; 

            % set regions between legs to nan
            if jleg > 1
                prev_intv_idx = s_t>tf(jleg-1) & s_t<ti(jleg);
                normAC(prev_intv_idx)= (z(prev_intv_idx)-z_CB)/(z_CT-z_CB);
            end
        end
        
        gen(ifile).normAC = normAC;
        
        for jleg = 1:nleg
            z_CB = fb.(fbvar)(ifile).z_CB(jleg);
            z_CT = fb.(fbvar)(ifile).z_CT(jleg);
            
            tleg_idx = s_t>=ti(jleg) & s_t<=tf(jleg);
            T_CB = nanmean(gen(ifile).s_ta(normAC(tleg_idx)<0.05 & normAC(tleg_idx)>-0.05))+273.15;
            r_CB = nanmean(gen(ifile).s_mr(normAC(tleg_idx)<0.05 & normAC(tleg_idx)>-0.05))/1000;
            p_CB = nanmean(gen(ifile).s_ps(normAC(tleg_idx)<0.05 & normAC(tleg_idx)>-0.05));
            
            
            % adiabatic q_l as a function of z
            ql_adb_lin = adiab_ql(z_CB,z_CT,T_CB,r_CB,p_CB)*1000; 
            z_lin = linspace(z_CB,z_CT,length(ql_adb_lin));
            % observed q_l as a function of time
            filtered_idx = gen(ifile).s_lwc_pdi<0.01;

            s_lwc_pdi = gen(ifile).s_lwc_pdi; %threshold for LWC
            rhoa = gen(ifile).s_rhoa;

            ql_obs = s_lwc_pdi./rhoa;
            % since different legs might have different cloud base
            % altitude
            if nleg == 1
                % convert the linear adiabatic q_l to the profile as a
                % function time
                ql_adb_prof = interp1(z_lin,ql_adb_lin,z);
            elseif nleg > 1

                if jleg == 1
                    ql_adb_prof = interp1(z_lin,ql_adb_lin,z(s_t<=tf(jleg)));
                end
                if jleg > 1 && jleg < nleg
                    ql_adb_prof = [ql_adb_prof;interp1(z_lin,...
                        ql_adb_lin,z(s_t>tf(jleg-1) & s_t<=tf(jleg)))];
                end
                if jleg == nleg
                    ql_adb_prof = [ql_adb_prof;interp1(z_lin,...
                        ql_adb_lin,z(s_t>tf(jleg-1)))];
                end
            end
        end
        
        
        
        
%         return
        
        % set the mean(a_ntot) below lb_cb and w/o LWC to be the
        % background a_ntot

        lb_cb = mean(fb.(fbvar)(ifile).z_CT);

        cloudlegs_i = fb.(fbvar)(ifile).ti;
        cloudlegs_f = fb.(fbvar)(ifile).tf;

        s_t = clouds.(camp)(ifile).s_t;
        a_t = clouds.(camp)(ifile).a_t;

        allinst_commontime = mintersect(floor(s_t), a_t);

        cmt_ipdi = find(ismember(floor(s_t), allinst_commontime));
        cmt_ipcasp = find(ismember(a_t, allinst_commontime));

%             cm_normAC = gen(ifile).normAC(idxpdi); % where normAC and PCASP data are measured at the same time
        cm_z = gen(ifile).s_ap(cmt_ipdi);
        cm_t = clouds.(camp)(ifile).a_t(cmt_ipcasp);
        cm_lwc = gen(ifile).s_lwc_pdi(cmt_ipdi);
        cm_a_ntot = clouds.(camp)(ifile).a_ntot(cmt_ipcasp);
        cm_s_ntot = clouds.(camp)(ifile).s_ntot_pdi(cmt_ipdi);
        cm_s_disp_pdi = gen(ifile).s_disp_pdi(cmt_ipdi);

        if ~isempty(cloudlegs_i)
            for ileg = 1:length(cloudlegs_i)
                ti = cloudlegs_i(ileg);
                tf = cloudlegs_f(ileg);
                ti_idx = findInSorted(cm_t, ti);
                tf_idx = findInSorted(cm_t, tf);

                if ti_idx<0
                    continue
                end
                
                % ti and tf used for cloud related properties (reldisp, s_ntot, actfrac) 
                % does not sample region below the cloud base
                ti_c = ti; 
                tf_c = tf;
                
                % to sample some extra distance below the cloud top, in case
                % the flight doesnt only go from low to high
                if cm_z(ti_idx) < cm_z(tf_idx)
                    ti = ti - 300;
                else
                    tf = tf + 300;
                end

                z_min = min(cm_z(cm_t < tf & cm_t > ti));
                z_max = max(cm_z(cm_t < tf & cm_t > ti));
                z_max_sampled = (z_min+z_max)/2;
   
                aerCMS = @(x) calcMeanSampsize(x, cm_t < tf & cm_t > ti & cm_z < z_max_sampled);
                cldCMS = @(x) calcMeanSampsize(x, cm_t < tf_c & cm_t > ti_c & cm_z < z_max_sampled & cm_s_ntot > 25);
                
                try
%                         drpCMS = @(x) calcMeanSampsize(x, cm_t < tf & cm_t > ti & gen(ifile).s_ap < z_max_sampled);

                    [gen(ifile).a_ntot_CB(ileg), gen(ifile).a_ntot_CB_sampsize(ileg)] = ...
                        aerCMS(cm_a_ntot);
                    [gen(ifile).s_ntot_CB(ileg), gen(ifile).s_ntot_CB_sampsize(ileg)] = ...
                        cldCMS(cm_s_ntot);
                    [gen(ifile).s_actfrac_CB(ileg), gen(ifile).s_actfrac_CB_sampsize(ileg)] = ...
                        cldCMS(cm_s_ntot./cm_a_ntot);
                    [gen(ifile).reldisp_CB(ileg), gen(ifile).reldisp_CB_sampsize(ileg)] = ...
                        cldCMS(cm_s_disp_pdi);
                catch
                    [gen(ifile).a_ntot_CB(ileg), gen(ifile).a_ntot_CB_sampsize(ileg), ...
                        gen(ifile).s_ntot_CB(ileg), gen(ifile).s_ntot_CB_sampsize(ileg), ...
                        gen(ifile).s_actfrac_CB(ileg), gen(ifile).s_actfrac_CB_sampsize(ileg), ...
                        gen(ifile).reldisp_CB(ileg), gen(ifile).reldisp_CB_sampsize(ileg)] = ...
                        deal(nan);
                end
            end
        else
            [gen(ifile).a_ntot_CB, gen(ifile).a_ntot_CB_sampsize, ...
                gen(ifile).s_ntot_CB, gen(ifile).s_ntot_CB_sampsize, ...
                gen(ifile).s_actfrac_CB, gen(ifile).s_actfrac_CB_sampsize, ...
                gen(ifile).reldisp_CB, gen(ifile).reldisp_CB_sampsize] = ...
                deal(nan);
        end
        
        
        
        
        
        
        
%         end
        
%         %% ------------------------ test -------------------------
%         if do_test == true
%             figure
%             plot(cm_t, cm_a_ntot)
%             hold on
%             plot(cm_t, (cm_z<lb_cb & cm_lwc<0.01)*1000)
%             yyaxis right
%             plot(cm_t, cm_z)
%         end
%         %% ---------------------- test end -----------------------
            
        AF = ql_obs./ql_adb_prof;
        AF(filtered_idx) = nan;
        AF(AF>2)=nan;
% %         AF = AF/max(AF); %normalize in case the cloud base was not 
%         % determined clearly
% 
        gen(ifile).ql_adb_prof = ql_adb_prof;
        gen(ifile).AF = AF;
        
        % calculate the mean droplet diameter
        gen(ifile).s_meand_pdi = sum(dlogD*gen(ifile).(nvar).*binmean',2)./gen(ifile).(Nvar);
        
        % calculate the deviation from gamma distribution
        a = 2*gen(ifile).s_reff_pdi./gen(ifile).s_meand_pdi;
        gen(ifile).calc_reldisp = real(sqrt((a-1)/2));
        gen(ifile).dev_gamma = abs(gen(ifile).calc_reldisp ./ gen(ifile).s_disp_pdi - 1);
        
    end
    
    fields=fieldnames(gen);
    n=0;
    
    for ifile=1:length(gen)
        display(ifile)
        %Find indices with non-NaN relative dispersion
        inds = find(~isnan(gen(ifile).(dvar)));
%         if length(inds)>100 %If the flight has at least 100 non-NaN reldisp values
            n=n+1;
            %Put the data for this flight into the master data structure 'clouds'
            for jfield=1:length(fields)
                try
                    clouds.(camp)(n).(fields{jfield})=gen(ifile).(fields{jfield})(:,:);
                catch
                    clouds.(camp)(n).(fields{jfield})=gen(ifile).(fields{jfield})(1,:);
                end
            end
            %Put in some extra info
            clouds.(camp)(n).file=[filedir,files(ifile).name];
            clouds.(camp)(n).filenum = ifile;
            clouds.(camp)(n).binlims = binlims;
            clouds.(camp)(n).binmean = binmean;
            clouds.(camp)(n).samplesize = length(inds);
            clouds.(camp)(n).maxN = max(gen(ifile).(Nvar));
%             clouds.(camp)(n).s_meand_pdi = ...
%                 sum(dlogD*clouds.(camp)(n).(nvar).*binmean',2)./clouds.(camp)(n).(Nvar);
%             clouds.(campaign)(n).M1 = ...
%                 sum(dlogD*clouds.(campaign)(n).(nvar).*binmean',2);
%             clouds.(campaign)(n).M1_5 = ...
%                 sum(dlogD*clouds.(campaign)(n).(nvar).*binmean'.^1.5,2);
%             clouds.(campaign)(n).M2 = ...
%                 sum(dlogD*clouds.(campaign)(n).(nvar).*binmean'.^2,2);
%             clouds.(campaign)(n).M3 = ...
%                 sum(dlogD*clouds.(campaign)(n).(nvar).*binmean'.^3,2);
%             clouds.(campaign)(n).M6 = ...
%                 sum(dlogD*clouds.(campaign)(n).(nvar).*binmean'.^6,2);
    end
end

% finishingTaskSound

% %% SEE boundary layers
% close all
% 
% campaigns={'vocalspdi','masepdi','postpdi'};
% c=3;
% camp = campaigns{c};
% 
% for iday = 1:length(clouds.(camp))
%     
%     s_t = clouds.(camp)(iday).s_t;
%     a_t = clouds.(camp)(iday).a_t;
%     
%     [pdi_pcasp_commontime, pdi_pcasp_commontime_ipdi, pdi_pcasp_commontime_ipcasp] = ...
%         intersect(floor(s_t), a_t);
%     
%     figure('Position',[1370 336 744 649])
%     
% %     subplot(3,1,1)
%     a_t = clouds.(camp)(iday).a_t;
%     s_thet = clouds.(camp)(iday).s_thet;
%     s_ap = clouds.(camp)(iday).s_ap;
%     a_ntot = clouds.(camp)(iday).a_ntot;
%     s_lwc_pdi = clouds.(camp)(iday).s_lwc_pdi;
%     line(a_ntot(pdi_pcasp_commontime_ipcasp), s_ap(pdi_pcasp_commontime_ipdi),...
%         'linestyle','none','marker','.','color',[0 0.4470 0.7410]);
%     xlabel('aerosol # conc cc^{-1}')
% %     xlim([0 max(ylim)])
%     ax1 = gca; % current axes
%     ax1.XColor = [0 0.4470 0.7410];
%     ax1.YColor = [0 0.4470 0.7410];
%     set(gca,'fontsize',18)
%     ylim([0 max(s_ap)])
%     ax1_pos = ax1.Position;
%     ax2 = axes('Position',ax1_pos,...
%         'XAxisLocation','top',...
%         'YAxisLocation','right',...
%         'Color','none');
%     line(s_thet(s_thet>0), s_ap(s_thet>0),...
%         'linestyle','none','marker','.','color','r');
%     ax2.XColor = 'r';
%     ax2.YColor = 'r';
%     xlabel('Potential temperature [K]')
%     ylim([0 max(s_ap)])
%     set(gca,'fontsize',18)
% 
% %     plot(hskp_t, hskp_z)
% %     yyaxis right
% %     plot(hskp_t, s_thet)
%     
% end

%% Save the data
save('clouds.mat','clouds', '-v7.3')
% finishingTaskSound