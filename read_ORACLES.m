clear
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
load clouds.mat
% cloud_backup=load('clouds_backup.mat');
%%
pdi_files = dir('ORACLES_raw/PDI_P3/*.ict');
pdi_dates = arrayfun(@(x) pdi_files(x).name(1:6), 1:length(pdi_files),...
    'UniformOutput', false)';

binedge = readmatrix('ORACLES_raw/PDI_P3/bin_edges.csv','Range',[1 2]);
binmean = (binedge(:,2)+binedge(:,1))/2;
dlogD = log10(binedge(1,2)/binedge(1,1));

pcasp_files = dir('ORACLES_raw/PCASP_P3/*.ict');
pcasp_dates = arrayfun(@(x) pcasp_files(x).name(1:6), 1:length(pcasp_files),...
    'UniformOutput', false)';

% fb_dates = arrayfun(@(x) ORACLES_flight_basics(x).date_str,1:length(ORACLES_flight_basics),...
%     'UniformOutput',false)';

hskp_files = dir('ORACLES_raw/Hskping_P3/*.ict');
hskp_dates = arrayfun(@(x) hskp_files(x).name(14:19), 1:length(hskp_files),...
    'UniformOutput', false)';

allinst_commondate = mintersect(pdi_dates, pcasp_dates, hskp_dates);
cmd_ipdi = find(ismember(pdi_dates, allinst_commondate));
cmd_ipcasp = find(ismember(pcasp_dates, allinst_commondate));
cmd_ihskp = find(ismember(hskp_dates, allinst_commondate));
% cmd_idxcoma = find(ismember(coma_dates, allinst_commondate));

% finishingTaskSound

%% 

Dbins_100_grtr = find(binmean>100,1,'first');
% 
% try
%     clouds = rmfield(clouds,'oraclespdi');
% catch
% end

for iday = 1:length(allinst_commondate)
    ipdi = cmd_ipdi(iday);
    ipcasp = cmd_ipcasp(iday);
    ihskp = cmd_ihskp(iday);
%     icoma = cmd_idxcoma(iday);
    

    %% reading pdi
    conc{iday,1}(:,:) = readmatrix([pdi_files(ipdi).folder '/' pdi_files(ipdi).name],...
        'FileType','text','NumHeaderLines',160);
    
    s_t = conc{iday}(:,1);
    s_conc_pdi = conc{iday}(:,2:end);
    s_ntot_pdi = sum(s_conc_pdi*dlogD,2);
    s_ntot_100_grtr_pdi = sum(s_conc_pdi(:,Dbins_100_grtr:end)*dlogD,2);
    s_meand_pdi = sum(s_conc_pdi.*binmean'*dlogD,2)./s_ntot_pdi;
    s_lwc_pdi = sum(s_conc_pdi.*binmean'.^3*dlogD,2)*pi/6*1e-6;
    s_std_pdi = arrayfun(@(x) std(binmean, s_conc_pdi(x,:)), 1:length(s_t))';
    s_disp_pdi = s_std_pdi./s_meand_pdi;
    
    % invalid index
    ividx = s_ntot_pdi<25 | s_meand_pdi < 3 | s_ntot_100_grtr_pdi>0; 
    % set the lower bound for datapoints statistical properties
    s_std_pdi(ividx) = nan;
    s_disp_pdi(ividx) = nan;
    
    
    clouds.oraclespdi(iday).file = pdi_files(iday).name;
    clouds.oraclespdi(iday).s_t = s_t;
    clouds.oraclespdi(iday).s_ntot_pdi = s_ntot_pdi;
%     clouds.oraclespdi(iday).s_ntot_100_grtr_pdi = s_ntot_100_grtr_pdi;
    clouds.oraclespdi(iday).s_conc_pdi = s_conc_pdi;
    clouds.oraclespdi(iday).s_meand_pdi = s_meand_pdi;
    clouds.oraclespdi(iday).s_lwc_pdi = s_lwc_pdi;
    clouds.oraclespdi(iday).s_std_pdi = s_std_pdi;
    clouds.oraclespdi(iday).s_disp_pdi = s_disp_pdi;
    
    
    %% reading pcasp
    aeros{iday,1}(:,:) = readmatrix([pcasp_files(ipcasp).folder '/' pcasp_files(ipcasp).name],...
        'FileType','text','Range',[71 1]);
    
    vidx = aeros{iday}(:,32)>0;
    
    a_t_raw = aeros{iday}(vidx,1);
    a_conc = aeros{iday}(vidx,2:31);
%     a_ntot = aeros{iday}(vidx,32);
%     pause
    a_ntot = sum(a_conc,2);
    a_meanD_ex = aeros{iday}(vidx,36);
%     a_supto = aeros{iday}(vidx,33);
%     a_totin = aeros{iday}(vidx,37);
%     a_timvc = aeros{iday}(vidx,38);
    
    [pdi_pcasp_cmt, pdi_pcasp_cmt_ipdi, pdi_pcasp_cmt_ipcasp] = ...
        intersect(s_t,a_t_raw);
    
    s_ntot_pdi(isnan(s_ntot_pdi)) = 0;
    
    clouds.oraclespdi(iday).a_date_str = pcasp_dates{ipcasp};
%     clouds.oraclespdi(iday).a_t = a_t;
    
%     clouds.oraclespdi(iday).a_conc = a_conc;
%     
%     clouds.oraclespdi(iday).a_ntot_ex = a_ntot;
%     a_ntot(pdi_pcasp_cmt_idx_in_pcasp) = a_ntot(pdi_pcasp_cmt_idx_in_pcasp) + ...
%         s_ntot_pdi(pdi_pcasp_cmt_idx_in_pdi);
%     clouds.oraclespdi(iday).a_ntot = a_ntot; % real a_ntot, including the clouds droplets
%     clouds.oraclespdi(iday).a_supto = a_supto;
%     clouds.oraclespdi(iday).a_totin = a_totin;
%     clouds.oraclespdi(iday).a_timvc = a_timvc;
    
    %% reading hskp
    hskp_header_line = csvread([hskp_files(ihskp).folder '/' hskp_files(ihskp).name],...
        0, 0, [0 0 0 0]);
    hskp{iday,1}(:,:) = readmatrix([hskp_files(ihskp).folder '/' hskp_files(ihskp).name],...
        'FileType','text','Range',[hskp_header_line 1]);
    hskp_t_raw = hskp{iday}(:,1);
    a_rh = hskp{iday}(:,end);
    a_rh(a_rh<0) = nan;
    
%     fb_t = ORACLES_flight_basics(ihskp).t_1Hz;
%     fb_z = ORACLES_flight_basics(ihskp).z_1Hz;
%     fb_T = ORACLES_flight_basics(ihskp).T_1Hz;
%     
    [hskp_pcasp_cmt, hskp_pcasp_cmt_ihskp, hskp_pcasp_cmt_ipcasp]=...
        intersect(hskp_t_raw, a_t_raw);
    
    %% saving pcasp and fb
    a_t = a_t_raw(hskp_pcasp_cmt_ipcasp);
    a_conc = a_conc(hskp_pcasp_cmt_ipcasp,:);
    a_ntot_ex = a_ntot(hskp_pcasp_cmt_ipcasp);
    a_meanD_ex = a_meanD_ex(hskp_pcasp_cmt_ipcasp);
    
    a_ntot(pdi_pcasp_cmt_ipcasp) = a_ntot(pdi_pcasp_cmt_ipcasp) + ...
        s_ntot_pdi(pdi_pcasp_cmt_ipdi);
    
    s_actfrac = nan(size(s_ntot_pdi));
    s_actfrac(pdi_pcasp_cmt_ipdi) = s_ntot_pdi(pdi_pcasp_cmt_ipdi)./a_ntot(pdi_pcasp_cmt_ipcasp);
    
    a_ntot = a_ntot(hskp_pcasp_cmt_ipcasp); % real a_ntot, including the clouds droplets
    
    
%     a_supto = a_supto(fb_pcasp_commontime_idxpcasp);
%     a_totin = a_totin(fb_pcasp_commontime_idxpcasp);
%     a_timvc = a_timvc(fb_pcasp_commontime_idxpcasp);
    
    
    a_t = hskp_t_raw(hskp_pcasp_cmt_ihskp);
    a_rh = a_rh(hskp_pcasp_cmt_ihskp);
    a_z = hskp{iday}(hskp_pcasp_cmt_ihskp,6)*0.3048;
    a_wz = hskp{iday}(hskp_pcasp_cmt_ihskp,12);
    hskp_thet = hskp{iday}(hskp_pcasp_cmt_ihskp,19);
    hskp_ta = hskp{iday}(hskp_pcasp_cmt_ihskp,18);
    hskp_p = hskp{iday}(hskp_pcasp_cmt_ihskp,23);
    
    %calculating air temperature using ideal gas law
    hskp_pres_raw = hskp{iday}(:,23)*100;
    hskp_ta_raw = hskp{iday}(:,18)+273.15;
    Rd = 287;
    hskp_rhoa = hskp_pres_raw./(Rd*hskp_ta_raw);
    
    hskp_mr_raw = hskp{iday}(:,31);
    hskp_qv = hskp_mr_raw/1000./(1+hskp_mr_raw/1000)*1000;
    
    
    clouds.oraclespdi(iday).a_t = a_t;
    clouds.oraclespdi(iday).a_conc = a_conc;
    clouds.oraclespdi(iday).a_ntot_ex = a_ntot_ex;
    clouds.oraclespdi(iday).a_meanD_ex = a_meanD_ex;
    clouds.oraclespdi(iday).a_ntot = a_ntot;
    clouds.oraclespdi(iday).s_actfrac = s_actfrac;
    clouds.oraclespdi(iday).a_z = a_z;
    clouds.oraclespdi(iday).a_wz = a_wz;
    clouds.oraclespdi(iday).a_thet = hskp_thet;
    clouds.oraclespdi(iday).a_T = hskp_ta;
    clouds.oraclespdi(iday).a_rh = a_rh;
    clouds.oraclespdi(iday).a_p = hskp_p;
    
    %% calculating q_total
    [hskp_pdi_cmt, hskp_pdi_cmt_ihskp, hskp_pdi_cmt_ipdi] = intersect(hskp_t_raw, s_t);
    s_ql = zeros(size(s_lwc_pdi));
    s_ql(hskp_pdi_cmt_ipdi) = s_lwc_pdi(hskp_pdi_cmt_ipdi)./hskp_rhoa(hskp_pdi_cmt_ihskp);
    
    s_mr = zeros(size(s_lwc_pdi));
    s_mr(hskp_pdi_cmt_ipdi) = hskp_mr_raw(hskp_pdi_cmt_ihskp);
    s_qv = s_mr/1000./(1+s_mr/1000)*1000;
    
    s_qt = s_ql + s_qv;
    clouds.oraclespdi(iday).a_mr = hskp_mr_raw;
    clouds.oraclespdi(iday).s_qv = s_qv;
    clouds.oraclespdi(iday).s_qt = s_qt;
    clouds.oraclespdi(iday).s_ql = s_ql;
    clouds.oraclespdi(iday).cmt = mintersect(a_t,s_t,hskp_t_raw);
end

% finishingTaskSound

% %% visualizing where the drizzle is
% for iday = 1:length(allinst_commondate)
%     figure
%     plot(clouds.oraclespdi(iday).s_lwc_pdi)
% end

% %% trying to determine where the BL is
% close all
% % % [fb_pdi_commondate, fb_pdi_commondate_idxpdi, fb_pdi_commondate_idxfb]=...
% % %     intersect(pdi_dates, fb_dates, 'sorted');
% % 
% for iday = 1:length(allinst_commondate)
%     
%     s_t = clouds.oraclespdi(iday).s_t;
%     hskp_t = clouds.oraclespdi(iday).hskp_t;
%     
%     [pdi_pf_commontime, pdi_pf_commontime_ipdi, pdi_pf_commontime_ipf] = ...
%         intersect(s_t, hskp_t);
%     
%     figure('Position',[1370 336 744 649])
%     
% %     subplot(3,1,1)
%     hskp_t = clouds.oraclespdi(iday).hskp_t;
%     hskp_thet = clouds.oraclespdi(iday).hskp_thet;
%     hskp_z = clouds.oraclespdi(iday).hskp_z;
%     a_ntot = clouds.oraclespdi(iday).a_ntot;
%     lwc = clouds.oraclespdi(iday).s_lwc_pdi;
% %     line(lwc(pdi_pf_commontime_ipdi), hskp_z(pdi_pf_commontime_ipf),...
% %         'linestyle','none','marker','.','color',[0 0.4470 0.7410]);
% %     xlabel('aerosol # conc cc^{-1}')
% %     xlim([0 max(ylim)])
% %     ax1 = gca; % current axes
% %     ax1.XColor = [0 0.4470 0.7410];
% %     ax1.YColor = [0 0.4470 0.7410];
% %     set(gca,'fontsize',18)
% %     ylim([0 max(hskp_z)])
% %     ax1_pos = ax1.Position;
% %     ax2 = axes('Position',ax1_pos,...
% %         'XAxisLocation','top',...
% %         'YAxisLocation','right',...
% %         'Color','none');
% %     line(hskp_thet(hskp_thet>0), hskp_z(hskp_thet>0),...
% %         'linestyle','none','marker','.','color','r');
% %     ax2.XColor = 'r';
% %     ax2.YColor = 'r';
% %     xlabel('Potential temperature [K]')
% %     ylim([0 max(hskp_z)])
%     set(gca,'fontsize',18)
% % 	subplot(3,1,1)
% %     plot(hskp_t, hskp_z)
% %     yyaxis right
% %     plot(hskp_t, hskp_thet,'LineWidth',2)
% %     ylabel('potential')
%     
%     subplot(2,1,1)
%     plot(hskp_t, hskp_z)
%     yyaxis right
%     plot(hskp_t, a_ntot,'LineWidth',2)
%     ylabel('aerosol conc')
%     
%     subplot(2,1,2)
%     plot(hskp_t, hskp_z)
%     yyaxis right
%     plot(hskp_t(pdi_pf_commontime_ipf), lwc(pdi_pf_commontime_ipdi),'LineWidth',2)
%     ylabel('LWC')
%     
% end
%%
% 
load oraclespdi_flight_basics.mat



for iday = 1:length(allinst_commondate)
%     clouds.oraclespdi(iday).hskp_blt = hskp_blt(iday);
%     clouds.oraclespdi(iday).hskp_cloudlegs_i = hskp_cloudlegs_i{iday};
%     clouds.oraclespdi(iday).hskp_cloudlegs_f = hskp_cloudlegs_f{iday};
    
    cloudlegs_i = oraclespdi_flight_basics(iday).ti;
    cloudlegs_f = oraclespdi_flight_basics(iday).tf;
    z_CB = oraclespdi_flight_basics(iday).z_CB;
    z_CT = oraclespdi_flight_basics(iday).z_CT;
    
    s_t = clouds.oraclespdi(iday).s_t;
    s_lwc_pdi = clouds.oraclespdi(iday).s_lwc_pdi;
    a_ntot = clouds.oraclespdi(iday).a_ntot;
    a_t = clouds.oraclespdi(iday).a_t;
    a_z = clouds.oraclespdi(iday).a_z;
    s_ntot = clouds.oraclespdi(iday).s_ntot_pdi;
    s_disp_pdi = clouds.oraclespdi(iday).s_disp_pdi;
    
    allinst_commontime = mintersect(s_t, a_t);
    cmt_ipdi = find(ismember(s_t, allinst_commontime));
    cmt_ipf = find(ismember(a_t, allinst_commontime));
    
    cm_t = s_t(cmt_ipdi);
    cm_z = a_z(cmt_ipf);
    cm_lwc = s_lwc_pdi(cmt_ipdi);
    cm_a_ntot = a_ntot(cmt_ipf);
    cm_s_ntot = s_ntot(cmt_ipdi);
    cm_disp = s_disp_pdi(cmt_ipdi);
    
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
                ti = ti - 600;
            else
                tf = tf + 600;
            end
            
            z_min = min(cm_z(cm_t < tf & cm_t > ti));
            z_max = max(cm_z(cm_t < tf & cm_t > ti));
            z_max_sampled = (z_min+z_max)/2; 
            
            % only sample the datapoints lower than mid cloud to prevent
            % sampling collision-coalescence
            
            aerCMS = @(x) calcMeanSampsize(x, cm_t < tf & cm_t > ti & cm_z < z_max_sampled);
            cldCMS = @(x) calcMeanSampsize(x, cm_t < tf_c & cm_t > ti_c & cm_z < z_max_sampled & cm_s_ntot > 25);
            
            try
                [clouds.oraclespdi(iday).a_ntot_CB(ileg), clouds.oraclespdi(iday).a_ntot_CB_sampsize(ileg)] = ...
                    aerCMS(cm_a_ntot);
                [clouds.oraclespdi(iday).s_ntot_CB(ileg), clouds.oraclespdi(iday).s_ntot_CB_sampsize(ileg)] = ...
                    cldCMS(cm_s_ntot);
                [clouds.oraclespdi(iday).s_actfrac_CB(ileg), clouds.oraclespdi(iday).s_actfrac_CB_sampsize(ileg)] = ...
                    cldCMS(cm_s_ntot./cm_a_ntot);
                [clouds.oraclespdi(iday).reldisp_CB(ileg), clouds.oraclespdi(iday).reldisp_CB_sampsize(ileg)] = ...
                    cldCMS(cm_disp);
            
            catch
                [clouds.oraclespdi(iday).a_ntot_CB(ileg), clouds.oraclespdi(iday).a_ntot_CB_sampsize(ileg), ...
                    clouds.oraclespdi(iday).s_ntot_CB(ileg), clouds.oraclespdi(iday).s_ntot_CB_sampsize(ileg), ...
                    clouds.oraclespdi(iday).s_actfrac_CB(ileg), clouds.oraclespdi(iday).s_actfrac_CB_sampsize(ileg), ...
                    clouds.oraclespdi(iday).reldisp_CB(ileg), clouds.oraclespdi(iday).reldisp_CB_sampsize(ileg)] = ...
                    deal(nan);
            end
        end
    else
        [clouds.oraclespdi(iday).a_ntot_CB, clouds.oraclespdi(iday).a_ntot_CB_sampsize, ...
            clouds.oraclespdi(iday).s_ntot_CB, clouds.oraclespdi(iday).s_ntot_CB_sampsize, ...
            clouds.oraclespdi(iday).s_actfrac_CB, clouds.oraclespdi(iday).s_actfrac_CB_sampsize, ...
            clouds.oraclespdi(iday).reldisp_CB, clouds.oraclespdi(iday).reldisp_CB_sampsize] = ...
            deal(nan);
    end
    
end

%% Save the data
save('clouds.mat','clouds', '-v7.3')

% finishingTaskSound