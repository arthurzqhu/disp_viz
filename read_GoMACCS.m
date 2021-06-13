clear
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
load clouds.mat

% tic 

close all
pdi_files = dir('GoMACCS_raw/06*/whole*');
cabin_files = dir('GoMACCS_raw/06*/TO_CABIN*');
binedges = readmatrix('GoMACCS_raw/bin_edges.txt');
binmean = (binedges(1:end-1) + binedges(2:end))/2;
dlogD = log10(binedges(2)/binedges(1));

pdi_dates = arrayfun(@(x) pdi_files(x).name(25:end-4), 1:length(pdi_files),...
    'UniformOutput', false)';
cabin_dates = arrayfun(@(x) cabin_files(x).name(12:end-4), 1:length(cabin_files),...
    'UniformOutput', false)';

% finishingTaskSound

%%

% try
%     clouds = rmfield(clouds,'gomaccspdi');
% catch
% end

for iday = 1:length(pdi_files)
    conc_dat{iday,1}(:,:) = readmatrix([pdi_files(iday).folder '/' pdi_files(iday).name],...
        'FileType','text');
    
    s_t = conc_dat{iday}(:,1);
    s_conc_pdi = conc_dat{iday}(:,2:end);
    s_ntot_pdi = sum(s_conc_pdi*dlogD, 2);
    s_meand_pdi = sum(s_conc_pdi.*binmean'*dlogD,2)./s_ntot_pdi;
    s_lwc_pdi = sum(s_conc_pdi.*binmean'.^3*dlogD,2)*pi/6*1e-6;
    s_std_pdi = arrayfun(@(x) std(binmean, s_conc_pdi(x,:)), 1:length(s_t))';
    s_disp_pdi = s_std_pdi./s_meand_pdi;
    
    % invalid index
    ividx = s_ntot_pdi<25 | s_meand_pdi < 3; 
    % set the lower bound for datapoints statistical properties
    s_std_pdi(ividx) = nan;
    s_disp_pdi(ividx) = nan;
    
    clouds.gomaccspdi(iday).file = pdi_files(iday).name(25:end-4);
    clouds.gomaccspdi(iday).s_t = s_t;
    clouds.gomaccspdi(iday).s_ntot_pdi = s_ntot_pdi;
    clouds.gomaccspdi(iday).s_conc_pdi = s_conc_pdi;
    clouds.gomaccspdi(iday).s_meand_pdi = s_meand_pdi;
    clouds.gomaccspdi(iday).s_lwc_pdi = s_lwc_pdi;
    clouds.gomaccspdi(iday).s_std_pdi = s_std_pdi;
    clouds.gomaccspdi(iday).s_disp_pdi = s_disp_pdi;
    %% read cabin
    
    cabin{iday,1}(:,:) = readmatrix([cabin_files(iday).folder '/' cabin_files(iday).name],...
        'FileType','text','Range',[2 1]);
    
    a_t = cabin{iday}(:,1);
%     aeros_conc = cabin{iday}(:,2:31);
    a_ntot = cabin{iday}(:,22);
    
    vidx = a_ntot>0;
    
    a_t = a_t(vidx,:);
    a_ntot = a_ntot(vidx,:);
    a_lat = cabin{iday}(vidx,3); a_lat(a_lat==0)=nan;
    a_lon = -cabin{iday}(vidx,4); a_lon(a_lon==0)=nan;
    
    % fix the nan if it's during the sampling
    nan_tlist=find(isnan(a_lon));
    
    for inan=1:sum(isnan(a_lon))
        prev_idx=nan_tlist(inan)-1;
        next_idx=nan_tlist(inan)+1;
        
        a_lat(nan_tlist(inan))=(a_lat(prev_idx)+a_lat(next_idx))/2;
        a_lon(nan_tlist(inan))=(a_lon(prev_idx)+a_lon(next_idx))/2;
    end
    
    
    a_Vtot_ex = cabin{iday}(vidx,23);
    a_wz = cabin{iday}(vidx,19);
    a_z = cabin{iday}(vidx,16);
    a_p = cabin{iday}(vidx,15);
    a_T = cabin{iday}(vidx,12);
    a_rh = cabin{iday}(vidx,14);
    a_lwc_gerb = cabin{iday}(vidx,33);
    a_lwc_wire = cabin{iday}(vidx,34);
    a_T(a_T==0)=nan;
    
    a_thet = (a_T+273.15).*(1000./a_p).^(287/1004);
    
    a_meanV_ex = a_Vtot_ex./a_ntot;
    a_meanD_ex = (a_meanV_ex/(pi/6)).^(1/3);
    
    [~, cmt_ipdi, cmt_icabin] = intersect(s_t,...
        a_t);
    
    %% calculating s_qt
    Lv = 2.5e6;
    Rv = 462;
    Rd = 287;
    
    a_Tk = a_T+273.15;
    a_es = 6.11*exp(Lv/Rv*(1/273.15-1./a_Tk));
    a_e = a_es.*a_rh/100;
    a_mr = a_e*Rd./(Rv*(a_p-a_e));
    a_qv = a_mr./(1+a_mr)*1000;
    
    a_rhoa = a_p*100./(Rd*a_Tk);
    
    s_ql = zeros(size(s_lwc_pdi));
    s_ql(cmt_ipdi) = s_lwc_pdi(cmt_ipdi)./a_rhoa(cmt_icabin);
    s_qv = zeros(size(s_lwc_pdi));
    s_qv(cmt_ipdi) = a_qv(cmt_icabin);
    s_qt = s_ql + s_qv;
    
    
    clouds.gomaccspdi(iday).s_qt = s_qt;
%     cm_z = a_zp(commontime_icabin);
%     cm_lwc = s_lwc_pdi(commontime_ipdi);
%     cm_a_ntot = a_ntot(commontime_icabin);
%     a_ntot_CB = nanmean(cm_a_ntot(cm_z < 700 & cm_lwc<0.01));
%     a_ntot_CB_sampsize = length(find(cm_z < 700 & cm_lwc<0.01));
%     reldisp_CB = nanmean(s_disp_pdi);
%     reldisp_CB_sampsize = length(find(~isnan(s_disp_pdi)));
    
    
    
    clouds.gomaccspdi(iday).a_date = cabin_files(iday).name(12:end-4);
    clouds.gomaccspdi(iday).a_t = a_t;
    clouds.gomaccspdi(iday).a_lat = a_lat;
    clouds.gomaccspdi(iday).a_lon = a_lon;
    clouds.gomaccspdi(iday).a_mr = a_mr*1000;
    clouds.gomaccspdi(iday).a_ntot_ex = a_ntot;
    clouds.gomaccspdi(iday).a_meanD_ex = a_meanD_ex;
    clouds.gomaccspdi(iday).a_rhoa = a_rhoa;
    clouds.gomaccspdi(iday).a_z = a_z;
    clouds.gomaccspdi(iday).a_wz = a_wz;
    clouds.gomaccspdi(iday).a_p = a_p;
    clouds.gomaccspdi(iday).a_T = a_T;
    clouds.gomaccspdi(iday).a_thet = a_thet;
    clouds.gomaccspdi(iday).a_rh = a_rh;
    clouds.gomaccspdi(iday).a_lwc_gerb = a_lwc_gerb;
    clouds.gomaccspdi(iday).a_lwc_wire = a_lwc_wire;
    
    a_ntot(cmt_icabin) = a_ntot(cmt_icabin) + ...
        s_ntot_pdi(cmt_ipdi);
    clouds.gomaccspdi(iday).a_ntot = a_ntot;
    clouds.gomaccspdi(iday).s_ql = s_ql;
    clouds.gomaccspdi(iday).s_actfrac = nan(size(s_ntot_pdi)); 
    clouds.gomaccspdi(iday).s_actfrac(cmt_ipdi) = s_ntot_pdi(cmt_ipdi)./a_ntot(cmt_icabin);
%     clouds.gomaccspdi(iday).a_ntot_CB = a_ntot_CB;
%     clouds.gomaccspdi(iday).a_ntot_CB_sampsize = a_ntot_CB_sampsize;
%     clouds.gomaccspdi(iday).reldisp_CB = reldisp_CB;
%     clouds.gomaccspdi(iday).reldisp_CB_sampsize = reldisp_CB_sampsize;
end

% finishingTaskSound

% %%
% close all
% for iday = 1:length(cabin_files)
%     
%     [common_time, commontime_ipdi, commontime_icabin] = ...
%         intersect(clouds.gomaccspdi(iday).s_t, clouds.gomaccspdi(iday).a_t);
%     
%     figure('Position',[1721 1 860 984])
%     
%     a_t = clouds.gomaccspdi(iday).a_t;
%     s_t = clouds.gomaccspdi(iday).s_t;
%     a_thet = clouds.gomaccspdi(iday).a_thet;
%     a_z = clouds.gomaccspdi(iday).a_z;
%     a_ntot = clouds.gomaccspdi(iday).a_ntot;
%     lwc = clouds.gomaccspdi(iday).s_lwc_pdi;
% %     line(lwc(commontime_ipdi), a_z(commontime_icabin),...
% %         'linestyle','none','marker','.','color',[0 0.4470 0.7410]);
% % %     xlabel('aerosol # conc cc^{-1}')
% %     
% % %     xlim([0 max(ylim)])
% %     ax1 = gca; % current axes
% %     ax1.XColor = [0 0.4470 0.7410];
% %     ax1.YColor = [0 0.4470 0.7410];
% %     set(gca,'fontsize',18)
% %     ylim([0 max(a_z)])
% %     ax1_pos = ax1.Position;
% %     ax2 = axes('Position',ax1_pos,...
% %         'XAxisLocation','top',...
% %         'YAxisLocation','right',...
% %         'Color','none');
% %     line(a_thet(a_thet>0), a_z(a_thet>0),...
% %         'linestyle','none','marker','.','color','r');
% %     ax2.XColor = 'r';
% %     ax2.YColor = 'r';
% %     xlabel('Potential temperature [K]')
% %     ylim([0 max(a_z)])
% %     set(gca,'fontsize',18)
% 
%     subplot(2,1,1)
%     plot(a_t, a_z)
%     yyaxis right
%     plot(a_t, a_ntot,'LineWidth',2)
%     ylabel('aerosol conc')
%     
%     subplot(2,1,2)
%     plot(a_t, a_z)
%     yyaxis right
%     plot(a_t(commontime_icabin), lwc(commontime_ipdi),'LineWidth',2)
%     ylabel('LWC')
% end


%%

% hskp_blt = [nan;1000;1330;1200;2300;nan;nan;520;1500;2300;2300;nan;2500;3800;2950;2300;...
%     3000;2900;nan;2500;2600];
% 
% hskp_cloudlegs_i = {74415; 59682; 67378; 53860; 57542; []; []; [65936 79344]; 55511; ...
%     69167; [65068 69123]; []; [74414 77439 79142]; 68860; [69346 75085]; 66660;...
%     [57008 59597 66336]; [52476 57630 62704]; []; 63896; 60387};
% hskp_cloudlegs_f = {74486; 63781; 69559; 58636; 62235; []; []; [68127 81925]; 64369; ...
%     74719; [67858 71429]; []; [75230 78015 79575]; 73248; [73278 77021]; 72176;...
%     [58776 63727 67527]; [54585 61441 64731]; []; 66207; 64443};

load gomaccspdi_flight_basics.mat

for iday = 1:length(pdi_files)
    
%     clouds.gomaccspdi(iday).hskp_blt = hskp_blt(iday);
%     clouds.gomaccspdi(iday).hskp_cloudlegs_i = hskp_cloudlegs_i{iday};
%     clouds.gomaccspdi(iday).hskp_cloudlegs_f = hskp_cloudlegs_f{iday};
    
    cloudlegs_i = gomaccspdi_flight_basics(iday).ti;
    cloudlegs_f = gomaccspdi_flight_basics(iday).tf;
    
    s_t = clouds.gomaccspdi(iday).s_t;
    s_lwc_pdi = clouds.gomaccspdi(iday).s_lwc_pdi;
    a_ntot = clouds.gomaccspdi(iday).a_ntot;
    a_t = clouds.gomaccspdi(iday).a_t;
    a_z = clouds.gomaccspdi(iday).a_z;
    s_disp_pdi = clouds.gomaccspdi(iday).s_disp_pdi;
    
    [~, cmt_ipdi,cmt_icabin]=...
        intersect(s_t, a_t);
    
    cm_t = a_t(cmt_icabin);
    cm_z = a_z(cmt_icabin);
    cm_lwc = s_lwc_pdi(cmt_ipdi);
    cm_s_ntot = clouds.gomaccspdi(iday).s_ntot_pdi(cmt_ipdi);
    cm_a_ntot = a_ntot(cmt_icabin);
    cm_disp = s_disp_pdi(cmt_ipdi);
    
    clouds.gomaccspdi(iday).cmt = cm_t;
    
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
            
            [clouds.gomaccspdi(iday).a_ntot_CB(ileg), clouds.gomaccspdi(iday).a_ntot_CB_sampsize(ileg)] = ...
                aerCMS(cm_a_ntot);
            [clouds.gomaccspdi(iday).s_ntot_CB(ileg), clouds.gomaccspdi(iday).s_ntot_CB_sampsize(ileg)] = ...
                cldCMS(cm_s_ntot);
            [clouds.gomaccspdi(iday).s_actfrac_CB(ileg), clouds.gomaccspdi(iday).s_actfrac_CB_sampsize(ileg)] = ...
                cldCMS(cm_s_ntot./cm_a_ntot);
            [clouds.gomaccspdi(iday).reldisp_CB(ileg), clouds.gomaccspdi(iday).reldisp_CB_sampsize(ileg)] = ...
                cldCMS(cm_disp);
        end
    else
        [clouds.gomaccspdi(iday).a_ntot_CB, clouds.gomaccspdi(iday).a_ntot_CB_sampsize, ...
            clouds.gomaccspdi(iday).s_ntot_CB, clouds.gomaccspdi(iday).s_ntot_CB_sampsize, ...
            clouds.gomaccspdi(iday).s_actfrac_CB, clouds.gomaccspdi(iday).s_actfrac_CB_sampsize, ...
            clouds.gomaccspdi(iday).reldisp_CB, clouds.gomaccspdi(iday).reldisp_CB_sampsize] = ...
            deal(nan);
    end
end

%%
save('clouds.mat','clouds', '-v7.3')
% toc
% finishingTaskSound