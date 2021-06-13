clear
load clouds.mat

%%
% close all
cas_files = dir('MASE_raw/senum-cas/*.txt');
% cip_files = dir('MASE_raw/senum-cip/*.txt');
% ccn_files = dir('MASE_raw/wang-ccn/*.txt');
% pcasp_files = dir('MASE_raw/hubbe-pcasp/*.txt');
summary_files = dir('MASE_raw/springston-summary/*.txt');

cas_binedge = readmatrix('MASE_raw/senum-cas/bin_edges.csv', 'Range', [1 2]);
cas_binmean = (cas_binedge(:,2)+cas_binedge(:,1))/2;

%%

clear summary_dat pcasp_dat cas_dat ccn_dat

for iday = 1:length(cas_files)
    %% read cas
    cas_opts = fixedWidthImportOptions('NumVariables', 25);
    cas_opts.VariableWidths = [11 12 repelem(18, 22) 17];
    cas_opts.VariableNames = {'Date','Time','CAS_Bin_1','CAS_Bin_2','CAS_Bin_3',...
        'CAS_Bin_4','CAS_Bin_5','CAS_Bin_6','CAS_Bin_7','CAS_Bin_8','CAS_Bin_9',...
        'CAS_Bin_10','CAS_Bin_11','CAS_Bin_12','CAS_Bin_13','CAS_Bin_14','CAS_Bin_15',...
        'CAS_Bin_16','CAS_Bin_17','CAS_Bin_18','CAS_Bin_19','CAS_Bin_20','CAS_tot_N',...
        'CAS_tot_A','CAS_tot_V'};
    cas_opts.DataLines = [40, Inf];
    cas_opts.VariableTypes = [{'char','char'} repelem({'double'},23)];
    cas_dat{iday,1}(:,:) = readtable([cas_files(iday).folder '/' cas_files(iday).name],...
        cas_opts);
    
    clouds.masecas(iday).s_date = cas_files(iday).name(1:6);
    
    s_time_raw = [cas_dat{iday}{:,'Date'} cas_dat{iday}{:,'Time'}];
    s_datetime_str = strcat(s_time_raw(:,1)," ",s_time_raw(:,2));
    s_timedtmum = datenum(s_datetime_str,'yyyy-mm-dd HH:MM:SS.FFF');
    
    s_t = (s_timedtmum - floor(s_timedtmum(1)))*86400;
    clouds.masecas(iday).s_t = s_t;
    
    s_conc_cas = cas_dat{iday}{:,3:22};
    s_ntot_cas = cas_dat{iday}{:,23};
    s_meand_cas = sum(s_conc_cas*cas_binmean,2)./s_ntot_cas;
    s_lwc_cas = sum(s_conc_cas*cas_binmean.^3,2)*pi/6*1e-6;
    s_std_cas = arrayfun(@(x) nanstd(cas_binmean, s_conc_cas(x,:)), 1:length(s_t))';
    s_disp_cas = s_std_cas./s_meand_cas;
    
    % filter out invalid index
    ividx = s_ntot_cas<25 | s_meand_cas<3;
    s_std_cas(ividx) = nan;
    s_disp_cas(ividx) = nan;
    
    clouds.masecas(iday).s_conc_cas = s_conc_cas;
    clouds.masecas(iday).s_ntot_cas = s_ntot_cas;
    clouds.masecas(iday).s_meand_cas = s_meand_cas;
    clouds.masecas(iday).s_lwc_cas = s_lwc_cas;
    clouds.masecas(iday).s_std_cas = s_std_cas;
    clouds.masecas(iday).s_disp_cas = s_disp_cas;
    
    %% summary
    summ_opts = delimitedTextImportOptions('NumVariables', 63);
    summ_opts.VariableNamesLine = 35;
    summ_opts.Delimiter = '\t';
    summ_opts.DataLines = 40;
    summ_opts.VariableTypes = [{'char','char'} repelem({'double'},61)];
    
    summary_dat{iday,1}(:,:) = readtable([summary_files(iday).folder '/' summary_files(iday).name],...
        summ_opts);
    
    clouds.masecas(iday).s_ap = summary_dat{iday}{:,7};
    clouds.masecas(iday).s_ta = summary_dat{iday}{:,8};
    clouds.masecas(iday).s_thet = summary_dat{iday}{:,9};
    clouds.masecas(iday).s_rh = summary_dat{iday}{:,11};
    clouds.masecas(iday).ccn_a = summary_dat{iday}{:,18};
    clouds.masecas(iday).ccn_b = summary_dat{iday}{:,17};
    clouds.masecas(iday).a_ntot = summary_dat{iday}{:,22} + s_ntot_cas; % <- includes the ones in the cloud droplets
    clouds.masecas(iday).s_ntot_casrpt = summary_dat{iday}{:,25};
    clouds.masecas(iday).s_lwc_gerb = summary_dat{iday}{:,31};
    clouds.masecas(iday).s_lwc_hotw = summary_dat{iday}{:,32};
    clouds.masecas(iday).s_actfrac = s_ntot_cas./clouds.masecas(iday).a_ntot;
    
%     %% pcasp
%     pcasp_opts = fixedWidthImportOptions('NumVariables', 26);
%     pcasp_opts.VariableWidths = [11 12 13 repelem(18, 22) 17];
%     pcasp_opts.VariableNames = {'Date','Time','Flow_Amb','PCASP_Bin_1','PCASP_Bin_2',...
%         'PCASP_Bin_3','PCASP_Bin_4','PCASP_Bin_5','PCASP_Bin_6','PCASP_Bin_7',...
%         'PCASP_Bin_8','PCASP_Bin_9','PCASP_Bin_10','PCASP_Bin_11','PCASP_Bin_12',...
%         'PCASP_Bin_13','PCASP_Bin_14','PCASP_Bin_15','PCASP_Bin_16','PCASP_Bin_17',...
%         'PCASP_Bin_18','PCASP_Bin_19','PCASP_Bin_20','PCASP_tot_N','PCASP_tot_A',...
%         'PCASP_tot_V'};
%     pcasp_opts.DataLines = [40, Inf];
%     pcasp_opts.VariableTypes = [{'char','char'} repelem({'double'},24)];
%     pcasp_dat{iday,1}(:,:) = readtable([pcasp_files(iday).folder '/' pcasp_files(iday).name],...
%         pcasp_opts);
%     
%     clouds.masecas(iday).a_ntot = pcasp_dat{iday}{:, 'PCASP_tot_N'};
%     
%     clouds.masecas(iday).a_date = pcasp_files(iday).name(1:6);
%     
%     a_time_raw = [pcasp_dat{iday}{:,'Date'} pcasp_dat{iday}{:,'Time'}];
%     a_datetime_str = strcat(a_time_raw(:,1)," ",a_time_raw(:,2));
%     a_timedtmum = datenum(a_datetime_str,'yyyy-mm-dd HH:MM:SS.FFF');
%     clouds.masecas(iday).a_t = (a_timedtmum - floor(a_timedtmum(1)))*86400;
%    
%         
%     %% read_cscn
%     ccn_opts = fixedWidthImportOptions('NumVariables', 13);
%     ccn_opts.VariableWidths = [11 repelem(12, 11) 11];
%     ccn_opts.VariableNamesLine=35;
%     ccn_opts.VariableNames = {'Date','Time','CCNs_1','CCNs_0_6','CCNs_0_3','CCNs_0_2',...
%         'CCNs_0_1','CCNs_0_08','CCNs_0_06','CCNs_0_04','CCNs_0_02','DMT_CCN1',...
%         'DMT_CCN2'};
% 	ccn_opts.DataLines = [40, Inf];
% 	ccn_opts.VariableTypes = [{'char','char'} repelem({'double'},11)];
%     
%     ccn_dat{iday,1}(:,:) = readtable([ccn_files(iday).folder '/' ccn_files(iday).name],...
%         ccn_opts);
%     
%     clouds.masecas(iday).ccn_date = ccn_files(iday).name(1:6);
%     
%     ccn_time_raw = [ccn_dat{iday}{:,'Date'} ccn_dat{iday}{:,'Time'}];
%     ccn_datetime_str = strcat(ccn_time_raw(:,1)," ",ccn_time_raw(:,2));
%     ccn_timedtmum = datenum(ccn_datetime_str,'yyyy-mm-dd HH:MM:SS.FFF');
%     clouds.masecas(iday).ccn_t = (ccn_timedtmum - floor(ccn_timedtmum(1)))*86400;
%     clouds.masecas(iday).ccn_tab = ccn_dat{iday}{:,3:end};
end

% %%
save('clouds.mat','clouds');

% finishingTaskSound