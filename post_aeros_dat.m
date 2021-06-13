clear
cd ~/MEGAsync/grad/research/aerosol_reldisp/datasets/
load clouds.mat
post_aer_files = dir('post/POST_pcasp/pcasp_1hz_*');
% post_ccn_files = dir('post/POST_CCN/CCN_1hz_*');
% post_caps_files = dir('post/POST_caps/CASF_1hz_RF*');


%%
pdi_date_start_char = regexp(clouds.postpdi(1).file,'/');

if isempty(pdi_date_start_char)
    pdi_date_start_char = 0;
end

pdi_dates = arrayfun(@(x) clouds.postpdi(x).file(pdi_date_start_char+1:pdi_date_start_char+6),...
    1:length(clouds.postpdi), 'UniformOutput', false)';

pcasp_dates = arrayfun(@(x) post_aer_files(x).name(18:23), 1:length(post_aer_files),...
    'UniformOutput', false)';

pcasp_binedges = readmatrix([post_aer_files(1).folder '/pcasp_bnds.txt'], 'Range',...
    [2 1 21 2]);

pcasp_binmean = readmatrix([post_aer_files(1).folder '/pcasp_bnds.txt'], 'Range',...
    [1 3 21 3]);

% ccn_dates = arrayfun(@(x) post_ccn_files(x).name(16:21), 1:length(post_ccn_files),...
%     'UniformOutput', false)';
% 
% caps_dates = arrayfun(@(x) post_caps_files(x).name(17:22), 1:length(post_caps_files),...
%     'UniformOutput', false)';
% cas_binmean = readmatrix([post_caps_files(1).folder '/bin_mean.txt'], 'FileType',...
%     'text');



allinst_commondates = mintersect(pdi_dates, pcasp_dates); %, ccn_dates);
cmd_ipdi = find(ismember(pdi_dates, allinst_commondates));
cmd_ipcasp = find(ismember(pcasp_dates, allinst_commondates));
% cmd_iccn = find(ismember(ccn_dates, allinst_commondates));
% cmd_icaps = find(ismember(caps_dates, allinst_commondates));
%%
for iday = 1:length(allinst_commondates)
    ipdi = cmd_ipdi(iday);
    ipcasp = cmd_ipcasp(iday);
%     iccn = cmd_iccn(iday);
%     icaps = cmd_icaps(iday);
    
    aeros = readmatrix([post_aer_files(ipcasp).folder '/' post_aer_files(ipcasp).name],...
            'FileType','text','Range',[13 2]);
    time_length = size(aeros,1);
    aeros(aeros<0)=nan;
    a_meanD_ex = arrayfun(@(x) wmean(pcasp_binmean,aeros(x,:)),1:time_length)';
    
    aeros_dt = readmatrix([post_aer_files(ipcasp).folder '/' post_aer_files(ipcasp).name],...
            'FileType','text','Range',[13 1 time_length+13 1], 'OutputType', 'string');

    aeros_dtnum = datenum(aeros_dt,'yymmddHHMMSS');
    
    clouds.postpdi(ipdi).a_date = pcasp_dates(ipdi);
    clouds.postpdi(ipdi).a_t = round((aeros_dtnum-floor(aeros_dtnum))*86400,3);
    clouds.postpdi(ipdi).a_ntot_ex = sum(aeros,2);
    clouds.postpdi(iday).a_meanD_ex = a_meanD_ex;
    
    [commontime, commontime_ipdi, commontime_ipcasp] = ...
        intersect(floor(clouds.postpdi(ipdi).s_t), clouds.postpdi(ipdi).a_t);
        
    s_ntot = clouds.postpdi(ipdi).s_ntot_pdi;
    s_ntot(isnan(s_ntot))=0;
    clouds.postpdi(ipdi).a_ntot = sum(aeros,2);
    clouds.postpdi(ipdi).a_ntot(commontime_ipcasp) = ...
        clouds.postpdi(ipdi).a_ntot(commontime_ipcasp) + s_ntot(commontime_ipdi);
    
    % calculating activation fraction
    
    clouds.postpdi(ipdi).s_actfrac = nan(size(s_ntot));
    clouds.postpdi(ipdi).s_actfrac(commontime_ipdi) = s_ntot(commontime_ipdi)./clouds.postpdi(ipdi).a_ntot(commontime_ipcasp);
    
%     ccn = readmatrix([post_ccn_files(iccn).folder '/' post_ccn_files(iccn).name],...
%         'FileType','text','Range',[13 2]);
%     time_length = size(ccn,1);
%     
%     ccn_dt = readmatrix([post_ccn_files(iccn).folder '/' post_ccn_files(iccn).name],...
%         'FileType','text','Range',[13 1 time_length+13 1], 'OutputType', 'string');
% 
%     ccn_dtnum = datenum(ccn_dt,'yymmddHHMMSS');
%     
%     ccn_a = ccn(:,9); ccn_a(ccn_a<=0) = nan;
%     ccn_b = ccn(:,10); ccn_b(ccn_b<=0) = nan;
%     
%     clouds.postpdi(ipdi).ccn_date = ccn_dates(iccn);
%     clouds.postpdi(ipdi).ccn_t = round((ccn_dtnum-floor(ccn_dtnum))*86400,3);
%     clouds.postpdi(ipdi).ccn_a = ccn_a;
%     clouds.postpdi(ipdi).ccn_b = ccn_b;
    %% caps
%     cas_conc = readmatrix([post_caps_files(icaps).folder '/' post_caps_files(icaps).name],...
%         'FileType','text','Range',[1 2]);
%     cas_conc(cas_conc<0) = nan;
%     time_length = size(cas_conc,1);
%     cas_dt = readmatrix([post_caps_files(icaps).folder '/' post_caps_files(icaps).name],...
%         'FileType','text','Range',[15 1 time_length+15 1],'OutputType','string');
%     wrong_dtype_idx = find(contains(cas_dt,'e+10'));
%     
%     % use for loop to make sure two consecutive data with wrong data type are treated
%     % correctly
%     for idx = wrong_dtype_idx'
%         prev_correct = cas_dt(idx-1);
%         correct_dtype = ['0' num2str(str2num(prev_correct)+1)];
%         cas_dt(idx) = correct_dtype;
%     end
%     
%     cas_dtnum = datenum(cas_dt,'yymmddHHMMSS');
%     
%     cas_t = round((cas_dtnum-floor(cas_dtnum))*86400);
%     
%     cas_ntot = sum(cas_conc,2);
%     cas_meand = sum(cas_conc.*cas_binmean',2)./cas_ntot;
%     cas_lwc = sum(cas_conc.*cas_binmean'.^3,2)*pi/6*1e-6;
%     cas_std = arrayfun(@(x) nanstd(cas_binmean, cas_conc(x,:)), 1:length(cas_t))';
%     cas_disp = cas_std./cas_meand;
%     
%     ividx = cas_ntot<25 | cas_meand<3;
%     cas_std(ividx) = nan;
%     cas_disp(ividx) = nan;
    
    
    s_t = floor(clouds.postpdi(iday).s_t);
    normAC = clouds.postpdi(iday).normAC;
    
    filt_crit = normAC>0 & normAC<1;
    
    s_t = s_t(filt_crit);
    
    s_ntot = clouds.postpdi(iday).s_ntot_pdi(filt_crit);
    s_disp = clouds.postpdi(iday).s_disp_pdi(filt_crit);
    s_meand = clouds.postpdi(iday).drpsz(filt_crit);
    
%     [~, cmt_icas, cmt_ipdi] = intersect(cas_t,s_t);
%     
%     mean_cas_disp(iday) = nanmean(cas_disp);
    mean_s_disp(iday) = nanmean(s_disp);
    
%     %% cas vs pdi evaluation
%     figure(1)
%     hold on
%     plot(s_t(cmt_ipdi), cas_disp(cmt_icas)./s_disp(cmt_ipdi),'.');
%     xlabel('seconds since 00Z');
%     ylabel('\epsilon_{cas}/\epsilon_{pdi}')
%     hold off
%     
%     figure(2)
%     hold on
%     plot(s_disp(cmt_ipdi), cas_disp(cmt_icas)./s_disp(cmt_ipdi),'.');
%     xlabel('\epsilon_{pdi}')
%     ylabel('\epsilon_{cas}/\epsilon_{pdi}')
%     hold off
%     
%     figure(3)
%     hold on
%     plot(s_ntot(cmt_ipdi), cas_disp(cmt_icas)./s_disp(cmt_ipdi),'.');
%     xlabel('N_{d,pdi}')
%     ylabel('\epsilon_{cas}/\epsilon_{pdi}')
%     hold off
%     
%     figure(4)
%     hold on
%     plot(s_meand(cmt_ipdi), cas_disp(cmt_icas)./s_disp(cmt_ipdi),'.');
%     xlabel('$$\overline{D}_{pdi}$$', 'Interpreter','latex')
%     ylabel('\epsilon_{cas}/\epsilon_{pdi}')
%     hold off
end

% %%
% scatter(mean_cas_disp, mean_s_disp)
% refline(1,0)
% xlabel('CAS')
% ylabel('PDI')
% grid
%%
save('clouds.mat','clouds', '-v7.3')