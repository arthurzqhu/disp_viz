clear
cd ~/MEGAsync/grad/research/aerosol_reldisp/datasets/
load clouds.mat
vocals_aer_files = dir('vocals/VOCALS_pcasp/PCASP_08*');
% vocals_ccn_files = dir('vocals/VOCALS_CCN/CCN_08*');
% vocals_caps_files = dir('vocals/VOCALS_caps/caps_08*');

%%
pdi_date_start_char = regexp(clouds.vocalspdi(1).file,'/');

if isempty(pdi_date_start_char)
    pdi_date_start_char = 0;
end

pdi_dates = arrayfun(@(x) clouds.vocalspdi(x).file(pdi_date_start_char+1:pdi_date_start_char+6),...
    1:length(clouds.vocalspdi), 'UniformOutput', false)';

pcasp_dates = arrayfun(@(x) vocals_aer_files(x).name(7:12), 1:length(vocals_aer_files),...
    'UniformOutput', false)';

pcasp_binedges = readmatrix([vocals_aer_files(1).folder '/pcasp_bnds.txt'], 'Range',...
    [2 1 21 2]);

pcasp_binmean = readmatrix([vocals_aer_files(1).folder '/pcasp_bnds.txt'], 'Range',...
    [1 3 21 3]);

% caps_dates = arrayfun(@(x) vocals_caps_files(x).name(6:11), 1:length(vocals_caps_files),...
%     'UniformOutput',false)';
% cas_binedge = readmatrix([vocals_caps_files(1).folder '/cas_bnds.txt'], 'Range',...
%     [2 1 21 2]);
% cas_binmean = readmatrix([vocals_caps_files(1).folder '/cas_bnds.txt'], 'Range',...
%     [1 3 21 3]);
% 
% ccn_dates = arrayfun(@(x) vocals_ccn_files(x).name(5:10), 1:length(vocals_ccn_files),...
%     'UniformOutput', false)';

allinst_commondates = mintersect(pdi_dates, pcasp_dates);
cmd_ipdi = find(ismember(pdi_dates, allinst_commondates));
cmd_ipcasp = find(ismember(pcasp_dates, allinst_commondates));
% cmd_iccn = find(ismember(ccn_dates, allinst_commondates));
% cmd_icaps = find(ismember(caps_dates, allinst_commondates));
%%
% tic     
for iday = 1:length(allinst_commondates)
    ipdi = cmd_ipdi(iday);
    ipcasp = cmd_ipcasp(iday);
%     iccn = cmd_iccn(iday);
%     icaps = cmd_icaps(iday);
%     
    aeros_dat = readmatrix([vocals_aer_files(ipcasp).folder '/' vocals_aer_files(ipcasp).name],...
            'FileType','text','Range',[2 2]);
    aeros = aeros_dat(:,1:20);
    
    time_length = size(aeros,1);
    aeros(aeros<0)=nan;
    a_meanD_ex = arrayfun(@(x) wmean(pcasp_binmean,aeros(x,:)),1:time_length)';

    aeros_dt = readmatrix([vocals_aer_files(ipcasp).folder '/' vocals_aer_files(ipcasp).name],...
            'FileType','text','Range',[2 39 time_length+2 39], 'OutputType', 'string');

    aeros_dtnum = datenum(aeros_dt,'HH:MM:SS');
    
    clouds.vocalspdi(iday).a_date = pcasp_dates(ipcasp);
    clouds.vocalspdi(iday).a_t = round((aeros_dtnum-floor(aeros_dtnum))*86400,3);
    clouds.vocalspdi(iday).a_ntot_ex = sum(aeros,2);
    clouds.vocalspdi(iday).a_meanD_ex = a_meanD_ex;
    
    [commontime, commontime_ipdi, commontime_ipcasp] = ...
        intersect(floor(clouds.vocalspdi(iday).s_t), clouds.vocalspdi(iday).a_t);
    
    
    % add s_ntot_pdi into a_ntot
    s_ntot = clouds.vocalspdi(iday).s_ntot_pdi;
    s_ntot(isnan(s_ntot))=0;
    clouds.vocalspdi(iday).a_ntot = sum(aeros,2);
    clouds.vocalspdi(iday).a_ntot(commontime_ipcasp) = ...
        clouds.vocalspdi(iday).a_ntot(commontime_ipcasp) + s_ntot(commontime_ipdi);
    
    % calculating activation fraction
    
    clouds.vocalspdi(iday).s_actfrac = zeros(size(s_ntot));
    clouds.vocalspdi(iday).s_actfrac(commontime_ipdi) = s_ntot(commontime_ipdi)./clouds.vocalspdi(iday).a_ntot(commontime_ipcasp);
    
    
%     ccn = readmatrix([vocals_ccn_files(iccn).folder '/' vocals_ccn_files(iccn).name],...
%         'FileType','text','Range',[1 2]);
%     time_length = size(ccn,1);
%     
%     ccn_dt = readmatrix([vocals_ccn_files(iccn).folder '/' vocals_ccn_files(iccn).name],...
%         'FileType','text','Range',[1 1 time_length 1], 'OutputType', 'string');
%     vidx = ccn_dt~='-9999';
%     ividx = ccn_dt=='-9999';
%     ccn_dt(ividx)=[];
%     
%     ccn_dtnum = datenum(ccn_dt,'HH:MM:SS');
%     
%     ccn_a = ccn(vidx,11); ccn_a(ccn_a<=0) = nan;
%     ccn_b = ccn(vidx,12); ccn_b(ccn_b<=0) = nan;
%     
%     clouds.vocalspdi(iday).ccn_date = ccn_dates(iccn);
%     clouds.vocalspdi(iday).ccn_t = round((ccn_dtnum-floor(ccn_dtnum))*86400,3);
%     clouds.vocalspdi(iday).ccn_a = ccn_a;
%     clouds.vocalspdi(iday).ccn_b = ccn_b;
%     
%     %%
%     cas_dat = readmatrix([vocals_caps_files(icaps).folder '/' vocals_caps_files(icaps).name],...
%         'FileType','text');
%     time_length = size(cas_dat,1);
%     cas_conc = cas_dat(1:end, 2:21);
%     cas_conc(cas_conc==-9999)=nan;
%     cas_dt = readmatrix([vocals_caps_files(icaps).folder '/' vocals_caps_files(icaps).name],...
%         'FileType','text', 'Range', [1 106 time_length 106],'OutputType','string');
%     cas_dtnum = datenum(cas_dt, 'HH:MM:SS');
%     
%     cas_t = round((cas_dtnum - floor(cas_dtnum))*86400);
%     cas_ntot = sum(cas_conc,2);
%     cas_meand = sum(cas_conc.*cas_binmean',2)./cas_ntot;
%     cas_lwc = sum(cas_conc.*cas_binmean'.^3,2)*pi/6*1e-6;
%     cas_std = arrayfun(@(x) nanstd(cas_binmean, cas_conc(x,:)), 1:length(cas_t))';
%     cas_disp = cas_std./cas_meand;
%     
%     ividx = cas_ntot<25 | cas_meand<3;
%     cas_std(ividx) = nan;
%     cas_disp(ividx) = nan;
    
    
    s_t = floor(clouds.vocalspdi(iday).s_t);
%     s_ap = clouds.vocalspdi(iday).s_ap;
    normAC = clouds.vocalspdi(iday).normAC;
%     load vocalspdi_flight_basics.mat
%     for ileg = 1:length(vocalspdi_flight_basics(iday).ti)
%         ti = vocalspdi_flight_basics(iday).ti(ileg);
%         tf = vocalspdi_flight_basics(iday).tf(ileg);
%         z_min = min(s_ap(s_t<tf & s_t>ti));
%         z_max = max(s_ap(s_t<tf & s_t>ti));
%         z_maxsampled = (z_min + z_max)/2;
%     end
    filt_crit = normAC>0 & normAC<1;
    s_t = s_t(filt_crit);
    
    s_ntot = clouds.vocalspdi(iday).s_ntot_pdi(filt_crit);
    s_disp = clouds.vocalspdi(iday).s_disp_pdi(filt_crit);
    s_meand = clouds.vocalspdi(iday).drpsz(filt_crit);
    
%     [~, cmt_icas, cmt_ipdi] = intersect(cas_t,s_t);
    
%     mean_cas_disp(iday) = nanmean(cas_disp);
    mean_s_disp(iday) = nanmean(s_disp);
    
%     %% cas vs pdi evaluation
% %     close all
%     figure(1)
%     hold on
%     plot(s_t(cmt_ipdi), cas_disp(cmt_icas)./s_disp(cmt_ipdi)/1.4087,'.');
%     xlabel('seconds since 00Z');
%     ylabel('\epsilon_{cas}/\epsilon_{pdi}')
%     hold off
%     
%     figure(2)
%     hold on
%     plot(s_disp(cmt_ipdi), cas_disp(cmt_icas)./s_disp(cmt_ipdi)/1.4087,'.');
%     xlabel('\epsilon_{pdi}')
%     ylabel('\epsilon_{cas}/\epsilon_{pdi}')
%     hold off
%     
%     figure(3)
%     hold on
%     plot(s_ntot(cmt_ipdi), cas_disp(cmt_icas)./s_disp(cmt_ipdi)/1.4087,'.');
%     xlabel('N_{d,pdi}')
%     ylabel('\epsilon_{cas}/\epsilon_{pdi}')
%     hold off
%     
%     figure(4)
%     hold on
%     plot(s_meand(cmt_ipdi), cas_disp(cmt_icas)./s_disp(cmt_ipdi)/1.4087,'.');
%     xlabel('$$\overline{D}_{pdi}$$', 'Interpreter','latex')
%     ylabel('\epsilon_{cas}/\epsilon_{pdi}')
%     hold off
%     
%     
% %     plot(drpsz(cmt_ipdi),cas_meand(cmt_icas),'.'); refline(1,0)
%     
% %     clouds
% %     clouds.vocalspdi(iday).cas_t = cas_t;
    
end

% toc

%%
% scatter(mean_cas_disp, mean_s_disp)
% refline(1,0)
% xlabel('CAS')
% ylabel('PDI')
% grid
%%
save('clouds.mat','clouds', '-v7.3')

% finishingTaskSound