clear
cd ~/MEGAsync/grad/research/aerosol_reldisp/datasets/
load clouds.mat
mase_aer_files = dir('mase/MASE_pcasp/PCASP_05*');
% pcasp_bin_edges = readmatrix('mase/MASE_pcasp/pcasp_bin_edge.csv');

%%
pdi_date_start_char = regexp(clouds.masepdi(1).file,'/');

if isempty(pdi_date_start_char)
    pdi_date_start_char = 0;
end

pdi_dates = arrayfun(@(x) clouds.masepdi(x).file(pdi_date_start_char+1:pdi_date_start_char+6),...
    1:length(clouds.masepdi), 'UniformOutput', false)';

pcasp_dates = arrayfun(@(x) mase_aer_files(x).name(7:12), 1:length(mase_aer_files),...
    'UniformOutput', false)';

pcasp_binedges = readmatrix([mase_aer_files(1).folder '/pcasp_bin_edge.csv'], 'Range',...
    [2 1 21 2]);

pcasp_binmean = readmatrix([mase_aer_files(1).folder '/pcasp_bin_edge.csv'], 'Range',...
    [1 3 21 3]);

allinst_commondates = mintersect(pdi_dates, pcasp_dates);
cmd_ipdi = find(ismember(pdi_dates, allinst_commondates));
cmd_ipcasp = find(ismember(pcasp_dates, allinst_commondates));
%%
% tic     
for iday = 1:length(allinst_commondates)
    ipdi = cmd_ipdi(iday);
    ipcasp = cmd_ipcasp(iday);
%     
    aeros_dat = readmatrix([mase_aer_files(ipcasp).folder '/' mase_aer_files(ipcasp).name],...
            'FileType','text','Range',[1 7]);
    aeros = aeros_dat;
    
    time_length = size(aeros,1);
    aeros(aeros<0)=nan;
    a_meanD_ex = arrayfun(@(x) wmean(pcasp_binmean,aeros(x,:)),1:time_length)';
    
    aeros_dt = readmatrix([mase_aer_files(ipcasp).folder '/' mase_aer_files(ipcasp).name],...
            'FileType','text','Range',[1 6 time_length+1 6], 'OutputType', 'string');

    aeros_dtnum = datenum(aeros_dt,'mm/dd/yyyy HH:MM:SS');
    
    a_ntot_ex = sum(aeros,2);
    a_t = round((aeros_dtnum-floor(aeros_dtnum))*86400,3);
%     clouds.masepdi(iday).a_date = pcasp_dates(ipcasp);
    clouds.masepdi(iday).a_t = a_t;
    clouds.masepdi(iday).a_ntot_ex = a_ntot_ex;
    clouds.masepdi(iday).a_meanD_ex = a_meanD_ex;
    
    [cmt, cmt_ipdi, cmt_ipcasp] = ...
        intersect(floor(clouds.masepdi(iday).s_t), a_t);
    
    
    % add s_ntot_pdi into a_ntot
    s_ntot = clouds.masepdi(iday).s_ntot_pdi;
    s_ntot(isnan(s_ntot))=0;
    a_ntot = a_ntot_ex; % make them the same length and share the same data
    a_ntot(cmt_ipcasp) = a_ntot_ex(cmt_ipcasp) + s_ntot(cmt_ipdi);
    
    clouds.masepdi(iday).a_ntot = a_ntot;
    
    % calculating activation fraction
    s_actfrac = zeros(size(s_ntot));
    s_actfrac(cmt_ipdi) = s_ntot(cmt_ipdi)./a_ntot(cmt_ipcasp);
    
    clouds.masepdi(iday).s_actfrac = s_actfrac;
    
    
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