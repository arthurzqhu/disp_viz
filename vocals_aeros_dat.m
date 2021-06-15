close all
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
clearvars -except clouds
if ~exist('clouds','var') load clouds.mat, end
vocals_aer_files = dir('vocals/VOCALS_pcasp/PCASP_08*');

%% get rid of the old fieldnames first
for iday = 1:length(clouds.vocalspdi)
   vocals_fn=fieldnames(clouds.vocalspdi(iday));
   ext_idx = startsWith(vocals_fn,'a_') | startsWith(vocals_fn,'ccn_');
   newstrct(iday)=rmfield(clouds.vocalspdi(iday),vocals_fn(ext_idx));
end
clouds.vocalspdi=newstrct;
%%
pdi_date_start_char = regexp(clouds.vocalspdi(1).file,'/');

if isempty(pdi_date_start_char)
   pdi_date_start_char = 0;
end

pdi_dates = arrayfun(@(x) ...
   clouds.vocalspdi(x).file(pdi_date_start_char+1:pdi_date_start_char+6),...
   1:length(clouds.vocalspdi), 'UniformOutput', false)';

pcasp_dates = arrayfun(@(x) vocals_aer_files(x).name(7:12),...
   1:length(vocals_aer_files),...
   'UniformOutput', false)';

pcasp_binedges = readmatrix([vocals_aer_files(1).folder '/pcasp_bnds.txt'],...
   'Range',...
   [2 1 21 2]);

pcasp_binmean = readmatrix([vocals_aer_files(1).folder '/pcasp_bnds.txt'],...
   'Range',...
   [2 3 21 3]);

allinst_commondates = mintersect(pdi_dates, pcasp_dates);
cmd_ipdi = find(ismember(pdi_dates, allinst_commondates));
cmd_ipcasp = find(ismember(pcasp_dates, allinst_commondates));
%%
for iday = 1:length(allinst_commondates)
   ipdi_d = cmd_ipdi(iday);
   ipcasp_d = cmd_ipcasp(iday);
   
   %% get the common time
   pcasp_dat = readmatrix([vocals_aer_files(ipcasp_d).folder '/' ...
      vocals_aer_files(ipcasp_d).name],...
      'FileType','text','Range',[2 2]);
   pcasp_t_len = size(pcasp_dat,1);
   pcasp_t_str = readmatrix([vocals_aer_files(ipcasp_d).folder '/' ...
      vocals_aer_files(ipcasp_d).name],'FileType','text','Range',...
      [2 39 pcasp_t_len+2 39], 'OutputType', 'string');
   pcasp_dtnum = datenum(pcasp_t_str,'HH:MM:SS');
   pcasp_t = round((pcasp_dtnum-floor(pcasp_dtnum))*86400,3);
   pdi_t = round(clouds.vocalspdi(iday).s_t);
   [s_t,ipdi_t,ipcasp_t] = intersect(pdi_t,pcasp_t);
   
   %% downsample the existing s_* variables first
   vocals_fn=fieldnames(clouds.vocalspdi(iday));
   if length(clouds.vocalspdi(iday).s_t)~=length(s_t)
      for ifn = 1:length(vocals_fn)
         if size(clouds.vocalspdi(iday).(vocals_fn{ifn}),1)==length(pdi_t)
            clouds.vocalspdi(iday).(vocals_fn{ifn})=...
               clouds.vocalspdi(iday).(vocals_fn{ifn})(ipdi_t,:);
         end
      end
   end
   % round the time anyway
   clouds.vocalspdi(iday).s_t=s_t;
   
   %%
   s_conc_pcasp = pcasp_dat(ipcasp_t,1:20);
   s_conc_pcasp(s_conc_pcasp<0)=nan;
   s_meand_pcasp = arrayfun(@(x) wmean(pcasp_binmean,s_conc_pcasp(x,:)),...
      1:length(s_t))';
   s_ntot_pcasp = sum(s_conc_pcasp,2);
   
   s_ntot_pdi = clouds.vocalspdi(iday).s_ntot_pdi;
   s_ntot_pdi(isnan(s_ntot_pdi))=0;
   s_ntot_aer = s_ntot_pcasp + s_ntot_pdi;
   
   s_actfac=s_ntot_pdi./s_ntot_aer;
   
   clouds.vocalspdi(iday).s_ntot_pcasp = s_ntot_pcasp;
   clouds.vocalspdi(iday).s_meand_pcasp = s_meand_pcasp;
   clouds.vocalspdi(iday).s_ntot_aer = s_ntot_aer;
   clouds.vocalspdi(iday).s_actfrac = s_actfac;

end
%%
% save('clouds.mat','clouds', '-v7.3')

% finishingTaskSound