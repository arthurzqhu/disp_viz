close all
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
clearvars -except clouds
if ~exist('clouds','var') load clouds.mat, end
mase_aer_files = dir('mase/MASE_pcasp/PCASP_05*');

% %% get rid of the old fieldnames first
% for iday = 1:length(clouds.masepdi)
%    mase_fn=fieldnames(clouds.masepdi(iday));
%    ext_idx = startsWith(mase_fn,'a_') | startsWith(mase_fn,'ccn_');
%    newstrct(iday)=rmfield(clouds.masepdi(iday),mase_fn(ext_idx));
% end
% clouds.masepdi=newstrct;
%%
pdi_date_start_char = regexp(clouds.masepdi(1).file,'/');

if isempty(pdi_date_start_char)
   pdi_date_start_char = 0;
end

pdi_dates = arrayfun(@(x) ...
   clouds.masepdi(x).file(pdi_date_start_char+1:pdi_date_start_char+6),...
   1:length(clouds.masepdi), 'UniformOutput', false)';
pcasp_dates = arrayfun(@(x) mase_aer_files(x).name(7:12), ...
   1:length(mase_aer_files),'UniformOutput', false)';
pcasp_binedges = readmatrix([mase_aer_files(1).folder '/pcasp_bin_edge.csv'],...
   'Range',[2 1 21 2]);
pcasp_binmean = readmatrix([mase_aer_files(1).folder '/pcasp_bin_edge.csv'],...
   'Range',[2 3 21 3]);

allinst_commondates = mintersect(pdi_dates, pcasp_dates);
cmd_ipdi = find(ismember(pdi_dates, allinst_commondates));
cmd_ipcasp = find(ismember(pcasp_dates, allinst_commondates));
%%

for iday = 1:length(allinst_commondates)
   ipdi_d = cmd_ipdi(iday);
   ipcasp_d = cmd_ipcasp(iday);
   %% get the common time
   pcasp = readmatrix([mase_aer_files(ipcasp_d).folder '/' ...
      mase_aer_files(ipcasp_d).name],'FileType','text','Range',[1 7]);
   pcasp_t_len = size(pcasp,1);
   pcasp(pcasp<0)=nan;
   pcasp_t_str = readmatrix([mase_aer_files(ipcasp_d).folder '/' ...
      mase_aer_files(ipcasp_d).name],'FileType','text','Range',...
      [1 6 pcasp_t_len+1 6], 'OutputType', 'string');
   aeros_dtnum = datenum(pcasp_t_str,'mm/dd/yyyy HH:MM:SS');
   pcasp_t = round((aeros_dtnum-floor(aeros_dtnum))*86400,3);
   pdi_t = round(clouds.masepdi(iday).s_t);
   [s_t,ipdi_t,ipcasp_t] = intersect(pdi_t,pcasp_t);
   %% downsample the existing s_* variables first
   mase_fn=fieldnames(clouds.masepdi(iday));
   if length(clouds.masepdi(iday).s_t)~=length(s_t)
      for ifn = 1:length(mase_fn)
         if size(clouds.masepdi(iday).(mase_fn{ifn}),1)==length(pdi_t)
            clouds.masepdi(iday).(mase_fn{ifn})=...
               clouds.masepdi(iday).(mase_fn{ifn})(ipdi_t,:);
         end
      end
   end
   % round the time anyway
   clouds.masepdi(iday).s_t=s_t;
   %%
   
   s_conc_pcasp = pcasp(ipcasp_t,:);
   s_conc_pcasp(s_conc_pcasp<0)=nan;
   s_meand_pcasp = arrayfun(@(x) wmean(pcasp_binmean,s_conc_pcasp(x,:)),...
      1:length(s_t))';
   s_ntot_pcasp = sum(s_conc_pcasp,2);
   
   s_ntot_pdi = clouds.masepdi(iday).s_ntot_pdi;
   s_ntot_pdi(isnan(s_ntot_pdi))=0;
   s_ntot_aer = s_ntot_pcasp + s_ntot_pdi;
   
   s_actfac=s_ntot_pdi./s_ntot_aer;
   
   clouds.masepdi(iday).s_ntot_pcasp = s_ntot_pcasp;
   clouds.masepdi(iday).s_meand_pcasp = s_meand_pcasp;
   clouds.masepdi(iday).s_ntot_aer = s_ntot_aer;
   clouds.masepdi(iday).s_actfrac = s_actfac;
   
   
end

%%
% save('clouds.mat','clouds', '-v7.3')