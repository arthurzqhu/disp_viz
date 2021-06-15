close all
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
clearvars -except clouds
if ~exist('clouds','var') load clouds.mat, end
post_aer_files = dir('post/POST_pcasp/pcasp_1hz_*');

%% get rid of the old fieldnames first
for iday = 1:length(clouds.postpdi)
   post_fn=fieldnames(clouds.postpdi(iday));
   ext_idx = startsWith(post_fn,'a_') | startsWith(post_fn,'ccn_');
   newstrct(iday)=rmfield(clouds.postpdi(iday),post_fn(ext_idx));
end
clouds.postpdi=newstrct;
%%
pdi_date_start_char = regexp(clouds.postpdi(1).file,'/');

if isempty(pdi_date_start_char)
   pdi_date_start_char = 0;
end

pdi_dates = arrayfun(@(x) ...
   clouds.postpdi(x).file(pdi_date_start_char+1:pdi_date_start_char+6),...
   1:length(clouds.postpdi), 'UniformOutput', false)';
pcasp_dates = arrayfun(@(x) post_aer_files(x).name(18:23), ...
   1:length(post_aer_files),'UniformOutput', false)';
pcasp_binedges = readmatrix([post_aer_files(1).folder '/pcasp_bnds.txt'],...
   'Range',[2 1 21 2]);
pcasp_binmean = readmatrix([post_aer_files(1).folder '/pcasp_bnds.txt'],...
   'Range',[2 3 21 3]);

allinst_commondates = mintersect(pdi_dates, pcasp_dates); 
cmd_ipdi = find(ismember(pdi_dates, allinst_commondates));
cmd_ipcasp = find(ismember(pcasp_dates, allinst_commondates));

%%
for iday = 1:length(allinst_commondates)
   ipdi_d = cmd_ipdi(iday);
   ipcasp_d = cmd_ipcasp(iday);
   
   %% get the common time
   pcasp = readmatrix([post_aer_files(ipcasp_d).folder '/' ...
      post_aer_files(ipcasp_d).name],'FileType','text','Range',[13 2]);
   pcasp_t_len = size(pcasp,1);
   pcasp(pcasp<0)=nan;
   pcasp_t_str = readmatrix([post_aer_files(ipcasp_d).folder '/' ...
      post_aer_files(ipcasp_d).name],'FileType','text','Range',...
      [13 1 pcasp_t_len+13 1], 'OutputType', 'string');
   pcasp_dtnum = datenum(pcasp_t_str,'yymmddHHMMSS');
   pcasp_t = round((pcasp_dtnum-floor(pcasp_dtnum))*86400,3);
   pdi_t = round(clouds.postpdi(iday).s_t);
   [s_t,ipdi_t,ipcasp_t] = intersect(pdi_t,pcasp_t);
   
   %% downsample the existing s_* variables first
   post_fn=fieldnames(clouds.postpdi(iday));
   if length(clouds.postpdi(iday).s_t)~=length(s_t)
      for ifn = 1:length(post_fn)
         if size(clouds.postpdi(iday).(post_fn{ifn}),1)==length(pdi_t)
            clouds.postpdi(iday).(post_fn{ifn})=...
               clouds.postpdi(iday).(post_fn{ifn})(ipdi_t,:);
         end
      end
   end
   % round the time anyway
   clouds.postpdi(iday).s_t=s_t;
   %%
   s_conc_pcasp = pcasp(ipcasp_t,:);
   s_conc_pcasp(s_conc_pcasp<0)=nan;
   s_meand_pcasp = arrayfun(@(x) wmean(pcasp_binmean,s_conc_pcasp(x,:)),...
      1:length(s_t))';
   s_ntot_pcasp = sum(s_conc_pcasp,2);
   
   s_ntot_pdi = clouds.postpdi(iday).s_ntot_pdi;
   s_ntot_pdi(isnan(s_ntot_pdi))=0;
   s_ntot_aer = s_ntot_pcasp + s_ntot_pdi;
   
   s_actfac=s_ntot_pdi./s_ntot_aer;
   
   clouds.postpdi(iday).s_ntot_pcasp = s_ntot_pcasp;
   clouds.postpdi(iday).s_meand_pcasp = s_meand_pcasp;
   clouds.postpdi(iday).s_ntot_aer = s_ntot_aer;
   clouds.postpdi(iday).s_actfrac = s_actfac;
   
end

%%
% save('clouds.mat','clouds', '-v7.3')