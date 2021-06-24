clear
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
load clouds.mat
%%
pdi_files=dir('ORACLES_raw/PDI_P3/*.ict');
pdi_dates=arrayfun(@(x) pdi_files(x).name(1:6), 1:length(pdi_files),...
   'UniformOutput', false)';

binedge=readmatrix('ORACLES_raw/PDI_P3/bin_edges.csv','Range',[1 2]);
binmean=(binedge(:,2)+binedge(:,1))/2;
dlogD=log10(binedge(1,2)/binedge(1,1));

pcasp_files=dir('ORACLES_raw/PCASP_P3/*.ict');
pcasp_dates=arrayfun(@(x) pcasp_files(x).name(1:6), 1:length(pcasp_files),...
   'UniformOutput', false)';

hskp_files=dir('ORACLES_raw/Hskping_P3/*.ict');
hskp_dates=arrayfun(@(x) hskp_files(x).name(14:19), 1:length(hskp_files),...
   'UniformOutput', false)';

allinst_commondate=mintersect(pdi_dates, pcasp_dates, hskp_dates);
cmd_ipdi=find(ismember(pdi_dates, allinst_commondate));
cmd_ipcasp=find(ismember(pcasp_dates, allinst_commondate));
cmd_ihskp=find(ismember(hskp_dates, allinst_commondate));
% cmd_idxcoma=find(ismember(coma_dates, allinst_commondate));

% finishingTaskSound

%%

Dbins_100_grtr=find(binmean>100,1,'first');

try
   clouds=rmfield(clouds,'oraclespdi');
catch
end

%%
for iday=1:length(allinst_commondate)
   %% get the common date indices
   ipdi_d=cmd_ipdi(iday);
   ipcasp_d=cmd_ipcasp(iday);
   ihskp_d=cmd_ihskp(iday);

   %% get the files
   pdi{iday,1}(:,:)=readmatrix([pdi_files(ipdi_d).folder '/' ...
      pdi_files(ipdi_d).name],'FileType','text','NumHeaderLines',160);
   pcasp{iday,1}(:,:)=readmatrix([pcasp_files(ipcasp_d).folder '/' ...
      pcasp_files(ipcasp_d).name],'FileType','text','Range',[71 1]);
   hskp_header_line=csvread([hskp_files(ihskp_d).folder '/' ...
      hskp_files(ihskp_d).name],0, 0, [0 0 0 0]);
   hskp{iday,1}(:,:)=readmatrix([hskp_files(ihskp_d).folder '/' ...
      hskp_files(ihskp_d).name],'FileType','text','Range',[hskp_header_line 1]);

   %% get the common time
   pdi_t=pdi{iday}(:,1);
   pcasp_vidx=pcasp{iday}(:,32)>=0;
   pcasp{iday}=pcasp{iday}(pcasp_vidx,:);
   pcasp_t=pcasp{iday}(:,1);
   hskp_t=hskp{iday}(:,1);
   s_t=mintersect(pdi_t,pcasp_t,hskp_t);
   cmt_ipdi=find(ismember(pdi_t, s_t));
   cmt_ipcasp=find(ismember(pcasp_t, s_t));
   cmt_ihskp=find(ismember(hskp_t, s_t));

   %% get pdi data
   s_conc_pdi=pdi{iday}(cmt_ipdi,2:end);
   s_ntot_pdi=sum(s_conc_pdi*dlogD,2);
   s_ntot_100_grtr_pdi=sum(s_conc_pdi(:,Dbins_100_grtr:end)*dlogD,2);
   s_meand_pdi=sum(s_conc_pdi.*binmean'*dlogD,2)./s_ntot_pdi;
   s_lwc_pdi=sum(s_conc_pdi.*binmean'.^3*dlogD,2)*pi/6*1e-6;
   s_std_pdi=arrayfun(@(x) std(binmean, s_conc_pdi(x,:)), 1:length(s_t))';
   s_disp_pdi=s_std_pdi./s_meand_pdi;

   % invalid index
   ividx=s_ntot_pdi<25 | s_meand_pdi<3 | s_ntot_100_grtr_pdi>0;
   % set the lower bound for datapoints statistical properties
   s_std_pdi(ividx)=nan;
   s_disp_pdi(ividx)=nan;
   s_ntot_pdi(isnan(s_ntot_pdi))=0;
   s_disp_pdi(s_disp_pdi==0)=nan;
   %% get pcasp
   s_conc_pcasp=pcasp{iday}(cmt_ipcasp,2:31);
%    s_ntot_pcasp=sum(s_conc_pcasp,2);
   s_ntot_pcasp=pcasp{iday}(cmt_ipcasp,32);
   s_meand_pcasp=pcasp{iday}(cmt_ipcasp,36);
   s_ntot_aer=s_ntot_pcasp + s_ntot_pdi;
   s_actfrac=s_ntot_pdi./s_ntot_aer;
   %% get hskp
   s_rh=hskp{iday}(cmt_ihskp,end);
   s_rh(s_rh<0)=nan;
   s_ap=hskp{iday}(cmt_ihskp,6)*0.3048;
   s_wz=hskp{iday}(cmt_ihskp,12);
   s_thet=hskp{iday}(cmt_ihskp,19);
   s_ta=hskp{iday}(cmt_ihskp,18);
   s_p=hskp{iday}(cmt_ihskp,23);
   
   %calculating air temperature using ideal gas law
   s_p_raw=s_p*100;
   s_ta_raw=s_ta+273.15;
   Rd=287;
   s_rhoa=s_p_raw./(Rd*s_ta_raw);
   
   s_mr=hskp{iday}(cmt_ihskp,31);
   s_qv=s_mr/1000./(1+s_mr/1000)*1000;
   s_ql=s_lwc_pdi./s_rhoa;
   s_qt=s_ql + s_qv;
   
   clouds.oraclespdi(iday).file=pdi_files(iday).name;
   clouds.oraclespdi(iday).s_t=s_t;
   clouds.oraclespdi(iday).s_ntot_pdi=s_ntot_pdi;
   clouds.oraclespdi(iday).s_conc_pdi=s_conc_pdi;
   clouds.oraclespdi(iday).s_meand_pdi=s_meand_pdi;
   clouds.oraclespdi(iday).s_lwc_pdi=s_lwc_pdi;
   clouds.oraclespdi(iday).s_std_pdi=s_std_pdi;
   clouds.oraclespdi(iday).s_disp_pdi=s_disp_pdi;
   clouds.oraclespdi(iday).s_conc_cpasp=s_conc_pcasp;
   clouds.oraclespdi(iday).s_ntot_pcasp=s_ntot_pcasp;
   clouds.oraclespdi(iday).s_meand_pcasp=s_meand_pcasp;
   clouds.oraclespdi(iday).s_ntot_aer=s_ntot_aer;
   clouds.oraclespdi(iday).s_actfrac=s_actfrac;
   clouds.oraclespdi(iday).s_ap=s_ap;
   clouds.oraclespdi(iday).s_wz=s_wz;
   clouds.oraclespdi(iday).s_thet=s_thet;
   clouds.oraclespdi(iday).s_ta=s_ta;
   clouds.oraclespdi(iday).s_rh=s_rh;
   clouds.oraclespdi(iday).s_p=s_p;
   clouds.oraclespdi(iday).s_rhoa=s_rhoa;
   clouds.oraclespdi(iday).s_mr=s_mr;
   clouds.oraclespdi(iday).s_qv=s_qv;
   clouds.oraclespdi(iday).s_ql=s_ql;
   clouds.oraclespdi(iday).s_qt=s_qt;
end

%%
%
load oraclespdi_flight_basics.mat



for iday=1:length(allinst_commondate)
   
   cloudlegs_i=oraclespdi_flight_basics(iday).ti;
   cloudlegs_f=oraclespdi_flight_basics(iday).tf;
   z_CB=oraclespdi_flight_basics(iday).z_CB;
   z_CT=oraclespdi_flight_basics(iday).z_CT;
   
   s_t=clouds.oraclespdi(iday).s_t;
   s_lwc_pdi=clouds.oraclespdi(iday).s_lwc_pdi;
   s_ntot_aer=clouds.oraclespdi(iday).s_ntot_aer;
   s_ap=clouds.oraclespdi(iday).s_ap;
   s_ntot_pdi=clouds.oraclespdi(iday).s_ntot_pdi;
   s_disp_pdi=clouds.oraclespdi(iday).s_disp_pdi;
   
   if ~isempty(cloudlegs_i)
      for ileg=1:length(cloudlegs_i)
         ti=cloudlegs_i(ileg);
         tf=cloudlegs_f(ileg);
         
         ti_idx=findInSorted(s_t, ti);
         tf_idx=findInSorted(s_t, tf);
         
         if ti_idx<0
            continue
         end
         
         % ti and tf used for cloud related properties (reldisp, s_ntot, actfrac)
         % does not sample region below the cloud base
         ti_c=ti;
         tf_c=tf;
         
         % to sample some extra distance below the cloud top, in case
         % the flight doesnt only go from low to high
         if s_ap(ti_idx) < s_ap(tf_idx)
            ti=ti - 600;
         else
            tf=tf + 600;
         end
         
         z_min=min(s_ap(s_t < tf & s_t > ti));
         z_max=max(s_ap(s_t < tf & s_t > ti));
         z_max_sampled=(z_min+z_max)/2;
         
         % only sample the datapoints lower than mid cloud to prevent
         % sampling collision-coalescence
         
         aerCMS=@(x) calcMeanSampsize(x, ...
            s_t < tf & s_t > ti & s_ap < z_max_sampled);
         cldCMS=@(x) calcMeanSampsize(x, ...
            s_t < tf_c & s_t > ti_c & s_ap < z_max_sampled & s_ntot_pdi > 25);
         
         try
            [clouds.oraclespdi(iday).a_ntot_CB(ileg), ...
               clouds.oraclespdi(iday).a_ntot_CB_sampsize(ileg)]=...
               aerCMS(s_ntot_aer);
            [clouds.oraclespdi(iday).s_ntot_CB(ileg), ...
               clouds.oraclespdi(iday).s_ntot_CB_sampsize(ileg)]=...
               cldCMS(s_ntot_pdi);
            [clouds.oraclespdi(iday).s_actfrac_CB(ileg), ...
               clouds.oraclespdi(iday).s_actfrac_CB_sampsize(ileg)]=...
               cldCMS(s_ntot_pdi./s_ntot_aer);
            [clouds.oraclespdi(iday).reldisp_CB(ileg), ...
               clouds.oraclespdi(iday).reldisp_CB_sampsize(ileg)]=...
               cldCMS(s_disp_pdi);
            
         catch
            [clouds.oraclespdi(iday).a_ntot_CB(ileg), ...
               clouds.oraclespdi(iday).a_ntot_CB_sampsize(ileg), ...
               clouds.oraclespdi(iday).s_ntot_CB(ileg), ...
               clouds.oraclespdi(iday).s_ntot_CB_sampsize(ileg), ...
               clouds.oraclespdi(iday).s_actfrac_CB(ileg), ...
               clouds.oraclespdi(iday).s_actfrac_CB_sampsize(ileg), ...
               clouds.oraclespdi(iday).reldisp_CB(ileg), ...
               clouds.oraclespdi(iday).reldisp_CB_sampsize(ileg)]=...
               deal(nan);
         end
      end
   else
      [clouds.oraclespdi(iday).a_ntot_CB, ...
         clouds.oraclespdi(iday).a_ntot_CB_sampsize, ...
         clouds.oraclespdi(iday).s_ntot_CB, ...
         clouds.oraclespdi(iday).s_ntot_CB_sampsize, ...
         clouds.oraclespdi(iday).s_actfrac_CB, ...
         clouds.oraclespdi(iday).s_actfrac_CB_sampsize, ...
         clouds.oraclespdi(iday).reldisp_CB, ...
         clouds.oraclespdi(iday).reldisp_CB_sampsize]=...
         deal(nan);
   end
   
end

%% Save the data
leg_inspection_oracles;
save('clouds.mat','clouds', '-v7.3')

% finishingTaskSound
