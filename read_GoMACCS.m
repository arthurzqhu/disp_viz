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

try
   clouds = rmfield(clouds,'gomaccspdi');
catch
end

for iday = 1:length(pdi_files)
   %% get the files location
   conc_dat{iday,1}(:,:) = readmatrix([pdi_files(iday).folder '/' pdi_files(iday).name],...
      'FileType','text');
   cabin{iday,1}(:,:) = readmatrix([cabin_files(iday).folder '/' cabin_files(iday).name],...
      'FileType','text','Range',[2 1]);
   
   %% get the common time
   
   pdi_t = conc_dat{iday}(:,1);
   cabin_t = cabin{iday}(:,1);
   [~, cmt_ipdi,cmt_icabin]=intersect(pdi_t, cabin_t);
   
   %% get pdi data
   s_t = pdi_t(cmt_ipdi);
   s_conc_pdi = conc_dat{iday}(cmt_ipdi,2:end);
   s_ntot_pdi = sum(s_conc_pdi*dlogD, 2);
   s_meand_pdi = sum(s_conc_pdi.*binmean'*dlogD,2)./s_ntot_pdi;
   s_lwc_pdi = sum(s_conc_pdi.*binmean'.^3*dlogD,2)*pi/6*1e-6;
   s_std_pdi = arrayfun(@(x) std(binmean, s_conc_pdi(x,:)), 1:size(s_conc_pdi,1))';
   s_disp_pdi = s_std_pdi./s_meand_pdi;
   
   % invalid index
   ividx = s_ntot_pdi<25 | s_meand_pdi < 3;
   % set the lower bound for datapoints statistical properties
   s_std_pdi(ividx) = nan;
   s_disp_pdi(ividx) = nan;
   
   clouds.gomaccspdi(iday).file = pdi_files(iday).name(25:end-4);
   %% get cabin data
   
   s_ntot_pcasp = cabin{iday}(cmt_icabin,22);
   s_lat = cabin{iday}(cmt_icabin,3); s_lat(s_lat==0)=nan;
   s_lon = -cabin{iday}(cmt_icabin,4); s_lon(s_lon==0)=nan;
   
   % fix the nan if it's during the sampling
   nan_tlist=find(isnan(s_lon));
   
   for inan=1:sum(isnan(s_lon))
      prev_idx=nan_tlist(inan)-1;
      next_idx=nan_tlist(inan)+1;
      
      s_lat(nan_tlist(inan))=(s_lat(prev_idx)+s_lat(next_idx))/2;
      s_lon(nan_tlist(inan))=(s_lon(prev_idx)+s_lon(next_idx))/2;
   end
   
   
   s_Vtot_ex = cabin{iday}(cmt_icabin,23);
   s_wz = cabin{iday}(cmt_icabin,19);
   s_ap = cabin{iday}(cmt_icabin,16);
   s_p = cabin{iday}(cmt_icabin,15);
   s_T = cabin{iday}(cmt_icabin,12);
   s_rh = cabin{iday}(cmt_icabin,14);
   s_lwc_gerb = cabin{iday}(cmt_icabin,33);
   s_lwc_wire = cabin{iday}(cmt_icabin,34);
   s_T(s_T==0)=nan;
   
   s_thet = (s_T+273.15).*(1000./s_p).^(287/1004);
   
   s_meanV_ex = s_Vtot_ex./s_ntot_pcasp;
   s_meanD_ex = (s_meanV_ex/(pi/6)).^(1/3);
   
   %% calculating s_qt
   Lv = 2.5e6;
   Rv = 462;
   Rd = 287;
   
   s_Tk = s_T+273.15;
   s_es = 6.11*exp(Lv/Rv*(1/273.15-1./s_Tk));
   s_e = s_es.*s_rh/100;
   s_mr = s_e*Rd./(Rv*(s_p-s_e));
   s_qv = s_mr./(1+s_mr)*1000;
   s_rhoa = s_p*100./(Rd*s_Tk);
   
   s_ql = s_lwc_pdi./s_rhoa;
   s_qt = s_ql + s_qv;
   
   
   clouds.gomaccspdi(iday).s_qt = s_qt;
   clouds.gomaccspdi(iday).s_t = s_t;
   clouds.gomaccspdi(iday).s_ntot_pdi = s_ntot_pdi;
   clouds.gomaccspdi(iday).s_conc_pdi = s_conc_pdi;
   clouds.gomaccspdi(iday).s_meand_pdi = s_meand_pdi;
   clouds.gomaccspdi(iday).s_lwc_pdi = s_lwc_pdi;
   clouds.gomaccspdi(iday).s_std_pdi = s_std_pdi;
   clouds.gomaccspdi(iday).s_disp_pdi = s_disp_pdi;
   clouds.gomaccspdi(iday).s_lat = s_lat;
   clouds.gomaccspdi(iday).s_lon = s_lon;
   clouds.gomaccspdi(iday).s_mr = s_mr*1000;
   clouds.gomaccspdi(iday).s_ntot_ex = s_ntot_pcasp;
   clouds.gomaccspdi(iday).s_meanD_ex = s_meanD_ex;
   clouds.gomaccspdi(iday).s_rhoa = s_rhoa;
   clouds.gomaccspdi(iday).s_ap = s_ap;
   clouds.gomaccspdi(iday).s_wz = s_wz;
   clouds.gomaccspdi(iday).s_p = s_p;
   clouds.gomaccspdi(iday).s_T = s_T;
   clouds.gomaccspdi(iday).s_thet = s_thet;
   clouds.gomaccspdi(iday).s_rh = s_rh;
   clouds.gomaccspdi(iday).s_lwc_gerb = s_lwc_gerb;
   clouds.gomaccspdi(iday).s_lwc_wire = s_lwc_wire;
   
   s_ntot_aer = s_ntot_pcasp + s_ntot_pdi;
   clouds.gomaccspdi(iday).s_ntot_aer = s_ntot_aer;
   clouds.gomaccspdi(iday).s_ql = s_ql;
   clouds.gomaccspdi(iday).s_actfrac = s_ntot_pdi./s_ntot_aer;
end
%%

load gomaccspdi_flight_basics.mat

for iday = 1:length(pdi_files)
   
   cloudlegs_i = gomaccspdi_flight_basics(iday).ti;
   cloudlegs_f = gomaccspdi_flight_basics(iday).tf;
   
   s_t = clouds.gomaccspdi(iday).s_t;
   s_lwc_pdi = clouds.gomaccspdi(iday).s_lwc_pdi;
   s_ntot_aer = clouds.gomaccspdi(iday).s_ntot_aer;
   s_ap = clouds.gomaccspdi(iday).s_ap;
   s_disp_pdi = clouds.gomaccspdi(iday).s_disp_pdi;
   s_ntot_pdi = clouds.gomaccspdi(iday).s_ntot_pdi;
   s_actfrac = clouds.gomaccspdi(iday).s_actfrac;
   
   if ~isempty(cloudlegs_i)
      for ileg = 1:length(cloudlegs_i)
         ti = cloudlegs_i(ileg);
         tf = cloudlegs_f(ileg);
         ti_idx = findInSorted(s_t, ti);
         tf_idx = findInSorted(s_t, tf);
         
         if ti_idx<0
            continue
         end
         
         % ti and tf used for cloud related properties (reldisp, s_ntot, actfrac)
         % does not sample region below the cloud base
         ti_c = ti;
         tf_c = tf;
         
         % to sample some extra distance below the cloud top, in case
         % the flight doesnt only go from low to high
         if s_ap(ti_idx) < s_ap(tf_idx)
            ti = ti - 300;
         else
            tf = tf + 300;
         end
         
         z_min = min(s_ap(s_t < tf & s_t > ti));
         z_max = max(s_ap(s_t < tf & s_t > ti));
         z_max_sampled = (z_min+z_max)/2;
         
         aerCMS = @(x) calcMeanSampsize(x, s_t < tf & s_t > ti & s_ap < z_max_sampled);
         cldCMS = @(x) calcMeanSampsize(x, s_t < tf_c & s_t > ti_c & s_ap < z_max_sampled & s_ntot_pdi > 25);
         
         [clouds.gomaccspdi(iday).a_ntot_CB(ileg), clouds.gomaccspdi(iday).a_ntot_CB_sampsize(ileg)] = ...
            aerCMS(s_ntot_aer);
         [clouds.gomaccspdi(iday).s_ntot_CB(ileg), clouds.gomaccspdi(iday).s_ntot_CB_sampsize(ileg)] = ...
            cldCMS(s_ntot_pdi);
         [clouds.gomaccspdi(iday).s_actfrac_CB(ileg), clouds.gomaccspdi(iday).s_actfrac_CB_sampsize(ileg)] = ...
            cldCMS(s_actfrac);
         [clouds.gomaccspdi(iday).reldisp_CB(ileg), clouds.gomaccspdi(iday).reldisp_CB_sampsize(ileg)] = ...
            cldCMS(s_disp_pdi);
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