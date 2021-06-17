close all
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
clearvars -except clouds
if ~exist('clouds','var') load clouds.mat, end
% ignoring the mase campaign because it has no vertical velocity records
campaigns={'vocalspdi','postpdi','oraclespdi','gomaccspdi'};
camp_proper = {'VOCALS','POST','ORACLES','GoMACCS'};
indvar_str = {'normAC', 's_actfrac', 's_ntot_aer', 's_lwc_pdi','AF'};
indvar_proper = {'NH', 'LAF', 'N_a', 'LWC', 'AF'};

colororder=colororder;
color_order={};
for i=1:size(colororder,1)
    color_order{i}=colororder(i,:);
end

%%
nc = length(campaigns);
nvar = length(indvar_str);

all_reldisp = cell(nc,1);
all_X = cell(nc,nvar);
all_wz = cell(nc,1);


for c = 1:nc
   indvar_flt = [];
   reldisp = [];
   wz = [];
   camp = campaigns{c};

   for ivar = 1:nvar
      for iday = 1:length(clouds.(camp))
         XX = clouds.(camp)(iday).(indvar_str{ivar});
         indvar_flt = [indvar_flt; XX];
      end
   end
   indvar_mat = reshape(indvar_flt,[],nvar);

   for iday = 1:length(clouds.(camp))
      reldisp = [reldisp;clouds.(camp)(iday).s_disp_pdi];
      wz = [wz;clouds.(camp)(iday).s_wz];
   end

   wz(wz<-20|wz>20)=nan;
   vidx = ~isnan(sum([indvar_mat reldisp],2));
   all_reldisp{c} = reldisp(vidx);
   all_wz{c} = wz(vidx);

   for ivar = 1:nvar
      all_X{c,ivar} = indvar_mat(vidx,ivar);
   end

end

%%
close all
rainbow=getPyPlot_cMap('rainbow',11);
coolwarm=getPyPlot_cMap('coolwarm',11);

w_str = {'strong updraft (>80th)', 'weak updraft (60-80th)',...
   'quiescent (40-60th)','weak downdraft (20-40th)','strong downdraft (<20th)'};
nw = length(w_str);

for iw=1:nw
   close all
   figure('Position',[1159 247 1422 730])
   tl=tiledlayout(nc,nvar,'TileSpacing','compact');
   for c=1:nc
      clear w_idx
      
      % get adjusted 20th, 40th, 60th, 80th percentile of wind velocity and 
      % divide into 5 groups
      w_dd_idx=all_wz{c}<=0; % downdraft index
      w_ud_idx=all_wz{c}>0; % updraft index
      
      % e.g., 40th percentile of downdraft = (adjusted) 20th percentile of all
      % vertical wind
      w_20_val=prctile(all_wz{c}(w_dd_idx),40); 
      w_40_val=prctile(all_wz{c}(w_dd_idx),80); 
      w_60_val=prctile(all_wz{c}(w_ud_idx),20); 
      w_80_val=prctile(all_wz{c}(w_ud_idx),60); 
      
      w_idx(1,:)=all_wz{c}>w_80_val; % strong updraft index
      w_idx(2,:)=all_wz{c}>w_60_val & all_wz{c}<=w_80_val; % weak updraft index
      w_idx(3,:)=all_wz{c}>w_40_val & all_wz{c}<=w_60_val; % quiescent index
      w_idx(4,:)=all_wz{c}>w_20_val & all_wz{c}<=w_40_val; % weak downdraft index
      w_idx(5,:)=all_wz{c}<=w_20_val; % strong downdraft index

      for ivar = 1:nvar
         nexttile
         
         XX=all_X{c,ivar}(w_idx(iw,:));
         YY=all_reldisp{c}(w_idx(iw,:));
         
         scatter(XX,YY,10,all_wz{c}(w_idx(iw,:)),'filled')
         
         if any(strcmp(indvar_str{ivar},{'NH','LAF','AF'}))
            xlim([0 1])
         else
            xlim([prctile(XX,3) prctile(XX,97)])
         end
         
         ylim([0 1])
         cbar=colorbar;
         cbar.Label.String='w [m/s]';
         colormap(coolwarm)
         wz_bound=prctile(abs(all_wz{c}),95);
         caxis([ -wz_bound wz_bound ])
         if c==nc xlabel(indvar_proper{ivar}), end
         if ivar==1 ylabel(camp_proper{c}), end
         
         % draw mean and errorbar
         xbins = linspace(min(xlim),max(xlim),11);
         mu_reldisp = zeros(length(xbins)-1,1);
         
         for ibin = 1:length(xbins)-1
            idx_at_ibin = XX>xbins(ibin) ...
               & XX<=xbins(ibin+1);
            
            if sum(idx_at_ibin)>25
               mu_reldisp(ibin) = nanmean(YY(idx_at_ibin));
               sigma_dv(ibin) = nanstd(YY(idx_at_ibin));
            else
               mu_reldisp(ibin) = nan;
               sigma_dv(ibin) = nan;
            end
         end
         
         hold on
         xbinmean = (xbins(2:end)+xbins(1:end-1))/2;
         errorbar(xbinmean,mu_reldisp,sigma_dv,'LineWidth',2,...
            'color','k')
         hold off
         
         set(gca,'fontsize',14)
      end
      title(tl,w_str{iw},'fontsize',20,'fontweight','bold')
   end
   saveas(gcf,['plots/' w_str{iw} '.png'])
   pause(0.5)
end