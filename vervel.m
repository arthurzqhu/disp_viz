close all
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
clearvars -except clouds
if ~exist('clouds','var') load clouds.mat, end
% ignoring the mase campaign because it has no vertical velocity records
campaigns={'vocalspdi','postpdi','oraclespdi','gomaccspdi'};
camp_proper = {'VOCALS','POST','ORACLES','GoMACCS'};
indvar_str = {'normAC', 's_actfrac', 's_ntot_aer', 's_ntot_pdi','AF'};
indvar_proper = {'NH', 'LAF', 'N_a', 'N_d', 'AF'};
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

   wz(wz<-10|wz>30)=nan;
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
w_str = {'strong updraft', 'weak updraft', 'quiescent',...
   'weak downdraft','strong downdraft'};
nw = length(w_str);

for iw=1:nw
   figure('Position',[1159 247 1422 730])
   tl=tiledlayout(nc,nvar,'TileSpacing','compact');
   for c=1:nc
      clear w_idx
      w_idx(1,:)=all_wz{c}>.5; % strong updraft index
      w_idx(2,:)=all_wz{c}<=.5 & all_wz{c}>=.1; % weak updraft index
      w_idx(3,:)=all_wz{c}<=.1 & all_wz{c}>=-.1; % quiescent index
      w_idx(4,:)=all_wz{c}>=-.5 & all_wz{c}<-.1; % weak downdraft index
      w_idx(5,:)=all_wz{c}<-.5; % strong downdraft index

      for ivar = 1:nvar
         nexttile
         scatter(all_X{c,ivar}(w_idx(iw,:)),all_reldisp{c}(w_idx(iw,:)),...
            10,all_wz{c}(w_idx(iw,:)),'filled','MarkerFaceAlpha',0.3)
         
         xlim([prctile(all_X{c,ivar}(w_idx(iw,:)),3) prctile(all_X{c,ivar}(w_idx(iw,:)),97)])
         ylim([0 1])
         colorbar
         colormap(rainbow)
         wz_bound=prctile(abs(all_wz{c}),95);
         caxis([ -wz_bound wz_bound ])
         if c==nc xlabel(indvar_proper{ivar}), end
         if ivar==1 ylabel(camp_proper{c}), end
         set(gca,'fontsize',14)
      end
      title(tl,w_str{iw},'fontsize',20,'fontweight','bold')
   end
   saveas(gcf,['plots/' w_str{iw} '.png'])
   pause(0.5)
end