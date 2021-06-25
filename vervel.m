close all
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
clearvars -except clouds
if ~exist('clouds','var') load clouds.mat, end
% ignoring the mase campaign because it has no vertical velocity records
campaigns={'vocalspdi','postpdi','oraclespdi','gomaccspdi'};
camp_proper={'VOCALS','POST','ORACLES','GoMACCS'};
indvar_str={'normAC', 's_actfrac', 's_ntot_aer', 's_ntot_pcasp',...
   's_lwc_pdi','AF'};
indvar_proper={'NH', 'LAF', 'N_a', 'N_{ia}','LWC', 'AF'};

colororder=colororder;
color_order={};
for i=1:size(colororder,1)
   color_order{i}=colororder(i,:);
end

%%
close all
rainbow=getPyPlot_cMap('rainbow',11);
coolwarm=getPyPlot_cMap('coolwarm',11);

nc=length(campaigns);
nvar=length(indvar_str);
figure('Position',[1324 247 981 730])

for c=1:nc
   camp=campaigns{c};
   datestrs=regexp(clouds.(camp)(1).file,'\d+');
   datestre=datestrs+5;

   for iday=1:length(clouds.(camp))
      tl=tiledlayout('flow');
      wz=clouds.(camp)(iday).s_wz;
      wz(wz<-20|wz>20)=nan;
      wz_bound=prctile(abs(wz),95);
      reldisp=clouds.(camp)(iday).s_disp_pdi;
      date_str=clouds.(camp)(iday).file(datestrs:datestre);
      for ivar=1:nvar
         XX=clouds.(camp)(iday).(indvar_str{ivar});

         nexttile
         scatter(XX,reldisp,10,wz,'filled','MarkerFaceAlpha',0.7)
         ylim([0 1])

         xlabel(indvar_proper{ivar})
         ylabel(tl,'\epsilon')
         caxis([-wz_bound wz_bound])
         cbar=colorbar;
         cbar.Label.String='w [m/s]';
         colormap(coolwarm)
         set(gca,'fontsize',16)

      end

      title(tl,[camp_proper{c} ' ' date_str],...
            'fontsize',20,'fontweight','bold')
      saveas(gcf,['plots/vervel/' camp_proper{c} ' ' date_str ...
            '.png'])
   end
end
