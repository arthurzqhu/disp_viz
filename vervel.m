close all
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
clearvars -except clouds
if ~exist('clouds','var') load clouds.mat, end
% ignoring the mase campaign because it has no vertical velocity records
campaigns={'vocalspdi','postpdi','oraclespdi','gomaccspdi'};
camp_proper={'VOCALS','POST','ORACLES','GoMACCS'};
indvar_str={'s_wz', 's_actfrac'};
indvar_proper={'w', 'LAF'};
notwellmixed={[6 8 9 10 11];[];[];[]};

colororder=colororder;
color_order={};
for i=1:size(colororder,1)
   color_order{i}=colororder(i,:);
end

%%
close all
rainbow=getPyPlot_cMap('rainbow',11);
coolwarm=getPyPlot_cMap('coolwarm',11);
coolwarm_r=getPyPlot_cMap('coolwarm_r',11);

nc=length(campaigns);
nvar=length(indvar_str);
figure('Position',[963 790 1169 277])

for c=1%:nc
   camp=campaigns{c};
   datestrs=regexp(clouds.(camp)(1).file,'\d+');
   datestre=datestrs+5;

   for iday=7%:length(clouds.(camp))
      tl=tiledlayout('flow');
      wz=clouds.(camp)(iday).s_wz;
      wz(wz<-20|wz>20)=nan;
      nh=clouds.(camp)(iday).normAC;
      rh=clouds.(camp)(iday).s_rh;
      laf=clouds.(camp)(iday).s_actfrac;

      vidx=~isnan(nh+wz+laf);
      wz_bound=prctile(abs(wz),95);
      reldisp=clouds.(camp)(iday).s_disp_pdi;
      date_str=clouds.(camp)(iday).file(datestrs:datestre);

      if ismember(iday,notwellmixed{c})
         date_title=[date_str '*'];
      else
         date_title=date_str;
      end

      for ivar=1:nvar
         XX=clouds.(camp)(iday).(indvar_str{ivar});
         [~,~,~,~,stats]=regress(reldisp,[ones(length(XX),1) XX]);
         rsqstr=['R^2=',sprintf('%0.3f',stats(1))];

         nexttile
         if strcmp(indvar_proper{ivar},'w')
            scatter(XX(vidx),reldisp(vidx),10,nh(vidx),'filled','MarkerFaceAlpha',0.7)
            caxis([0 1])
            cbar=colorbar;
            cbar.Label.String='NH';
         else
            scatter(XX(vidx),reldisp(vidx),10,wz(vidx),'filled','MarkerFaceAlpha',0.7)
            caxis([-wz_bound wz_bound])
            cbar=colorbar;
            cbar.Label.String='w';
         end
         ylim([0 1])

         axisdim=axis;
         text((axisdim(2)-axisdim(1))*0.7+axisdim(1),...
                 (axisdim(4)-axisdim(3))*0.6+axisdim(3),...
                 rsqstr,'fontsize',16)

         xlabel(indvar_proper{ivar})
         ylabel(tl,'\epsilon','fontsize',20)
         
         
         set(gca,'fontsize',20)

      end
      
      nexttile
      scatter(XX(vidx),reldisp(vidx),10,nh(vidx),'filled','MarkerFaceAlpha',0.7)
      caxis([0 1])
      cbar=colorbar;
      cbar.Label.String='NH';
      colormap(coolwarm_r)
      ylim([0 1])
      axisdim=axis;
      text((axisdim(2)-axisdim(1))*0.7+axisdim(1),...
              (axisdim(4)-axisdim(3))*0.6+axisdim(3),...
              rsqstr,'fontsize',16)

      xlabel(indvar_proper{ivar})
      
      set(gca,'fontsize',20)
      
      title(tl,[camp_proper{c} ' ' date_title],...
            'fontsize',24,'fontweight','bold')
      saveas(gcf,['plots/' camp_proper{c} ' ' date_str ...
            '.png'])
   end
end
