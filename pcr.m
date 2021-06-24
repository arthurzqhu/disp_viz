close all
clearvars -except clouds
if ~exist('clouds','var') load clouds.mat, end
campaigns={'vocalspdi','masepdi','postpdi','oraclespdi','gomaccspdi'};
camp_proper = {'VOCALS','MASE','POST','ORACLES','GoMACCS'};
nc=length(campaigns);

dp_color='#FFBE89';%'#A8CCCC';

mosaic_rat=[1 5 10 20];
%% PCR
close all


demean = @(x) (x - nanmean(x));
zscore = @(x) (x - nanmean(x))./nanstd(x);

for ir=1:length(mosaic_rat)
   figure('Position',[1010 198 1415 779])
   tl1=tiledlayout('flow');
   figure('Position',[1010 198 1415 779])
   tl2=tiledlayout('flow');
   for c=1:nc
      %%
      clear XX YY

      camp=campaigns{c};
      [reldisp,lwc,actfrac,ntot_aer,ntot_pdi,ntot_pcasp,AF,NH,ETqt,ETT,w,...
         ]=...
         deal([]);

      for iday=1:length(clouds.(camp))

         reldisp_raw=clouds.(camp)(iday).s_disp_pdi;
         vidx=~isnan(reldisp_raw);
         reldisp=[reldisp;reldisp_raw(vidx)];

         lwc=[lwc;clouds.(camp)(iday).s_lwc_pdi(vidx)];
         actfrac=[actfrac;clouds.(camp)(iday).s_actfrac(vidx)];
         ntot_pdi=[ntot_pdi;clouds.(camp)(iday).s_ntot_pdi(vidx)];
         ntot_aer=[ntot_aer;clouds.(camp)(iday).s_ntot_aer(vidx)];
         ntot_pcasp=[ntot_pcasp;clouds.(camp)(iday).s_ntot_pcasp(vidx)];
         AF=[AF;clouds.(camp)(iday).AF(vidx)];
         NH=[NH;clouds.(camp)(iday).normAC(vidx)];
         ETqt=[ETqt;clouds.(camp)(iday).ent_ratio_qt(vidx)];
         ETT=[ETT;clouds.(camp)(iday).ent_ratio_T(vidx)];
         %w=[w;clouds.(camp)(iday).s_wz(vidx)];

      end

      YY_raw=reldisp;
      XX_raw=[lwc actfrac AF NH];

      YY=mosaicify(YY_raw,mosaic_rat(ir));
      for ip=1:size(XX_raw,2)
         XX(:,ip)=mosaicify(XX_raw(:,ip),mosaic_rat(ir));
      end

      XX_pos=XX;
      XX_pos(XX<=0)=nan;
      YYlog=log10(YY);
      XXlog=log10(XX_pos);

      XX_demean=zscore(XX);
      XXlog_demean=zscore(XXlog);

      [PCAcoeff{ir,c},~,PCAlatt]=pca(XX_demean);
      PCAexp=PCAlatt/sum(PCAlatt);

      [PCAcoeff_log{ir,c},~,PCAlatt_log]=pca(XXlog_demean);
      PCAexp_log=PCAlatt_log/sum(PCAlatt_log);
      % nfeat=3;
      nfeat=size(XX,2);


      XX_new=[ones(size(XX,1),1) XX*PCAcoeff{ir,c}(:,1:nfeat)];
      XXlog_new=[ones(size(XX,1),1) XXlog*PCAcoeff_log{ir,c}(:,1:nfeat)];


      [blog{ir,c},blogint,~,~,statslog]=regress(YYlog,XXlog_new);
      Yhatlog=10.^(XXlog_new*blog{ir,c});

      [b{ir,c},bint,~,~,stats]=regress(YY,XX_new);
      Yhat=XX_new*b{ir,c};
      
      figure(1)
      nexttile
      scatter(Yhatlog,YY,5,'filled','MarkerFaceColor',dp_color,...
         'MarkerFaceAlpha',0.3)
      pred_str = ['R^2','=', sprintf('%0.3f',...
                           statslog(1))];
      text(.7,.6,pred_str,'FontSize',16)
      xlim([0 1])
      ylim([0 1])
      refline(1,0)
      xlabel(tl1,'predicted \epsilon','fontsize',24)
      ylabel(tl1,'observed \epsilon','fontsize',24)
      set(gca,'fontsize',20)
      if ir==length(mosaic_rat) xlabel([camp_proper{c} ' \epsilon']), end
      if c==1 ylabel(['ratio=' num2str(mosaic_rat(ir))]), end

      figure(2)
      nexttile
      scatter(Yhat,YY,5,'filled','MarkerFaceColor',dp_color,...
         'MarkerFaceAlpha',0.3)
      pred_str = ['R^2','=', sprintf('%0.3f',...
                           stats(1))];
      text(.7,.6,pred_str,'FontSize',16)
      xlim([0 1])
      ylim([0 1])
      refline(1,0)
      xlabel(tl2,'predicted \epsilon','fontsize',24)
      ylabel(tl2,'observed \epsilon','fontsize',24)
      set(gca,'fontsize',20)
      if ir==length(mosaic_rat) xlabel([camp_proper{c} ' \epsilon']), end
      if c==1 ylabel(['ratio=' num2str(mosaic_rat(ir))]), end


   end
end

%%
saveas(figure(1),'plots/PCR scale no w power.png')
saveas(figure(2),'plots/PCR scale no w linear.png')
