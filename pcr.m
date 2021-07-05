close all
clearvars -except clouds
if ~exist('clouds','var') load clouds.mat, end
campaigns={'vocalspdi','masepdi','postpdi','oraclespdi','gomaccspdi'};
camp_proper={'VOCALS','MASE','POST','ORACLES','GoMACCS'};
nc=length(campaigns);

dp_color='#FFBE89';%'#A8CCCC';
psz=20;
falp=[.1 .4 .7 1];
arat=[1 10 100 1000];
downsample_rat=arat;
train_ratio=1;
%% PCR
close all

findpred=@(yhat,yy) 1-nansum((yhat-yy).^2)/nansum((yy-nanmean(yy)).^2);
demean=@(x) (x - nanmean(x));
zscore=@(x) (x - nanmean(x))./nanstd(x);

% figure('Position',[1010 198 1415 779])
% tl1=tiledlayout('flow');
figure('Position',[1010 198 1415 779])
tl2=tiledlayout('flow');

for ir=4%1:length(arat)
   
   for c=3%1:nc
      
      clear XX YY mlwc mpdi mAF mpcasp maer mactfrac mreldisp
      
      binedges=load(['bin_edges_' campaigns{c}(1:end-3) '.csv']);
      binmean=(binedges(:,2)+binedges(:,3))/2;
      dlogD=log10(binedges(2,3)/binedges(2,2));

      camp=campaigns{c};
      ndays=length(clouds.(camp));
      days_analyzed=1:ndays;
      %     remove the days that have incomplete flights in vocals
      if strcmp(camp,'vocalspdi')
         %                 days_analyzed(ismember(days_analyzed,[8,10,11,13]))=[];
      elseif strcmp(camp,'oraclespdi')
         days_analyzed(ismember(days_analyzed, [2,3,5,14,15,16,17,18]))=[];
      end

      [s_lwc_pdi,normAC,s_ntot_aer,s_ntot_pdi,s_ntot_pcasp,ql_adb,rhoa,...
         s_actfrac_imp,AF_imp,reldisp,reldisp_imp]=deal([]);

      [m_lwc_pdi,m_ntot_pdi,m_ntot_aer,m_actfrac,m_reldisp,m_AF,m_rhoa,...
         m_ql_adb]=deal(nan(length(days_analyzed),1));
      
      for iday=days_analyzed
         clear reldisp_sd
         % indvars
         s_lwc_pdi=[s_lwc_pdi;mosaicify(clouds.(camp)(iday).s_lwc_pdi,arat(ir))];
         normAC=[normAC;mosaicify(clouds.(camp)(iday).normAC,arat(ir))];
         s_ntot_pdi=[s_ntot_pdi;mosaicify(clouds.(camp)(iday).s_ntot_pdi,arat(ir))];
         s_ntot_aer=[s_ntot_aer;mosaicify(clouds.(camp)(iday).s_ntot_aer,arat(ir))];
         s_ntot_pcasp=[s_ntot_pcasp;mosaicify(clouds.(camp)(iday).s_ntot_pcasp,arat(ir))];
         s_actfrac_imp=[s_actfrac_imp;mosaicify(clouds.(camp)(iday).s_actfrac,arat(ir))];
         ql_adb=[ql_adb;mosaicify(clouds.(camp)(iday).ql_adb_prof,arat(ir))];
         rhoa=[rhoa;mosaicify(clouds.(camp)(iday).s_rhoa,arat(ir))];
         AF_imp=[AF_imp;mosaicify(clouds.(camp)(iday).AF,arat(ir))];
         
         vidx=clouds.(camp)(iday).s_ntot_pdi>10;
         m_lwc_pdi(iday)=nanmean(clouds.(camp)(iday).s_lwc_pdi(vidx));
         m_ntot_pdi(iday)=nanmean(clouds.(camp)(iday).s_ntot_pdi(vidx));
         m_ntot_aer(iday)=nanmean(clouds.(camp)(iday).s_ntot_aer(vidx));
         m_actfrac(iday)=m_ntot_pdi(iday)/m_ntot_aer(iday);
         m_rhoa(iday)=nanmean(clouds.(camp)(iday).s_rhoa(vidx));
         m_ql_adb(iday)=nanmean(clouds.(camp)(iday).ql_adb_prof(vidx));
%          m_AF_imp=nanmean(clouds.(camp)(iday).AF(vidx));
         
         m_conc_pdi=nanmean(clouds.(camp)(iday).s_conc_pdi(vidx,:));
         m_reldisp(iday)=std(binmean,m_conc_pdi,'omitnan')./...
               wmean(binmean,m_conc_pdi);
         m_AF(iday)=m_lwc_pdi(iday)./m_rhoa(iday)./m_ql_adb(iday);
         
         
         % depvar
         s_conc_pdi=mosaicify(clouds.(camp)(iday).s_conc_pdi,arat(ir));
         for itime=1:size(s_conc_pdi,1)
            ntot=sum(s_conc_pdi(itime,:))*dlogD;
            if ntot<10 
               reldisp_sd(itime)=nan;
            else
               reldisp_sd(itime)=std(binmean,s_conc_pdi(itime,:),'omitnan')./...
               wmean(binmean,s_conc_pdi(itime,:));
            end
         end
         
         reldisp=[reldisp;reldisp_sd'];
         reldisp_imp=[reldisp_imp;...
            mosaicify(clouds.(camp)(iday).s_disp_pdi,arat(ir))];
      end
      
      s_actfrac=s_ntot_pdi./s_ntot_aer;
      ql_obs=s_lwc_pdi./rhoa;
      AF=ql_obs./ql_adb;
      
      YY=reldisp_imp;
      XX=[s_lwc_pdi s_actfrac_imp normAC AF_imp];
      
      mYY=m_reldisp;
      mXX=[m_lwc_pdi m_actfrac];
      
%       ndat=length(reldisp);
%       downsample_idx=randsample(ndat,floor(ndat/downsample_rat(ir)));
      
      tot_sampsz=size(XX,1);
      tr_sampsz=ceil(tot_sampsz*train_ratio); % sample size for training
      ts_sampsz=tot_sampsz-tr_sampsz; % sample size for testing
      tr_idx=randsample(tot_sampsz,tr_sampsz);
      ts_idx=randsample(tot_sampsz,ts_sampsz);
      
      XX_tr=XX(tr_idx,:);
      YY_tr=YY(tr_idx);
      
      XX_ts=XX(ts_idx,:);
      YY_ts=YY(ts_idx);
      
%       XX_pos=XX;
%       XX_pos(XX<=0)=nan;
%       YYlog=log10(YY);
%       XXlog=log10(XX_pos);

      XX_demean_tr=zscore(XX_tr);
      XX_demean_ts=zscore(XX_ts);
%       XXlog_demean=zscore(XXlog);
      mXX_demean=zscore(mXX);

      [PCAcoeff{ir,c},~,PCAlatt]=pca(XX_demean_tr);
      PCAexp=PCAlatt/sum(PCAlatt);

%       [PCAcoeff_log{ir,c},~,PCAlatt_log]=pca(XXlog_demean);
%       PCAexp_log=PCAlatt_log/sum(PCAlatt_log);
%       
      [PCAcoeffm{c},~,PCAlattm]=pca(mXX_demean);
      PCAexpm=PCAlattm/sum(PCAlattm);
      
      % nfeat=3;
      nfeat=size(XX,2);


      XX_new=[ones(size(XX_tr,1),1) XX_tr*PCAcoeff{ir,c}(:,1:nfeat)];
      XX_new_ts=[ones(size(XX_ts,1),1) XX_ts*PCAcoeff{ir,c}(:,1:nfeat)];
%       XXlog_new=[ones(size(XX,1),1) XXlog*PCAcoeff_log{ir,c}(:,1:nfeat)];
      mXX_new=[ones(size(mXX,1),1) mXX*PCAcoeffm{c}];

%       [blog{ir,c},blogint,~,~,~]=regress(YYlog,XXlog_new);
%       Yhatlog=10.^(XXlog_new*blog{ir,c});

      [b{ir,c},bint,~,~,~]=regress(YY_tr,XX_new);
      Yhat_tr=XX_new*b{ir,c};
      vidx=~isnan(Yhat_tr+YY_tr);
      tr_rsq=1-sum((Yhat_tr(vidx)-YY_tr(vidx)).^2)/...
         sum((YY_tr(vidx)-mean(YY_tr(vidx))).^2);
      
      Yhat_ts=XX_new_ts*b{ir,c};
      vidx=~isnan(Yhat_ts+YY_ts);
      ts_rsq=1-sum((Yhat_ts(vidx)-YY_ts(vidx)).^2)/...
         sum((YY_ts(vidx)-mean(YY_ts(vidx))).^2);
      
      [bm{ir,c},bintm,~,~,~]=regress(mYY,mXX_new);
      Yhatm=mXX_new*bm{ir,c};
      
      nexttile
      scatter(Yhatm,mYY,psz,'filled','MarkerFaceColor',dp_color,...
         'MarkerFaceAlpha',1)
      vidx=~isnan(Yhatm+mYY);
      PiC=findpred(Yhatm(vidx),mYY(vidx));
      pred_str=['R^2','=', sprintf('%0.3f',PiC)];
      xlim([0 1])
      ylim([0 1])
      
      axisdim=axis;
      textx=(axisdim(2)-axisdim(1))*0.7+axisdim(1);
      texty=(axisdim(4)-axisdim(3))*0.6+axisdim(3);
      text(textx,texty,pred_str,'FontSize',16)
      
      
      refline(1,0)
      set(gca,'fontsize',20)
      xlabel([camp_proper{c} ' \epsilon'])
      
%       figure(1)
%       nexttile
%       scatter(Yhatlog,YY,psz,'filled','MarkerFaceColor',dp_color,...
%          'MarkerFaceAlpha',falp(ir))
%       vidx=~isnan(Yhatlog+YY);
%       PiClog=findpred(Yhatlog(vidx),YY(vidx));
%       pred_str=['R^2','=', sprintf('%0.3f',...
%                            PiClog)];
%       text(.7,.6,pred_str,'FontSize',16)
%       xlim([0 1])
%       ylim([0 1])
%       refline(1,0)
%       xlabel(tl1,'predicted \epsilon','fontsize',24)
%       ylabel(tl1,'observed \epsilon','fontsize',24)
%       set(gca,'fontsize',20)
%       if ir==length(arat) xlabel([camp_proper{c} ' \epsilon']), end
%       if c==1 ylabel(['ratio=' num2str(arat(ir))]), end

%       figure(2)
%       nexttile
%       scatter(Yhat_tr,YY_tr,psz,'filled','MarkerFaceColor',dp_color,...
%          'MarkerFaceAlpha',falp(ir))
%       vidx=~isnan(Yhat_tr+YY_tr);
% %       PiC=findpred(Yhat_tr(vidx),YY_tr(vidx));
%       pred_str=['R^2','=', sprintf('%0.3f',...
%                            tr_rsq)];
%       text(.7,.6,pred_str,'FontSize',16)
%       xlim([0 1])
%       ylim([0 1])
%       refline(1,0)
%       xlabel(tl2,'predicted \epsilon','fontsize',24)
%       ylabel(tl2,'observed \epsilon','fontsize',24)
%       set(gca,'fontsize',20)
%       if ir==length(arat) xlabel([camp_proper{c} ' \epsilon']), end
%       if c==1 ylabel(['ratio=' num2str(arat(ir))]), end

   end
end

%%
% saveas(figure(1),'plots/PCR randsample power.png')
% saveas(figure(2),'plots/PCR randsample linear.png')
