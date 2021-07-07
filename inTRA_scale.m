close all
cd '~/MEGASync/grad/research/aerosol_reldisp/datasets/'
clearvars -except clouds
if ~exist('clouds','var') load clouds.mat, end

campaigns={'vocals','mase','post','oracles','gomaccs'};
ctitle={'VOCALS','MASE','POST','ORACLES','GoMACCS'};
drp_instr={'pdi','pdi','pdi','pdi','pdi'};
dpc='#FB9966'; % datapoint color
% axis_color='#336C8A';
nc=length(campaigns);
colororder=colororder;
color_order={};
for i=1:size(colororder,1)
    color_order{i}=colororder(i,:);
end
%% inflight averaging and cloud averaging (without subsampling)
close all
clear X avgcoeff avgtest_rsq avgtrain_rsq

arat=round(2.^(.5:.5:10));%[1:9 10:10:90 100:100:1000]; % averaging ratio
falp=(erf(log10(arat)-1)+1)/2;
ntrial=20;
[Sctest_rsq,Sctrain_rsq]=deal(zeros(length(arat),ntrial));
Sccoeff={};
max_ratio=[1 2 3 5 10 20 50];
iratio=2;
nlyrs=20;

tic
for itrial=1:ntrial
   itrial
   for ir=1:length(arat)

%       figure('Position',[1032 199 1549 778]);
%       tl=tiledlayout('flow');
      [allSc_YYss,allSc_XXss]=deal([]);

      for c=1:3
         [s_lwc_pdi,normAC,s_ntot_aer,s_ntot_pdi,s_ntot_pcasp,ql_adb,rhoa,...
            s_actfrac_imp,reldisp,reldisp_imp]=deal([]);

         binedges=load(['bin_edges_' campaigns{c} '.csv']);
         binmean=(binedges(:,2)+binedges(:,3))/2;
         dlogD=log10(binedges(2,3)/binedges(2,2));

         camp=[campaigns{c} drp_instr{c}];

         ndays=length(clouds.(camp));
         days_analyzed=1:ndays;
         %     remove the days that have incomplete flights in vocals
         if strcmp(camp,'vocalspdi')
            %                 days_analyzed(ismember(days_analyzed,[8,10,11,13]))=[];
         elseif strcmp(camp,'oraclespdi')
            days_analyzed(ismember(days_analyzed, [2,3,5,14,15,16,17,18]))=[];
         end

         [m_lwc_pdi,m_ntot_pdi,m_ntot_aer,m_actfrac,m_reldisp,m_AF,m_rhoa,...
            m_ql_adb]=...
            deal(zeros(length(days_analyzed),1));

         for iday=days_analyzed
            reldisp_sd=[];
            % indvars
            s_lwc_pdi=[s_lwc_pdi;mosaicify(clouds.(camp)(iday).s_lwc_pdi,arat(ir))];
            normAC=[normAC;mosaicify(clouds.(camp)(iday).normAC,arat(ir))];
            s_ntot_pdi=[s_ntot_pdi;mosaicify(clouds.(camp)(iday).s_ntot_pdi,arat(ir))];
            s_ntot_aer=[s_ntot_aer;mosaicify(clouds.(camp)(iday).s_ntot_aer,arat(ir))];
            s_ntot_pcasp=[s_ntot_pcasp;mosaicify(clouds.(camp)(iday).s_ntot_pcasp,arat(ir))];
            s_actfrac_imp=[s_actfrac_imp;mosaicify(clouds.(camp)(iday).s_actfrac,arat(ir))];
            ql_adb=[ql_adb;mosaicify(clouds.(camp)(iday).ql_adb_prof,arat(ir))];
            rhoa=[rhoa;mosaicify(clouds.(camp)(iday).s_rhoa,arat(ir))];

            vidx=clouds.(camp)(iday).s_ntot_pdi>10;
%             m_lwc_pdi(iday)=nanmean(clouds.(camp)(iday).s_lwc_pdi(vidx));
%             m_ntot_pdi(iday)=nanmean(clouds.(camp)(iday).s_ntot_pdi(vidx));
%             m_ntot_aer(iday)=nanmean(clouds.(camp)(iday).s_ntot_aer(vidx));
%             m_actfrac(iday)=m_ntot_pdi(iday)/m_ntot_aer(iday);
%             m_rhoa(iday)=nanmean(clouds.(camp)(iday).s_rhoa(vidx));
%             m_ql_adb(iday)=nanmean(clouds.(camp)(iday).ql_adb_prof(vidx));
% 
%             m_conc_pdi=nanmean(clouds.(camp)(iday).s_conc_pdi(vidx,:));
%             m_reldisp(iday)=std(binmean,m_conc_pdi,'omitnan')./...
%                   wmean(binmean,m_conc_pdi);
%             m_AF(iday)=m_lwc_pdi(iday)./m_rhoa(iday)./m_ql_adb(iday);


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
         
         if strcmp(camp,'postpdi')
            % unbias the vertical sampling
            normAC_incloud=normAC>=0&normAC<=1;
            [N_raw,edges,lyrs]=histcounts(normAC(normAC_incloud),nlyrs);
            minNval=min(N_raw(N_raw>0));
            N_raw(N_raw>max_ratio(iratio)*minNval)=...
                           max_ratio(iratio)*minNval;
            rdsmp_idx=[];
            for ilyr=1:nlyrs
               rdsmp_idx_lyr=randsample(find(lyrs==ilyr),N_raw(ilyr));
               rdsmp_idx=[rdsmp_idx;rdsmp_idx_lyr];
            end
            
            s_lwc_pdi=s_lwc_pdi(rdsmp_idx);
            normAC=normAC(rdsmp_idx);
            rhoa=rhoa(rdsmp_idx);
            ql_adb=ql_adb(rdsmp_idx);
            reldisp=reldisp(rdsmp_idx);
            s_ntot_pdi=s_ntot_pdi(rdsmp_idx);
            s_ntot_aer=s_ntot_aer(rdsmp_idx);
            
         end

         s_actfrac=s_ntot_pdi./s_ntot_aer;
         ql_obs=s_lwc_pdi./rhoa;
         AF=ql_obs./ql_adb;

         x1=s_lwc_pdi;
         x2=s_ntot_aer;
         x3=s_actfrac;
         x4=normAC;
         x5=AF;

         if strcmp(camp,'gomaccspdi')
            X_raw=[ones(size(x1)) x3 x5];
         elseif strcmp(camp,'oraclespdi')
            X_raw=[ones(size(x1)) x1 x2];
         else
            X_raw=[ones(size(x1)) x3 x4];
         end
         vidx=~isnan(sum(X_raw,2));
         X=X_raw(vidx,:);

         yy=reldisp(vidx);

%          [avgcoeff{ir}(c,:), avgtrain_rsq{ir}(c), avgtest_rsq{ir}(c)]=...
%             regress_tt(yy, X, .7);

%          yhat=X*avgcoeff{ir}(c,:)';

%          nexttile
%          scatter(yhat,yy,'filled',...
%             'MarkerFaceAlpha',falp(ir),'MarkerFaceColor',dpc)
%          title(ctitle{c},'FontSize',16)
%          refline(1,0)
%          ylim([0 1])
%          xlim([0 1])
         
%          vidx=~isnan(yhat+yy);
   %       rsq=1-sum((yhat(vidx)-yy(vidx)).^2)/sum((yy(vidx)-mean(yy(vidx))).^2);
   %       rsq=avgtest_rsq{ir}(c);
         
%          pred_str=['\Pi_C','=', sprintf('%0.3f', avgtest_rsq{ir}(c))];
%          axisdim=axis;
%          textx=(axisdim(2)-axisdim(1))*0.7+axisdim(1);
%          texty=(axisdim(4)-axisdim(3))*0.6+axisdim(3);
%          text(textx,texty,pred_str,'FontSize',16)

         if ismember(c,1:3)
            XXSc_cell{c}=X;
            YYSc_cell{c}=yy;
         end

      end

      Sc_ss=min(cellfun(@length,YYSc_cell));

      for c=1:3
         ss_idx=randsample(length(YYSc_cell{c}),Sc_ss);
         allSc_XXss=[allSc_XXss; XXSc_cell{c}(ss_idx,:)];
         allSc_YYss=[allSc_YYss; YYSc_cell{c}(ss_idx)];
   %       allSc_XXss=[allSc_XXss; XXSc_cell{c}];
   %       allSc_YYss=[allSc_YYss; YYSc_cell{c}];
      end



      [Sccoeff{ir,itrial},Sctrain_rsq(ir,itrial),Sctest_rsq(ir,itrial)]=...
            regress_tt(allSc_YYss, allSc_XXss, .7);

%       yhatSc=allSc_XXss*Sccoeff{ir,itrial};
%       nexttile
%       scatter(yhatSc,allSc_YYss,'filled',...
%          'MarkerFaceAlpha',falp(ir),'MarkerFaceColor',dpc)
%       title('All Sc','FontSize',16)
%       refline(1,0)
%       ylim([0 1])
%       xlim([0 1])
      
%       vidx=~isnan(yhatSc+allSc_YYss);
   %    rsq=1-sum((yhatSc(vidx)-allSc_YYss(vidx)).^2)/sum((allSc_YYss(vidx)-mean(allSc_YYss(vidx))).^2);
%       rsq=Sctest_rsq(ir);
%       pred_str=['\Pi_C','=', sprintf('%0.3f',rsq)];
%       axisdim=axis;
%       textx=(axisdim(2)-axisdim(1))*0.7+axisdim(1);
%       texty=(axisdim(4)-axisdim(3))*0.6+axisdim(3);
%       text(textx,texty,pred_str,'FontSize',16)
%       
%       title(tl,['ratio=' num2str(arat(ir))],'fontsize',20,'fontweight','bold')
%       saveas(gcf,['plots/scale dep test rsq ratio=' num2str(arat(ir)) ' proper.png'])

   end
end
toc
%%
close all
load scale_analysis
figure('Position',[1282 951 1142 386])
tl=tiledlayout('flow');

rsq_mean=mean(Sctest_rsq,2);
rsq_std=std(Sctest_rsq')';

% errorbar(arat,rsq_mean,rsq_std,'LineWidth',2)
nexttile
hold on
scalex=arat'*55/1e3;

fill([scalex;flipud(scalex)],...
   [rsq_mean-rsq_std;flipud(rsq_mean+rsq_std)],...
   color_order{1},'FaceAlpha',0.2,'LineStyle','none')
      
plot(scalex,rsq_mean,'LineWidth',2,'Color',color_order{1})
grid
ylim([-.2 .6])
set(gca,'XScale','log')
ylabel('\Pi_{Sc}')
xticks([1e-2 1e-1 1e0 1e1 1e2])
set(gca,'FontSize',24)
% saveas(gcf,'plots/pi vs scale.png')

coeff_itr=zeros(ntrial,size(Sccoeff{1,1},1));
coeff_mu=zeros(length(arat),size(Sccoeff{1,1},1));
coeff_std=zeros(length(arat),size(Sccoeff{1,1},1));

for ir=1:length(arat)

   for itrial=1:ntrial
   coeff_itr(itrial,:)=Sccoeff{ir,itrial};
   end

   coeff_mu(ir,:)=mean(coeff_itr);
   coeff_std(ir,:)=std(coeff_itr);
end
nexttile
hold on
plot(scalex,coeff_mu,'LineWidth',2)
for i=1:3
   fill([scalex;flipud(scalex)],...
      [coeff_mu(:,i)-coeff_std(:,i);flipud(coeff_mu(:,i)+coeff_std(:,i))],...
      color_order{i},'FaceAlpha',0.2,'LineStyle','none')
end



grid
set(gca,'XScale','log')
legend('Int','m_{LAF}','m_{NH}','Location','best')
xlabel(tl,'Scale [km]','fontsize',24)
ylabel('Coefficients')
xticks([1.e-2 1e-1 1e0 1e1 1e2])
set(gca,'FontSize',24)
saveas(gcf,'plots/pi & coeff vs scale.png')