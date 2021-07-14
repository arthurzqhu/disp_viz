close all
cd '~/MEGASync/grad/research/aerosol_reldisp/datasets/'
clearvars -except clouds
if ~exist('clouds','var') load clouds.mat, end

campaigns={'vocalspdi','masepdi','postpdi','oraclespdi','gomaccspdi'};
ctitle={'VOCALS','MASE','POST','ORACLES','GoMACCS'};
dpc='#FFBE89'; % datapoint color
% axis_color='#336C8A';
nc=length(campaigns);
colororder=colororder;
color_order={};
for i=1:size(colororder,1)
   color_order{i}=colororder(i,:);
end

%%
close all
max_ratio=[1 2 3 5 10 20 50];
ntrial=1;
trial_allSc_rsq=zeros(ntrial,length(max_ratio));
doplot=true;
iratio=1;
figure('Position',[1668 63 913 914])
tl1=tiledlayout('flow');

figure('Position',[1392 287 1244 593])
tl2=tiledlayout('flow');

for c=1:nc
   camp=campaigns{c};
   fb=load([camp,'_flight_basics.mat']);
   fbvar=[camp,'_flight_basics'];
   
   nlyrs=20;
   samp_frac=0.5;
   
   ndays=length(clouds.(camp));
   days_analyzed=1:ndays;
   %     remove the days that have incomplete flights in vocals
   if strcmp(camp,'vocalspdi')
      %                 days_analyzed(ismember(days_analyzed,[8,10,11,13]))=[];
   elseif strcmp(camp,'oraclespdi')
      days_analyzed(ismember(days_analyzed, [2,3,5,14,15,16,17,18]))=[];
   end
   
   [s_ntot_aer,s_ap,thet,s_ntot_pcasp,s_ntot_pdi,s_actfrac,s_lwc_pdi,normAC,...
      AF,reldisp]=deal([]);
   
   for iday=days_analyzed
      
      % get the unfiltered time first for later use
      s_t=clouds.(camp)(iday).s_t;
      
      cloudlegs_i=fb.(fbvar)(iday).ti;
      cloudlegs_f=fb.(fbvar)(iday).tf;
      
      if ~isempty(cloudlegs_i)
         for ileg=1:length(cloudlegs_i)
            % get the initial and final time for each cloud
            ti=cloudlegs_i(ileg);
            tf=cloudlegs_f(ileg);
            
            % filtering criteria
            filt_crit=find(s_t>ti & s_t<tf & ...
               clouds.(camp)(iday).s_ntot_pdi>25 & ...
               clouds.(camp)(iday).normAC>0 & ...
               clouds.(camp)(iday).normAC<1.1);
            
            s_t_leg=s_t(filt_crit);
            
            reldisp=[reldisp;clouds.(camp)(iday).s_disp_pdi(filt_crit)];
            normAC=[normAC;clouds.(camp)(iday).normAC(filt_crit)];
            s_ntot_pdi=[s_ntot_pdi;clouds.(camp)(iday).s_ntot_pdi(filt_crit)];
            s_ap=[s_ap;clouds.(camp)(iday).s_ap(filt_crit)];
            thet=[thet;clouds.(camp)(iday).s_thet(filt_crit)];
            s_ntot_aer=[s_ntot_aer;clouds.(camp)(iday).s_ntot_aer(filt_crit)];
            s_ntot_pcasp=[s_ntot_pcasp;clouds.(camp)(iday).s_ntot_pcasp(filt_crit)];
            s_actfrac=[s_actfrac;clouds.(camp)(iday).s_actfrac(filt_crit)];
            s_lwc_pdi=[s_lwc_pdi;clouds.(camp)(iday).s_lwc_pdi(filt_crit)];
            AF=[AF;clouds.(camp)(iday).AF(filt_crit)];
            
         end
      end
   end
   
   AF(AF==1)=nan;
   
   switch c
      case 3
         s_ntot_aer(s_ntot_aer>1500)=nan;
      case 4
         s_ntot_aer(s_ntot_aer>3000)=nan;
      case 5
         s_actfrac(s_actfrac==1)=nan;
         s_ntot_aer(s_ntot_aer>3500)=nan;
         s_ntot_pdi(s_ntot_pdi>3000)=nan;
   end
   
   figure(1)
   nexttile
   scatter(s_ntot_aer,reldisp,20,'filled','MarkerFaceColor',dpc,...
      'MarkerFaceAlpha',0.1)
   xbins = linspace(min(xlim),max(xlim),11);
   mu_reldisp = zeros(length(xbins)-1,1);
   for ibin = 1:length(xbins)-1
      idx_at_ibin = s_ntot_aer>xbins(ibin) ...
          & s_ntot_aer<=xbins(ibin+1);

      if sum(idx_at_ibin)>25
          mu_reldisp(ibin) = nanmean(reldisp(idx_at_ibin));
          sigma_dv(ibin) = nanstd(reldisp(idx_at_ibin));
      else
          mu_reldisp(ibin) = nan;
          sigma_dv(ibin) = nan;
      end
   end

   xbinmean = (xbins(2:end)+xbins(1:end-1))/2;
   hold on
   errorbar(xbinmean,mu_reldisp,sigma_dv,'LineWidth',2)
   ylim([0 1])
   ylabel([ctitle{c}, ' \epsilon'])
   set(gca, 'fontsize', 18)
   if c==1, title('(a)'), end
   if c==5, xlabel('Aerosol Concentration (#/cc)'), end
   
   nexttile
   scatter(s_ntot_pdi,reldisp,20,'filled','MarkerFaceColor',dpc,...
      'MarkerFaceAlpha',0.1)
   xbins = linspace(min(xlim),max(xlim),11);
   mu_reldisp = zeros(length(xbins)-1,1);
   for ibin = 1:length(xbins)-1
      idx_at_ibin = s_ntot_pdi>xbins(ibin) ...
          & s_ntot_pdi<=xbins(ibin+1);

      if sum(idx_at_ibin)>25
          mu_reldisp(ibin) = nanmean(reldisp(idx_at_ibin));
          sigma_dv(ibin) = nanstd(reldisp(idx_at_ibin));
      else
          mu_reldisp(ibin) = nan;
          sigma_dv(ibin) = nan;
      end
   end

   xbinmean = (xbins(2:end)+xbins(1:end-1))/2;
   hold on
   errorbar(xbinmean,mu_reldisp,sigma_dv,'LineWidth',2)
   ylim([0 1])
   set(gca,'YTick',[], 'fontsize', 18)
   if c==1, title('(b)'), end
   if c==5, xlabel('Droplet Concentration (#/cc)'), end
   
   nexttile
   scatter(s_actfrac,reldisp,20,'filled','MarkerFaceColor',dpc,...
      'MarkerFaceAlpha',0.1)
   xbins = linspace(min(xlim),max(xlim),11);
   mu_reldisp = zeros(length(xbins)-1,1);
   for ibin = 1:length(xbins)-1
      idx_at_ibin = s_actfrac>xbins(ibin) ...
          & s_actfrac<=xbins(ibin+1);

      if sum(idx_at_ibin)>25
          mu_reldisp(ibin) = nanmean(reldisp(idx_at_ibin));
          sigma_dv(ibin) = nanstd(reldisp(idx_at_ibin));
      else
          mu_reldisp(ibin) = nan;
          sigma_dv(ibin) = nan;
      end
   end

   xbinmean = (xbins(2:end)+xbins(1:end-1))/2;
   hold on
   errorbar(xbinmean,mu_reldisp,sigma_dv,'LineWidth',2)
   ylim([0 1])
   set(gca,'YTick',[], 'fontsize', 18)
   if c==1, title('(c)'), end
   if c==5, xlabel('Local Activation Fraction'), end
   
   figure(2)
   nexttile
   scatter(s_actfrac,AF,20,'filled',...
            'MarkerFaceAlpha',0.3,'MarkerFaceColor',dpc)
   xlim([0 1])
   ylim([0 1])
   
   xbins = linspace(min(xlim),max(xlim),11);
   mu_AF = zeros(length(xbins)-1,1);
   for ibin = 1:length(xbins)-1
      idx_at_ibin = s_actfrac>xbins(ibin) ...
          & s_actfrac<=xbins(ibin+1);

      if sum(idx_at_ibin)>25
          mu_AF(ibin) = nanmean(AF(idx_at_ibin));
          sigma_dv(ibin) = nanstd(AF(idx_at_ibin));
      else
          mu_AF(ibin) = nan;
          sigma_dv(ibin) = nan;
      end
   end
   xbinmean = (xbins(2:end)+xbins(1:end-1))/2;
   hold on
   errorbar(xbinmean,mu_AF,sigma_dv,'LineWidth',2)
   
   title(ctitle(c),'FontWeight','normal')
   set(gca,'fontsize',20)
   xlabel(tl2,'Local Activation Fraction','fontsize',24,'fontweight','bold')
   ylabel(tl2,'Activation Fraction','fontsize',24,'fontweight','bold')
end

exportgraphics(figure(1),'../paper in progress/figs/fig2.jpg','Resolution',300)
exportgraphics(figure(2),'../paper in progress/figs/fig3.jpg','Resolution',300)
