clearvars -except clouds
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
if ~exist('clouds','var') load clouds.mat, end
campaigns={'vocalspdi','postpdi','oraclespdi','gomaccspdi'};
camp_proper={'VOCALS','POST','ORACLES','GoMACCS'};
typevar={'s_ntot_aer', 's_ntot_pdi', 's_actfrac', 's_rh', 'AF','s_lwc_pdi'};
typename={'N_a (#/cc)','N_d (#/cc)', 'LAF','RH',...
   'Adiabatic Fraction','LWC'};
typetitle={'(a)','(b)','(c)','RH','AF','LWC'};
dp_color='#FFBE89';%'#A8CCCC';
axis_color='#638FAA';
nc=length(campaigns);
colororder=colororder;
color_order={};
for i=1:size(colororder,1)
    color_order{i}=colororder(i,:);
end
%%

all_Y=cell(nc,1);
all_X=cell(nc,length(typevar));
moist_t=cell(nc,1);
all_wz = cell(nc,1);

for c=1:nc
   %%
   camp=campaigns{c};
   fb=load([camp,'_flight_basics.mat']);
   fbvar=[camp,'_flight_basics'];
   
   fldnm=fieldnames(clouds.(camp)(1));
   
   % commontitles='INTRA-cloud correlation \epsilon vs ';
   set(groot, 'DefaultFigurePosition', [1399 1 913 984])
   
   for iday=1:length(clouds.(camp))
      
      depvar=clouds.(camp)(iday).s_disp_pdi;
      vidx=~isnan(depvar);
      all_Y{c}=[all_Y{c}; depvar(vidx)];
      
      for ivar=1:length(typevar)
         indvar=clouds.(camp)(iday).(typevar{ivar});
         indvar(isinf(indvar))=nan;
         all_X{c, ivar}=[all_X{c, ivar}; indvar(vidx)];
      end
      
      wz=clouds.(camp)(iday).s_wz;
      wz(wz<-20|wz>20)=nan;
      all_wz{c}=[all_wz{c}; wz(vidx)];

   end
end

%% shaded errorbar and compare w

vars=[1 3 5 6];
nvar=length(vars);
all_Xlimmin=[0 0 0 90 0 0
   0 0 0 90 0 0
   0 0 0 90 0 0
   0 0 0 90 0 0];

all_Xlimmax=[800 600 1 110 1 1
% 1200 300 1 110 1 1
1500 600 1 110 1 1
12000 600 1 110 1 1
3000 2500 1 110 1 1];

w_str = {'strong updraft','quiescent','strong downdraft'};
cord=[2 1 5];
nw = length(w_str);

close all
tl=tiledlayout(nc,nvar,'TileSpacing','compact');
for c=1:nc
   for icol=1:nvar
      %%
      clear w_idx
      
      nexttile
%       subplot(nc,nvar,nvar*(c-1)+icol)
      ivar=vars(icol);
      
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
%       w_idx(2,:)=all_wz{c}>w_60_val & all_wz{c}<=w_80_val; % weak updraft index
      w_idx(2,:)=all_wz{c}>w_40_val & all_wz{c}<=w_60_val; % quiescent index
%       w_idx(2,:)=~isnan(all_wz{c}); % quiescent index
%       w_idx(4,:)=all_wz{c}>w_20_val & all_wz{c}<=w_40_val; % weak downdraft index
      w_idx(3,:)=all_wz{c}<=w_20_val; % strong downdraft index
      
      
%       xlim([0 all_Xlimmax(c,ivar)])
      
      xx=all_X{c, ivar};
      yy=all_Y{c};
      
      xbins=linspace(all_Xlimmin(c,ivar),all_Xlimmax(c,ivar),21)';
      xbinmean_def=(xbins(2:end)+xbins(1:end-1))/2;
      
      
      for iw=1:nw
         clear xbinmeanw mu_reldisp_w sigma_dv_w
         j=1;
         for ibin=1:length(xbins)-1
            idx_allw=xx>xbins(ibin) ...
               & xx<=xbins(ibin+1);
            
            idx_w=intersect(find(idx_allw),find(w_idx(iw,:)));
            
            if length(idx_w)>10
               mu_reldisp_w(j,1)=nanmean(yy(idx_w));
               sigma_dv_w(j,1)=nanstd(yy(idx_w));
               xbinmeanw(j,1)=xbinmean_def(ibin);
               j=j+1;
            end
         end
         
         hold on
         fill([xbinmeanw;flipud(xbinmeanw)],...
         [mu_reldisp_w-sigma_dv_w;flipud(mu_reldisp_w+sigma_dv_w)],...
         color_order{cord(iw)},'linestyle','none','FaceAlpha',0.2,...
         'HandleVisibility','off')
         plot(xbinmeanw,mu_reldisp_w,'Color',color_order{cord(iw)},...
            'LineWidth',1,'DisplayName',w_str{iw})
         hold off
         ylim([0 .7])
      end
%       legend('strong updraft','quiescent','strong downdraft')
      
      if c==nc
         xlabel(typename{ivar})
      end
%       
      if icol>1
         set(gca,'YTick',[])
      else
         ylabel([camp_proper{c}, ' \epsilon'])
      end
      
      set(gca, 'fontsize', 16)
   end
end

lg=legend('show','Location','northeastoutside');
lg.Layout.Tile = 'north';

% saveas(gcf,'plots/INTRA-cloud_allcamps_w.png')

% %% plot LAF vs var
% 
% close all
% 
% for icol=4
%    figure('Position',[1392 287 1244 593])
%    tl=tiledlayout('flow');
%    
%    
%    for ic=1:5
%       nexttile
%       
%       indvar=all_X{ic,3};
%       depvar=all_X{ic,icol};
%       
%       if strcmp(typevar{icol},'AF')
%          depvar(depvar>1)=nan;
%          depvar(depvar<0)=0;
%       end
%       
%       scatter(indvar,depvar,20,'filled',...
%          'MarkerFaceAlpha',0.3,'MarkerFaceColor',dp_color)
%       if strcmp(typevar{icol},'AF'); ylim([0 1]); end
%       if strcmp(typetitle{icol},'ER'); ylim([-.1 .3]); end
%       xlim([0 1])
%       
%       hold on
%       
%       %         [comparison_mtx,hist_centers]=hist3([indvar,...
%       %             depvar],'Nbins',[30,30]);
%       %         contour(hist_centers{1},hist_centers{2},comparison_mtx','LineWidth',2)
%       
%       xbins=linspace(min(xlim),max(xlim),11);
%       mu_dv=zeros(length(xbins)-1,1);
%       sigma_dv=zeros(length(xbins)-1,1);
%       
%       for ibin=1:length(xbins)-1
%          idx_allw=indvar>xbins(ibin) ...
%             & indvar<=xbins(ibin+1);
%          mu_dv(ibin)=nanmean(depvar(idx_allw));
%          sigma_dv(ibin)=nanstd(depvar(idx_allw));
%          
%          if sum(idx_allw)>100
%             mu_dv(ibin)=nanmean(depvar(idx_allw));
%             sigma_dv(ibin)=nanstd(depvar(idx_allw));
%          else
%             mu_dv(ibin)=nan;
%             sigma_dv(ibin)=nan;
%          end
%       end
%       
%       
%       xbinmean=(xbins(2:end)+xbins(1:end-1))/2;
%       
%       errorbar(xbinmean,mu_dv,sigma_dv,'LineWidth',2)
%       
%       title(gca,camp_proper{ic},'FontWeight','normal')
%       set(gca,'fontsize',20)
%    end
%    
%    xlabel(tl,'Local Activation Fraction','fontsize',24,'fontweight','bold')
%    ylabel(tl,typename{icol},'fontsize',24,'fontweight','bold')
% end
% 
% saveas(gcf,'plots/AF_vs_LAF.png')
% 
% %%
% clear lt25d
% for iday=1:21
%    lt25=sum(clouds.gomaccspdi(iday).s_conc_pdi(:,1:70),2)./sum(clouds.gomaccspdi(iday).s_conc_pdi,2);
%    lt25d(iday)=nanmean(lt25);
% end
