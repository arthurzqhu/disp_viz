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
%%

close all

max_ratio=[2 3 5 10 20 50];
ntrial=1;
trial_allSc_rsq=zeros(ntrial,length(max_ratio));
doplot=true;
disp_thres=.4;

allSc_Xss=[];
allSc_Yss=[];

if doplot; set(0, 'DefaultFigurePosition',[1032 199 1549 778]); end

tic
for iratio=2%:length(max_ratio)
   
   for itrial=1:ntrial
      
      if mod(itrial,10)==0; disp(itrial); end
      
      for c=1:nc
         
         
         clear b X y y_hat pred ac_arr bincenters
         
         clear slp_raw intcpt_raw rsq_raw
         camp=[campaigns{c} drp_instr{c}];
         fb=load([camp,'_flight_basics.mat']);
         fbvar=[camp,'_flight_basics'];
         isOorG=any(strcmp({'oraclespdi';'gomaccspdi'},camp));
         
         % maxN=arrayfun(@(x) clouds.(camppdi)(x).maxN, 1:length(clouds.(camppdi)));
         min_dtpt=100;
         
         % exclude the date where not enough data points are present
         dtpt_pdi=arrayfun(@(x) length(clouds.(camp)(x).s_t),1:length(clouds.(camp)));
         vdate_pdi=find(dtpt_pdi>min_dtpt)';
         
         
         % close all
         epsvar=['s_disp_' drp_instr{c}];
         Nvar=['s_ntot_' drp_instr{c}];
         
         ndays=length(vdate_pdi);
         
         nbins=20;
         samp_frac=0.5;
         upper_limit=2e4;
         %     layers=3;
         layer_seg=[0 0.2 0.8 1.2];
         allb=[];
         all_pred_dailyreg=[];
         allX=[];
         ally=[];
         allnormAC=[];
         all_y_hat_dailyreg=[];
         reglr_allX=[];
         
         ndays=length(clouds.(camp));
         days_analyzed=1:ndays;
         %     remove the days that have incomplete flights in vocals
         if strcmp(camp,'vocalspdi')
            %                 days_analyzed(ismember(days_analyzed,[8,10,11,13]))=[];
         elseif strcmp(camp,'oraclespdi')
            days_analyzed(ismember(days_analyzed, [2,3,5,14,15,16,17,18]))=[];
         end
         
         
         % close all
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
                     clouds.(camp)(iday).(Nvar)>25);
                  
                  s_t_leg=s_t(filt_crit);
                  
                  reldisp=clouds.(camp)(iday).(epsvar)(filt_crit);
                  s_ntot_pdi=clouds.(camp)(iday).(Nvar)(filt_crit);
                  s_ap=clouds.(camp)(iday).s_ap(filt_crit);
                  
                  normAC=clouds.(camp)(iday).normAC(filt_crit);
                  thet=clouds.(camp)(iday).s_thet(filt_crit);
                  
                  s_ntot_aer=clouds.(camp)(iday).s_ntot_aer(filt_crit);
                  s_ntot_pcasp=clouds.(camp)(iday).s_ntot_pcasp(filt_crit);
                  s_actfrac=clouds.(camp)(iday).s_actfrac(filt_crit);
                  s_lwc_pdi=clouds.(camp)(iday).s_lwc_pdi(filt_crit);
                  ent_ratio_qt=clouds.(camp)(iday).ent_ratio_qt(filt_crit);
                  AF=clouds.(camp)(iday).AF(filt_crit);
                  
                  indvar_raw=normAC;
                  depvar_raw=reldisp;
                  
                  % make sure indvar and depvar are a pair of non-NaNs
                  vidx=~isnan(indvar_raw+depvar_raw);
                  indvar=indvar_raw(vidx);
                  depvar=depvar_raw(vidx);
                  % color=normAC(vidx);
                  
                  if sum(~isnan(depvar))>100 && ~strcmp(camp,'gomaccspdi')
                     % set the maximum subsample size of each bin
                     max_per_bin=floor(length(indvar)/nbins);
                     [N_raw,edges,bin]=histcounts(indvar,nbins);
                     
                     % in case it's still too biased
                     if max(N_raw)/min(N_raw(N_raw>0))>max_ratio(iratio) && c==3
                        minNval=min(N_raw(N_raw>0));
                        N_raw(N_raw>max_ratio(iratio)*minNval)=...
                           max_ratio(iratio)*minNval;
                     end
                     
                     N=ceil(samp_frac*N_raw);
                     N(N>max_per_bin)=max_per_bin;
                     N(N<5)=0;
                     
                     rand_samp_idx=[];
                     for ibin=1:nbins
                        rand_samp_idx_bin=randsample(find(bin==ibin),N(ibin));
                        rand_samp_idx=[rand_samp_idx; rand_samp_idx_bin];
                     end
                     
                     % SubSampled data
                     indvar_ss=indvar(rand_samp_idx);
                     depvar_ss=depvar(rand_samp_idx);
                     bin_ss=bin(rand_samp_idx);
                     
                     % to prevent matrix dimension error below
                     if length(bin_ss) ~= length(N)
                        ac_arr_raw=accumarray(bin_ss,depvar_ss);
                        ac_arr_raw(ac_arr_raw==0)=[];
                        ac_arr(N~=0)=ac_arr_raw;
                        ac_arr(N==0)=0;
                     end
                     
                     bindep=ac_arr./N;
                     bincenters=(edges(1:end-1)+edges(2:end))/2;
                  elseif strcmp(camp,'gomaccspdi')
                     rand_samp_idx=true(size(s_ntot_pdi));
                  end
                  
                  x1=s_lwc_pdi(rand_samp_idx);
                  x2=s_ntot_pdi(rand_samp_idx);
                  x3=s_actfrac(rand_samp_idx);
                  x4=normAC(rand_samp_idx);
                  x5=AF(rand_samp_idx);

                  y=reldisp(rand_samp_idx);

                  if strcmp(camp,'gomaccspdi')
                     X=[ones(size(x1)) x3 x5];
                  elseif strcmp(camp,'oraclespdi')
                     X=[ones(size(x1)) x1 x2];
                  else
                     X=[ones(size(x1)) x3 x4];
                  end

                  b=regress(y,X);
                  y_hat=X*b;
                  
                  reglr_X=[ones(size(x1)) x1 x2 x3 x4 x5];
                  pred=1-nansum((y_hat-y).^2)/nansum((y-nanmean(y)).^2);
                  allX=[allX;X];
                  allnormAC=[allnormAC;normAC];
                  reglr_allX=[reglr_allX; reglr_X];
                  ally=[ally;y];
                  all_y_hat_dailyreg=[all_y_hat_dailyreg;y_hat];
                  allb=[allb;b'];
                  all_pred_dailyreg=[all_pred_dailyreg;pred];
                  warning('off')
               end
            end
         end
         
         if ~doplot && itrial==1 && c==3
            length(ally)
         end
         
         
         %             disp_ge_thres=ally<=disp_thres;
         %
         %             ally=ally(disp_ge_thres);
         %             allX=allX(disp_ge_thres);
         %             all_y_hat_dailyreg=all_y_hat_dailyreg(disp_ge_thres);
         %             reglr_allX=reglr_allX(disp_ge_thres,:);
         
         meanB{c}=regress(ally, allX);
         %             meanB{c}
         [coeff(c,:), train_rsq(c), test_rsq(c)]=regress_tt(ally, allX, 0.7);
         all_y_hat_campmean=allX*meanB{c};
         notnan_idx_cm=~isnan(all_y_hat_campmean);
         notnan_idx_dr=~isnan(all_y_hat_dailyreg);
         pred_campmean=1-nansum((all_y_hat_campmean-ally).^2)/nansum((ally(notnan_idx_cm)-nanmean(ally(notnan_idx_cm))).^2);
         pred_dailyreg=1-nansum((all_y_hat_dailyreg-ally).^2)/nansum((ally(notnan_idx_dr)-nanmean(ally(notnan_idx_dr))).^2);

         allcamp_X{c}=[reglr_allX(:,1) reglr_allX(:,4:5)];
         allcamp_Y{c}=ally;
         %
         y_hat_combined=[all_y_hat_dailyreg all_y_hat_campmean];
         pred_combined=[pred_dailyreg test_rsq(c)];
         title_string={'Specialized prediction', 'Campaign prediction'};
         y_hat_string={'dailyreg', 'campmean'};
         
         for itype=2
            if doplot
               fh(itype)=figure(itype);
               sp_loc=c;
               if c>3; sp_loc=c+1; end
               subplot(2,4,sp_loc)
               hold on
               plot([0 1], [0 1])
               %         close all
               %                     gclr=dpc.*allnormAC;
               %                     gclr(gclr>1)=1;
               sp=scatter(y_hat_combined(:,itype),ally,...
                  'filled','MarkerFaceAlpha',0.1, 'MarkerFaceColor',dpc);
               
               cont_nbins=floor((length(ally))^(1/3));
               [cmprsn_mtx,hist_cens]=hist3([y_hat_combined(:,itype),ally],...
                  'nbins',[cont_nbins,cont_nbins]);
               contour(hist_cens{1},hist_cens{2},cmprsn_mtx','LineWidth',2)
               %         cbar=colorbar;
               %         cbar.Label.String='Frequency of occurrence';
               ax_lim=[0 1];
               xlim(ax_lim)
               ylim(ax_lim)
               %     set(gca,'XTick',[0.1 0.3 0.5])
               %     set(gca,'YTick',[.2 .4 .6 .8])
               set(gca,'fontsize',24)
               
               
               %% drawing the quartile plot
               y_quartl=num2cell(prctile(ally, 0:25:100));
               [y_min, y_25, y_med, y_75, y_max]=deal(y_quartl{:});
               
               x_quartl=num2cell(prctile(y_hat_combined(:,itype), 0:25:100));
               [x_min, x_25, x_med, x_75, x_max]=deal(x_quartl{:});
               
               
               line([x_min, x_25],[ax_lim(1) ax_lim(1)],'linewidth',9,'color',dpc)
               plot(x_med,ax_lim(1),'.','MarkerSize',20,'color',dpc)
               line([x_75, x_max],[ax_lim(1) ax_lim(1)],'linewidth',9,'color',dpc)
               
               line([ax_lim(1) ax_lim(1)],[y_min, y_25],'linewidth',9,'color',dpc)
               plot(ax_lim(1),y_med,'.','MarkerSize',20,'color',dpc)
               line([ax_lim(1) ax_lim(1)],[y_75, y_max],'linewidth',9,'color',dpc)
               
               %         xlabel('Predicted relative dispersion')
               %         ylabel('Measured relative dispersion')
               title(ctitle{c})%,'Color',axis_color)
               
               
               hold off
               %% text for predictibility
               pred_str=['\Pi_',title_string{itype}(1),'=', sprintf('%0.3f',...
                  pred_combined(itype))];
               text(.7,.6,pred_str,'FontSize',24)%,'Color',axis_color)
            end
         end
      end
      %
      reglr_set=[reglr_allX ally];
      
      
      
      camp_ss_size=min(cellfun(@length,allcamp_Y));
      
      for c=1:3
         ss_idx=randsample(length(allcamp_Y{c}),camp_ss_size);
         allSc_Xss=[allSc_Xss; allcamp_X{c}(ss_idx,:)];
         allSc_Yss=[allSc_Yss; allcamp_Y{c}(ss_idx)];
      end
      
      
      allSc_meanB=regress(allSc_Yss, allSc_Xss);
      allSc_tbl=fitlm(allSc_Xss(:,2:end),allSc_Yss);
      [allSc_coeff, allSc_train_rsq, allSc_test_rsq]=regress_tt(allSc_Yss, allSc_Xss, 0.7);
      allSc_Ysshat=allSc_Xss*allSc_meanB;
      
      trial_allSc_rsq(itrial,iratio)=allSc_test_rsq;
   end
   
   
end

toc


% %%
%
% if ~doplot
%     figure('Position',[1658 404 812 533])
%     mean_rsq=mean(trial_allSc_rsq);
%     std_rsq=std(trial_allSc_rsq);
%
%     errorbar(max_ratio,mean_rsq,std_rsq,'LineWidth',2)
%     xlabel('max/min # datapoints')
%     ylabel('\Pi_{Sc}')
%     set(gca,'FontSize',18)
%     % xticks(1:length(max_ratio))
%     % xlim([0 7])
%     % xticklabels(cellfun(@num2str,num2cell(max_ratio),'UniformOutput',false))
%     saveas(gcf,'plots/pi vs bias.png')
% end

% %% Allst plot
if doplot
   for itype=2
      fh(itype)=figure(itype);
      sgt=sgtitle(title_string{itype});
      sgt.FontSize=36;
      sgt.FontWeight='bold';
      %     sgt.Color=axis_color;
      
      %         whitebg([.95 .95 .95])
      
      
      subplot(2,4,7)
      hold on
      set(gca,'Color',[.95 .95 .95])
      set(gca,'XColor','k')
      set(gca,'YColor','k')
      plot([0 1], [0 1])
      %         close all
      sp=scatter(allSc_Ysshat,allSc_Yss,'filled',...
         'MarkerFaceAlpha',0.1,'MarkerFaceColor',dpc);
      
      cont_nbins=floor((length(allSc_Yss))^(1/3));
      [cmprsn_mtx,hist_cens]=hist3([allSc_Ysshat,allSc_Yss],'nbins',...
         [cont_nbins,cont_nbins]);
      contour(hist_cens{1},hist_cens{2},cmprsn_mtx','LineWidth',2)
      ax_lim=[0 1];
      xlim(ax_lim)
      ylim(ax_lim)
      set(gca,'fontsize',24)
      
      
      % drawing the quartile plot
      y_quartl=num2cell(prctile(allSc_Yss, 0:25:100));
      [y_min, y_25, y_med, y_75, y_max]=deal(y_quartl{:});
      
      x_quartl=num2cell(prctile(allSc_Ysshat, 0:25:100));
      [x_min, x_25, x_med, x_75, x_max]=deal(x_quartl{:});
      
      
      line([x_min, x_25],[ax_lim(1) ax_lim(1)],'linewidth',9,'color',dpc)
      plot(x_med,ax_lim(1),'.','MarkerSize',20,'color',dpc)
      line([x_75, x_max],[ax_lim(1) ax_lim(1)],'linewidth',9,'color',dpc)
      
      line([ax_lim(1) ax_lim(1)],[y_min, y_25],'linewidth',9,'color',dpc)
      plot(ax_lim(1),y_med,'.','MarkerSize',20,'color',dpc)
      line([ax_lim(1) ax_lim(1)],[y_75, y_max],'linewidth',9,'color',dpc)
      
      title('All Stratocumulus')
      pred_str=['\Pi_{Sc}=', sprintf('%0.3f',allSc_test_rsq)];
      text(.7,.6,pred_str,'FontSize',24)%,'Color',axis_color)
      hold off
      %% instruction for reading the quartile plot
      subplot('Position',[0.72 0.3 0.2 0.35])
      set(gca,'Color','none')
      set(gca,'XColor','none')
      set(gca,'YColor','none')
      
      hold on
      
      
      
      
      
      box_minx=(ax_lim(2)-ax_lim(1))*.2+ax_lim(1);
      box_maxx=(ax_lim(2)-ax_lim(1))*.8+ax_lim(1);
      box_miny=(ax_lim(2)-ax_lim(1))*.1+ax_lim(1);
      box_maxy=(ax_lim(2)-ax_lim(1))*.9+ax_lim(1);
      box_wid=box_maxx-box_minx;
      box_len=box_maxy-box_miny;
      
      rectangle('Position',[box_minx box_miny box_wid box_len],...
         'EdgeColor',dpc,'LineWidth',1,'FaceColor',[.95 .95 .95])
      
      
      
      line([box_minx+box_wid*.2, box_minx+box_wid*.2],[box_miny+box_len*.15 box_miny+box_len*.4],...
         'linewidth',9,'color',dpc)
      plot(box_minx+box_wid*.2,box_miny+box_len*.5,'.','MarkerSize',20,'color',dpc)
      line([box_minx+box_wid*.2, box_minx+box_wid*.2],[box_miny+box_len*.6 box_miny+box_len*.85],...
         'linewidth',9,'color',dpc)
      
      text(box_minx+box_wid*.3,box_miny+box_len*.15, 'min','color',dpc,'FontSize',24)
      text(box_minx+box_wid*.3,box_miny+box_len*.4, '1st quartile','color',dpc,'FontSize',24)
      text(box_minx+box_wid*.3,box_miny+box_len*.5, 'median','color',dpc,'FontSize',24)
      text(box_minx+box_wid*.3,box_miny+box_len*.6, '3rd quartile','color',dpc,'FontSize',24)
      text(box_minx+box_wid*.3,box_miny+box_len*.85, 'max','color',dpc,'FontSize',24)
      
      set(gca,'XLim',ax_lim)
      set(gca,'YLim',ax_lim)
      
      hold off
      %         else
      
      %         end
   end
   
   for itype=2
      for c=1:5
         fh(itype)=figure(itype);
         sp_loc=c;
         if c>3; sp_loc=c+1; end
         subplot(2,4,sp_loc)
         
         %         set(gca,'Xcolor',axis_color)
         %         set(gca,'Ycolor',axis_color)
      end
      figure(itype)
      h_AxTop=createOverlayAxis(get(subplot(2,4,5)),get(subplot(2,4,3)));
      xlabel('Predicted relative dispersion')
      ylabel('Measured relative dispersion')
      set(gca,'FontSize',28)
   end
   
   
   
   
   
   
   %    saveas(figure(1),'plots/sp_meas_vs_pred_eps_u.png') %not using figure anymore
%        saveas(figure(2),'plots/gen_meas_vs_pred_eps_u.png')
end