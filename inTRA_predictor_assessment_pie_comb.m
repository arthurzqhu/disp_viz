close all
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
clearvars -except clouds
if ~exist('clouds','var') load clouds.mat, end

campaigns={'vocals','mase','post','oracles','gomaccs'};
ctitle={'VOCALS','MASE','POST','ORACLES','GoMACCS'};
drp_instr={'pdi','pdi','pdi','pdi','pdi'};

color_order=colororder;
piecmap_hex=['#d7191c'; '#fdae61'; '#ffffbf'; '#abd9e9'; '#2c7bb6'];
piecmap=sscanf(piecmap_hex','#%2x%2x%2x',[3,size(piecmap_hex,1)]).'/255;
linestyle_order={'-','--',':','-.',':*',};

nc=length(campaigns);

%%
close all

max_ratio=[2 3 5 10 20 50];
iratio=3;

figure('Position',[1459 89 1122 896])
for c=1:nc
   subplot(nc,1,c)
   %     set(gca,'Color',[.8 .8 .8])
end


clear camp_B_reldisp

for c=1:nc
   
   clear b X y y_hat pred
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
   
   %%
   % close all
   epsvar=['s_disp_' drp_instr{c}];
   Nvar=['s_ntot_' drp_instr{c}];
   
   ndays=length(vdate_pdi);
   
   nbins=20;
   samp_frac=0.5;
   upper_limit=2e4;
   %     layers=3;
   layer_seg=[0 0.2 0.8 1.2];
   allb=cell(31,1);
   all_pred_dailyreg=[];
   all_coefcons_dailyreg=[]; % the consistency/reldisp of the coefficients
   allX=[];
   ally=[];
   all_y_hat_dailyreg=[];
   icld=0;
   sampsize_cld=[];
   pred_diag=[];
   ndays=length(clouds.(camp));
   days_analyzed=1:ndays;
   
   %     remove the days that have incomplete flights in vocals
   if contains(camp,'vocals')
%       days_analyzed(ismember(days_analyzed,[8,10,11,13]))=[];
   elseif contains(camp,'oracles')
      days_analyzed(ismember(days_analyzed, [2,3,5,14,15,16,17,18]))=[];
   end
   %%
   
   % close all
   for iday=days_analyzed
      %%
      % get the unfiltered time first for later use
      s_t=floor(clouds.(camp)(iday).s_t);
      
      cleg_i=fb.(fbvar)(iday).ti;
      cleg_f=fb.(fbvar)(iday).tf;
      %         T_BL=fb.(fbvar)(iday).T_BL;
      %         T_FB=fb.(fbvar)(iday).T_FB;
      %%
      if ~isempty(cleg_i)
         for ileg=1:length(cleg_i)
            icld=icld + 1;
            
            % get the initial and final time for each cloud
            ti=cleg_i(ileg);
            tf=cleg_f(ileg);
            
            % filtering criteria
            filt_crit=find(s_t > ti & s_t < tf &...
               clouds.(camp)(iday).(Nvar) > 25);
            
            s_t_leg=s_t(filt_crit);
            reldisp=clouds.(camp)(iday).(epsvar)(filt_crit);
            s_ntot_pdi=clouds.(camp)(iday).(Nvar)(filt_crit);
            s_ap=clouds.(camp)(iday).s_ap(filt_crit);
            normAC=clouds.(camp)(iday).normAC(filt_crit);
            thet=clouds.(camp)(iday).s_thet(filt_crit);
            AF=clouds.(camp)(iday).AF(filt_crit);
            ent_ratio_T=clouds.(camp)(iday).ent_ratio_T(filt_crit);
            
            s_ntot_aer=clouds.(camp)(iday).s_ntot_aer(filt_crit);
            s_ntot_pcasp=clouds.(camp)(iday).s_ntot_pcasp(filt_crit);
            s_actfrac=clouds.(camp)(iday).s_actfrac(filt_crit);
            s_lwc=clouds.(camp)(iday).s_lwc_pdi(filt_crit);
            s_meand=clouds.(camp)(iday).s_meand_pdi(filt_crit);
            s_qt=clouds.(camp)(iday).s_qt(filt_crit);
            ent_ratio_qt=clouds.(camp)(iday).ent_ratio_qt(filt_crit);
            

            % setting the indvar for subsampling
            indvar_raw=normAC;
            depvar_raw=reldisp;
            
            % make sure indvar and depvar are a pair of non-NaNs
            vidx=~isnan(indvar_raw+depvar_raw);
            indvar=indvar_raw(vidx);
            depvar=depvar_raw(vidx);
            color=normAC(vidx);
            
            if sum(~isnan(depvar))>100
               
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
            end

            x1=s_lwc(rand_samp_idx);
            x2=s_ntot_pdi(rand_samp_idx);
            x3=s_actfrac(rand_samp_idx);
            x4=normAC(rand_samp_idx);
            x5=AF(rand_samp_idx);

            y=reldisp(rand_samp_idx);

            X=[ones(size(x1)) x5 x4 x3 x2 x1];
            
            % all possible combination of features expressed in binary
            feat_comb_bin_raw=logical(dec2bin(0:2^size(X,2)-1) - '0');
            
            % given that we need at least a constant (the bias term) and one 
            % feature to predict, so everything starts with 0 or '100000' is 
            % not needed for our purposes
            feat_comb_bin_unsorted=feat_comb_bin_raw(2^size(X,2)/2+2:end,:);
            
            % sorting feature combination binary matrix in ascending order of number
            % of features
            [~, fcb_idx]=sort(sum(feat_comb_bin_unsorted,2));
            feat_comb_bin=feat_comb_bin_unsorted(fcb_idx,:);
            
            % create a temporary data structure for y_hat_dailyreg and b (coeff)
            y_hat_tmp=[];
            
            for icombo=1:size(feat_comb_bin,1)

               % selected the predictor set
               X_sel=X(:,feat_comb_bin(icombo,:));
               [b, ~, pred]=regress_tt(y, X_sel, 1);
               y_hat=X_sel*b;
               all_pred_dailyreg(icombo,icld)=pred;
               
               if icombo == 3
                  pred_diag=[pred_diag pred];
               end
               
               allb{icombo}=vertcat(allb{icombo}, b');
               allb{icombo}(allb{icombo}==0)=nan; % set 0 to nan so that it wont affect the averaging
               y_hat_tmp=[y_hat_tmp; y_hat'];
            end
            
            %                 allb=horzcat(allb, b_tmp);
            all_y_hat_dailyreg=horzcat(all_y_hat_dailyreg, y_hat_tmp);
            allX=[allX;X];
            ally=[ally;y];
            
            sampsize_cld=[sampsize_cld;sum(~isnan(x3))];
            
            
         end % ileg
      end % if ~isempty(cleg_i)
   end % iday
   
   %     meanB=mean(allb)';
   %     camp_B_reldisp(c,:)=std(allb)./abs(meanB)';
   %     all_y_hat_campmean=allX*meanB;
   %     pred_campmean=1-nansum((all_y_hat_campmean-ally).^2)/nansum((ally-nanmean(ally)).^2);
   %     pred_dailyreg=1-nansum((all_y_hat_dailyreg-ally).^2)/nansum((ally-nanmean(ally)).^2);
   %     nanmean(allpred)
   %     nanmedian(allpred)
   %     abs(nanstd(allb)./nanmean(allb))
   %%
   mean_allb=cellfun(@nanmean, allb, 'UniformOutput', false);
   std_allb=cellfun(@nanstd, allb, 'UniformOutput', false);
   cons_allb=cellfun(@(x,y) abs(y)./abs(x), std_allb, mean_allb, 'UniformOutput',false);
   mean_cons_allb=cellfun(@mean, cons_allb);
   
   %     mean_pred_dailyreg=(all_pred_dailyreg*sampsize_cld)/sum(sampsize_cld);
   % to ignore nans
   for icombo=1:size(feat_comb_bin,1)
      notnan_idx_dr=~isnan(all_y_hat_dailyreg(icombo,:));
      mean_pred_dailyreg(icombo)=1-nansum((all_y_hat_dailyreg(icombo,:)'-ally).^2)/nansum((ally(notnan_idx_dr)-nanmean(ally(notnan_idx_dr))).^2);
   end
   %
   %     figure('Position',[92 148 1325 650])
   %     plot(mean_cons_allb,'LineWidth', 2)
   %     ylabel('mean consistency of coefficients (1/\epsilon_{coef}; averaged over all predictors)')
   %     yyaxis right
   %     plot(mean_pred_dailyreg, 'LineWidth', 2)
   %     ylabel('mean predictability given the combination of features (averaged over all days)')
   %     grid
   %     ylim([-1 1])
   %
   %     xlabel('Combination #')
   %     set(gca,'fontsize', 16)
   %     saveas(gcf, ['plots/abstract/why i chose certain predictors ', camp '.png'])
   
   %% testing how predictable X is using mean coefficients
   %     all_y_hat_campmean=[];
   %     all_pred_campmean=[];
   
   for icombo=1:size(feat_comb_bin,1)
      allX_sel=allX(:,feat_comb_bin(icombo,:));
      [coeff{icombo,c},train_rsq(icombo,c), test_rsq(icombo,c)]=...
         regress_tt(ally, allX_sel, .7);
      
      
      y_hat_camp_mean=allX_sel*coeff{icombo,c};
      pred_campmean=1-nansum((y_hat_camp_mean-ally).^2)/nansum((ally-nanmean(ally)).^2);
      %         all_y_hat_campmean=[all_y_hat_campmean; y_hat_camp_mean'];
      %         all_pred_campmean=[all_pred_campmean; pred_campmean];
   end
   
   %% plot to see which combination of predictors has the highest consistency
   
   ax_pos=get(gca,'Position');
   
   subplot(5,1,c)
   ax_numb(c)=gca;
   yyaxis left
   set(gca,'YColor',color_order(6,:))
   ylim([0 6])
   xlim([1 31])
   yticks(0:3:6)
   set(gca, 'XTick',1:31)
   set(gca, 'XColor', 'none')
   ytk=get(gca,'ytick').';  % get the tick values as column vector
   set(gca,'yticklabel',horzcat(num2str(ytk), repelem(' ',length(ytk),1)))
   
   if strcmp('oraclespdi',camp)
      hold on
      pt2=patch([1 31 31 1], [max(ylim) max(ylim) 0 0],...
         [.8 .8 .8],'edgecolor','none');
      pt=patch([5.45 6.55 6.55 5.45], [max(ylim) max(ylim) 0 0],...
         [.99 .99 .99],'edgecolor','none');
      rectangle('Position',[5.5 0 1 6],'LineWidth',3,...
         'EdgeColor',piecmap(4,:))
      rectangle('Position',[5.5 0 1 6],'LineWidth',3,...
         'LineStyle','--','EdgeColor',piecmap(5,:))

   elseif strcmp('gomaccspdi',camp)
      hold on
      pt2=patch([1 31 31 1], [max(ylim) max(ylim) 0 0],...
         [.8 .8 .8],'edgecolor','none');
      pt=patch([13.45 14.55 14.55 13.45], [max(ylim) max(ylim) 0 0],...
         [.99 .99 .99],'edgecolor','none');
      rectangle('Position',[13.5 0 1 6],'LineWidth',3,...
         'EdgeColor',piecmap(1,:))
      rectangle('Position',[13.5 0 1 6],'LineWidth',3,...
         'LineStyle','--','EdgeColor',piecmap(3,:))

   else
      hold on
      pt2=patch([1 31 31 1], [max(ylim) max(ylim) 0 0],...
         [.8 .8 .8],'edgecolor','none');
      pt=patch([10.45 11.55 11.55 10.45], [max(ylim) max(ylim) 0 0],...
         [.99 .99 .99],'edgecolor','none');
      rectangle('Position',[10.5 0 1 6],'LineWidth',3,...
         'EdgeColor',piecmap(2,:))
      rectangle('Position',[10.5 0 1 6],'LineWidth',3,...
         'LineStyle','--','EdgeColor',piecmap(3,:))

   end
   
   set(get(get(pt2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
   set(get(get(pt,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
   
   plot(mean_cons_allb, '-*', 'LineWidth', 2,'HandleVisibility','off','Color',...
      color_order(6,:))
   
   yyaxis right
   %     hold on
   plot(mean_pred_dailyreg, ':o', 'LineWidth', 2, 'DisplayName',...
      '\Pi_s','Color',color_order(7,:))
   plot(test_rsq(:,c), '-o', 'LineWidth', 2,'DisplayName', '\Pi_c','Color',...
      color_order(7,:))
   yticks(0:.5:1)
   ytk=get(gca,'ytick').';  % get the tick values as column vector
   set(gca,'yticklabel',horzcat(repelem(' ',length(ytk),1),num2str(ytk)))
   %     ylabel('Predictability of a predictor set')
   set(gca, 'YGrid', 'off', 'XGrid', 'on')
   ylim([0 1])
   hleg=legend('hide');
   set(gca, 'YColor',color_order(7,:))
   set(gca, 'fontsize', 18)
   title(ctitle{c})
   cbar=colorbar;
   cbar.Ticks=.1:.2:.9;
   cbar.TickLabels={'AF', 'NH', 'LAF', 'N_d', 'LWC'};
   cbar.Position=getCbarPos(subplot(5,1,1),subplot(5,1,5));
   subplot(5,1,c)
   colormap(piecmap)
   bubblepie(1:31,zeros(31,1)-.08,ones(31,1),double(feat_comb_bin(:,2:end)))
   
   hold off
   
   
   
end

hleg.Visible='on';
set(0,'DefaultLegendAutoUpdate','off')
set(hleg,'Position',[0.8512, 0.9145, 0.0539, 0.0642]);
% text(-0.1,0.5,'Consistency of the reg. coeff. of a predictor set',...
%     'rotation',90, 'horizontalalignment','center',...
%     'verticalalignment','bottom');


% make up for the mysteriously lost bubblepie
subplot(5,1,1)
colormap(piecmap)
bubblepie(1:31,zeros(31,1)-.08,ones(31,1),double(feat_comb_bin(:,2:end)))

h_AxTop=createOverlayAxis(ax_numb(5),ax_numb(1));
xlabel('Combination of predictors')
ylabel('Consistency of the reg. coeff. of a predictor set')
set(gca,'FontSize', 18)
set(gca,'YColor',color_order(6,:))

% saveas(gcf, 'plots/(pie) why i chose certain predictors_u.png')
exportgraphics(gcf,'../paper in progress/figs/fig4.jpg','Resolution',300)