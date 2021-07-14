close all
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
clearvars -except clouds
if ~exist('clouds','var') load clouds.mat, end

title_combined = {'aerosol concentration', 'droplet concentration', 'activation fraction'};
subtitle_combined = {'(a)', '(b)', '(c)'};
Xlabel_combined = {'Aerosol concentration (#/cc)', 'Droplet concentration (#/cc)',...
   'Local activation fraction'};
camp_name = {'VOCALS','MASE','POST','ORACLES','GoMACCS'};
filename_combined = {'aer','drp','actfrac'};
campaigns={'vocalspdi','masepdi','postpdi','oraclespdi','gomaccspdi'};
% %
% % finishingTaskSound

% pi_a_int_ntot = [2 10 36 372 22000];
% pi_s_ntot = [21.3 76.9 201.2 564.6 1944.3];
% pi_a_ntot = pi_a_int_ntot + pi_s_ntot;
% pi_actfrac = pi_s_ntot./pi_a_ntot;
% pi_reldisp = [.4 .38 .35 .28 .28];
% pi_dbar = [16.6 15.1 12.7 8.6 7.6];
% pi_sigma = [6.7 5.7 4.5 2.4 2.1];

%%
color_order_hex = {'#B47157', '#D75455', '#6A4C9C', '#2060AF', '#4F726C'};
color_order = arrayfun(@(x) hex2rgb(color_order_hex{x}),1:length(color_order_hex),...
   'UniformOutput', false);
axis_color = '#638FAA';
close all

% c=2;
indvar_lumped = {[] [] []};
reldisp_lumped = {[] [] []};
set(groot, 'DefaultFigurePosition', [1441,-123,1302,350])
sgt = sgtitle('Inter-cloud correlation of relative dispersion with different environmental variables');
sgt.FontSize = 22;
sgt.FontWeight = 'bold';

% set(groot,'DefaultLegendAutoUpdate','off')
for c = 1:5
   camp=campaigns{c};
   
   % for iday = 1:length(clouds.(camp))
   %
   %     ipcasp = commondate_idx_in_pcasp(iday);
   %     ipdi = commondate_idx_in_pdi(iday);
   %
   %     [commontime, commontime_idx_in_pdi, commontime_idx_in_pcasp] = ...
   %         intersect(floor(clouds.vocalspdi(ipdi).s_t), clouds.vocalspcasp(ipcasp).a_t);
   %
   %     figure
   %     plot(clouds.vocalspdi(ipdi).s_ap(commontime_idx_in_pdi),...
   %         clouds.vocalspcasp(ipcasp).a_ntot(commontime_idx_in_pcasp),'.')
   %
   %     yyaxis right
   %     plot(clouds.vocalspdi(ipdi).s_ap(commontime_idx_in_pdi),...
   %         clouds.vocalspdi(ipdi).s_lwc_pdi(commontime_idx_in_pdi),'.')
   % end
   
   
   if strcmp(camp, 'postpdi')
      a_ntot_CB = arrayfun(@(x) wmean(clouds.(camp)(x).a_ntot_CB,...
         clouds.(camp)(x).a_ntot_CB_sampsize), 1:length(clouds.(camp)))';
      s_ntot_CB = arrayfun(@(x) wmean(clouds.(camp)(x).s_ntot_CB,...
         clouds.(camp)(x).s_ntot_CB_sampsize), 1:length(clouds.(camp)))';
      s_actfrac_CB = arrayfun(@(x) wmean(clouds.(camp)(x).s_actfrac_CB,...
         clouds.(camp)(x).s_actfrac_CB_sampsize), 1:length(clouds.(camp)))';
      reldisp_CB = arrayfun(@(x) wmean(clouds.(camp)(x).reldisp_CB,...
         clouds.(camp)(x).reldisp_CB_sampsize), 1:length(clouds.(camp)))';
      a_ntot_spsz = arrayfun(@(x) sum(clouds.(camp)(x).a_ntot_CB_sampsize),...
         1:length(clouds.(camp)))';
      s_ntot_spsz = arrayfun(@(x) sum(clouds.(camp)(x).s_ntot_CB_sampsize),...
         1:length(clouds.(camp)))';
      s_actfrac_spsz = arrayfun(@(x) sum(clouds.(camp)(x).s_actfrac_CB_sampsize),...
         1:length(clouds.(camp)))';
      reldisp_spsz = arrayfun(@(x) sum(clouds.(camp)(x).reldisp_CB_sampsize),...
         1:length(clouds.(camp)))';
   else
      a_ntot_CB = arrayfun(@(x) clouds.(camp)(x).a_ntot_CB, ...
         1:length(clouds.(camp)),'UniformOutput',false);
      a_ntot_spsz = arrayfun(@(x) clouds.(camp)(x).a_ntot_CB_sampsize,...
         1:length(clouds.(camp)),'UniformOutput',false);
      s_ntot_CB = arrayfun(@(x) clouds.(camp)(x).s_ntot_CB, ...
         1:length(clouds.(camp)),'UniformOutput',false);
      s_ntot_spsz = arrayfun(@(x) clouds.(camp)(x).s_ntot_CB_sampsize,...
         1:length(clouds.(camp)),'UniformOutput',false);
      s_actfrac_CB = arrayfun(@(x) clouds.(camp)(x).s_actfrac_CB, ...
         1:length(clouds.(camp)),'UniformOutput',false);
      s_actfrac_spsz = arrayfun(@(x) clouds.(camp)(x).s_actfrac_CB_sampsize,...
         1:length(clouds.(camp)),'UniformOutput',false);
      reldisp_CB = arrayfun(@(x) clouds.(camp)(x).reldisp_CB,...
         1:length(clouds.(camp)), 'UniformOutput',false);
      reldisp_spsz = arrayfun(@(x) clouds.(camp)(x).reldisp_CB_sampsize,...
         1:length(clouds.(camp)),'UniformOutput',false);
   end
   
   
   
   
   
   %%
   if ~strcmp(camp, 'postpdi')
      a_ntot_CB = cell2mat(a_ntot_CB)';
      a_ntot_spsz = cell2mat(a_ntot_spsz)';
      
      s_ntot_CB = cell2mat(s_ntot_CB)';
      s_ntot_spsz = cell2mat(s_ntot_spsz)';
      
      s_actfrac_CB = cell2mat(s_actfrac_CB)';
      s_actfrac_spsz = cell2mat(s_actfrac_spsz)';
      
      reldisp_CB = cell2mat(reldisp_CB)';
      reldisp_spsz = cell2mat(reldisp_spsz)';
   end
   
   %%
   % presumably Chandrakar's injected aerosol is mostly, if not all, CCNs?
   indvar_combined = {a_ntot_CB s_ntot_CB s_actfrac_CB};
   indvar_spsz_combined = {a_ntot_spsz s_ntot_spsz s_actfrac_spsz};
   
   
   for itype = 1:length(indvar_combined)
      cf = subplot(1,3,itype);
      
      
      %         cf = figure(itype); % current figure
      vidx = find(indvar_spsz_combined{itype}>25 & reldisp_spsz>25);
      x = indvar_combined{itype}(vidx);
      y = reldisp_CB(vidx);
      [xmax_v, xmax_i] = max(x);
      
      hold on
      hSc = scatter(x, y, 'filled','DisplayName',camp_name{c},'CData',...
         color_order{c});
      
      % ylim([.1 .4])
      
      
      if itype == 3
         p = polyfit(x, y, 1);
         y_hat = polyval(p, x);
         plot(x, y_hat, '-', 'Color', hSc.CData,...
            'HandleVisibility','off')
         stats = regstats(y, x, 'linear');
         
      else
         p = polyfit(log10(x), y, 1);
         y_hat = polyval(p, log10(x));
         plot(x, y_hat, '-', 'Color', ...
            hSc.CData, 'HandleVisibility','off')
         set(gca,'Xscale','log')
         stats = regstats(y, log10(x), 'linear');
         
      end
      pval = stats.tstat.pval;
      rsq = stats.rsquare;
      
      textx=xmax_v*1.05;
      texty=y_hat(xmax_i)*1.05;
      
      set(gca,'fontsize',20)
      xlabel(Xlabel_combined{itype})
      title(subtitle_combined{itype})%,'Color',axis_color)
      
      if itype == 1
         ylabel('Relative dispersion \epsilon')
      else
         set(gca,'YTick',[])
      end
      
      
      if textx>max(x)
         textx=max(x)*.95;
      elseif textx<min(x)
         textx=min(x)*1.05;
      end
      
      if texty>max(y_hat)
         texty=max(y_hat)*.95;
      elseif texty<min(y_hat)
         texty=min(y_hat)*1.05;
      end

      text(textx, texty,sprintf('%.4f', rsq),...
         'Color',hSc.CData,'FontSize',16, 'FontWeight', 'bold')
      
      
      %         indvar_lumped{itype} = [indvar_lumped{itype};indvar_combined{itype}(vidx)];
      %         reldisp_lumped{itype} = [reldisp_lumped{itype};reldisp_CB(vidx)];
      %
      %         [indvar_sorted{itype}, indvar_sort_idx{itype}] = sort(indvar_lumped{itype});
      %         reldisp_sorted{itype} = reldisp_lumped{itype}(indvar_sort_idx{itype});
   end
   
   %     sprintf('p-value is %f, R^2 is %f', pval(2), rsq)
   % finishingTaskSound
end



% plot(pi_a_ntot, pi_reldisp,'--', 'LineWidth',2)
% %%
% for itype = 1:length(indvar_combined)
%     subplot(1,3,itype)
%     p = polyfit(log(indvar_sorted{itype}), reldisp_sorted{itype},1);
%     log_fit_line = polyval(p, log(indvar_sorted{itype}));
%     plot(indvar_sorted{itype}, p(2)+p(1)*log(indvar_sorted{itype}), ':', 'LineWidth',2)
%
%
%     %     log_fit_line_eqn = sprintf('$y = %0.4f { } log(x) + %0.4f$ \n$p$-value $=%0.4e$ \n$R^2=%0.4f$',...
%     %         p(1), p(2), pval(2), rsq);
%
%     log_fit_line_eqn = sprintf('$R^2=%0.4f$', rsq);
%
%     text(indvar_sorted{itype}(end)*0.3,log_fit_line(end)*1.05,log_fit_line_eqn,'Interpreter',...
%         'latex','FontSize',16)
%
% end

hleg = legend('show');
hleg.Position=[0.9094, 0.7057, 0.0864, 0.2671];
hleg.FontSize=16;

% hleg.EdgeColor=axis_color;
% hleg.TextColor=axis_color;

% sgt.Color = axis_color;
% set(0,'defaultfigurecolor','#D8D8D8')
% whitebg([.95 .95 .95])
% fig = gcf;
% fig.Color = '#F3F3F3';
% fig.InvertHardcopy = 'off';

% for itype = 1:length(indvar_combined)
%     subplot(1,3,itype)
%     set(gca,'Xcolor',axis_color)
%     set(gca,'Ycolor',axis_color)
% end

% set(gca,'Xcolor',axis_color)
% set(gca,'Ycolor',axis_color)
% saveas(gcf,'plots/INTER-cloud_3.jpg')

exportgraphics(gcf,'../paper in progress/figs/fig1.jpg','Resolution',300)

% %%
% % pause
% close all
% scatter(indvar_lumped,reldisp_lumped); hold on
% plot(pi_a_ntot, pi_reldisp,'--');
%
% set(gca,'Xscale','log')
%
% stats = regstats(reldisp_lumped, log(indvar_lumped), 'linear');
% pval = stats.tstat.pval;
% rsq = stats.rsquare;
% sprintf('p-value is %f, R^2 is %f', pval(2), rsq)
