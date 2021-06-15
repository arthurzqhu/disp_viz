clearvars -except clouds
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
if ~exist('clouds','var') load clouds.mat, end
campaigns={'vocalspdi','masepdi','postpdi','oraclespdi','gomaccspdi'};
camp_proper = {'VOCALS','MASE','POST','ORACLES','GoMACCS'};
typevar = {'s_ntot_aer', 's_ntot_pdi', 's_actfrac', 'AF', 'a_meanD_ex',...
    'ent_ratio_T'};
typename = {'N_{ia} (#/cc)',...
    'N_d (#/cc)', 'LAF',...
    'Adiabatic Fraction','Mean interstitial aerosol diameter [\mum]',...
    'Entrainment Ratio'};
typetitle = {'(a)','(b)','(c)','AF','D_a','ER'};
dp_color = '#FFBE89';%'#A8CCCC';
axis_color = '#638FAA';

%%

all_reldisp = cell(length(campaigns),1);
all_X = cell(length(campaigns),length(typevar));
normAC = cell(length(campaigns),1);
moist_t = cell(length(campaigns),1);

for c = 1:5
    %%
    camp = campaigns{c};
    fb = load([camp,'_flight_basics.mat']);
    fbvar = [camp,'_flight_basics'];
    
    fldnm = fieldnames(clouds.(camp)(1));
    
    % commontitles = 'INTRA-cloud correlation \epsilon vs ';
    set(groot, 'DefaultFigurePosition', [1399 1 913 984])
    
    for iday = 1:length(clouds.(camp))
        
        % get the unfiltered time first for later use
        s_t_raw = floor(clouds.(camp)(iday).s_t);
        a_t_raw = clouds.(camp)(iday).a_t;
        
        isOorG=any(strcmp({'oraclespdi';'gomaccspdi'},campaigns{c}));
        
        if isOorG
            cmt_raw = clouds.(camp)(iday).cmt;
            try typevar{contains(typevar,'normAC')}='a_normAC'; end
        else
            try typevar{contains(typevar,'normAC')}='normAC'; end
        end
        
        cloud_i = fb.(fbvar)(iday).ti;
        cloud_f = fb.(fbvar)(iday).tf;
        
        for icld = 1:length(cloud_i)
            ti = cloud_i(icld);
            tf = cloud_f(icld);
            
            % filtering criteria
            s_inleg_idx = s_t_raw > ti & s_t_raw < tf;
            a_inleg_idx = a_t_raw > ti & a_t_raw < tf;
            if isOorG; c_inleg_idx=cmt_raw>ti & cmt_raw<tf; end
            min_ntot_idx = clouds.(camp)(iday).s_ntot_pdi > 25;
            
            s_filt_crit = find(s_inleg_idx & min_ntot_idx);
            a_filt_crit = find(a_inleg_idx);
            if isOorG; c_filt_crit = find(c_inleg_idx); end
            
            s_t_leg = s_t_raw(s_filt_crit);
            a_t_leg = a_t_raw(a_filt_crit);
            if isOorG; c_t_leg = cmt_raw(c_filt_crit); end
            
            if ~isOorG
                [~, cmt_ipdi{iday,icld}, cmt_ipcasp{iday,icld}] = ...
                    intersect(s_t_leg,a_t_leg);
            else
                cmt=mintersect(s_t_leg,a_t_leg,c_t_leg);
                cmt_ipdi{iday,icld}=find(ismember(s_t_leg,cmt));
                cmt_ipcasp{iday,icld}=find(ismember(a_t_leg,cmt));
                cmt_ic{iday,icld}=find(ismember(c_t_leg,cmt));
            end
            
            s_all_filt_crit = s_filt_crit(cmt_ipdi{iday,icld});
            a_all_filt_crit = a_filt_crit(cmt_ipcasp{iday,icld});
            if isOorG
                c_all_filt_crit=c_filt_crit(cmt_ic{iday,icld});
                normAC{c} = [normAC{c}; clouds.(camp)(iday).a_normAC(a_all_filt_crit)];
            else
                normAC{c} = [normAC{c}; clouds.(camp)(iday).normAC(s_all_filt_crit)];
                moist_t{c} = [moist_t{c}; clouds.(camp)(iday).s_qt(s_all_filt_crit)];
            end
            
            for ivar = 1:length(typevar)
                if contains(typevar{ivar},'a_')
                    indvar = clouds.(camp)(iday).(typevar{ivar})(a_all_filt_crit);
                elseif isOorG && any(strcmp(typevar{ivar},{'AF','ent_ratio_T'}))
                    indvar = clouds.(camp)(iday).(typevar{ivar})(c_all_filt_crit);
                else
                    indvar = clouds.(camp)(iday).(typevar{ivar})(s_all_filt_crit);
                end
                indvar(isinf(indvar))=nan;
                all_X{c, ivar} = [all_X{c, ivar}; indvar];
            end
            depvar = clouds.(camp)(iday).s_disp_pdi(s_all_filt_crit);
            all_reldisp{c} = [all_reldisp{c}; depvar];
        end
    end
end

%% regression
c=1;
laf=all_X{c, 3};
af=all_X{c, 4};
y=all_reldisp{c};
nh=normAC{c};
s_qt=moist_t{c};

X1=[ones(size(laf)) laf];
X2=[ones(size(laf)) laf af];
X3=[ones(size(laf)) af];
X4=[ones(size(laf)) laf af nh];
X5=[ones(size(laf)) laf nh];
X6=[ones(size(laf)) laf nh s_qt af];

[b1, ~, res1, ~, stats1] = regress(y, X1);
[b2, ~, res2, ~, stats2] = regress(y, X2);
[b3, ~, res3, ~, stats3] = regress(y, X3);
[b4, ~, res4, ~, stats4] = regress(y, X4);
[b5, ~, res5, ~, stats5] = regress(y, X5);
[b6, ~, res6, ~, stats6] = regress(y, X6);
%%
vars=[1 2 3];

close all
for c = 1:5
    for icol = 1:3
        subplot(5,3,3*(c-1)+icol)
        ivar=vars(icol);
        
        xx=all_X{c, ivar};
        yy=all_reldisp{c};
        
        if any(strcmp(typevar{ivar},{'a_ntot_ex','s_ntot_pdi'}))
            xx=xx(xx<prctile(all_X{c, ivar},98));
            yy=yy(xx<prctile(all_X{c, ivar},98));
        end
        
        XX=[ones(size(xx)) xx];
        [bb, ~, rres, ~, sstats] = regress(yy, XX);
        
        scatter(xx,yy,20,'filled',...
            'MarkerFaceAlpha',0.1,'MarkerFaceColor',dp_color)
        hold on
        
        %         [comparison_mtx,hist_centers] = hist3([all_X{c, ivar},...
        %             all_reldisp{c}],'Nbins',[40,40]);
        %         contour(hist_centers{1},hist_centers{2},comparison_mtx','LineWidth',2)
        %         cbar = colorbar;
        %         cbar.Label.String = 'Frequency of occurrence';
        
        %         xlabel(typename{ivar})
%         ylim([0 1])

        
        
        if ivar == 3
            xlim([0 1])
        end
        
%         if c == 1
%             title(typetitle{ivar})%,'Color',axis_color)
%         end
        title(sprintf('R^2=%0.3f', sstats(1)))
        
        if c < 5
        else
            xlabel(typename{ivar})
        end
        
        if icol > 1
            set(gca,'YTick',[])
        else
            ylabel([camp_proper{c}, ' \epsilon'])
        end
        
        if strcmp(typevar{ivar},'AF')
            xlim([0 1])
        end
        
        xbins = linspace(min(xlim),max(xlim),11);
        mu_reldisp = zeros(length(xbins)-1,1);
        
        for ibin = 1:length(xbins)-1
            idx_at_ibin = xx>xbins(ibin) ...
                & xx<=xbins(ibin+1);
            
            if sum(idx_at_ibin)>25
                mu_reldisp(ibin) = nanmean(yy(idx_at_ibin));
                sigma_dv(ibin) = nanstd(yy(idx_at_ibin));
            else
                mu_reldisp(ibin) = nan;
                sigma_dv(ibin) = nan;
            end
        end
        
        xbinmean = (xbins(2:end)+xbins(1:end-1))/2;
        %
        errorbar(xbinmean,mu_reldisp,sigma_dv,'LineWidth',2)
        
        %         if ivar == 1
        %         %         title([commontitles,typename{ivar}]);
        %             ylabel('Relative dispersion \epsilon')
        %         else
        %             set(gca,'YTick',[])
        %         %         title(['vs ', typename{ivar}])
        %         end
        set(gca, 'fontsize', 16)
        
    end
    
    
    %     sgt = sgtitle('Intra-cloud correlation of relative dispersion with different environmental variables - VOCALS');
    %     sgt.FontSize = 20;
    %     sgt.FontWeight = 'bold';
end

% whitebg([.95 .95 .95])
% fig = gcf;
% fig.Color = '#F3F3F3';
% fig.InvertHardcopy = 'off';

for isubp = 1:15
    subplot(5,3,isubp)
    %     set(gca,'Xcolor',axis_color)
    %     set(gca,'Ycolor',axis_color)
end

% 
%  saveas(gcf, 'plots/INTRA-cloud_allcamps.png')

%% plot LAF vs var

close all

for icol = 4
    figure('Position',[1392 287 1244 593])
    tl = tiledlayout('flow');
    
    
    for ic = 1:5
        nexttile
        
        indvar = all_X{ic,3};
        depvar = all_X{ic,icol};
        
        if strcmp(typevar{icol},'AF')
            depvar(depvar>1)=nan;
            depvar(depvar<0)=0;
        end

        scatter(indvar,depvar,20,'filled',...
            'MarkerFaceAlpha',0.3,'MarkerFaceColor',dp_color)
        if strcmp(typevar{icol},'AF'); ylim([0 1]); end
        if strcmp(typetitle{icol},'ER'); ylim([-.1 .3]); end
        xlim([0 1])
        
        hold on
        
%         [comparison_mtx,hist_centers] = hist3([indvar,...
%             depvar],'Nbins',[30,30]);
%         contour(hist_centers{1},hist_centers{2},comparison_mtx','LineWidth',2)

        xbins = linspace(min(xlim),max(xlim),11);
        mu_dv = zeros(length(xbins)-1,1);
        sigma_dv = zeros(length(xbins)-1,1);

        for ibin = 1:length(xbins)-1
            idx_at_ibin = indvar>xbins(ibin) ...
                & indvar<=xbins(ibin+1);
            mu_dv(ibin) = nanmean(depvar(idx_at_ibin));
            sigma_dv(ibin) = nanstd(depvar(idx_at_ibin));
            
            if sum(idx_at_ibin)>100
                mu_dv(ibin) = nanmean(depvar(idx_at_ibin));
                sigma_dv(ibin) = nanstd(depvar(idx_at_ibin));
            else
                mu_dv(ibin) = nan;
                sigma_dv(ibin) = nan;
            end
        end
        

        xbinmean = (xbins(2:end)+xbins(1:end-1))/2;

        errorbar(xbinmean,mu_dv,sigma_dv,'LineWidth',2)
        
        title(gca,camp_proper{ic},'FontWeight','normal')
        set(gca,'fontsize',20)
    end
    
    xlabel(tl,'Local Activation Fraction','fontsize',24,'fontweight','bold')
    ylabel(tl,typename{icol},'fontsize',24,'fontweight','bold')
end

saveas(gcf,'plots/AF_vs_LAF.png')

%%
clear lt25d
for iday=1:21
    lt25=sum(clouds.gomaccspdi(iday).s_conc_pdi(:,1:70),2)./sum(clouds.gomaccspdi(iday).s_conc_pdi,2); 
    lt25d(iday)=nanmean(lt25);
end
