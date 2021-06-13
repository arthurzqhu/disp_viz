clear
close all
cd '~/Box/grad/research/aerosol_reldisp/datasets/'
load clouds.mat

campaigns={'vocals','mase','post','oracles','gomaccs'};
ctitle={'VOCALS','MASE','POST','ORACLES','GoMACCS'};
drp_instr={'pdi','pdi','pdi','pdi','pdi'};

%%
close all
set(0, 'DefaultFigurePosition',[763 65 1707 872])
for c = 1:5
     
    clear b X y y_hat pred
    
    clear slp_raw intcpt_raw rsq_raw
    camp = [campaigns{c} drp_instr{c}];
    fb = load([camp,'_flight_basics.mat']);
    fbvar = [camp,'_flight_basics'];
    
    % maxN = arrayfun(@(x) clouds.(camppdi)(x).maxN, 1:length(clouds.(camppdi)));
    min_dtpt = 100;
    
    % exclude the date where not enough data points are present
    dtpt_pdi = arrayfun(@(x) length(clouds.(camp)(x).s_t),1:length(clouds.(camp)));
    vdate_pdi = find(dtpt_pdi>min_dtpt)';
    

    % close all
    epsvar = ['s_disp_' drp_instr{c}];
    Nvar = ['s_ntot_' drp_instr{c}];

    ndays = length(vdate_pdi);
    
    nbins = 20;
    samp_frac = 0.5;
    upper_limit = 2e4;
%     layers = 3;
    layer_seg = [0 0.2 0.8 1.2];
    allb = [];
    all_pred_dailyreg = [];
    allX = [];
    ally = [];
    all_y_hat_dailyreg = [];
    
    ndays = length(clouds.(camp));
    days_analyzed = 1:ndays;
    
%     remove the days that have incomplete flights in vocals
    if strcmp(camp,'vocalspdi')
%         days_analyzed(ismember(days_analyzed,[8,10,11,13]))=[];
    elseif strcmp(camp,'oraclespdi')
        days_analyzed(ismember(days_analyzed, [2,3,5,14,15,16,17,18]))=[];
    end
    

% close all
    for iday = days_analyzed
        
        % get the unfiltered time first for later use
        s_t_unfilt = floor(clouds.(camp)(iday).s_t);
        a_t_unfilt = clouds.(camp)(iday).a_t;
        
        cloudlegs_i = fb.(fbvar)(iday).ti;
        cloudlegs_f = fb.(fbvar)(iday).tf;
%         T_BL = fb.(fbvar)(iday).T_BL;
%         T_FB = fb.(fbvar)(iday).T_FB;
        
        if ~isempty(cloudlegs_i)
            for ileg = 1:length(cloudlegs_i)
                
                % get the initial and final time for each cloud
                ti = cloudlegs_i(ileg);
                tf = cloudlegs_f(ileg);
                
                % filtering criteria
                s_filt_crit = find(s_t_unfilt > ti & s_t_unfilt < tf &...
                    clouds.(camp)(iday).(Nvar) > 25);
                a_filt_crit = find(a_t_unfilt > ti & a_t_unfilt < tf);
                
                s_t_leg = floor(clouds.(camp)(iday).s_t(s_filt_crit));
                
                if c == 5 
                    a_t_leg = s_t_leg;
                else
                    a_t_leg = clouds.(camp)(iday).a_t(a_filt_crit);
                end
                
                [cmt{c,iday,ileg}, cmt_ipdi{c,iday,ileg}, cmt_ipcasp{c,iday,ileg}] = ...
                    intersect(s_t_leg,a_t_leg);
                
                s_all_filt_crit = s_filt_crit(cmt_ipdi{c,iday,ileg});
                a_all_filt_crit = a_filt_crit(cmt_ipcasp{c,iday,ileg});
                
                reldisp = clouds.(camp)(iday).(epsvar)(s_all_filt_crit);
                s_ntot = clouds.(camp)(iday).(Nvar)(s_all_filt_crit);
%                 s_ap = clouds.(camp)(iday).s_ap(s_all_filt_crit);
                if strcmp(camp,'oraclespdi') || strcmp(camp,'gomaccspdi')
                    normAC = clouds.(camp)(iday).a_normAC(a_all_filt_crit);
                    try %#ok<TRYNC>
                        thet = clouds.(camp)(iday).a_thet(a_all_filt_crit);
                    end
                elseif ~strcmp(camp,'gomaccspdi')
                    normAC = clouds.(camp)(iday).normAC(s_all_filt_crit);
                    thet = clouds.(camp)(iday).s_thet(s_all_filt_crit);
                end
                a_ntot = clouds.(camp)(iday).a_ntot(a_all_filt_crit);
                a_ntot_ex = clouds.(camp)(iday).a_ntot_ex(a_all_filt_crit);
                s_actfrac = clouds.(camp)(iday).s_actfrac(s_all_filt_crit);
                s_lwc = clouds.(camp)(iday).s_lwc_pdi(s_all_filt_crit);
                
%                 s_qt = clouds.(camp)(iday).s_qt(s_all_filt_crit);
%                 ent_ratio_T = clouds.(camp)(iday).ent_ratio_T(s_all_filt_crit);
%                 ent_ratio_qt = clouds.(camp)(iday).ent_ratio_qt(s_all_filt_crit);
                
        
        
%                 if ismember(c,[1 2 5])
%                     s_rh = clouds.(camp)(iday).s_rh(s_all_filt_crit);
%                 elseif ismember(c, [3 4])
%                     s_rh = clouds.(camp)(iday).a_rh(a_all_filt_crit);
%                 end
%                 s_ss = (s_rh-100)/100;
%                 s_ss(s_ss<0) = nan;

                indvar_raw = normAC;
                depvar_raw = reldisp;

                % make sure indvar and depvar are a pair of non-NaNs
                indvar = indvar_raw(~isnan(indvar_raw) & ~isnan(depvar_raw));
                depvar = depvar_raw(~isnan(indvar_raw) & ~isnan(depvar_raw)); 
%                 color = normAC(~isnan(indvar_raw) & ~isnan(depvar_raw));

                if sum(~isnan(depvar))>100 && ~strcmp(camp,'gomaccspdi')
                    
                    % set the maximum subsample size of each bin
                    max_per_bin = floor(length(indvar)/nbins);
                    [N_raw,edges,bin] = histcounts(indvar,nbins);
                    
                    N = ceil(samp_frac*N_raw);
                    N(N>max_per_bin) = max_per_bin;
                    N(N<5) = 0;
                    
                    rand_samp_idx = [];
                    for ibin = 1:nbins
                        rand_samp_idx_bin = randsample(find(bin==ibin),N(ibin));
                        rand_samp_idx = [rand_samp_idx; rand_samp_idx_bin];
                    end
                    
                    % SubSampled data
                    indvar_ss = indvar(rand_samp_idx);
                    depvar_ss = depvar(rand_samp_idx);
                    bin_ss = bin(rand_samp_idx);
                    
                    % to prevent matrix dimension error below
                    if length(bin_ss) ~= length(N)
                        ac_arr_raw = accumarray(bin_ss,depvar_ss);
                        ac_arr_raw(ac_arr_raw==0) = [];
                        ac_arr(N~=0) = ac_arr_raw;
                        ac_arr(N==0) = 0;
                    end
                    
                    bindep = ac_arr./N;
                    bincenters=(edges(1:end-1)+edges(2:end))/2;
                elseif strcmp(camp,'gomaccspdi')
                    rand_samp_idx = true(size(s_ntot));
                end
                
%                 figure
%                 scatter(indvar(rand_samp_idx), depvar(rand_samp_idx), [], color(rand_samp_idx), '.')
                x1 = s_ntot(rand_samp_idx);
                x2 = normAC(rand_samp_idx);
                x3 = thet(rand_samp_idx);
                x4 = s_actfrac(rand_samp_idx);
%                 x5 = ent_ratio_T(rand_samp_idx);
%                 x6 = ent_ratio_qt(rand_samp_idx);
                
                y = reldisp(rand_samp_idx);
                
                X = [ones(size(x1)) x4 x2];
%                 X(:,2:end) = X(:,2:end) - mean(X(:,2:end));
                
                b = regress(y,X);
                y_hat = X*b;

                
                pred = 1-nansum((y_hat-y).^2)/nansum((y-nanmean(y)).^2);
                allX = [allX;X];
                ally = [ally;y];
                all_y_hat_dailyreg = [all_y_hat_dailyreg;y_hat];
                allb = [allb;b'];
                all_pred_dailyreg = [all_pred_dailyreg;pred];
                
            end
        end
    end
    
    meanB = regress(ally, allX);
    [coeff(c,:), train_rsq(c), test_rsq(c)] = regress_tt(ally, allX, 0.7);
    all_y_hat_campmean = allX*coeff(c,:)';
    pred_campmean = 1-nansum((all_y_hat_campmean-ally).^2)/nansum((ally-nanmean(ally)).^2);
    pred_dailyreg = 1-nansum((all_y_hat_dailyreg-ally).^2)/nansum((ally-nanmean(ally)).^2);
%     nanmean(allpred)
%     nanmedian(allpred)
%     abs(nanstd(allb)./nanmean(allb))


% 
    y_hat_combined = [all_y_hat_dailyreg all_y_hat_campmean];
    pred_combined = [pred_dailyreg pred_campmean];
    title_string = {'Self predictability', 'Cross predictability'};
    y_hat_string = {'dailyreg', 'campmean'};
    
    for itype = 1:2
        figure(itype)
        subplot(2,3,c)
        hold on
        plot([0 1], [0 1])
        dp_color = '#E16B8C';
%         close all
        sp = scatter(y_hat_combined(:,itype),ally,'filled','MarkerFaceAlpha',0.1,'MarkerFaceColor',dp_color);
        [comparison_mtx,hist_centers] = hist3([y_hat_combined(:,itype),ally],'nbins',[30,30]);
        contour(hist_centers{1},hist_centers{2},comparison_mtx','LineWidth',2)
        cbar = colorbar;
        cbar.Label.String = 'Frequency of occurrence';
        
        ax_lim = [0 1];
        xlim(ax_lim)
        ylim(ax_lim)
    %     set(gca,'XTick',[0.1 0.3 0.5])
    %     set(gca,'YTick',[.2 .4 .6 .8])
        set(gca,'fontsize',16)
        
        %% drawing the quartile plot
        y_quartl = num2cell(prctile(ally, 0:25:100));
        [y_min, y_25, y_med, y_75, y_max] = deal(y_quartl{:});

        x_quartl = num2cell(prctile(y_hat_combined(:,itype), 0:25:100));
        [x_min, x_25, x_med, x_75, x_max] = deal(x_quartl{:});


        line([x_min, x_25],[ax_lim(1) ax_lim(1)],'linewidth',3,'color',dp_color)
        plot(x_med,ax_lim(1),'.','MarkerSize',20,'color',dp_color)
        line([x_75, x_max],[ax_lim(1) ax_lim(1)],'linewidth',3,'color',dp_color)

        line([ax_lim(1) ax_lim(1)],[y_min, y_25],'linewidth',3,'color',dp_color)
        plot(ax_lim(1),y_med,'.','MarkerSize',20,'color',dp_color)
        line([ax_lim(1) ax_lim(1)],[y_75, y_max],'linewidth',3,'color',dp_color)

        xlabel('Predicted relative dispersion')
        ylabel('Measured relative dispersion')
        title(ctitle{c})
        
        %% instruction for reading the quartile plot
        box_minx = (ax_lim(2)-ax_lim(1))*0.725+ax_lim(1);
        box_maxx = (ax_lim(2)-ax_lim(1))*1+ax_lim(1);
        box_miny = (ax_lim(2)-ax_lim(1))*0.6+ax_lim(1);
        box_maxy = (ax_lim(2)-ax_lim(1))*1+ax_lim(1);
        box_wid = box_maxx-box_minx;
        box_len = box_maxy-box_miny;
        
        rectangle('Position',[box_minx box_miny box_wid box_len], 'FaceColor','#FFFFFF',...
            'EdgeColor',dp_color,'LineWidth',1)
        line([box_minx+box_wid*.2, box_minx+box_wid*.2],[box_miny+box_len*.15 box_miny+box_len*.4],...
            'linewidth',3,'color',dp_color)
        plot(box_minx+box_wid*.2,box_miny+box_len*.5,'.','MarkerSize',20,'color',dp_color)
        line([box_minx+box_wid*.2, box_minx+box_wid*.2],[box_miny+box_len*.6 box_miny+box_len*.85],...
            'linewidth',3,'color',dp_color)

        text(box_minx+box_wid*.3,box_miny+box_len*.15, 'min','color',dp_color)
        text(box_minx+box_wid*.3,box_miny+box_len*.4, '1st quartile','color',dp_color)
        text(box_minx+box_wid*.3,box_miny+box_len*.5, 'median','color',dp_color)
        text(box_minx+box_wid*.3,box_miny+box_len*.6, '3rd quartile','color',dp_color)
        text(box_minx+box_wid*.3,box_miny+box_len*.85, 'max','color',dp_color)

        
        hold off
        
        %% text for predictibility
        pred_str = ['\Pi = ', sprintf('%0.4f \n(predictability of fit)',pred_combined(itype))];
        text(box_minx,box_miny-box_len*0.2,pred_str)
    end
    
end

for itype = 1:2
    figure(itype)
    sgt = sgtitle(title_string{itype});
    sgt.FontSize = 30;
    sgt.FontWeight = 'bold';
end

% saveas(figure(1),'plots/abstract/self_meas_vs_pred_eps.png')
% saveas(figure(2),'plots/abstract/cross_meas_vs_pred_eps.png')