clear
close all
cd '~/Box/grad/research/aerosol_reldisp/datasets/'
load clouds.mat

campaigns={'vocals','mase','post','oracles','gomaccs'};
ctitle={'VOCALS','MASE','POST','ORACLES','GoMACCS'};
drp_instr={'pdi','pdi','pdi','pdi','pdi'};

%%
close all

clear camp_B_reldisp

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
    allb = cell(31,1);
    all_pred_dailyreg = [];
    all_coefcons_dailyreg = []; % the consistency/reldisp of the coefficients
    allX = [];
    ally = [];
    all_y_hat_dailyreg = [];
    icld = 0;
    ndays = length(clouds.(camp));
    days_analyzed = 1:ndays;
    
%     remove the days that have incomplete flights in vocals
    if c==1
%         days_analyzed(ismember(days_analyzed,[8,10,11,13]))=[];
    elseif c==4
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
                icld = icld + 1;
                
                % get the initial and final time for each cloud
                ti = cloudlegs_i(ileg);
                tf = cloudlegs_f(ileg);
                
                % filtering criteria
                s_filt_crit = find(s_t_unfilt > ti & s_t_unfilt < tf &...
                    clouds.(camp)(iday).(Nvar) > 25);
                a_filt_crit = find(a_t_unfilt > ti & a_t_unfilt < tf);
                
                s_t_leg = floor(clouds.(camp)(iday).s_t(s_filt_crit));
                
                if c == ~strcmp(camp,'masepdi')
                    a_t_leg = s_t_leg;
                else
                    a_t_leg = clouds.(camp)(iday).a_t(a_filt_crit);
                end
                
                [~, cmt_ipdi{c,iday,ileg}, cmt_ipcasp{c,iday,ileg}] = ...
                    intersect(s_t_leg,a_t_leg);
                
                s_all_filt_crit = s_filt_crit(cmt_ipdi{c,iday,ileg});
                a_all_filt_crit = a_filt_crit(cmt_ipcasp{c,iday,ileg});
                
                reldisp = clouds.(camp)(iday).(epsvar)(s_all_filt_crit);
                s_ntot = clouds.(camp)(iday).(Nvar)(s_all_filt_crit);
%                 s_ap = clouds.(camp)(iday).s_ap(s_all_filt_crit);
                if strcmp(camp,'oraclespdi') || strcmp(camp,'gomaccspdi')
                    normAC = clouds.(camp)(iday).a_normAC(a_all_filt_crit);
                    try
                        thet = clouds.(camp)(iday).a_thet(a_all_filt_crit);
                    catch
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
                ent_ratio_qt = clouds.(camp)(iday).ent_ratio_qt(s_all_filt_crit);
                
        
        
%                 if ismember(c,[1 2 5])
%                     s_rh = clouds.(camp)(iday).s_rh(s_all_filt_crit);
%                 elseif ismember(c, [3 4])
%                     s_rh = clouds.(camp)(iday).a_rh(a_all_filt_crit);
%                 end
%                 s_ss = (s_rh-100)/100;
%                 s_ss(s_ss<0) = nan;

                % setting the indvar for subsampling
                indvar_raw = normAC; 
                depvar_raw = reldisp;

                % make sure indvar and depvar are a pair of non-NaNs
                indvar = indvar_raw(~isnan(indvar_raw) & ~isnan(depvar_raw));
                depvar = depvar_raw(~isnan(indvar_raw) & ~isnan(depvar_raw)); 
                color = normAC(~isnan(indvar_raw) & ~isnan(depvar_raw));

                if sum(~isnan(depvar))>100
                    
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
                end
                
%                 figure
%                 scatter(indvar(rand_samp_idx), depvar(rand_samp_idx), [], color(rand_samp_idx), '.')
                x1 = a_ntot(rand_samp_idx);
                x2 = s_ntot(rand_samp_idx);
                x3 = s_actfrac(rand_samp_idx);
                x4 = normAC(rand_samp_idx);
                x5 = ent_ratio_qt(rand_samp_idx);
                
                y = reldisp(rand_samp_idx);
                
                X = [ones(size(x1)) x5 x4 x3 x2 x1];
%                 X(:,2:end) = X(:,2:end) - mean(X(:,2:end));
                
                % all possible combination of features expressed in binary
                feat_comb_bin_raw = logical(dec2bin(0:2^size(X,2)-1) - '0');
                
                % given that we need at least a constant (the bias term) and one feature
                % to predict, so everything starts with 0 or '100000' is not needed for
                % our purposes
                feat_comb_bin_unsorted = feat_comb_bin_raw(2^size(X,2)/2+2:end,:);
                
                % sorting feature combination binary matrix in ascending order of number
                % of features
                [~, fcb_idx] = sort(sum(feat_comb_bin_unsorted,2));
                feat_comb_bin = feat_comb_bin_unsorted(fcb_idx,:);
                
                % create a temporary data structure for y_hat_dailyreg and b (coef)
                y_hat_tmp = [];
%                 b_tmp = cell(31,1);
                for icombo = 1:size(feat_comb_bin,1)
                    % selected X:
                    X_sel = X(:,feat_comb_bin(icombo,:));
                    b = regress(y,X_sel);
                    y_hat = X_sel*b;
                    pred = 1-nansum((y_hat-y).^2)/nansum((y-nanmean(y)).^2);
                    all_pred_dailyreg(icombo,icld) = pred;
                    allb{icombo} = vertcat(allb{icombo}, b');
                    allb{icombo}(allb{icombo}==0)=nan; % set 0 to nan so that it wont affect the averaging
                    y_hat_tmp = [y_hat_tmp; y_hat'];
                end

%                 allb = horzcat(allb, b_tmp);
                all_y_hat_dailyreg = horzcat(all_y_hat_dailyreg, y_hat_tmp);
                allX = [allX;X];
                ally = [ally;y];

                
                
            end
        end
    end
    
%     meanB = mean(allb)';
%     camp_B_reldisp(c,:) = std(allb)./abs(meanB)';
%     all_y_hat_campmean = allX*meanB;
%     pred_campmean = 1-nansum((all_y_hat_campmean-ally).^2)/nansum((ally-nanmean(ally)).^2);
%     pred_dailyreg = 1-nansum((all_y_hat_dailyreg-ally).^2)/nansum((ally-nanmean(ally)).^2);
%     nanmean(allpred)
%     nanmedian(allpred)
%     abs(nanstd(allb)./nanmean(allb))
    %%
    mean_allb = cellfun(@nanmean, allb, 'UniformOutput', false);
    std_allb = cellfun(@nanstd, allb, 'UniformOutput', false);
    cons_allb = cellfun(@(x,y) abs(y)./abs(x), std_allb, mean_allb, 'UniformOutput',false);
    mean_cons_allb = cellfun(@mean, cons_allb);
    mean_pred_dailyreg = nanmean(all_pred_dailyreg,2);
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
%     all_y_hat_campmean = [];
%     all_pred_campmean = [];

    for icombo = 1:size(feat_comb_bin,1)
        allX_sel = allX(:,feat_comb_bin(icombo,:));
        [coeff{icombo,c},train_rsq(icombo,c), test_rsq(icombo,c)] = ...
            regress_tt(ally, allX_sel, .7);
        
        
        y_hat_camp_mean = allX_sel*coeff{icombo,c};
        pred_campmean = 1-nansum((y_hat_camp_mean-ally).^2)/nansum((ally-nanmean(ally)).^2);
%         all_y_hat_campmean = [all_y_hat_campmean; y_hat_camp_mean'];
%         all_pred_campmean = [all_pred_campmean; pred_campmean];
    end

    %% plot to see which combination of predictors has the highest consistency
%     close all
    
    all_xlabels = {'N_a            ', 'N_d            ', 'LAF            ',...
        'NH            ', 'ER            '};
    figure('Position',[1761 227 820 650]); %[92 148 1325 650])
    
    nsp = subplot(10,1,1:5);
    
    yyaxis left
    plot(mean_cons_allb,'LineWidth', 2,'HandleVisibility','off')
    ylabel('mean consistency of coefficients')
    color_order = get(gca,'ColorOrder');
    set(gca,'YColor',color_order(1,:))
    ylim([0 6])
    
    yyaxis right 
    
    hold on
    plot(mean_pred_dailyreg, 'LineWidth', 2, 'LineStyle',':','DisplayName',...
        'cloud R^2')
    plot(train_rsq(:,c), 'LineWidth', 2, 'LineStyle','-', 'DisplayName', 'campaign R^2 - training')
    plot(test_rsq(:,c), 'LineWidth', 1, 'LineStyle','--', 'DisplayName', 'campaign R^2 - test')
    
    ylabel('mean predictability of coefficients')
    set(gca, 'YGrid', 'off', 'XGrid', 'on')
    ylim([0 1])
    legend('show')
    xlim([1 31])
    set(gca, 'XTick',1:31)
    set(gca, 'XColor', 'none')
    

%     xlabel('Combination #')
    set(gca,'fontsize', 12)
    set(0,'DefaultLegendAutoUpdate','off')
    
    for isubp = 1:5
        sp(isubp) = subplot(10,1,isubp+5);
        scatter(find(feat_comb_bin(:,7-isubp)), 0*find(feat_comb_bin(:,7-isubp)), 100,...
            'filled')
%         set(gca,'YColor','none')
        set(gca,'YTick',[])
        set(gca,'Color','none')
        ylim([0 1])
        xlim([1 31])
        set(gca, 'XTick',1:31)
        set(gca, 'XTick',1:31)
        set(gca, 'XTickLabel', [])
        ylabel(all_xlabels{isubp})
        set(gca,'YColor',color_order(1,:))
        set(gca, 'YGrid', 'off', 'XGrid', 'on')
        set(gca,'fontsize', 12)
    end
    
    dy = sp(1).Position(2) - sp(2).Position(2);
    
    for isubp = 1:5
%         sp(isubp).Position(2) = sp(isubp).Position(2)+dy;
        sp(isubp).Position(4) = sp(isubp).Position(4)+dy;
    end
    

% saveas(gcf, ['plots/abstract/(new) why i chose certain predictors ', camp '.png'])

end