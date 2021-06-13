clear
cd '~/Box/grad/research/datasets/PDI Data/'
load obs_clouds_wholeday.mat

%%
campaigns={'vocalspdi','masepdi','postpdi'};
ctitle={'VOCALS','MASE','POST'};
instr={'pdi','pdi','pdi'};

thresN=75;
thresRH=80;
thresT=3; %minimum temperature, just to make sure we avoid mixed-phase
thresLWC=0.05;
thresPTS=100; %minimum number of data samples

nlayers = 76;
layer_left_edges = linspace(0,0.75,nlayers);
layer_thickness = 0.3;
layer_edges = [layer_left_edges' layer_left_edges'+layer_thickness];
layer_mean=(layer_edges(1:end-1)+layer_edges(2:end))/2; % layer center diameter

do_plot = true; % whether to plot the slopes figure
do_animate = true;
do_save = true;

%%
clear wghts wghts_day
% for c = 1:length(campaigns)

% -------------------------- tweakable stuff for trial runs --------------------------
c=3; 

% nlayers=1;
% days_analyzed=1;
% ------------------------------------------------------------------------------------

% for c = 1:length(campaigns)
    camp=campaigns{c}; ndays = length(clouds.(camp)); days_analyzed = 1:ndays;
    %Get variable names
    Nvar=['s_ntot_',instr{c}];
    nvar=['s_conc_',instr{c}];
    Lvar=['s_lwc_',instr{c}];
    Dvar='drpsz';
    epsvar=['s_disp_',instr{c}]; %relative dispersion
    stdvar=['s_std_',instr{c}];
    Tvar='s_ta';
    RHvar='s_rh';
    ACvar='normAC';
    Rvar=['s_R_',instr{c}];

    if c==1
        days_analyzed(days_analyzed==10 | days_analyzed==13 | days_analyzed==11)=[]; 
        % remove the days that have incomplete flights in vocals
    end

    for icl = 65 %1:nlayers % iterator for different each layer of the cloud
        for iday = 1%days_analyzed
            binlims = clouds.(camp)(iday).binlims;
            binmean = clouds.(camp)(iday).binmean;
            dlogD = log10(binlims(2)/binlims(1));

            t = clouds.(camp)(iday).s_t;
            AC = clouds.(camp)(iday).(ACvar);
            t_cl = t(AC>layer_edges(icl,1) & AC<layer_edges(icl,2));
            t_cl_idx = ismember(t,t_cl);

            [alleps,allstd,alll,alln,alld,allR,allac]=deal([]);

            %Find all data that meet the thresholds
            %Remember clouds contains all data with non-NaN relative dispersion
            cldpts_all_cp=find(clouds.(camp)(iday).(Lvar)>thresLWC & ...
                clouds.(camp)(iday).(Nvar)>thresN & ...
                clouds.(camp)(iday).(Tvar)>thresT & ...
                clouds.(camp)(iday).(RHvar)>thresRH);

            % find the cloud droplets that are both valid by the standard above
            % and in the specific cloud layer
            cldpts = intersect(cldpts_all_cp, find(t_cl_idx));
            ncldpts(c,iday,icl) = length(cldpts);
            
            if ncldpts(c,iday,icl)>=thresPTS
                alln=clouds.(camp)(iday).(Nvar)(cldpts,:);
                alld=clouds.(camp)(iday).(Dvar)(cldpts,:);
                alll=clouds.(camp)(iday).(Lvar)(cldpts,:);
                allnd=clouds.(camp)(iday).(nvar)(cldpts,:);
                allac=clouds.(camp)(iday).(ACvar)(cldpts,:);
                meand(c,iday,icl) = mean(alld);
                allstd=clouds.(camp)(iday).(stdvar)(cldpts);
                alleps=clouds.(camp)(iday).(epsvar)(cldpts);
                alldev_gamma = clouds.(camp)(iday).dev_gamma(cldpts);

                allm0=sum(dlogD*allnd.*binmean'.^0,2);
                allm1=sum(dlogD*allnd.*binmean'.^1,2);
                allm2=sum(dlogD*allnd.*binmean'.^2,2);
                allm3=sum(dlogD*allnd.*binmean'.^3,2);
                allm4=sum(dlogD*allnd.*binmean'.^4,2);
                allm5=sum(dlogD*allnd.*binmean'.^5,2);
                allm6=sum(dlogD*allnd.*binmean'.^6,2);
                allR=clouds.(camp)(iday).(Rvar)(cldpts,:);
            else
                ncldpts(c,iday,icl) = 0;
            end
%             if ~isempty(alleps)
%                 AC_D_relation{c}(iday,icl,:) = polyfit(allac, alld, 1);
%                 X_raw = [ones(size(allm0)) allm0 allm1 allm3]; % allm1 == N * D_bar
%                 [X_norm, X_mu, X_sigma] = featureNormalize(X_raw);
% 
%                 % wghts(c,iday,icl,:) = normalEqn(X_norm, alleps);
%                 initial_wghts = zeros(size(X_norm,2),1);
%                 lambda = 1;
% 
%                 options = optimset('GradObj', 'on', 'MaxIter', 50, 'Display', 'off');
%                 wghts{c}(:,iday,icl) = ...
%                     fminunc(@(t)(costFunction(t, X_norm, alleps, lambda)), initial_wghts, options);
%                 eps_imp = X_norm*wghts{c}(:,iday,icl);
%                 m = size(X_raw,1);
%                 cost(c,iday,icl) = 1/(2*m)*sum((eps_imp-alleps).^2);
%             end
        end
    end
    
%     try
%         wghts{c}(wghts{c}==0)=nan;
%         meand(meand==0) = nan;
%     catch
%     end
    
%     for icl = 1:nlayers
%         wghts_dailyavg(c,:,icl) = nanmean(wghts{c}(1:4,:,icl), 2);
%         meand_dailyavg(c,icl) = nanmean(meand(c,:,icl));
%         cost_dailyavg(c,icl) = nanmean(cost(c,:,icl));
%     end
%     
%     for iday = days_analyzed
%         wghts_clavg(c,:,iday) = nanmean(wghts{c}(1:4,iday,:), 3);
%         meand_clavg(c,iday) = nanmean(meand(c,iday,:));
%         cost_clavg(c,iday) = nanmean(cost(c,iday,:));
%     end
% end

%% animating
close all
if do_animate==1
    f1 = figure('position',[183 505 1800 480]);
    for i=1:length(campaigns)
        a(i)=subplot(1,length(campaigns),i);
    end

    for icl = 1:nlayers
        cld_ax = axes('Position',[0.05, 0.16, 0.02, 0.7]);
        rectangle('Position',[0 layer_left_edges(icl) 1 layer_thickness],...
            'FaceColor',[.8,.8,.8],'EdgeColor','none')
        set(cld_ax,'XColor','none')
        ylim([0 max(layer_edges(:))])
        ylabel('Cloud')
        % prevent subsequent annotations from overlapping
        try
            delete(ant.Parent.Children)
        catch
        end
%         c = 2;
        for c=1:length(campaigns)
            axes(a(c))
            title(ctitle{c})
            pos = get(a, 'Position');
            pos{c}(1) = pos{c}(1)+1/15;
            pos{c}(2) = pos{c}(2);
            str = ['$\overline {D}$ = ', sprintf('%0.1f', meand_dailyavg(c,icl)),...
                ' $\mu$m, ', 'cost=',sprintf('%0.2e', cost_dailyavg(c,icl))];
            ant = annotation(f1,'textbox', pos{c}, 'String',str,...
                'FitBoxToText','on','Interpreter','latex','FontSize',18);
            ant.LineStyle='none';
            bar(wghts_dailyavg(c,:,icl))
            ylim([min(wghts_dailyavg(:)) max(wghts_dailyavg(:))])
            xticklabels({'Bias', '0th', '1st', '3rd'})
            xlim([0 5])
            xticks(1:4)
            xlabel('k in M_k (moment of the distribution)')
            ylabel('Weights value')
            set(gca,'FontSize',16)
            grid
        end
        
        F(icl) = getframe(gcf);
    end
    
    %%
    f2 = figure('position',[183 505 1800 480]);
    for i=1:length(campaigns)
        a(i)=subplot(1,length(campaigns),i);
    end
    
    for iday = days_analyzed
        try
            delete(ant.Parent.Children)
        catch
        end
        for c=1:length(campaigns)
            axes(a(c))
            title(ctitle{c})
            pos = get(a, 'Position');
            pos{c}(1) = pos{c}(1)+1/15;
            pos{c}(2) = pos{c}(2);
            str = ['$\overline {D}$ = ', sprintf('%0.1f', meand_clavg(c,iday)),' $\mu$m'];
            ant = annotation(f2,'textbox', pos{c}, 'String',str,...
                'FitBoxToText','on','Interpreter','latex','FontSize',18);
            ant.LineStyle='none';
            bar(wghts_clavg(c,:,iday))
            ylim([min(wghts_clavg(:)) max(wghts_clavg(:))])
            xticklabels({'Bias', '0th', '1st', '3rd'})
            xlim([0 5])
            xticks(1:4)
            xlabel('k in M_k (moment of the distribution)')
            ylabel('Weights value')
            set(gca,'FontSize',16)
            grid
        end
    end
    
    
end

if do_save==1
    v = VideoWriter('weights_vs_layer.mp4','MPEG-4');
    v.FrameRate=10;
    open(v)
    writeVideo(v,F)
    close(v)
end

%% linreg
close all
train_ratio = 0.7;
% p = polyfit(allac, alld, 1);
% fakeallm1_ac = polyval(p,allac).*allm0;
% fakealld = polyval(p,allac);

% fakeallm1 = (allm3./allm0).^(1/3).*allm0; 
% delta_dist = 1 - mean(allm1./fakeallm1);
% fakeallm1_m3 = fakeallm1 - delta_dist*allm1;

% X_dev_train = gen_rgrsnfeat([allm0 allm3], 1, 0);
% wghts_dev = normalEqn(X_dev_train.normmtx, allm1(X_dev_train.idx));
% 
% fakedev = X_dev_train.normmtx*wghts_dev;

% test_m1 = allm1.*(1+randn(length(allm1),1)/75);

% [X_train, ~] = gen_rgrsnfeat([allm0 allm1 allm3], train_ratio, 1);
% wghts = normalEqn(X_train.normmtx, alleps(X_train.idx));


[X_imp, X_cv] = gen_rgrsnfeat([allm0 allm1 allm3], 1, 0);
wghts = normalEqn(X_imp.normmtx, alleps(X_imp.idx));

% initial_wghts = zeros(size(X_norm,2),1);
% lambda = 0;
% 
% options = optimset('GradObj', 'on', 'MaxIter', 50, 'Display', 'off');
% wghts = ...
%     fminunc(@(t)(costFunction(t, X_norm, alleps, lambda)), initial_wghts, options);
% 
% X_raw2 = [ones(size(allm0)) allm0 (allm3./allm0).^(1/3).*allm0 allm3]; % allm1 == N * D_bar
% [X_norm2, X_mu2, X_sigma2] = featureNormalize(X_raw2);
eps_imp = X_imp.normmtx*wghts;
% eps_cv = X_cv.normmtx*wghts;
m_imp = length(X_imp.idx);
cost_imp = 1/(2*m_imp)*sum((eps_imp-alleps(X_imp.idx)).^2);
R_imp = corrcoef(eps_imp,alleps(X_imp.idx));
R_imp_sq = R_imp(1,2)^2;
% m_cv = length(X_cv.idx);
% cost_cv = 1/(2*m_cv)*sum((eps_cv-alleps(X_cv.idx)).^2) 
% calculated from the cross validation set to check how good the prediction
% is at estimating 
% eps_hat2 = X_norm2*wghts;

figure()
plot(eps_imp,alleps(X_imp.idx),'.')
rl = refline(1,0);
text(0.1*1.05, 0.4*0.97, sprintf('R^2 = %0.3f', R_imp_sq), 'FontSize',16)
% axis('equal')
xlabel('predicted reldisp')
ylabel('measured reldisp')
set(gca,'FontSize',16)
ylim([0.1 0.4])
xlim([0.1 0.4])
legend(rl, 'y = x')
% title([ctitle{c},' - Predicting w/ M_0, fake M_1, & M_3'])

% saveas(gcf,['groupmeeting/20190529/',ctitle{c},'_pred_meas_m0m3fm1.png'])

% figure()
% hold on; plot(alleps(X_imp.idx),'.'); plot(eps_imp); hold off
% xlabel('time series')
% ylabel('reldisp')
% legend('measured', 'predicted')
% set(gca,'fontsize',16)

% calculating cost

%%
scatter3(allm0(alldev_gamma<0.02), allm3(alldev_gamma<0.02), alleps(alldev_gamma<0.02), [], alldev_gamma(alldev_gamma<0.02))
c = colorbar;
c.Label.String = 'deviation from gamma';
xlabel('M_0')
ylabel('M_3')
zlabel('reldisp')