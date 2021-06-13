clear
cd '~/Box/grad/research/datasets/PDI Data/'
load obs_clouds_wholeday.mat

%%
clear meanvals meanvalsn

camps={'vocalspdi','masepdi','postpdi'};
ctitle={'VOCALS','MASE','POST'};
instr={'pdi','pdi','pdi'};


thresN=75;
thresRH=80;
thresT=3; %minimum temperature, just to make sure we avoid mixed-phase
thresLWC=0.01;
thresPTS=100; %minimum number of data samples
doresamp = false;

nlayers = 3;
z_l = -0.1;
z_h = 1.1;
layer_edges = [0 0.5 0.75 1.05]; %[-.05, .3, .7, 1.05]; %
layer_mean=(layer_edges(1:end-1)+layer_edges(2:end))/2; % layer center diameter

Moments = 0:0.25:2.5;

do_plot = true; % whether to plot the slopes figure
do_plot_corr = true;
do_plot_slope = false;
do_save = false;



%% separate out the cloud top mid bottom

clear slp slpn intcpt intcptn stdcam_intcpt stdcam_slp fitsampsize
clear binstd binstdn fitsampsizen mstdcam_slp_day mstdcam_slpn_day
clear mstdcam_intcpt_day mstdcam_intcptn_day mslp_day mslpn_day
clear mintcpt_day mintcptn_day meand

c=2;

close all

% for c=1:length(camps)
    
    campaign=camps{c};
    nc = length(clouds.(campaign));

    %Get variable names
    Nvar=['s_ntot_',instr{c}];
    nvar=['s_conc_',instr{c}];
    Hvar='s_ap';
    Vvar=['s_lwc_',instr{c}];
    Dvar='drpsz';
    Pvar='s_lwc_hvps';
    epsvar=['s_disp_',instr{c}]; %relative dispersion
    stdvar=['s_std_',instr{c}];
    Tvar='s_ta';
    RHvar='s_rh';
    distvar=['s_conc_',instr{c}];
    ACvar='normAC';
    Lvar=['s_lwc_',instr{c}];
    Thetvar='s_thet';
    Reffvar=['s_reff_',instr{c}];
    ndays = nc;

    days_analyzed = 1:ndays;
    if c==1
        days_analyzed(days_analyzed==10 | days_analyzed==13 | days_analyzed==11)=[]; % remove the 
    % days that have incomplete flights in vocals
    end
    for im = 1:length(Moments)
        
        for icl = 1:nlayers % iterator for different each layer of the cloud
            
            for iday = days_analyzed
            %     c = 1; 
                binlims = clouds.(campaign)(iday).binlims;
                binmean = clouds.(campaign)(iday).binmean;
                dlogD = log10(binlims(2)/binlims(1));


                t = clouds.(campaign)(iday).s_t;
                AC = clouds.(campaign)(iday).(ACvar);
                t_cl = t(AC>layer_edges(icl) & AC<layer_edges(icl+1));
                t_cl_idx = ismember(t,t_cl);

                [alleps,allstd,alln,allnormN,maxn,alld,maxd,mind,allnormD,...
                    allh,maxh,minh,allRH,maxRH,minRH,allnormRH,allsupsat,...
                    maxsupsat,alll,allnormL,maxl,allthet,allnormThet,maxthet,...
                    minthet,allreff,allnormReff,maxreff,minreff,allm1,allm2,...
                    allm3,allm6,allm1_5]=deal([]);
                %Find all data that meet the thresholds
                %Remember clouds contains all data with non-NaN relative dispersion
                cldpts_all_cp=find(clouds.(campaign)(iday).(Vvar)>thresLWC & ...
                    clouds.(campaign)(iday).(Nvar)>thresN & ...
                    clouds.(campaign)(iday).(Tvar)>thresT & ...
                    clouds.(campaign)(iday).(RHvar)>thresRH);

                % find the cloud droplets that are both valid by the standard above
                % and in the specific cloud layer
                cldpts = intersect(cldpts_all_cp, find(t_cl_idx));

                ncldpts(iday,icl) = length(cldpts);

                %If there are enough data points, give data a more convenient name,
                %and collect data in cells (the cell data isn't actually used).

                if ncldpts(iday,icl)>=thresPTS
                    alld=clouds.(campaign)(iday).(Dvar)(cldpts,:);
                    allnd=clouds.(campaign)(iday).(nvar)(cldpts,:);
                    allmx=sum(dlogD*allnd.*binmean'.^Moments(im),2);

                    meand(iday,icl) = mean(alld);
                    allstd=clouds.(campaign)(iday).(stdvar)(cldpts);
                    alleps=clouds.(campaign)(iday).(epsvar)(cldpts);
                else
                    ncldpts(iday,icl) = 0;
                end        

                if ~isempty(allstd) %If data exists
    %             try
                    nbins=20; %Number of cloud droplet concentration bins
                    indvar = allmx; 
%                     indnormVar = indvar/max(indvar);
                    %try
                    %Do some binning for cloud local data
                    [N,edges,bin]=histcounts(indvar,nbins);
                    %Find average relative dispersion in each bin
                    binstd = accumarray(bin,allstd)./N';
                    %Set to NaN if there are fewer than 5 points in the bin
                    binstd(N<5)=NaN;
                    meanvals{icl}(:,iday)=binstd;
                    centers=(edges(1:end-1)+edges(2:end))/2;
                    ref_binstd = sqrt((centers-1)/2);

                    %Also fit a line to the data
                    inds=~isnan(binstd);
                    x=centers(inds);
                    y=binstd(inds);


                    if length(x)>1
                        fitline = fit(x',y,'poly1');
                        slp(im,iday,icl) = fitline.p1;
                        intcpt(im,iday,icl) = fitline.p2;
                        fitsampsize(im,iday,icl) = sum(inds);
                        fitsampsize(fitsampsize(im,iday,icl)<5)=0;
                    else
                        slp(im,iday,icl) = nan;
                        intcpt(im,iday,icl) = nan;
                        fitsampsize(im,iday,icl) = 0;
                    end

                    %Do some binning now for reldisp
                    [N,edges,bin]=histcounts(indvar,nbins);
                    binstdn = accumarray(bin,alleps)./N';
                    binstdn(N<5)=NaN;
                    meanvalsn{icl}(:,iday)=binstdn;
                    centersn=(edges(1:end-1)+edges(2:end))/2;

                    %Also fit a line to the data
                    inds=~isnan(binstdn);
                    x=centersn(inds);
                    y=binstdn(inds);

                    if length(x)>1
                        fitline = fit(x',y,'poly1');
                        slpn(im,iday,icl) = fitline.p1;
                        intcptn(im,iday,icl) = fitline.p2;
                        fitsampsizen(im,iday,icl) = sum(inds);
                        fitsampsizen(fitsampsize(im,iday,icl)<5)=0;
                    else
                        slpn(im,iday,icl) = nan;
                        intcptn(im,iday,icl) = nan;
                        fitsampsizen(im,iday,icl) = 0;
                    end
        %             daymaxn{icl}(iday) = clouds.(campaign)(iday).maxN;
        
% %                 else
% %                     slp(im,iday,icl) = nan;
% %                     intcpt(im,iday,icl) = 0;
% %                     slpn(im,iday,icl) = 0;
% %                     intcptn(im,iday,icl) = 0;
% %                     fitsampsize(im,iday,icl) = 0;
% %                     fitsampsizen(im,iday,icl) = 0;
%                 end 
%             end
%             if do_plot == 1
%                 axes(a(icl))
%                 plot(Moments,nanmean(slp(im,:,icl)))
%                 
%                 axes(a(icl+nlayers))
%                 plot(Moments,nanmean(slpn(im,:,icl)))
                end
            end
        end   
        
        for icl = 1:nlayers
            slp(slp==0) = nan;
            slpn(slpn==0) = nan;

            mslp_day(im,icl)=nanmean(slp(im,:,icl));
            mslpn_day(im,icl)=nanmean(slpn(im,:,icl));

            intcpt(intcpt==0) = nan;
            intcptn(intcptn==0) = nan;

            mintcpt_day(im,icl) = nanmean(intcpt(im,:,icl));
            mintcptn_day(im,icl) = nanmean(intcptn(im,:,icl));
            
            meand(meand==0) = nan;
            meand_day(icl) = nanmean(meand(:,icl));

            for iday = days_analyzed

                stdcam_slp(im,iday,icl)=nanstd(slp(im,:,icl)/mslp_day(im,...
                    icl),fitsampsize(im,:,icl));
                
                stdcam_slpn(im,iday,icl)=nanstd(slpn(im,:,icl)/mslpn_day(im,...
                    icl),fitsampsizen(im,:,icl));

                stdcam_intcpt(im,iday,icl)=nanstd(intcpt(im,:,...
                    icl)/mintcpt_day(im,icl),fitsampsize(im,:,icl));
                
                stdcam_intcptn(im,iday,icl)=nanstd(intcptn(im,:,...
                    icl)/mintcptn_day(im,icl),fitsampsizen(im,:,icl));
            end
            
            mstdcam_slp_day(im,icl)=mean(stdcam_slp(im,:,icl));
            mstdcam_slpn_day(im,icl)=mean(stdcam_slpn(im,:,icl));
            mstdcam_intcpt_day(im,icl)=mean(stdcam_intcpt(im,:,icl));
            mstdcam_intcptn_day(im,icl)=mean(stdcam_intcptn(im,:,icl));
        end
    end

    if do_plot==1 && do_plot_corr==1
        f = figure('position',[183 374 1750 650]);
        for i=1:2*nlayers
            a(i)=subplot(2,nlayers,i);
            hold on
            set(gca,'linestyleorder','-')
        end
        set(gca,'linestyleorder','-')
        for icl = 1:nlayers
            axes(a(icl+nlayers))

            pos = get(a, 'position');
            pos{icl}(1) = pos{icl}(1)+1/15;
            str = ['$\overline {D}$ = ', num2str(meand_day(icl)),' $\mu$m'];
            annotation(f,'textbox', pos{icl}, 'String',str,...
                'FitBoxToText','on','Interpreter','latex','FontSize',20);
    %         set(ant,'Interpreter','latex')

            ylim([max([min(mstdcam_slp_day(:)) 0]) min([max(mstdcam_slp_day(:)) 1]) ])
            xlabel('k in M_k (moment of the distribution)')
            plot(Moments,mstdcam_slp_day(:,icl),'LineWidth',2)
            ylabel('Relative \sigma of slopes of \sigma(k)')

            yyaxis right
            ylim([max([min(mstdcam_intcpt_day(:)) 0]) min([max(mstdcam_intcpt_day(:)) 1]) ])
            plot(Moments,mstdcam_intcpt_day(:,icl),'LineWidth',2)
            ylabel('Relative \sigma of intercepts of \sigma(k)')
            set(gca,'FontSize',16)



            axes(a(icl))
            cptitle = [num2str(layer_edges(icl)),' - ',...
                            num2str(layer_edges(icl+1))];
            title(cptitle)
            ylim([max([min(mstdcam_slpn_day(:)) 0]) min([max(mstdcam_slpn_day(:)) 1]) ])
            xlabel('k in M_k (moment of the distribution)')
            plot(Moments,mstdcam_slpn_day(:,icl),'LineWidth',2)
            ylabel('Relative \sigma of slopes of \epsilon(k)')

            yyaxis right
            ylim([max([min(mstdcam_intcptn_day(:)) 0]) min([max(mstdcam_intcptn_day(:)) 1]) ])
            plot(Moments,mstdcam_intcptn_day(:,icl),'LineWidth',2)
            ylabel('Relative \sigma of intercepts of \epsilon(k)')

            set(gca,'FontSize',16)

        end

        if do_save==1
        saveas(gcf,['plots/momentsofDSD/',campaign, '_std_mx, # of moments=',...
            num2str(length(Moments)),', layers=',...
            num2str(nlayers),', z_l=', num2str(z_l), ', z_h=', num2str(z_h),...
            ', thresN=', num2str(thresN), ', thresPTS=', num2str(thresPTS),...
            ', thresRH=', num2str(thresRH),', thresT=', num2str(thresT), '.png'])
        end
    end
    
    
    
    
    if do_plot == 1 && do_plot_slope==1
        f = figure('position',[183 374 1850 650]);
        for i=1:2*nlayers
            a(i)=subplot(2,nlayers,i);
            hold on
            set(gca,'linestyleorder','-')
        end
        set(gca,'linestyleorder','-')
        
        for icl = 1:nlayers
            axes(a(icl+nlayers))

            pos = get(a, 'position');
            pos{icl}(1) = pos{icl}(1)+1/15;
            str = ['$\overline {D}$ = ', num2str(meand_day(icl)),' $\mu$m'];
            annotation(f,'textbox', pos{icl}, 'String',str,...
                'FitBoxToText','on','Interpreter','latex','FontSize',20);
    %         set(ant,'Interpreter','latex')

            xlabel('k in M_k (moment of the distribution)')
            ylim([min(mslp_day,[],[1 2]) max(mslp_day,[],[1 2])])
            semilogy(Moments(mslp_day(:,icl)<0),mslp_day(mslp_day(:,icl)<0,icl),'LineWidth',2)
            ylabel('Campaign-averaged Slopes of \sigma(k)')
            set(gca,'YScale','log')

            yyaxis right
            ylim([min(mintcpt_day,[],[1 2]) max(mintcpt_day,[],[1 2])])
            plot(Moments,mintcpt_day(:,icl),'LineWidth',2)
            ylabel('Campaign-averaged Intercepts of \sigma(k)')
            set(gca,'FontSize',16)



            axes(a(icl))
            cptitle = [num2str(layer_edges(icl)),' - ',...
                            num2str(layer_edges(icl+1))];
            title(cptitle)
            ylim([min(mslpn_day,[],[1 2]) max(mslpn_day,[],[1 2])])
            xlabel('k in M_k (moment of the distribution)')
            plot(Moments,mslpn_day(:,icl),'LineWidth',2)
            ylabel('Campaign-averaged Slopes of \epsilon(k)')
            set(gca,'YScale','log')

            yyaxis right
            ylim([min(mintcptn_day,[],[1 2]) max(mintcptn_day,[],[1 2])])
            plot(Moments,mintcptn_day(:,icl),'LineWidth',2)
            ylabel('Campaign-averaged Intercepts of \epsilon(k)')

            set(gca,'FontSize',16)
        end

        if do_save==1
        saveas(gcf,['plots/momentsofDSD/',campaign, '_slp_int_std_mx, # of moments=',...
            num2str(length(Moments)),', layers=',...
            num2str(nlayers),', z_l=', num2str(z_l), ', z_h=', num2str(z_h),...
            ', thresN=', num2str(thresN), ', thresPTS=', num2str(thresPTS),...
            ', thresRH=', num2str(thresRH),', thresT=', num2str(thresT), '.png'])
        end
    end
%         figure('Position',[800 600 800 400])
%         subplot(1,2,1)
%         bar(layer_mean,stdcam_slp)
%         ylabel('Standard Deviation/Mean of Fitted Slopes')
%         xlim([layer_edges(1) layer_edges(end)])
%         ylim([0 2])
%         legend('N','N/Nmax','Location','northwest')
%         % set(gca,'xticklabel',cptitle)
%         set(gca,'FontSize',18)
% 
%         subplot(1,2,2)
%         bar(layer_mean-0.05,mslp,0.4)
%         ylabel('Mean of Fitted Slopes')
%         yyaxis right
%         bar(layer_mean+0.05,mslpn,0.4);
%         ylabel('Mean of Fitted Normalized Slopes')
% 
%         hold on
%         yyaxis left
%         errorbar(layer_mean-0.05,mslp,stdcam_slp(:,1).*mslp','.',...
%             'color','red','LineWidth',2)
%         ax1_ylim = ylim;
% 
%         yyaxis right
%         errorbar(layer_mean+0.05,mslpn,stdcam_slp(:,2).*mslpn','.',...
%             'color','blue','LineWidth',2)
%         ax2_ylim = ylim;
% 
%         lim_ratio = ax2_ylim./ax1_ylim;
%         which_lim = find(lim_ratio~=inf & ~isnan(lim_ratio) & lim_ratio~=0);
%         lim_ratio = lim_ratio(which_lim);
% 
% 
%         if ismember(0,ax2_ylim) && ~ismember(0,ax1_ylim)
%             ax2_ylim(1+(2-which_lim)) = ax1_ylim(1+(2-which_lim))*lim_ratio;
%         elseif ismember(0,ax1_ylim) && ~ismember(0,ax2_ylim)
%             ax1_ylim(1+(2-which_lim)) = ax2_ylim(1+(2-which_lim))/lim_ratio;
%         else
%             [lim_ratio,which_lim] = max(lim_ratio);
%             ax2_ylim(1+(2-which_lim)) = ax1_ylim(1+(2-which_lim))*lim_ratio;
%         %     ax1_ylim(1+(2-which_lim)) = ax2_ylim(1+(2-which_lim))/lim_ratio;
%         end
%         % ax1_ylim = [ax1_ylim(1)*lim_ratio ax1_ylim(2)];
% 
%         yyaxis left
%         ylim(ax1_ylim)
%         yyaxis right
%         ylim(ax2_ylim)
%         xlim([layer_edges(1) layer_edges(end)])
%         xlabel('Cloud Layer')
% 
%     %     ylim([-0.6 0.1])
%         % set(gca,'xticklabel',cptitle)
%         set(gca,'FontSize',18)
%         hold off

    %     subplot(1,3,3)
    %     bar(layer_mean,mslpn)
    %     hold on
    %     errorbar(layer_mean,mslpn,stdcam_slp(:,2).*mslpn','.','LineWidth',2)
    %     xlim([layer_edges(1) layer_edges(end)])
    %     xlabel('Cloud Layer')
    %     ylabel('Mean of Fitted Slopes')
    %     xlim([layer_edges(1) layer_edges(end)])
    % %     ylim([-0.6 0.1])
    %     % set(gca,'xticklabel',cptitle)
    %     set(gca,'FontSize',18)
    %     hold off

        % ylim([0 0.6])
%         if do_save==1
%             saveas(gcf,['plots/slope/test/',campaign, '_reldisp_m',...
%                 num2str(Moments(im)),', layers=',num2str(nlayers), ', z_l=',...
%                 num2str(z_l), ', z_h=', num2str(z_h),', thresN=',...
%                 num2str(thresN), ', thresPTS=', num2str(thresPTS), ...
%                 ', thresRH=', num2str(thresRH), ', thresT=', ...
%                 num2str(thresT), '.png'])
%         end
%     end
% end
%%


% figure
% bar(stdcam_intcpt)
% ylabel('Standard Deviation/Mean of Fitted Intercepts')
% legend('N','N/Nmax','Location','best')
% set(gca,'xticklabel',cptitle)
% set(gca,'FontSize',18)
% saveas(gcf,'VOCALS_cloudtmb_intercept_wgtd_std.png')