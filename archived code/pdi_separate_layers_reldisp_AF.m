clear
cd '~/Box/grad/research/datasets/PDI Data/'
load obs_clouds_wholeday.mat

%%
clear meanvals meanvalsn
camps={'vocalspdi','masepdi','postpdi'};
ctitle={'VOCALS','MASE','POST'};
instr={'pdi','pdi','pdi'};


thresN=50;
thresRH=70;
thresT=5; %minimum temperature, just to make sure we avoid mixed-phase
thresLWC=0.01;
thresPTS=100; %minimum number of data samples
doresamp = false;

% z_l = -0.1;
% z_h = 1.1;
% layer_edges = linspace(z_l,z_h,nlayers+1); %[-.05, .3, .7, 1.05]; %
% layer_mean=(layer_edges(1:end-1)+layer_edges(2:end))/2; % layer center diameter
AF_edges = [0 0.25 0.5 0.9];
nlayers = length(AF_edges)-1;
AF_labels = {sprintf('%0.2f<AF<%0.2f', [AF_edges(1) AF_edges(2)]),...
    sprintf('%0.2f<AF<%0.2f', [AF_edges(2) AF_edges(3)]),...
    sprintf('%0.2f<AF<%0.2f', [AF_edges(3) AF_edges(4)])};
AF_mean = (AF_edges(1:end-1) + AF_edges(2:end))/2;

Moments = 0:0.25:3;

do_plot = true; % whether to plot the slopes figure
do_plot_corr = true;
do_plot_stdslp = true;
do_save = false;

%Initialize cells - but not used.
ALLN=cell(1,nlayers);
ALLSTD=cell(1,nlayers);
ALLH=cell(1,nlayers);
ALLEPS=cell(1,nlayers);
ALLD=cell(1,nlayers);
ALLL=cell(1,nlayers);
ALLThet=cell(1,nlayers);
ALLReff=cell(1,nlayers);
% ALLM1=cell(1,nlayers);
% ALLM2=cell(1,nlayers);
% ALLM3=cell(1,nlayers);
% ALLM6=cell(1,nlayers);
ALLPCP=cell(1,nlayers);



%% separate out the cloud top mid bottom

clear slp slpn intcpt intcptn stdcam_intcpt stdcam_slp fitsampsize
clear bineps binepsn fitsampsizen
c=3;

close all

% for c=1:length(camps)
    
    if do_plot == 1 && do_plot_corr == 1
        figure('position',[183 374 1328 667])
        for i=1:2*nlayers
            a(i)=subplot(2,nlayers,i);
            hold on
            set(gca,'linestyleorder','-')
        end
        set(gca,'linestyleorder','-')
    end
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
    afvar='AF';
    ndays = nc;

    days_analyzed = 1:ndays;
    if c==1
        days_analyzed(days_analyzed==10 | days_analyzed==13 | days_analyzed==11)=[]; % remove the 
    % days that have incomplete flights in vocals
    end

    for iday = days_analyzed
        binlims = clouds.(campaign)(iday).binlims;
        binmean = clouds.(campaign)(iday).binmean;
        dlogD = log10(binlims(2)/binlims(1));
        
        for iaf = 1:length(AF_labels) % iterator for different each layer of the cloud
        %     c = 1; 
        
            binlims = clouds.(campaign)(iday).binlims;
            binmean = clouds.(campaign)(iday).binmean;
            dlogD = log10(binlims(2)/binlims(1));
            t = clouds.(campaign)(iday).s_t;
            AF = clouds.(campaign)(iday).(afvar);
            t_cl = t(AF>AF_edges(iaf) & AF<AF_edges(iaf+1));
            t_cl_idx = ismember(t,t_cl);

            [alleps,allstd_log,alln,allnormN,maxn,alld,maxd,mind,allnormD,...
                allh,maxh,minh,allRH,maxRH,minRH,allnormRH,allsupsat,...
                maxsupsat,alll,allnormL,maxl,allthet,allnormThet,maxthet,...
                minthet,allreff,allnormReff,maxreff,minreff,allaf]=deal([]);
            %Find all data that meet the thresholds
            %Remember clouds contains all data with non-NaN relative dispersion
            cldpts_all_cp=find(clouds.(campaign)(iday).(Vvar)>thresLWC & ...
                clouds.(campaign)(iday).(Nvar)>thresN & ...
                clouds.(campaign)(iday).(Tvar)>thresT & ...
                clouds.(campaign)(iday).(RHvar)>thresRH);

            % find the cloud droplets that are both valid by the standard above
            % and in the specific cloud layer
            cldpts = intersect(cldpts_all_cp, find(t_cl_idx));

            ncldpts(iday,iaf) = length(cldpts);

            %If there are enough data points, give data a more convenient name,
            %and collect data in cells (the cell data isn't actually used).
            if ncldpts(iday,iaf)>=thresPTS
                alld=clouds.(campaign)(iday).(Dvar)(cldpts);
                alln=clouds.(campaign)(iday).(Nvar)(cldpts);
                allnd=clouds.(campaign)(iday).(nvar)(cldpts,:);
                allh=clouds.(campaign)(iday).(Hvar)(cldpts);
                allRH = clouds.(campaign)(iday).(RHvar)(cldpts);
                allsupsat = clouds.(campaign)(iday).(RHvar)(cldpts)-100;
                alll=clouds.(campaign)(iday).(Lvar)(cldpts);
                allthet=clouds.(campaign)(iday).(Thetvar)(cldpts);
                allreff=clouds.(campaign)(iday).(Reffvar)(cldpts);
                allaf = clouds.(campaign)(iday).(afvar)(cldpts);
                
                maxn(iday,iaf)=max(clouds.(campaign)(iday).(Nvar));
                allnormN=alln/maxn(iday,iaf);

                maxh(iday,iaf)=max(allh);
                minh(iday,iaf)=min(allh);
                allnormH=(allh-minh(iday,iaf))/(maxh(iday,iaf)-minh(iday,iaf));

                maxRH(iday,iaf)=max(allRH);
                minRH(iday,iaf)=min(allRH);
                allnormRH=(allRH-minRH(iday,iaf))/(maxRH(iday,iaf)-minRH(iday,iaf));

                maxsupsat(iday,iaf)=max(allsupsat);
                allnormsupsat=allsupsat/maxsupsat(iday,iaf);

                alleps=(clouds.(campaign)(iday).(epsvar)(cldpts));
%                 allstd_log=(clouds.(campaign)(iday).(stdvar)(cldpts));

                maxd(iday,iaf)=max(alld);
                mind(iday,iaf)=min(alld)-0.01;
                allnormD=(alld-mind(iday,iaf))/(maxd(iday,iaf)-mind(iday,iaf));
                
                maxl(iday,iaf)=max(clouds.(campaign)(iday).(Lvar));
                allnormL=alll/maxl(iday,iaf);
                
                maxthet(iday,iaf)=max(allthet);
                minthet(iday,iaf)=min(allthet);
                allnormThet=(allthet-minthet(iday,iaf))/(maxthet(iday,...
                    iaf)-minthet(iday,iaf));
                
                maxreff(iday,iaf)=max(allreff);
                minreff(iday,iaf)=min(allreff)-0.01;
                allnormReff=(allreff-minreff(iday,...
                    iaf))/(maxreff(iday,iaf)-minreff(iday,iaf));
                
                allmx=sum(dlogD*allnd.*binmean'.^2,2);

%                 ALLN{icl}(end+1:end+ncldpts(iday,icl))=alln;
%                 ALLEPS{icl}(end+1:end+ncldpts(iday,icl))=alleps_log;
%                 ALLSTD{icl}(end+1:end+ncldpts(iday,icl))=allstd_log;
%                 ALLH{icl}(end+1:end+ncldpts(iday,icl))=allh;
%                 ALLD{icl}(end+1:end+ncldpts(iday,icl))=alld;
%                 ALLL{icl}(end+1:end+ncldpts(iday,icl))=alll;
%                 ALLThet{icl}(end+1:end+ncldpts(iday,icl))=allthet;
%                 ALLReff{icl}(end+1:end+ncldpts(iday,icl))=allreff;
%                 ALLM1{icl}(end+1:end+ncldpts(iday,icl))=allm1;
%                 ALLM2{icl}(end+1:end+ncldpts(iday,icl))=allm2;
%                 ALLM3{icl}(end+1:end+ncldpts(iday,icl))=allm3;
%                 ALLM6{icl}(end+1:end+ncldpts(iday,icl))=allm6;
                %ALLPCP{c}(end+1:end+ncldpts)=clouds.(campaign)(iday).(Pvar)(cldpts);
                indvar = allmx; % sqrt(alln.*alll);%log(sum(dlogD*allnd.*binmean'.^Moments(im),2)); % 2*allreff./alld; %; allm2./alln %independent variable
%                 depvar_fit = sqrt((indvar-1)/2);
                indvar = indvar(~isnan(indvar));
                alleps = alleps(~isnan(indvar));
        %                 indvar = abs(depvar_fit-alleps); % let the independent
        %                 % variable be the deviation from gamma-like PDF
                indnormVar = indvar/max(indvar); % sqrt(allnormN.*allnormL); %allnormReff./allnormD;
            else
                ncldpts(iday,iaf) = 0;
            end   
            %try
            
            if ~isempty(alleps) %If data exists
%             try
                nbins=20; %Number of cloud droplet concentration bins
                
                
                %Do some binning for cloud local data
                
                [N,edges,bin]=histcounts(indvar,nbins);
                %Find average relative dispersion in each bin
                bineps = accumarray(bin,alleps)./N';
                
                %Set to NaN if there are fewer than 5 points in the bin
                bineps(N<5)=NaN;
                meanvals{iaf}(:,iday)=bineps;
                centers=(edges(1:end-1)+edges(2:end))/2;
                ref_bineps = sqrt((centers-1)/2);

                %Also fit a line to the data
                inds=~isnan(bineps);
                x=centers(inds);
                y=bineps(inds);
                
                
                if length(x)>1
                    fitline = fit(x',y,'poly1');
                    slp(iday,iaf) = fitline.p1;
                    intcpt(iday,iaf) = fitline.p2;
                    fitsampsize(iday,iaf) = sum(inds);
                    fitsampsize(fitsampsize(iday,iaf)<5)=0;
                else
                    slp(iday,iaf) = 0;
                    intcpt(iday,iaf) = 0;
                    fitsampsize(iday,iaf) = 0;
                end

                %Do some binning now for normalized droplet concentration
                [N,edges,bin]=histcounts(indnormVar,nbins);
                binepsn = accumarray(bin,alleps)./N';
                binepsn(N<5)=NaN;
                meanvalsn{iaf}(:,iday)=binepsn;
                centersn=(edges(1:end-1)+edges(2:end))/2;

                %Also fit a line to the data
                inds=~isnan(binepsn);
                x=centersn(inds);
                y=binepsn(inds);
                
                if length(x)>1
                    fitline = fit(x',y,'poly1');
                    slpn(iday,iaf) = fitline.p1;
                    intcptn(iday,iaf) = fitline.p2;
                    fitsampsizen(iday,iaf) = sum(inds);
                    fitsampsizen(fitsampsize(iday,iaf)<5)=0;
                else
                    slpn(iday,iaf) = 0;
                    intcptn(iday,iaf) = 0;
                    fitsampsizen(iday,iaf) = 0;
                end
    %             daymaxn{icl}(iday) = clouds.(campaign)(iday).maxN;


                if do_plot == 1 && do_plot_corr == 1
                    colors=colormap(parula(ndays));
                    axes(a(iaf))
                    cptitle = AF_labels{iaf};
                    title(cptitle)
                    set(gca,'colororder',colors)
                    plot(centers,bineps,'LineWidth',1.5)
%                     hold on
%                     plot(centers,ref_bineps,'color',[0.3 0.3 0.3])
                    xlabel('M_k')
                    ylabel('Relative Dispersion of the DSD')
        %             legend(num2str(iday))
                    set(gca,'ylim',[0.1 0.6])
                    set(gca,'FontSize',16)

                    axes(a(iaf+nlayers))
                    set(gca,'colororder',colors)
                    plot(centersn,binepsn,'LineWidth',1.5)
                    xlabel('M_k')
                    ylabel('Standard deviation of the DSD')
        %             legend(num2str(iday))
                    set(gca,'xlim',[0 1],'ylim',[0.1 0.6])
                    set(gca,'FontSize',16)
                end
            else
                slp(iday,iaf) = 0;
                intcpt(iday,iaf) = 0;
                slpn(iday,iaf) = 0;
                intcptn(iday,iaf) = 0;
                fitsampsize(iday,iaf) = 0;
                fitsampsizen(iday,iaf) = 0;
            end 
        end
    end
    
    if do_save == 1 && do_plot_corr == 1
        saveas(gcf,['plots/corr/',campaign, '_reldisp_Ntot_corr, layers=',...
            num2str(nlayers), ', z_l=', num2str(z_l), ', z_h=',...
            num2str(z_h), ', thresN=', num2str(thresN), ...
            ', thresPTS=', num2str(thresPTS), ', thresRH=', ...
            num2str(thresRH), ', thresT=', num2str(thresT), '.png'])
    end
    
    slp(slp==0) = nan;
    slpn(slpn==0) = nan;

    mslp=nanmean(slp);
    mslpn=nanmean(slpn);

    intcpt(intcpt==0) = nan;
    intcptn(intcptn==0) = nan;

    mintcpt = nanmean(intcpt);
    mintcptn = nanmean(intcptn);

    for iaf=1:nlayers
        stdcam_slp(iaf,1)=nanstd(slp(:,iaf)/mslp(iaf),fitsampsize(:,iaf));
        stdcam_slp(iaf,2)=nanstd(slpn(:,iaf)/mslpn(iaf),fitsampsizen(:,iaf));

        stdcam_intcpt(iaf,1)=nanstd(intcpt(:,iaf)/mintcpt(iaf),fitsampsize(:,iaf));
        stdcam_intcpt(iaf,2)=nanstd(intcptn(:,iaf)/mintcptn(iaf),fitsampsizen(:,iaf));
    end
    
    if do_plot == 1 && do_plot_stdslp == 1
        figure('Position',[800 600 800 400])
        subplot(1,2,1)
        bar(AF_mean,stdcam_slp)
        ylabel('Standard Deviation/Mean of Fitted Slopes')
        xlim([AF_edges(1) AF_edges(end)])
        ylim([0 2])
        legend('N','N/Nmax','Location','northwest')
        % set(gca,'xticklabel',cptitle)
        set(gca,'FontSize',18)

        subplot(1,2,2)
        bar(AF_mean-0.05,mslp,0.4)
        ylabel('Mean of Fitted Slopes')
        yyaxis right
        bar(AF_mean+0.05,mslpn,0.4);
        ylabel('Mean of Fitted Normalized Slopes')

        hold on
        yyaxis left
        errorbar(AF_mean-0.05,mslp,stdcam_slp(:,1).*mslp','.',...
            'color','red','LineWidth',2)
        ax1_ylim = ylim;

        yyaxis right
        errorbar(AF_mean+0.05,mslpn,stdcam_slp(:,2).*mslpn','.',...
            'color','blue','LineWidth',2)
        ax2_ylim = ylim;

        lim_ratio = ax2_ylim./ax1_ylim;
        which_lim = find(lim_ratio~=inf & ~isnan(lim_ratio) & lim_ratio~=0);
        lim_ratio = lim_ratio(which_lim);


        if ismember(0,ax2_ylim) && ~ismember(0,ax1_ylim)
            ax2_ylim(1+(2-which_lim)) = ax1_ylim(1+(2-which_lim))*lim_ratio;
        elseif ismember(0,ax1_ylim) && ~ismember(0,ax2_ylim)
            ax1_ylim(1+(2-which_lim)) = ax2_ylim(1+(2-which_lim))/lim_ratio;
        else
            [lim_ratio,which_lim] = max(lim_ratio);
            ax2_ylim(1+(2-which_lim)) = ax1_ylim(1+(2-which_lim))*lim_ratio;
        %     ax1_ylim(1+(2-which_lim)) = ax2_ylim(1+(2-which_lim))/lim_ratio;
        end
        % ax1_ylim = [ax1_ylim(1)*lim_ratio ax1_ylim(2)];

        yyaxis left
        ylim(ax1_ylim)
        yyaxis right
        ylim(ax2_ylim)
        xlim([AF_edges(1) AF_edges(end)])
        xlabel('Cloud Layer')

    %     ylim([-0.6 0.1])
        % set(gca,'xticklabel',cptitle)
        set(gca,'FontSize',18)
        hold off
    end
    
%     subplot(1,3,3)
%     bar(AF_mean,mslpn)
%     hold on
%     errorbar(AF_mean,mslpn,stdcam_slp(:,2).*mslpn','.','LineWidth',2)
%     xlim([AF_edges(1) AF_edges(end)])
%     xlabel('Cloud Layer')
%     ylabel('Mean of Fitted Slopes')
%     xlim([AF_edges(1) AF_edges(end)])
% %     ylim([-0.6 0.1])
%     % set(gca,'xticklabel',cptitle)
%     set(gca,'FontSize',18)
%     hold off
    
    % ylim([0 0.6])
    if do_save == 1 && do_plot_stdslp == 1
        saveas(gcf,['plots/slope/',campaign, '_reldisp_Ntot_slope, layers=',...
            num2str(nlayers), ', z_l=', num2str(z_l), ', z_h=', num2str(z_h),...
            ', thresN=', num2str(thresN), ', thresPTS=', num2str(thresPTS), ...
            ', thresRH=', num2str(thresRH), ', thresT=', num2str(thresT), '.png'])
    end
% end
%%


% figure
% bar(stdcam_intcpt)
% ylabel('Standard Deviation/Mean of Fitted Intercepts')
% legend('N','N/Nmax','Location','best')
% set(gca,'xticklabel',cptitle)
% set(gca,'FontSize',18)
% saveas(gcf,'VOCALS_cloudtmb_intercept_wgtd_std.png')