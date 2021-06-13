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

AF_edges = [0 0.25 0.5 0.9];
nlayers = length(AF_edges)-1;
AF_labels = {sprintf('%0.2f<AF<%0.2f', [AF_edges(1) AF_edges(2)]),...
    sprintf('%0.2f<AF<%0.2f', [AF_edges(2) AF_edges(3)]),...
    sprintf('%0.2f<AF<%0.2f', [AF_edges(3) AF_edges(4)])};
AF_mean = (AF_edges(1:end-1) + AF_edges(2:end))/2;

Moments = 0:0.25:2;

do_plot = true; % whether to plot the slopes figure
do_plot_corr = true;
do_plot_slope = false;
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
clear bineps binepsn fitsampsizen mstdcam_slp_day mstdcam_slpn_day
clear mstdcam_intcpt_day mstdcam_intcptn_day mslp_day mslpn_day
clear mintcpt_day mintcptn_day

c = 1;

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
    afvar='AF';
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
        
        for iaf = 1:nlayers % iterator for different each layer of the cloud
            
            for iday = days_analyzed
            %     c = 1; 
                binlims = clouds.(campaign)(iday).binlims;
                binmean = clouds.(campaign)(iday).binmean;
                dlogD = log10(binlims(2)/binlims(1));


                t = clouds.(campaign)(iday).s_t;
                AC = clouds.(campaign)(iday).(ACvar);
                AF = clouds.(campaign)(iday).(afvar);
                t_cl = t(AF>AF_edges(iaf) & AF<AF_edges(iaf+1));
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

                ncldpts(iday,iaf) = length(cldpts);

                %If there are enough data points, give data a more convenient name,
                %and collect data in cells (the cell data isn't actually used).

                if ncldpts(iday,iaf)>=thresPTS
                    alld=clouds.(campaign)(iday).(Dvar)(cldpts,:);
                    allnd=clouds.(campaign)(iday).(nvar)(cldpts,:);
                    allmx=sum(dlogD*allnd.*binmean'.^Moments(im),2);

                    meand(iday,iaf) = mean(alld);
    %                 maxn(iday,icl)=max(alln);
    %                 allnormN=alln/maxn(iday,icl);
    % 
    %                 maxh(iday,icl)=max(allh);
    %                 minh(iday,icl)=min(allh);
    %                 allnormH=(allh-minh(iday,icl))/(maxh(iday,icl)-minh(iday,icl));
    % 
    %                 maxRH(iday,icl)=max(allRH);
    %                 minRH(iday,icl)=min(allRH);
    %                 allnormRH=(allRH-minRH(iday,icl))/(maxRH(iday,icl)-minRH(iday,icl));
    % 
    %                 maxsupsat(iday,icl)=max(allsupsat);
    %                 allnormsupsat=allsupsat/maxsupsat(iday,icl);
    % 
                    alleps=clouds.(campaign)(iday).(epsvar)(cldpts);
                    allstd=clouds.(campaign)(iday).(stdvar)(cldpts);
    % 
    %                 maxd(iday,icl)=max(alld);
    %                 mind(iday,icl)=min(alld)-0.01;
    %                 allnormD=(alld-mind(iday,icl))/(maxd(iday,icl)-mind(iday,icl));
    %                 
    %                 maxl(iday,icl)=max(alll);
    %                 allnormL=alll/maxl(iday,icl);
    %                 
    %                 maxthet(iday,icl)=max(allthet);
    %                 minthet(iday,icl)=min(allthet);
    %                 allnormThet=(allthet-minthet(iday,icl))/(maxthet(iday,...
    %                     icl)-minthet(iday,icl));
    %                 
    %                 maxreff(iday,icl)=max(allreff);
    %                 minreff(iday,icl)=min(allreff)-0.01;
    %                 allnormReff=(allreff-minreff(iday,...
    %                     icl))/(maxreff(iday,icl)-minreff(iday,icl));


    %                 ALLN{icl}(end+1:end+ncldpts(iday,icl))=alln;
    %                 ALLEPS{icl}(end+1:end+ncldpts(iday,icl))=alleps;
    %                 ALLSTD{icl}(end+1:end+ncldpts(iday,icl))=allstd;
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
                else
                    ncldpts(iday,iaf) = 0;
                end        

                if ~isempty(alleps) %If data exists
    %             try
                    nbins=20; %Number of cloud droplet concentration bins
                    indvar = allmx; %; allm2./alln %independent variable
    %                 depvar_fit = sqrt((indvar-1)/2);
    %                 
    %                 indvar = abs(depvar_fit-alleps); % let the independent
    %                 % variable be the deviation from gamma-like PDF
                    indnormVar = indvar/max(indvar); %allnormReff./allnormD;
                    %try
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
                        slp(im,iday,iaf) = fitline.p1;
                        intcpt(im,iday,iaf) = fitline.p2;
                        fitsampsize(im,iday,iaf) = sum(inds);
                        fitsampsize(fitsampsize(im,iday,iaf)<5)=0;
                    else
                        slp(im,iday,iaf) = nan;
                        intcpt(im,iday,iaf) = nan;
                        fitsampsize(im,iday,iaf) = 0;
                    end

                    %Do some binning now for normalized droplet concentration
                    [N,edges,bin]=histcounts(indvar,nbins);
                    binstd = accumarray(bin,allstd)./N';
                    binstd(N<5)=NaN;
                    meanvalsn{iaf}(:,iday)=binstd;
                    centersn=(edges(1:end-1)+edges(2:end))/2;

                    %Also fit a line to the data
                    inds=~isnan(binstd);
                    x=centersn(inds);
                    y=binstd(inds);

                    if length(x)>1
                        fitline = fit(x',y,'poly1');
                        slpn(im,iday,iaf) = fitline.p1;
                        intcptn(im,iday,iaf) = fitline.p2;
                        fitsampsizen(im,iday,iaf) = sum(inds);
                        fitsampsizen(fitsampsize(im,iday,iaf)<5)=0;
                    else
                        slpn(im,iday,iaf) = nan;
                        intcptn(im,iday,iaf) = nan;
                        fitsampsizen(im,iday,iaf) = 0;
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
        
        for iaf = 1:nlayers
            slp(slp==0) = nan;
            slpn(slpn==0) = nan;

            mslp_day(im,iaf)=nanmean(slp(im,:,iaf));
            mslpn_day(im,iaf)=nanmean(slpn(im,:,iaf));

            intcpt(intcpt==0) = nan;
            intcptn(intcptn==0) = nan;

            mintcpt_day(im,iaf) = nanmean(intcpt(im,:,iaf));
            mintcptn_day(im,iaf) = nanmean(intcptn(im,:,iaf));
            
            meand(meand==0) = nan;
            meand_day(iaf) = nanmean(meand(:,iaf));

            for iday = days_analyzed

                stdcam_slp(im,iday,iaf)=nanstd(slp(im,:,iaf)/mslp_day(im,...
                    iaf),fitsampsize(im,:,iaf));
                
                stdcam_slpn(im,iday,iaf)=nanstd(slpn(im,:,iaf)/mslpn_day(im,...
                    iaf),fitsampsizen(im,:,iaf));

                stdcam_intcpt(im,iday,iaf)=nanstd(intcpt(im,:,...
                    iaf)/mintcpt_day(im,iaf),fitsampsize(im,:,iaf));
                
                stdcam_intcptn(im,iday,iaf)=nanstd(intcptn(im,:,...
                    iaf)/mintcptn_day(im,iaf),fitsampsizen(im,:,iaf));
            end
            
            mstdcam_slp_day(im,iaf)=mean(stdcam_slp(im,:,iaf));
            mstdcam_slpn_day(im,iaf)=mean(stdcam_slpn(im,:,iaf));
            mstdcam_intcpt_day(im,iaf)=mean(stdcam_intcpt(im,:,iaf));
            mstdcam_intcptn_day(im,iaf)=mean(stdcam_intcptn(im,:,iaf));
        end
    end

    if do_plot == 1 && do_plot_corr==1
        f = figure('position',[183 374 1800 640]);
        for i=1:2*nlayers
            a(i)=subplot(2,nlayers,i);
            hold on
            set(gca,'linestyleorder','-')
        end
        set(gca,'linestyleorder','-')
        for iaf = 1:nlayers
            axes(a(iaf))

            pos = get(a, 'position');
            pos{iaf}(1) = pos{iaf}(1)+1/12;
            str = ['$\overline {D}$ = ', num2str(meand_day(iaf)),' $\mu$m'];
            annotation(f,'textbox', pos{iaf}, 'String',str,...
                'FitBoxToText','on','Interpreter','latex');
    %         set(ant,'Interpreter','latex')

            cptitle = AF_labels{iaf};
            title(cptitle)
            ylim([max([min(mstdcam_slp_day(:)) 0]) min([max(mstdcam_slp_day(:)) 1]) ])
            xlabel('k in M_k (moment of the distribution)')
            plot(Moments,mstdcam_slp_day(:,iaf),'LineWidth',2)
            ylabel('Relative std of slopes')
            yyaxis right
            ylim([max([min(mstdcam_intcpt_day(:)) 0]) min([max(mstdcam_intcpt_day(:)) 1]) ])
            plot(Moments,mstdcam_intcpt_day(:,iaf),'LineWidth',2)
            ylabel('Relative std of intercepts')
            set(gca,'FontSize',16)

            axes(a(iaf+nlayers))
            ylim([max([min(mstdcam_slpn_day(:)) 0]) min([max(mstdcam_slpn_day(:)) 2]) ])
            xlabel('k in M_k (moment of the distribution)')
            plot(Moments,mstdcam_slpn_day(:,iaf),'LineWidth',2)
            ylabel('Std of slopes')
            yyaxis right
            ylim([max([min(mstdcam_intcptn_day(:)) 0]) min([max(mstdcam_intcptn_day(:)) 2]) ])
            plot(Moments,mstdcam_intcptn_day(:,iaf),'LineWidth',2)
            ylabel('Std of intercepts')

            set(gca,'FontSize',16)
        end

        if do_save==1
            saveas(gcf,['plots/slope/test/',campaign, '_std_mx, # of moments=',...
                num2str(length(Moments)),', layers=',...
                num2str(nlayers),', z_l=', num2str(z_l), ', z_h=', num2str(z_h),...
                ', thresN=', num2str(thresN), ', thresPTS=', num2str(thresPTS),...
                ', thresRH=', num2str(thresRH),', thresT=', num2str(thresT), '.png'])
        end
    end
    
    
    
    
    if do_plot == 1 && do_plot_slope==1
        f = figure('position',[183 374 1800 640]);
        for i=1:2*nlayers
            a(i)=subplot(2,nlayers,i);
            hold on
            set(gca,'linestyleorder','-')
        end
        set(gca,'linestyleorder','-')
        for iaf = 1:nlayers
            axes(a(iaf))

            pos = get(a, 'position');
            pos{iaf}(1) = pos{iaf}(1)+1/12;
            str = ['$\overline {D}$ = ', num2str(meand_day(iaf)),' $\mu$m'];
            annotation(f,'textbox', pos{iaf}, 'String',str,...
                'FitBoxToText','on','Interpreter','latex');
    %         set(ant,'Interpreter','latex')

            cptitle = AF_labels{iaf};
            title(cptitle)

            xlabel('k in M_k (moment of the distribution)')
            ylim([min(mslp_day(:)) max(mslp_day(:))])
            semilogy(Moments,mslp_day(:,iaf),'LineWidth',2)
            ylabel('Slopes')
            set(gca,'YScale','log')
            yyaxis right
            ylim([min(mintcpt_day(:)) max(mintcpt_day(:))])
            semilogy(Moments,mintcpt_day(:,iaf),'LineWidth',2)
            ylabel('Intercepts')
            set(gca,'FontSize',16)

            axes(a(iaf+nlayers))
            ylim([min(mslpn_day(:)) max(mslpn_day(:))])
            xlabel('k in M_k (moment of the distribution)')
            plot(Moments,mslpn_day(:,iaf),'LineWidth',2)
            ylabel('Campaign-averaged Slopes of \sigma(k)')
            set(gca,'YScale','log', 'FontSize',16)
            yyaxis right
            ylim([min(mintcptn_day(:)) max(mintcptn_day(:))])
            plot(Moments,mintcptn_day(:,iaf),'LineWidth',2)
            ylabel('Campaign-averaged Intercepts of \sigma(k)')
%             set(gca,'YScale','log', 'FontSize',16)
        end

        if do_save==1
            saveas(gcf,['plots/slope/test/',campaign, '_slp_int_std_mx, # of moments=',...
                num2str(length(Moments)),', layers=',...
                num2str(nlayers),', z_l=', num2str(z_l), ', z_h=', num2str(z_h),...
                ', thresN=', num2str(thresN), ', thresPTS=', num2str(thresPTS),...
                ', thresRH=', num2str(thresRH),', thresT=', num2str(thresT), '.png'])
        end
    end
% end