clear
cd '~/Box/grad/research/datasets/PDI Data/'
load obs_clouds_wholeday.mat

%%
clear meanvals meanvalsn

campaigns={'vocalspdi','masepdi','postpdi'};
ctitle={'VOCALS','MASE','POST'};
instr={'pdi','pdi','pdi'};


thresN=75;
thresRH=80;
thresT=3; %minimum temperature, just to make sure we avoid mixed-phase
thresLWC=0.01;
thresPTS=100; %minimum number of data samples
doresamp = false;

nlayers = 76;
layer_left_edges = linspace(0,0.75,nlayers);
layer_thickness = 0.3;
layer_edges = [layer_left_edges' layer_left_edges'+layer_thickness];
layer_mean=(layer_edges(1:end-1)+layer_edges(2:end))/2; % layer center diameter

Moments = 0:0.5:8;

do_plot = false; % whether to plot the slopes figure
do_plot_corr = true;
do_plot_slope = false;
do_animate = true;
do_save = true;

%% separate out the cloud top mid bottom

clear slp slpn intcpt intcptn stdcam_intcpt stdcam_slp fitsampsize
clear binstd binstdn fitsampsizen mstdcam_slp_day mstdcam_slpn_day
clear mstdcam_intcpt_day mstdcam_intcptn_day mslp_day mslpn_day
clear mintcpt_day mintcptn_day

% c=2;

close all

for c=1:length(campaigns)
    
    camp=campaigns{c};
    nc = length(clouds.(camp));

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
%     Thetvar='s_thet';
%     Reffvar=['s_reff_',instr{c}];

ndays = nc;

    days_analyzed = 1:ndays;
    
    if c==1
        days_analyzed(days_analyzed==10 | days_analyzed==13 | days_analyzed==11)=[]; 
        % remove the days that have incomplete flights in vocals
    end
    for im = 1:length(Moments)
        disp(im)
        for icl = 1:nlayers % iterator for different each layer of the cloud
            for iday = days_analyzed
                binlims = clouds.(camp)(iday).binlims;
                binmean = clouds.(camp)(iday).binmean;
                dlogD = log10(binlims(2)/binlims(1));

                t = clouds.(camp)(iday).s_t;
                AC = clouds.(camp)(iday).(ACvar);
                t_cl = t(AC>layer_edges(icl,1) & AC<layer_edges(icl,2));
                t_cl_idx = ismember(t,t_cl);

                [alleps,allstd,alln,alld]=deal([]);
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

                %If there are enough data points, give data a more convenient name,
                %and collect data in cells (the cell data isn't actually used).
                if ncldpts(c,iday,icl)>=thresPTS
                    alld=clouds.(camp)(iday).(Dvar)(cldpts,:);
                    allnd=clouds.(camp)(iday).(nvar)(cldpts,:);
                    allmx=sum(dlogD*allnd.*binmean'.^Moments(im),2);

                    meand(c,iday,icl) = mean(alld);
                    allstd=clouds.(camp)(iday).(stdvar)(cldpts);
                    alleps=clouds.(camp)(iday).(epsvar)(cldpts);
                    
                else
                    ncldpts(c,iday,icl) = 0;
                end        

                if ~isempty(allstd) %If data exists
                    
                    nbins=20; %Number of cloud droplet concentration bins
                    indvar = allmx; 
                    
                    %Do some binning for cloud local data
                    [N,edges,bin]=histcounts(indvar,nbins);
                    %Find average relative dispersion in each bin
                    binstd = accumarray(bin,allstd)./N';
                    %Set to NaN if there are fewer than 5 points in the bin
                    binstd(N<5)=NaN;
                    meanvals{c,icl}(:,iday)=binstd;
                    centers=(edges(1:end-1)+edges(2:end))/2;
                    ref_binstd = sqrt((centers-1)/2);

                    %Also fit a line to the data
                    inds=~isnan(binstd);
                    x=centers(inds);
                    y=binstd(inds);


                    if length(x)>1
                        fitline = fit(x',y,'poly1');
                        slp(c,im,iday,icl) = fitline.p1;
                        intcpt(c,im,iday,icl) = fitline.p2;
                        fitsampsize(c,im,iday,icl) = sum(inds);
                        fitsampsize(fitsampsize(c,im,iday,icl)<5)=0;
                    else
                        slp(c,im,iday,icl) = nan;
                        intcpt(c,im,iday,icl) = nan;
                        fitsampsize(c,im,iday,icl) = 0;
                    end

                    %Do some binning now for reldisp
                    [N,edges,bin]=histcounts(indvar,nbins);
                    binstdn = accumarray(bin,alleps)./N';
                    binstdn(N<5)=NaN;
                    meanvalsn{c,icl}(:,iday)=binstdn;
                    centersn=(edges(1:end-1)+edges(2:end))/2;

                    %Also fit a line to the data
                    inds=~isnan(binstdn);
                    x=centersn(inds);
                    y=binstdn(inds);

                    if length(x)>1
                        fitline = fit(x',y,'poly1');
                        slpn(c,im,iday,icl) = fitline.p1;
                        intcptn(c,im,iday,icl) = fitline.p2;
                        fitsampsizen(c,im,iday,icl) = sum(inds);
                        fitsampsizen(fitsampsizen(c,im,iday,icl)<5)=0;
                    else
                        slpn(c,im,iday,icl) = nan;
                        intcptn(c,im,iday,icl) = nan;
                        fitsampsizen(c,im,iday,icl) = 0;
                    end
                end
            end
        end   
%         toc
        for icl = 1:nlayers
            slp(slp==0) = nan;
            slpn(slpn==0) = nan;

            mslp_day(c,im,icl)=nanmean(slp(c,im,:,icl));
            mslpn_day(c,im,icl)=nanmean(slpn(c,im,:,icl));

            intcpt(intcpt==0) = nan;
            intcptn(intcptn==0) = nan;

            mintcpt_day(c,im,icl) = nanmean(intcpt(c,im,:,icl));
            mintcptn_day(c,im,icl) = nanmean(intcptn(c,im,:,icl));
            
            meand(meand==0) = nan;
            meand_day(c,icl) = nanmean(meand(c,:,icl));

            for iday = days_analyzed
                
                fitsampsize_vec = reshape(fitsampsize(c,im,:,icl),[],1);
                fitsampsizen_vec = reshape(fitsampsizen(c,im,:,icl),[],1);
                stdcam_slp(c,im,iday,icl)=nanstd(slp(c,im,:,icl)/mslp_day(c,...
                    im,icl),fitsampsize_vec);
                
                stdcam_slpn(c,im,iday,icl)=nanstd(slpn(c,im,:,icl)/mslpn_day(c,im,...
                    icl),fitsampsizen_vec);

                stdcam_intcpt(c,im,iday,icl)=nanstd(intcpt(c,im,:,...
                    icl)/mintcpt_day(c,im,icl),fitsampsize_vec);
                
                stdcam_intcptn(c,im,iday,icl)=nanstd(intcptn(c,im,:,...
                    icl)/mintcptn_day(c,im,icl),fitsampsizen_vec);
            end
            
            mstdcam_slp_day(c,im,icl)=mean(stdcam_slp(c,im,:,icl));
            mstdcam_slpn_day(c,im,icl)=mean(stdcam_slpn(c,im,:,icl));
            mstdcam_intcpt_day(c,im,icl)=mean(stdcam_intcpt(c,im,:,icl));
            mstdcam_intcptn_day(c,im,icl)=mean(stdcam_intcptn(c,im,:,icl));
        end
    end
end

%% plotting
clear F

close all
if do_animate==1 && do_plot_corr==1
    f = figure('position',[183 374 1800 800]);
    for i=1:2*length(campaigns)
        a(i)=subplot(2,length(campaigns),i);
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
            pos{c}(2) = pos{c}(2)-6.9/17;
            str = ['$\overline {D}$ = ', sprintf('%0.1f', meand_day(c,icl)),' $\mu$m'];
            ant = annotation(f,'textbox', pos{c}, 'String',str,...
                'FitBoxToText','on','Interpreter','latex','FontSize',18);
            ant.LineStyle='none';
            yyaxis left
            plot(Moments,mstdcam_slpn_day(c,:,icl),'LineWidth',2)
            ylim([max([min(mstdcam_slpn_day(c,:)) 0]) min([max(mstdcam_slpn_day(c,:)) 2]) ])
            xlabel('k in M_k (moment of the distribution)')
            ylabel('Relative \sigma of slopes of \epsilon_{DSD}(k)')
            yyaxis right
            plot(Moments,mstdcam_intcptn_day(c,:,icl),'LineWidth',2)
            ylim([max([min(mstdcam_intcptn_day(c,:)) 0]) min([max(mstdcam_intcptn_day(c,:)) 2]) ])
            ylabel('Relative \sigma of intercepts of \epsilon_{DSD}(k)')
            set(gca,'FontSize',16)
            
            
            axes(a(c+length(campaigns)))
            yyaxis left
            cptitle = sprintf('%0.2f - %0.2f', layer_edges(icl,1),...
                layer_edges(icl,2));
            title(cptitle)
            plot(Moments,mstdcam_slp_day(c,:,icl),'LineWidth',2)
            ylim([max([min(mstdcam_slp_day(c,:)) 0]) min([max(mstdcam_slp_day(c,:)) 2]) ])
            xlabel('k in M_k (moment of the distribution)')
            ylabel('Relative \sigma of slopes of \sigma_{DSD}(k)')

            yyaxis right
            plot(Moments,mstdcam_intcpt_day(c,:,icl),'LineWidth',2)
            ylim([max([min(mstdcam_intcpt_day(c,:)) 0]) min([max(mstdcam_intcpt_day(c,:)) 2]) ])
            ylabel('Relative \sigma of intercepts of \sigma_{DSD}(k)')
            set(gca,'FontSize',16)
            
        end
        
        F(icl) = getframe(gcf);
    end

end
%% save video
if do_save==1
    v = VideoWriter('reldisp_vs_layer.mp4','MPEG-4');
    v.FrameRate=5;
    open(v)
    writeVideo(v,F)
    close(v)
end