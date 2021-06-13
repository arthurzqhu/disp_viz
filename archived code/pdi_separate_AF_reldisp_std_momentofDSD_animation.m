clear
% cd '~/Box/grad/research/datasets/PDI Data/'
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

AF_left_edges = 0:0.01:0.70;
naf = length(AF_left_edges);
AF_edges = [AF_left_edges' AF_left_edges'+0.3];
AF_mean=(AF_edges(1:end-1)+AF_edges(2:end))/2; % layer center diameter

Moments = 0:0.5:10;

do_plot = false; % whether to plot the slopes figure
do_plot_corr = true;
do_plot_slope = false;
do_animate = true;
do_save = false;

%% separate out the cloud top mid bottom

clear slp slpn intcpt intcptn stdcam_intcpt stdcam_slp fitsampsize
clear binstd binstdn fitsampsizen mstdcam_slp_day mstdcam_slpn_day
clear mstdcam_intcpt_day mstdcam_intcptn_day mslp_day mslpn_day
clear mintcpt_day mintcptn_day

% c=2;

close all

for c=1:length(camps)
    
    campaign=camps{c};
    nc = length(clouds.(campaign));

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
    AFvar='AF';

ndays = nc;

    days_analyzed = 1:ndays;
    
    if c==1
        days_analyzed(days_analyzed==10 | days_analyzed==13 | days_analyzed==11)=[]; % remove the 
    % days that have incomplete flights in vocals
    
    end
    for im = 1:length(Moments)
        disp(im)
        for iaf = 1:naf % iterator for different each layer of the cloud
            for iday = days_analyzed
                binlims = clouds.(campaign)(iday).binlims;
                binmean = clouds.(campaign)(iday).binmean;
                dlogD = log10(binlims(2)/binlims(1));

                t = clouds.(campaign)(iday).s_t;
                AF = clouds.(campaign)(iday).(AFvar);
                t_af = t(AF>AF_edges(iaf,1) & AF<AF_edges(iaf,2));
                t_af_idx = ismember(t,t_af);

                [alleps,allstd,alln,alld]=deal([]);
                %Find all data that meet the thresholds
                %Remember clouds contains all data with non-NaN relative dispersion
                cldpts_all_cp=find(clouds.(campaign)(iday).(Lvar)>thresLWC & ...
                    clouds.(campaign)(iday).(Nvar)>thresN & ...
                    clouds.(campaign)(iday).(Tvar)>thresT & ...
                    clouds.(campaign)(iday).(RHvar)>thresRH);

                % find the cloud droplets that are both valid by the standard above
                % and in the specific cloud layer
                cldpts = intersect(cldpts_all_cp, find(t_af_idx));

                ncldpts(c,iday,iaf) = length(cldpts);

                %If there are enough data points, give data a more convenient name,
                %and collect data in cells (the cell data isn't actually used).
                if ncldpts(c,iday,iaf)>=thresPTS
                    alld=clouds.(campaign)(iday).(Dvar)(cldpts,:);
                    allnd=clouds.(campaign)(iday).(nvar)(cldpts,:);
                    allmx=sum(dlogD*allnd.*binmean'.^Moments(im),2);

                    meand(c,iday,iaf) = mean(alld);
                    allstd=clouds.(campaign)(iday).(stdvar)(cldpts);
                    alleps=clouds.(campaign)(iday).(epsvar)(cldpts);
                    
                else
                    ncldpts(c,iday,iaf) = 0;
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
                    meanvals{c,iaf}(:,iday)=binstd;
                    centers=(edges(1:end-1)+edges(2:end))/2;
                    ref_binstd = sqrt((centers-1)/2);

                    %Also fit a line to the data
                    inds=~isnan(binstd);
                    x=centers(inds);
                    y=binstd(inds);


                    if length(x)>1
                        fitline = fit(x',y,'poly1');
                        slp(c,im,iday,iaf) = fitline.p1;
                        intcpt(c,im,iday,iaf) = fitline.p2;
                        fitsampsize(c,im,iday,iaf) = sum(inds);
                        fitsampsize(fitsampsize(c,im,iday,iaf)<5)=0;
                    else
                        slp(c,im,iday,iaf) = nan;
                        intcpt(c,im,iday,iaf) = nan;
                        fitsampsize(c,im,iday,iaf) = 0;
                    end

                    %Do some binning now for reldisp
                    [N,edges,bin]=histcounts(indvar,nbins);
                    binstdn = accumarray(bin,alleps)./N';
                    binstdn(N<5)=NaN;
                    meanvalsn{c,iaf}(:,iday)=binstdn;
                    centersn=(edges(1:end-1)+edges(2:end))/2;

                    %Also fit a line to the data
                    inds=~isnan(binstdn);
                    x=centersn(inds);
                    y=binstdn(inds);

                    if length(x)>1
                        fitline = fit(x',y,'poly1');
                        slpn(c,im,iday,iaf) = fitline.p1;
                        intcptn(c,im,iday,iaf) = fitline.p2;
                        fitsampsizen(c,im,iday,iaf) = sum(inds);
                        fitsampsizen(fitsampsizen(c,im,iday,iaf)<5)=0;
                    else
                        slpn(c,im,iday,iaf) = nan;
                        intcptn(c,im,iday,iaf) = nan;
                        fitsampsizen(c,im,iday,iaf) = 0;
                    end
                end
            end
        end   
%         toc
        for iaf = 1:naf
            slp(slp==0) = nan;
            slpn(slpn==0) = nan;

            mslp_day(c,im,iaf)=nanmean(slp(c,im,:,iaf));
            mslpn_day(c,im,iaf)=nanmean(slpn(c,im,:,iaf));

            intcpt(intcpt==0) = nan;
            intcptn(intcptn==0) = nan;

            mintcpt_day(c,im,iaf) = nanmean(intcpt(c,im,:,iaf));
            mintcptn_day(c,im,iaf) = nanmean(intcptn(c,im,:,iaf));
            
            meand(meand==0) = nan;
            meand_day(c,iaf) = nanmean(meand(c,:,iaf));

            for iday = days_analyzed
                
                fitsampsize_vec = reshape(fitsampsize(c,im,:,iaf),[],1);
                fitsampsizen_vec = reshape(fitsampsizen(c,im,:,iaf),[],1);
                stdcam_slp(c,im,iday,iaf)=nanstd(slp(c,im,:,iaf)/mslp_day(c,...
                    im,iaf),fitsampsize_vec);
                
                stdcam_slpn(c,im,iday,iaf)=nanstd(slpn(c,im,:,iaf)/mslpn_day(c,im,...
                    iaf),fitsampsizen_vec);

                stdcam_intcpt(c,im,iday,iaf)=nanstd(intcpt(c,im,:,...
                    iaf)/mintcpt_day(c,im,iaf),fitsampsize_vec);
                
                stdcam_intcptn(c,im,iday,iaf)=nanstd(intcptn(c,im,:,...
                    iaf)/mintcptn_day(c,im,iaf),fitsampsizen_vec);
            end
            
            mstdcam_slp_day(c,im,iaf)=mean(stdcam_slp(c,im,:,iaf));
            mstdcam_slpn_day(c,im,iaf)=mean(stdcam_slpn(c,im,:,iaf));
            mstdcam_intcpt_day(c,im,iaf)=mean(stdcam_intcpt(c,im,:,iaf));
            mstdcam_intcptn_day(c,im,iaf)=mean(stdcam_intcptn(c,im,:,iaf));
        end
    end
end

%% plotting
clear F

close all
if do_animate==1 && do_plot_corr==1
    f = figure('position',[183 374 1800 800]);
    for i=1:2*length(camps)
        a(i)=subplot(2,length(camps),i);
    end
    
    for iaf = 1:naf
        
        % prevent subsequent annotations from overlapping
        try
            delete(ant.Parent.Children)
        catch
        end
%         c = 2;
        for c=1:length(camps)
            axes(a(c))
            title(ctitle{c})
            pos = get(a, 'Position');
            pos{c}(1) = pos{c}(1)+1/15;
            pos{c}(2) = pos{c}(2)-6.9/17;
            str = ['$\overline {D}$ = ', sprintf('%0.1f',meand_day(c,iaf)),' $\mu$m'];
            ant = annotation(f,'textbox', pos{c}, 'String',str,...
                'FitBoxToText','on','Interpreter','latex','FontSize',18);
            ant.LineStyle='none';
            yyaxis left
            plot(Moments,mstdcam_slpn_day(c,:,iaf),'LineWidth',2)
            ylim([max([min(mstdcam_slpn_day(c,:)) 0]) min([max(mstdcam_slpn_day(c,:)) 2]) ])
            xlabel('k in M_k (moment of the distribution)')
            ylabel('Relative \sigma of slopes of \epsilon(k)')
            yyaxis right
            plot(Moments,mstdcam_intcptn_day(c,:,iaf),'LineWidth',2)
            ylim([max([min(mstdcam_intcptn_day(c,:)) 0]) min([max(mstdcam_intcptn_day(c,:)) 2]) ])
            ylabel('Relative \sigma of intercepts of \epsilon(k)')
            set(gca,'FontSize',16)
            
            
            axes(a(c+length(camps)))
            yyaxis left
            cptitle = sprintf('%0.2f - %0.2f', AF_edges(iaf,1),...
                AF_edges(iaf,2));
            title(cptitle)
            plot(Moments,mstdcam_slp_day(c,:,iaf),'LineWidth',2)
            ylim([max([min(mstdcam_slp_day(c,:)) 0]) min([max(mstdcam_slp_day(c,:)) 2]) ])
            xlabel('k in M_k (moment of the distribution)')
            ylabel('Relative \sigma of slopes of \sigma(k)')

            yyaxis right
            plot(Moments,mstdcam_intcpt_day(c,:,iaf),'LineWidth',2)
            ylim([max([min(mstdcam_intcpt_day(c,:)) 0]) min([max(mstdcam_intcpt_day(c,:)) 2]) ])
            ylabel('Relative \sigma of intercepts of \sigma(k)')
            set(gca,'FontSize',16)
            
        end
        
        F(naf-iaf+1) = getframe(gcf);
    end

end
%% save video
if do_save==1
    v = VideoWriter('reldisp_vs_AF_M10.mp4','MPEG-4');
    v.FrameRate=10;
    open(v)
    writeVideo(v,F)
    close(v)
end