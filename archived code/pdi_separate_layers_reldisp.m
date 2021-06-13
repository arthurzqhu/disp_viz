%%
clear
cd '~/Box Sync/grad/research/datasets/PDI Data/'

load obs_clouds_wholeday.mat

clear meanvals meanvalsn
camps={'vocalspdi','masepdi','postpdi'};
ctitle={'VOCALS','MASE','POST'};
instr={'pdi','pdi','pdi'};
cptitle={'cloud base', 'midcloud', 'cloud top'};
bar_ticks = {'base', 'mid', 'top'};

thresN=100;
thresRH=70;
thresT=1; %minimum temperature, just to make sure we avoid mixed-phase
thresLWC=0.05;
thresPTS=100; %minimum number of data samples
doresamp = false;

%Initialize cells - but not used.
ALLN=cell(1,3);
ALLSTD=cell(1,3);
ALLH=cell(1,3);
ALLEPS=cell(1,3);
ALLD=cell(1,3);
ALLPCP=cell(1,3);

cl_name = {'CB';'CM';'CT'}; % cloud layer name

%% separate out the cloud top mid bottom

clear slp slpn intcpt intcptn

figure('position',[1183 374 1228 767])
for i=1:6
    a(i)=subplot(2,3,i);
    hold on
    set(gca,'linestyleorder','-')
end
set(gca,'linestyleorder','-')
%f4=figure('position',[1500 550 1000 700]);

for icl = 1:3 % iterator for different cloud base, midcloud, and cloud top
    c = 1; % only for vocals now
    campaign=camps{c};
    nc = length(clouds.(campaign));
    
    %Get variable names
    Nvar=['s_ntot_',instr{c}];
    Hvar='s_ap';
    Vvar=['s_lwc_',instr{c}];
    Dvar='drpsz';
    Pvar='s_lwc_hvps';
    epsvar=['s_disp_',instr{c}]; %relative dispersion
    stdvar=['s_std_',instr{c}];
    Tvar='s_ta';
    RHvar='s_rh';
    distvar=['s_conc_',instr{c}];
    
    ndays = nc;
    colors=colormap(parula(ndays));
    
    days_analyzed = 1:ndays;
    days_analyzed(days_analyzed==10 | days_analyzed==13 | days_analyzed==11)=[]; % remove the 
    % days that have incomplete flights
    
    for iday = 1:ndays
        
        t = clouds.(campaign)(iday).s_t;
        t_cl = clouds.(campaign)(iday).(['t_',cl_name{icl}]);
        t_cl_idx = ismember(t,t_cl);
        
        [alleps,allstd,allnormN,alln,alld,allnormD,maxn,maxsupsat,allh,...
            allRH,maxRH,minRH,allnormRH,allsupsat,maxh,minh]=deal([]);
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
            alld=clouds.(campaign)(iday).(Dvar)(cldpts);
            alln=clouds.(campaign)(iday).(Nvar)(cldpts);
            allh=clouds.(campaign)(iday).(Hvar)(cldpts);
            allRH = clouds.(campaign)(iday).(RHvar)(cldpts);
            allsupsat = clouds.(campaign)(iday).(RHvar)(cldpts)-100;
            
            maxn(iday,icl)=max(alln);
            allnormN=alln/maxn(iday,icl);
            
            maxh(iday,icl)=max(allh);
            minh(iday,icl)=min(allh);
            allnormH=(allh-minh(iday,icl))/(maxh(iday,icl)-minh(iday,icl));
            
            maxRH(iday,icl)=max(allRH);
            minRH(iday,icl)=min(allRH);
            allnormRH=(allRH-minRH(iday,icl))/(maxRH(iday,icl)-minRH(iday,icl));
            
            maxsupsat(iday,icl)=max(allsupsat);
            allnormsupsat=allsupsat/maxsupsat(iday,icl);
            
            alleps=clouds.(campaign)(iday).(epsvar)(cldpts);
            allstd=clouds.(campaign)(iday).(stdvar)(cldpts);
            
            allnormD=alld/max(alld);
            
            ALLN{icl}(end+1:end+ncldpts(iday,icl))=alln;
            ALLEPS{icl}(end+1:end+ncldpts(iday,icl))=alleps;
            ALLSTD{icl}(end+1:end+ncldpts(iday,icl))=allstd;
            ALLH{icl}(end+1:end+ncldpts(iday,icl))=allh;
            ALLD{icl}(end+1:end+ncldpts(iday,icl))=alld;
            %ALLPCP{c}(end+1:end+ncldpts)=clouds.(campaign)(iday).(Pvar)(cldpts);
        else
            ncldpts(iday,icl) = 0;
        end        
        
        if ~isempty(alleps) %If data exists
            nbins=20; %Number of cloud droplet concentration bins
            indvar = alln; %independent variable
            indnormVar = allnormN;
            %try
            %Do some binning for cloud local data
            [N,edges,bin]=histcounts(indvar,nbins);
            %Find average relative dispersion in each bin
            bineps = accumarray(bin,alleps)./N';
            %Set to NaN if there are fewer than 5 points in the bin
            bineps(N<5)=NaN;
            meanvals{icl}(:,iday)=bineps;
            centers=(edges(1:end-1)+edges(2:end))/2;
            
            %Also fit a line to the data
            inds=~isnan(bineps);
            x=centers(inds);
            y=bineps(inds);
            
            fitline = fit(x',y,'poly1');
            slp(iday,icl) = fitline.p1;
            intcpt(iday,icl) = fitline.p2;
            fitsampsize(iday,icl) = sum(inds);
            fitsampsize(fitsampsize(iday,icl)<5)=0;
            
            
            axes(a(icl))
            title(cptitle{icl})
            set(gca,'colororder',colors)
            plot(centers,bineps,'LineWidth',1.5)
            xlabel('Droplet Conc. (#/cm^3)')
            ylabel('Relative Dispersion')
            set(gca,'xlim',[0 700],'ylim',[0.1 0.6])
            set(gca,'FontSize',16)
            
            %Do some binning now for normalized droplet concentration
            [N,edges,bin]=histcounts(indnormVar,nbins);
            bineps = accumarray(bin,alleps)./N';
            bineps(N<5)=NaN;
            meanvalsn{icl}(:,iday)=bineps;
            centers=(edges(1:end-1)+edges(2:end))/2;
            
            %Also fit a line to the data
            inds=~isnan(bineps);
            x=centers(inds);
            y=bineps(inds);
            fitline = fit(x',y,'poly1');
            slpn(iday,icl) = fitline.p1;
            intcptn(iday,icl) = fitline.p2;
%             daymaxn{icl}(iday) = clouds.(campaign)(iday).maxN;
            
            axes(a(icl+3))
            set(gca,'colororder',colors)
            plot(centers,bineps,'LineWidth',1.5)
            xlabel('Normalized Droplet Conc.')
            ylabel('Relative Dispersion')
            set(gca,'xlim',[0 1],'ylim',[0.1 0.6])
            set(gca,'FontSize',16)
        end
        

    end
    
end
% saveas(gcf,'VOCALS_cloudtmb_reldisp.png')

slp(slp==0) = nan;
slpn(slpn==0) = nan;

mslp=nanmean(slp);
mslpn=nanmean(slpn);

intcpt(intcpt==0) = nan;
intcptn(intcptn==0) = nan;

mintcpt = nanmean(intcpt);
mintcptn = nanmean(intcptn);

for icl=1:3
    stdcam_slp(icl,1)=nanstd(slp(:,icl)/mslp(icl),fitsampsize(:,icl));
    stdcam_slp(icl,2)=nanstd(slpn(:,icl)/mslpn(icl),fitsampsize(:,icl));
    
    stdcam_intcpt(icl,1)=nanstd(intcpt(:,icl)/mintcpt(icl),fitsampsize(:,icl));
    stdcam_intcpt(icl,2)=nanstd(intcptn(:,icl)/mintcptn(icl),fitsampsize(:,icl));
end

figure('Position',[800 600 800 400])
subplot(1,2,1)
bar(stdcam_slp)
ylabel('Standard Deviation/Mean of Fitted Slopes')
legend('N','N/Nmax','Location','best')
% set(gca,'xticklabel',cptitle)
set(gca,'FontSize',18)

subplot(1,2,2)
bar(1:3,mslpn)
hold on
errorbar(1:3,mslpn,stdcam_slp(:,2).*mslpn','.','LineWidth',2)
xticklabels(bar_ticks)
ylabel('Mean of Fitted Slopes')
% set(gca,'xticklabel',cptitle)
set(gca,'FontSize',18)

% saveas(gcf,'VOCALS_cloudtmb_slope_wgtd_std.png')

% figure
% bar(stdcam_intcpt)
% ylabel('Standard Deviation/Mean of Fitted Intercepts')
% legend('N','N/Nmax','Location','best')
% set(gca,'xticklabel',cptitle)
% set(gca,'FontSize',18)
% saveas(gcf,'VOCALS_cloudtmb_intercept_wgtd_std.png')