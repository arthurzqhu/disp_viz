%%
clear
cd '~/Box Sync/grad/research/datasets/PDI Data/'

load obs_clouds_wholeday.mat

clear meanvals meanvalsn
camps={'vocalspdi','masepdi','postpdi'};
ctitle={'VOCALS','MASE','POST'};
instr={'pdi','pdi','pdi'};
cptitle={'cloud base', 'midcloud', 'cloud top'};

minN=50;
minRH=70;
minT=1; %minimum temperature, just to make sure we avoid mixed-phase
minLWC=0.01;
minPTS=100; %minimum number of data samples
doresamp = false;

%Initialize cells - but not used.
ALLN=cell(1,3);
ALLSTD=cell(1,3);
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
    
    for iday = 1:ndays
        
        t = clouds.(campaign)(iday).s_t;
        t_cl = clouds.(campaign)(iday).(['t_',cl_name{icl}]);
        t_cl_idx = ismember(t,t_cl);
        
        [alleps,allstd,allnormN,alln,alld,allnormD]=deal([]);
        %Find all data that meet the thresholds
        %Remember clouds contains all data with non-NaN relative dispersion
        cldpts_all_cp=find(clouds.(campaign)(iday).(Vvar)>minLWC & ...
            clouds.(campaign)(iday).(Nvar)>minN & ...
            clouds.(campaign)(iday).(Tvar)>minT & ...
            clouds.(campaign)(iday).(RHvar)>minRH);
        
        % find the cloud droplets that are both valid by the standard above
        % and in the specific cloud layer
        cldpts = intersect(cldpts_all_cp, find(t_cl_idx));
        
        ncldpts(iday,icl) = length(cldpts);
        
        %If there are enough data points, give data a more convenient name,
        %and collect data in cells (the cell data isn't actually used).
        if ncldpts(iday,icl)>=minPTS
            maxn(iday,icl)=clouds.(campaign)(iday).maxN_pcl(icl);
            allnormN=clouds.(campaign)(iday).(Nvar)(cldpts)/maxn(iday,icl);
            alln=clouds.(campaign)(iday).(Nvar)(cldpts);
            alleps=clouds.(campaign)(iday).(epsvar)(cldpts);
            allstd=clouds.(campaign)(iday).(stdvar)(cldpts);
            alld=clouds.(campaign)(iday).(Dvar)(cldpts);
            allnormD=clouds.(campaign)(iday).(Dvar)(cldpts)/max(alld);
            
            ALLN{icl}(end+1:end+ncldpts(iday,icl))=alln;
            ALLEPS{icl}(end+1:end+ncldpts(iday,icl))=alleps;
            ALLSTD{icl}(end+1:end+ncldpts(iday,icl))=allstd;
            ALLD{icl}(end+1:end+ncldpts(iday,icl))=alld;
            %ALLPCP{c}(end+1:end+ncldpts)=clouds.(campaign)(iday).(Pvar)(cldpts);
        else
            ncldpts(iday,icl) = 0;
        end        
        
        if ~isempty(allstd) %If data exists
            nbins=20; %Number of cloud droplet concentration bins
            %try
            %Do some binning for cloud local data
            [N,edges,bin]=histcounts(alln,nbins);
            %Find average relative dispersion in each bin
            binstd = accumarray(bin,allstd)./N';
            %Set to NaN if there are fewer than 5 points in the bin
            binstd(N<5)=NaN;
            meanvals{icl}(:,iday)=binstd;
            centers=(edges(1:end-1)+edges(2:end))/2;
            
            %Also fit a line to the data
            inds=~isnan(binstd);
            x=centers(inds);
            y=binstd(inds);
            
            fitline = fit(x',y,'poly1');
            slp(iday,icl) = fitline.p1;
            intcpt(iday,icl) = fitline.p2;
            fitsampsize(iday,icl) = sum(inds);
            fitsampsize(fitsampsize(iday,icl)<5)=0;
            
            
            axes(a(icl))
            title(cptitle{icl})
            set(gca,'colororder',colors)
            plot(centers,binstd,'LineWidth',1.5)
            xlabel('Droplet Conc. (#/cm^3)')
            ylabel('Relative Dispersion')
            set(gca,'xlim',[0 600],'ylim',[1 6])
            set(gca,'FontSize',16)
            
            %Do some binning now for normalized droplet concentration
            [N,edges,bin]=histcounts(allnormN,nbins);
            binstd = accumarray(bin,allstd)./N';
            binstd(N<5)=NaN;
            meanvalsn{icl}(:,iday)=binstd;
            centers=(edges(1:end-1)+edges(2:end))/2;
            
            %Also fit a line to the data
            inds=~isnan(binstd);
            x=centers(inds);
            y=binstd(inds);
            fitline = fit(x',y,'poly1');
            slpn(iday,icl) = fitline.p1;
            intcptn(iday,icl) = fitline.p2;
%             daymaxn{icl}(iday) = clouds.(campaign)(iday).maxN;
            
            axes(a(icl+3))
            set(gca,'colororder',colors)
            plot(centers,binstd,'LineWidth',1.5)
            xlabel('Normalized Droplet Conc.')
            ylabel('Relative Dispersion')
            set(gca,'xlim',[0 1],'ylim',[1 6])
            set(gca,'FontSize',16)
        end
    end
    
end
% saveas(gcf,'VOCALS_cloudtmb_reldisp.png')

%%
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

figure
bar(stdcam_slp)
ylabel('Standard Deviation/Mean of Fitted Slopes')
legend('N','N/Nmax','Location','best')
set(gca,'xticklabel',cptitle)
set(gca,'FontSize',18)
ylim([0 0.6])
% saveas(gcf,'VOCALS_cloudtmb_slope_wgtd_std.png')

figure
bar(stdcam_intcpt)
ylabel('Standard Deviation/Mean of Fitted Intercepts')
legend('N','N/Nmax','Location','best')
set(gca,'xticklabel',cptitle)
set(gca,'FontSize',18)
% saveas(gcf,'VOCALS_cloudtmb_intercept_wgtd_std.png')