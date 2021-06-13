%%
clear
cd '~/Box/grad/research/datasets/PDI Data/'

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
ALLPCP=cell(1,3);

cl_name = {'CB';'CM';'CT'}; % cloud layer name

%%
clear slp slpn intcpt intcptn

figure('position',[1183 374 1228 767])
for i=1:6
    a(i)=subplot(2,3,i);
    hold on
    set(gca,'linestyleorder','-')
end
set(gca,'linestyleorder','-')
%f4=figure('position',[1500 550 1000 700]);

%For each campaign ...
for c=1:length(camps)
    campaign=camps{c};
    
    nc = length(clouds.(campaign));
    
    %Get variable names
    Nvar=['s_ntot_',instr{c}];
    Vvar=['s_lwc_',instr{c}];
    Pvar='s_lwc_hvps';
    epsvar=['s_disp_',instr{c}]'; %relative dispersion
    Tvar='s_ta';
    RHvar='s_rh';
    distvar=['s_conc_',instr{c}];
    
    ndays = nc;
    colors=colormap(parula(ndays));
    
    %For each flight (day) ...
    for iday = 1:ndays
        [alleps,allnormN,alln,maxn]=deal([]);
        %Find all data that meet the thresholds
        %Remember clouds contains all data with non-NaN relative dispersion
        cldpts=find(clouds.(campaign)(iday).(Vvar)>minLWC & ...
            clouds.(campaign)(iday).(Nvar)>minN & ...
            clouds.(campaign)(iday).(Tvar)>minT & ...
            clouds.(campaign)(iday).(RHvar)>minRH);
        
        ncldpts = length(cldpts);
        %If there are enough data points, give data a more convenient name,
        %and collect data in cells (the cell data isn't actually used).
        if ncldpts>=minPTS
            maxn(iday)=clouds.(campaign)(iday).maxN;
            allnormN=clouds.(campaign)(iday).(Nvar)(cldpts)/maxn(iday);
            alln=clouds.(campaign)(iday).(Nvar)(cldpts);
            alleps=clouds.(campaign)(iday).(epsvar)(cldpts);
            
            ALLN{c}(end+1:end+ncldpts)=alln;
            ALLEPS{c}(end+1:end+ncldpts)=alleps;
            %ALLPCP{c}(end+1:end+ncldpts)=clouds.(campaign)(iday).(Pvar)(cldpts);
        end        
        
        if ~isempty(alleps) %If data exists
            nbins=20; %Number of cloud droplet concentration bins
            %try
            %Do some binning for cloud local data
            [N,edges,bin]=histcounts(alln,nbins);
            %Find average relative dispersion in each bin
            bineps = accumarray(bin,alleps)./N';
            %Set to NaN if there are fewer than 5 points in the bin
            bineps(N<5)=NaN;
            meanvals{c}(:,iday)=bineps;
            centers=(edges(1:end-1)+edges(2:end))/2;
            
            %Also fit a line to the data
            inds=~isnan(bineps);
            x=centers(inds);
            y=bineps(inds);
            
            fitline = fit(x',y,'poly1');
            slp{c}(iday) = fitline.p1;
            intcpt{c}(iday) = fitline.p2;
            
            axes(a(c))
            title(ctitle{c})
            set(gca,'colororder',colors)
            plot(centers,bineps,'LineWidth',1)
            xlabel('Droplet Conc. (#/cm^3)')
            ylabel('Relative Dispersion')
            set(gca,'xlim',[0 650],'ylim',[0.1 0.7],'FontSize',18)
            %Do some binning now for normalized droplet concentration
            [N,edges,bin]=histcounts(allnormN,nbins);
            bineps = accumarray(bin,alleps)./N';
            bineps(N<5)=NaN;
            meanvalsn{c}(:,iday)=bineps;
            centers=(edges(1:end-1)+edges(2:end))/2;
            
            %Also fit a line to the data
            inds=~isnan(bineps);
            x=centers(inds);
            y=bineps(inds);
            fitline = fit(x',y,'poly1');
            slpn{c}(iday) = fitline.p1;
            intcptn{c}(iday) = fitline.p2;
            daymaxn{c}(iday) = clouds.(campaign)(iday).maxN;
            
            axes(a(c+3))
            set(gca,'colororder',colors)
            plot(centers,bineps,'LineWidth',1)
            xlabel('Normalized Droplet Conc.')
            ylabel('Relative Dispersion')
            set(gca,'xlim',[0 1],'ylim',[0.1 0.7],'FontSize',18)
        end
    end
end
saveas(gcf,'camps_reldisp.png')
%%
mslp=nanmean([slp{:}]);
mslpn=nanmean([slpn{:}]);
for c=1:3
    stdcam_slp(c,1)=nanstd(slp{c}/mslp);
    stdcam_slp(c,2)=nanstd(slpn{c}/mslpn);
end
figure
bar(stdcam_slp)
ylabel('Relative std of the Fitted Slopes')
legend('N','N/Nmax','Location','northwest')
set(gca,'xticklabel',ctitle,'FontSize',18)
saveas(gcf,'camps_slope_std.png')
% figure
% hold on
% mslp=mean([slp{:}]);
% for c=1:3
%     scatter(daymaxn{c},slp{c}/mslp,'MarkerFaceColor','flat')
% end
% mslp=mean([slpn{:}]);
% set(gca,'colororderindex',1)
% for c=1:3
%     scatter(daymaxn{c},slpn{c}/mslp)
% end
% 
% 
% figure
% hold on
% for c=1:3
%     scatter(daymaxn{c},intcpt{c},'MarkerFaceColor','flat')
% end