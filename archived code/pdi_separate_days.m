%%

clear

cd '~/Box Sync/grad/research/datasets/PDI Data/'

load obs_clouds_wholeday.mat

clear meanvals meanvalsn
camps={'vocalspdi','masepdi','postpdi'};
ctitle={'VOCALS','MASE','POST'};
instr={'pdi','pdi','pdi'};
cptitle={'cloud base', 'midcloud', 'cloud top'};

minN=100;
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


%% compare among cloud top mid bottom 

%Data thresholds. Can play with these.

clear slp slpn intcpt intcptn

figure('position',[1183 374 1228 767])
for i=1:15
    a(i)=subplot(3,5,i);
    hold on
    set(gca,'linestyleorder','-')
end
set(gca,'linestyleorder','-')

for icl = 1:3 % iterator for different cloud base, midcloud, and cloud top
    c = 1; % only for vocals now
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
    colors=colormap(hsv(3));
    
    for iday = 1:ndays
        
        t = clouds.(campaign)(iday).s_t;
        t_cp = clouds.(campaign)(iday).(t_cp_name{icl});
        t_cp_idx = ismember(t,t_cp);
        
        [alleps,allnormN,alln,maxn]=deal([]);
        %Find all data that meet the thresholds
        %Remember clouds contains all data with non-NaN relative dispersion
        cldpts_all_cp=find(clouds.(campaign)(iday).(Vvar)>minLWC & ...
            clouds.(campaign)(iday).(Nvar)>minN & ...
            clouds.(campaign)(iday).(Tvar)>minT & ...
            clouds.(campaign)(iday).(RHvar)>minRH);
        
        % find the cloud droplets that are both valid by the standard above
        % and in the specific cloud portion
        cldpts = intersect(cldpts_all_cp, find(t_cp_idx));
        
        ncldpts = length(cldpts);
        %If there are enough data points, give data a more convenient name,
        %and collect data in cells (the cell data isn't actually used).
        if ncldpts>=minPTS
            maxn(iday,icl)=clouds.(campaign)(iday).maxN_pcl(icl);
            allnormN=clouds.(campaign)(iday).(Nvar)(cldpts)/maxn(iday,icl);
            alln=clouds.(campaign)(iday).(Nvar)(cldpts);
            alleps=clouds.(campaign)(iday).(epsvar)(cldpts);
            
            ALLN{icl}(end+1:end+ncldpts)=alln;
            ALLEPS{icl}(end+1:end+ncldpts)=alleps;
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
            meanvals{icl}(:,iday)=bineps;
            centers=(edges(1:end-1)+edges(2:end))/2;
            
            %Also fit a line to the data
            inds=~isnan(bineps);
            x=centers(inds);
            y=bineps(inds);
            
            fitline = fit(x',y,'poly1');
            slp{icl}(iday) = fitline.p1;
            intcpt{icl}(iday) = fitline.p2;
            
%             axes(a(iday))
% %             title(cptitle{iday})
%             set(gca,'colororder',colors)
%             plot(centers,bineps,'LineWidth',1)
%             xlabel('Droplet Conc. (#/cm^3)')
%             ylabel('Relative Dispersion')
%             set(gca,'xlim',[0 700],'ylim',[0.1 0.6])
%             set(gca,'FontSize',16)

            %Do some binning now for normalized droplet concentration
            [N,edges,bin]=histcounts(allnormN,nbins);
            bineps = accumarray(bin,alleps)./N';
            bineps(N<5)=NaN;
            meanvalsn{icl}(:,iday)=bineps;
            centers=(edges(1:end-1)+edges(2:end))/2;
            
            %Also fit a line to the data
            inds=~isnan(bineps);
            x=centers(inds);
            y=bineps(inds);
            fitline = fit(x',y,'poly1');
            slpn{icl}(iday) = fitline.p1;
            intcptn{icl}(iday) = fitline.p2;
            daymaxn{icl}(iday) = clouds.(campaign)(iday).maxN;
            
            axes(a(iday))
            set(gca,'colororder',colors)
            plot(centers,bineps,'LineWidth',1)
            xlabel('Normalized Droplet Conc.')
            ylabel('Relative Dispersion')
            set(gca,'xlim',[0 1],'ylim',[0.1 0.6])
            set(gca,'FontSize',13)
        end
    end
    
end
legend('Cloud base', 'Mid cloud', 'Cloud top')

saveas(gcf, 'VOCALS_cloudtmb_reldisp_eachday_norm.png')