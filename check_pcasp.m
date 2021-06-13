clear
close all

load clouds.mat
campaigns={'vocalspdi','masepdi','postpdi','oraclespdi','gomaccspdi'};
camp_proper = {'VOCALS','MASE','POST','ORACLES','GoMACCS'};
colororder=colororder;
%%
close all

depvar_comb=[];
indvar_comb=[];

pcasp_edges=readmatrix('ORACLES_raw/pcasp_bin_edge.csv');
pcasp_mean = (pcasp_edges(:,2)+pcasp_edges(:,1))/2;

sl_bound=20;

iday=3;
for c = 4
    %%
    %     c=1;
    camp = campaigns{c};
    isOorG=any(strcmp({'oraclespdi';'gomaccspdi'},campaigns{c}));
    
    ndays = length(clouds.(camp));
    days_analyzed = 1:ndays;
    
    %     remove the days that have incomplete flights in vocals
    if strcmp(camp,'vocalspdi')
        days_analyzed(ismember(days_analyzed,[8,10,11,13]))=[];
    elseif strcmp(camp,'oraclespdi')
        days_analyzed(ismember(days_analyzed, [2,3,5,14,15,16,17,18]))=[];
    end
    
    for iday = 23%days_analyzed
        
        a_ntot_raw=clouds.(camp)(iday).a_ntot;
        a_ntot_ex=clouds.(camp)(iday).a_ntot_ex;
        a_conc=clouds.(camp)(iday).a_conc;
        s_ntot_raw=clouds.(camp)(iday).s_ntot_pdi;
        a_t=clouds.(camp)(iday).a_t;
        s_t=round(clouds.(camp)(iday).s_t);
        s_meand=clouds.(camp)(iday).s_meand_pdi;
        
        [cmt,cmt_ipdi,cmt_ipcasp] = intersect(s_t,a_t);
        
        cmat = a_ntot_raw(cmt_ipcasp);
        cmai = a_ntot_ex(cmt_ipcasp);
        cmac = a_conc(cmt_ipcasp,:);
        cmst = s_ntot_raw(cmt_ipdi);
        
        if isOorG
            z_raw = clouds.(camp)(iday).a_z;
            cmz = z_raw(cmt_ipcasp);
        else
            z_raw = clouds.(camp)(iday).s_ap;
            cmz = z_raw(cmt_ipdi);
        end
        
        cmmeand = s_meand(cmt_ipdi);
        
        
%         indvar_comb=[indvar_comb;cmmeand];
%         depvar_comb=[depvar_comb;cmai];
        
        sml_idx=cmmeand<=sl_bound;% | isnan(cmmeand);
        lg_idx=cmmeand>sl_bound;
        no_idx=isnan(cmmeand);
        
%         cmac(sml_idx,:);
        figure('position',[1722 557 859 420])
        
        hold on
        
        errorbar(pcasp_mean,mean(cmac(sml_idx,:)),std(cmac(sml_idx,:)),...
            'LineWidth',1)
        errorbar(pcasp_mean,mean(cmac(lg_idx,:)),std(cmac(lg_idx,:)),...
            'LineWidth',1)
        errorbar(pcasp_mean,mean(cmac(no_idx,:)),std(cmac(no_idx,:)),...
            'LineWidth',1)
        
        set(gca,'YScale','log')
        set(gca,'XScale','log')
        set(gca,'FontSize',16)
        legend(['mean N_a where 0<D<=' num2str(sl_bound) '\mum'],...
            ['mean N_a where D>' num2str(sl_bound) '\mum'],...
            'mean N_a where no droplet')
        
%         ylim([1 1e4])
        xlabel('aerosol size [\mum]')
        ylabel('mean pcasp concentration [1/cc]')
        xticks([.1:.1:1 2 3])
        hold off
        grid
        saveas(gcf,['plots/3 a_conc comp. bound=' num2str(sl_bound) 'um.png'])
        
%         figure
%         plot(cmmeand,cmat,'.')
%         xlabel('mean diameter [\mum]')
%         ylabel('pcasp+pdi concentration [1/cc]')
%         xlim([0 30])
%         ylim([0 5000])
%         set(gca,'FontSize',16)
        
%     
% %         saveas(gcf,['plots/check_pcasp/' camp ' Na vs Dbar day ',...
% %             num2str(iday) '.png'])
% %         pause(0.5)
% %         %%
%         
%         figure
%         
%         hold on
%         plot(cmat,cmz,'.')
%         
%         plot(cmst,cmz,'.','Color',colororder(2,:))
%         
% %         if max(cmst)>0 && max(cmst)*5<max(cmat)
% %             xlim([0 max(cmst)*5])
% %         elseif max(cmst)>0 && max(cmst)*5>=max(cmat)
%         xlim([0 8000])
% %         end
%         
%         xlabel('concentration [1/cc]')
%         ylabel('altitude [m]')
%         set(gca,'FontSize',16)
%         legend('total aerosol','droplet')
%         title(camp_proper{c})
%         ylim([0 5000])
        
% % %         saveas(gcf,['plots/check_pcasp/' camp ' day ' num2str(iday) '.png'])
% %         hold off
% % %         pause(0.5)
    end
    
    
    
end