clear
close all
cd '/Users/arthurhu/Box/grad/research/aerosol_reldisp/datasets/'
vocals_files = dir('vocals/*.mat');
%%
for ifile=1:length(vocals_files)
    clear s_*
    vocals_date = vocals_files(ifile).name(1:6);
    vocals_folder_date = ['20',vocals_date(1:2),'-',vocals_date(3:4),'-',vocals_date(5:6)];
    load(['vocals/', vocals_files(ifile).name])

    orig_lwc=s_lwc_pdi;
    orig_conc=s_conc_pdi;
    orig_ntot=s_ntot_pdi;
    orig_reff=s_reff_pdi;
    orig_std=s_std_pdi;
    orig_disp=s_disp_pdi;
    orig_R=s_R_pdi;
    %%
    diams=load('dp_pdi.mat');
    binwidth=log10(diams.pdi_bins(1,2)/diams.pdi_bins(1,1));
    upd_lwc_all=readmatrix(['vocals-updated-raw/',vocals_folder_date,...
        '/whole_flight_lwc_pdi_20',vocals_date,'.csv']);
    % upd_wit=readmatrix('vocals-updated-raw/2008-10-19/whole_flight_lwc_pdi_20081019.csv');
    upd_lwc_val=upd_lwc_all(:,2);
    upd_t=upd_lwc_all(:,1);

    % sum(upd_val(4445,2:end))*binwidth
    %% check time interval consistency
    close all

    plot(upd_t(2:end)-upd_t(1:end-1)); hold on
    plot(s_t(2:end)-s_t(1:end-1)); hold off
    
    figure
    plot(s_t,s_lwc_pdi); hold on; 
    plot(upd_t,upd_lwc_val,'--'); hold off
    pause
    %% find initial and final time; they might not start or end at the same time stamp
%     close all
    orig_t=round(s_t);
    
    if orig_t(1)<upd_t(1)
        orig_ti_idx=find(orig_t==upd_t(1),1,'first');
        upd_ti_idx=1;
    else
        orig_ti_idx=1;
        upd_ti_idx=find(upd_t==orig_t(1),1,'first');
    end
    
    if orig_t(end)<upd_t(end)
        orig_tf_idx=length(orig_t);
        upd_tf_idx=find(upd_t==orig_t(end),1,'first');
    else
        orig_tf_idx=find(orig_t==upd_t(end),1,'first');
        upd_tf_idx=length(upd_t);
    end

    %% update lwc & conc & ntot & reff & std & disp & R?
%     close all
    s_lwc_pdi(orig_ti_idx:orig_tf_idx)=upd_lwc_val(upd_ti_idx:upd_tf_idx);
    
    upd_conc_all=readmatrix(['vocals-updated-raw/',vocals_folder_date,...
        '/whole_flight_conc_pdi_20',vocals_date,'.csv']);
    upd_conc_val=upd_conc_all(:,2:end);

    s_conc_pdi(orig_ti_idx:orig_tf_idx,:)=upd_conc_val(upd_ti_idx:upd_tf_idx,:);
    s_ntot_pdi(orig_ti_idx:orig_tf_idx)=sum(upd_conc_val(upd_ti_idx:upd_tf_idx,:),2)*binwidth;
    s_reff_pdi(orig_ti_idx:orig_tf_idx)=(upd_conc_val(upd_ti_idx:upd_tf_idx,:)*diams.dp_pdi.^3)./...
        (upd_conc_val(upd_ti_idx:upd_tf_idx,:)*diams.dp_pdi.^2)/2;
    s_std_pdi(orig_ti_idx:orig_tf_idx) = arrayfun(@(x) std(diams.dp_pdi, s_conc_pdi(x,:)),...
        orig_ti_idx:orig_tf_idx)';

    s_meand_pdi = sum(s_conc_pdi(orig_ti_idx:orig_tf_idx,:).*diams.dp_pdi'*binwidth,2)./s_ntot_pdi(orig_ti_idx:orig_tf_idx);
    s_disp_pdi(orig_ti_idx:orig_tf_idx) = s_std_pdi(orig_ti_idx:orig_tf_idx)./s_meand_pdi;

    clear s_meand_pdi
    save(['vocals/',vocals_files(ifile).name],'s_*')
end
%%
% imagesc(diams.dp_pdi,s_t,log(s_conc_pdi))
% colorbar
% caxis([3 9])
% ylim([4e4 5.5e4])
% plot(upd_conc_t,upd_wit(:,2),'.'); hold on
% plot(s_t, s_ntot_pdi,'.')