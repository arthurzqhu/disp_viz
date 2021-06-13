clear 
close all
load clouds.mat
load gomaccspdi_flight_basics.mat

z_CB = {nan; 550; 1090; 750; [1200 1000]; nan; nan; nan; [750 1200]; 1790; 1320; nan;...
    nan; 1630; [1990 1515]; 1820; [320 770 820]; [1200 1220]; nan; 770; 610};
z_CT = {nan; 925; 1330; 950; [2000 2300]; nan; nan; nan; [1060 1360]; 2300; 2300; nan;...
    nan; 3800; [2590 3190]; 2340; [1070 3240 1390]; [1980 2300]; nan; 1990; 2600};

%%
close all
iday=5;
s_t=clouds.gomaccspdi(iday).s_t;
a_z=clouds.gomaccspdi(iday).a_z;
a_t=clouds.gomaccspdi(iday).a_t;
a_lat=clouds.gomaccspdi(iday).a_lat;
a_lon=-clouds.gomaccspdi(iday).a_lon;
lwc=clouds.gomaccspdi(iday).s_lwc_pdi;

[t, cmt_ipdi, cmt_icabin] = intersect(s_t, a_t);

% plot(a_lon,a_lat,'.'); hold on
% imagesc(a_lon(cmt_icabin),a_lat(cmt_icabin),lwc(cmt_ipdi))
plot(t,lwc(cmt_ipdi),'.')
yyaxis right
plot(t,a_z(cmt_icabin))

%%
close all
for iday = 1:length(clouds.gomaccspdi)
    cloudleg_i = gomaccspdi_flight_basics(iday).ti;
    cloudleg_f = gomaccspdi_flight_basics(iday).tf; 
    a_t = clouds.gomaccspdi(iday).a_t;
    a_z = clouds.gomaccspdi(iday).a_z;
    a_thet = clouds.gomaccspdi(iday).a_thet;
    s_t = clouds.gomaccspdi(iday).s_t;
    s_qt = clouds.gomaccspdi(iday).s_qt;
    
    
    nleg=length(cloudleg_i);
    
    [t, cmt_ipdi, cmt_icabin] = intersect(s_t, a_t);
    
    ti_idx{iday} = arrayfun(@(x) findInSorted(t, cloudleg_i(x)), 1:nleg);
    tf_idx{iday} = arrayfun(@(x) findInSorted(t, cloudleg_f(x)), 1:nleg);
    
    z = clouds.gomaccspdi(iday).a_z(cmt_icabin);
    lwc = clouds.gomaccspdi(iday).s_lwc_pdi(cmt_ipdi);
    meand = clouds.gomaccspdi(iday).s_meand_pdi(cmt_ipdi);
    Tk = clouds.gomaccspdi(iday).a_T(cmt_icabin)+273.15;
    mr = clouds.gomaccspdi(iday).a_mr(cmt_icabin)/1e3;
    p = clouds.gomaccspdi(iday).a_p(cmt_icabin);
    ql = clouds.gomaccspdi(iday).s_ql(cmt_ipdi);
    
    qt_BL(iday,1) = mean(s_qt(z<mean(z_CB{iday})));
    qt_FT(iday,1) = mean(s_qt(z<mean(z_CT{iday})*1.3 & z>mean(z_CT{iday})*1.2));
    
    ent_ratio_qt = (s_qt-qt_BL(iday))./(qt_FT(iday)-qt_BL(iday));
    
    low_lwc_idx = lwc<0.01;
    
    normAC = zeros(size(z));
    for ileg = 1:nleg
        leg_ti = cloudleg_i(ileg);
        leg_tf = cloudleg_f(ileg);
        leg_CB = z_CB{iday}(ileg);
        leg_CT = z_CT{iday}(ileg);
        
        
        if ileg == 1
            be4_1st_cld = a_t<leg_ti;
            normAC(be4_1st_cld) = (a_z(be4_1st_cld)-leg_CB)/(leg_CT-leg_CB);
        end
        
        incloud_idx = a_t>=leg_ti & a_t<=leg_tf;
        normAC(incloud_idx) = (a_z(incloud_idx)-leg_CB)/(leg_CT-leg_CB);
        
        if ileg > 1
            prev_int_idx = a_t>cloudleg_f(ileg-1) & a_t<leg_ti;
            normAC(prev_int_idx) = (a_z(prev_int_idx)-leg_CB)/(leg_CT-leg_CB);
        end
        
        tleg_idx = t>=leg_ti & t<=leg_tf;
        T_CB = nanmean(Tk(normAC(tleg_idx)<.05 & normAC(tleg_idx) > -0.05));
        r_CB = nanmean(mr(normAC(tleg_idx)<.05 & normAC(tleg_idx) > -0.05));
        p_CB = nanmean(p(normAC(tleg_idx)<.05 & normAC(tleg_idx) > -0.05));
        
        ql_adb_lin = adiab_ql(z_CB{iday}(ileg),z_CT{iday}(ileg),T_CB,...
            r_CB,p_CB)*1000; 
        z_lin = linspace(z_CB{iday}(ileg),z_CT{iday}(ileg),...
            length(ql_adb_lin));
        
        if nleg == 1
            ql_adb_prof = interp1(z_lin,ql_adb_lin,z);
        elseif nleg > 1
            if ileg == 1
                ql_adb_prof = interp1(z_lin,ql_adb_lin,z(t<=cloudleg_f(ileg)));
            end
            
            if ileg > 1 && ileg < nleg
                ql_adb_prof = [ql_adb_prof;interp1(z_lin,...
                    ql_adb_lin,z(t>cloudleg_f(ileg-1) & t<=cloudleg_f(ileg)))];
            end
            
            if ileg == nleg
                ql_adb_prof = [ql_adb_prof;interp1(z_lin,...
                        ql_adb_lin,z(t>cloudleg_f(ileg-1)))];
            end
        end
        
        
        
%         if ~isnan(ti_idx{iday}(ileg))
%             z_leg = z(ti_idx{iday}(ileg):tf_idx{iday}(ileg));
%             lwc_leg = lwc(ti_idx{iday}(ileg):tf_idx{iday}(ileg));
%             t_leg = t(ti_idx{iday}(ileg):tf_idx{iday}(ileg));
%             meand_leg = meand(ti_idx{iday}(ileg):tf_idx{iday}(ileg));
%             figure
%             plot(lwc_leg,z_leg,'.')
%             title(['day ' num2str(iday) ', leg' num2str(ileg)])
%         end
    end
    
    Th_BL(iday) = nanmean(a_thet(normAC<0.05));
    Th_FT(iday) = nanmean(a_thet(normAC>1.01));
    
    ent_ratio_T = (a_thet(cmt_icabin)-Th_BL(iday))./(Th_FT(iday)-Th_BL(iday));
    
    gomaccspdi_flight_basics(iday).ti_idx = ti_idx{iday};
    gomaccspdi_flight_basics(iday).tf_idx = tf_idx{iday};
    gomaccspdi_flight_basics(iday).z_CB = z_CB{iday};
    gomaccspdi_flight_basics(iday).z_CT = z_CT{iday};
    
    if nleg>0
        clouds.gomaccspdi(iday).a_normAC = normAC;
        clouds.gomaccspdi(iday).ent_ratio_qt = ent_ratio_qt;
        clouds.gomaccspdi(iday).AF = ql./ql_adb_prof;
        clouds.gomaccspdi(iday).ent_ratio_T = ent_ratio_T;
    end
    %% pcasp validity testing
    a_t = clouds.gomaccspdi(iday).a_t;
    a_ntot = clouds.gomaccspdi(iday).a_ntot;
    ai_ntot = clouds.gomaccspdi(iday).a_ntot_ex;
    
    [cmt, cmt_ipdi, cmt_ipcasp] = intersect(t, a_t);
%     figure
%     plot(lwc(cmt_ipdi),ai_ntot(cmt_ipcasp),'.')
end

%%
% save('gomaccspdi_flight_basics.mat','gomaccspdi_flight_basics')
% save('clouds.mat','clouds', '-v7.3')