clear
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
load oraclespdi_flight_basics.mat
load clouds.mat

%%

close all

% eyeball the ti, tf, cloud top and bottom first
% then check lwc vs z to make sure
ti = {45925; []; [30580 35470 45990]; []; []; []; [45810 51070]; 36150;
    [35840 42210]; 47600; 49230; []; 43393; []; 46996; [];
    28136; 41333; 49364; 51653; []; 37258; [];
    [47482 49241]; []};

tf = {46300; []; [31220 36210 46710]; []; []; []; [46580 51360]; 37180;
    [36950 42810]; 48750; 49891; []; 44079; []; 47776; [];
    28989; 41590; 50324; 52805; []; 38073; [];
    [48102 49541]; []};

z_CB = {1000; []; [2118 1427 200]; []; []; []; [1363 1550]; 2780;
    [400 1100]; 2280; 1580; []; 2216; []; 500; []; 2120;
    2369; 1376; 1125; []; 3409; []; [2812 1095]; []};
z_CB=cellfun(@(x) x*0.3048,z_CB,'un',0); % convert ft to meters

z_CT = {2700; []; [2937 2343 1450]; []; []; []; [2584 2220]; 3500;
    [1885 2000]; 3538; 2993; []; 2818; []; 2500; []; 2500;
    2729; 2112; 4213; []; 4043; []; [3160 7487]; []};
z_CT=cellfun(@(x) x*0.3048,z_CT,'un',0);

for iday = 1:length(clouds.oraclespdi)
    %     ti{iday} = oraclespdi_flight_basics(iday).ti;
    %     tf{iday} = oraclespdi_flight_basics(iday).tf;
    %     z_CB{iday} = oraclespdi_flight_basics(iday).z_CB;
    %     z_CT{iday} = oraclespdi_flight_basics(iday).z_CT;
    
    s_t = floor(clouds.oraclespdi(iday).s_t);
    s_lwc = clouds.oraclespdi(iday).s_lwc_pdi;
    a_z = clouds.oraclespdi(iday).a_z;
    s_ntot = clouds.oraclespdi(iday).s_ntot_pdi;
    a_t = clouds.oraclespdi(iday).a_t;
    a_thet = clouds.oraclespdi(iday).a_thet;
    s_qt = clouds.oraclespdi(iday).s_qt;
    
    [~, cmt_ipdi, cmt_ihskp] = intersect(s_t, a_t);
    
    nleg = length(ti{iday});
    
    t = s_t(cmt_ipdi);
    lwc = s_lwc(cmt_ipdi);
    z = a_z(cmt_ihskp);
    thet = a_thet(cmt_ihskp);
    Tk = clouds.oraclespdi(iday).a_T(cmt_ihskp)+273.15;
    mr = clouds.oraclespdi(iday).a_mr(cmt_ihskp)/1e3;
    p = clouds.oraclespdi(iday).a_p(cmt_ihskp);
    ql = clouds.oraclespdi(iday).s_ql(cmt_ipdi);
    
%     thet_BL(iday,1) = mean(thet(z<mean(z_CB{iday})));
%     thet_FT(iday,1) = mean(thet(z<mean(z_CT{iday})*1.3 & z>mean(z_CT{iday})*1.2));
    
    
    qt_BL(iday,1) = mean(s_qt(z<mean(z_CB{iday})));
    qt_FB(iday,1) = mean(s_qt(z<mean(z_CT{iday})*1.3 & z>mean(z_CT{iday})*1.2));
    
    ent_ratio_qt = (s_qt-qt_BL(iday))./(qt_FB(iday)-qt_BL(iday));
    
    ti_idx{iday} = arrayfun(@(x) findInSorted(t, ti{iday}(x)), 1:nleg);
    tf_idx{iday} = arrayfun(@(x) findInSorted(t, tf{iday}(x)), 1:nleg);
%     
%     for ileg=1:length(ti{iday})
%         figure('Position',[1773 440 560 420])
%         plot(lwc(ti_idx{iday}(ileg):tf_idx{iday}(ileg)),...
%             z(ti_idx{iday}(ileg):tf_idx{iday}(ileg)),'.')
%         
%         try
%             yline(z_CB{iday}(ileg),'--','linewidth',3)
%             yline(z_CT{iday}(ileg),'--','linewidth',3)
%         catch
%         end
%         
%         %     plot(t,lwc,'.')
%         %     yyaxis right
%         %     plot(t,z)
%     end
    
    
%     dcm = datacursormode(gcf);
%     
%     set(dcm, 'UpdateFcn', @customDataCursorUpdateFcn, 'Enable', 'On');
    % Here's the function that specifies 5 decimal places
    %
    
    
    %
    low_lwc_idx = lwc<0.01;
    
    normAC = zeros(size(a_z));
%
    for ileg = 1:nleg
        leg_ti = ti{iday}(ileg);
        leg_tf = tf{iday}(ileg);
        leg_CB = z_CB{iday}(ileg);
        leg_CT = z_CT{iday}(ileg);

        if ileg == 1
            a_be4_1st_cld = a_t<leg_ti;
            normAC(a_be4_1st_cld) = (a_z(a_be4_1st_cld)-leg_CB)/(leg_CT-leg_CB);
        end

        a_incloud_idx = a_t>=leg_ti & a_t<=leg_tf;
        normAC(a_incloud_idx) = (a_z(a_incloud_idx)-leg_CB)/(leg_CT-leg_CB);

        if ileg > 1
            prev_int_idx = a_t>tf{iday}(ileg-1) & a_t<leg_ti;
            normAC(prev_int_idx) = (a_z(prev_int_idx)-leg_CB)/(leg_CT-leg_CB);
        end
        
        tleg_idx = t>=leg_ti & t<=leg_tf;
        T_CB = mean(Tk(normAC(tleg_idx)<.05 & normAC(tleg_idx) > -0.05));
        r_CB = mean(mr(normAC(tleg_idx)<.05 & normAC(tleg_idx) > -0.05));
        p_CB = mean(p(normAC(tleg_idx)<.05 & normAC(tleg_idx) > -0.05));
        
        ql_adb_lin = adiab_ql(leg_CB,leg_CT,T_CB,r_CB,p_CB)*1000; 
        z_lin = linspace(leg_CB,leg_CT,length(ql_adb_lin));
        
        if nleg == 1
            ql_adb_prof = interp1(z_lin,ql_adb_lin,z);
        elseif nleg > 1
            if ileg == 1
                ql_adb_prof = interp1(z_lin,ql_adb_lin,z(t<=leg_tf));
            end
            
            if ileg > 1 && ileg < nleg
                ql_adb_prof = [ql_adb_prof;interp1(z_lin,...
                    ql_adb_lin,z(t>tf{iday}(ileg-1) & t<=tf{iday}(ileg)))];
            end
            
            if ileg == nleg
                ql_adb_prof = [ql_adb_prof;interp1(z_lin,...
                        ql_adb_lin,z(t>tf{iday}(ileg-1)))];
            end
        end
        
        
    end
    
    thet_BL(iday) = nanmean(a_thet(normAC<0.05));
    thet_FT(iday) = nanmean(a_thet(normAC>1.01));
    
    ent_ratio_T = (thet-thet_BL(iday))./(thet_FT(iday)-thet_BL(iday));

    oraclespdi_flight_basics(iday).ti = ti{iday};
    oraclespdi_flight_basics(iday).tf = tf{iday};
    oraclespdi_flight_basics(iday).ti_idx = ti_idx{iday};
    oraclespdi_flight_basics(iday).tf_idx = tf_idx{iday};
    oraclespdi_flight_basics(iday).z_CB = z_CB{iday};
    oraclespdi_flight_basics(iday).z_CT = z_CT{iday};
    oraclespdi_flight_basics(iday).thet_BL = thet_BL(iday);
    oraclespdi_flight_basics(iday).thet_FT = thet_FT(iday);
    oraclespdi_flight_basics(iday).qt_BL = qt_BL(iday);
    oraclespdi_flight_basics(iday).qt_FB = qt_FB(iday);
%
    if nleg>0
        clouds.oraclespdi(iday).a_normAC = normAC;
        clouds.oraclespdi(iday).ent_ratio_T = ent_ratio_T;
        clouds.oraclespdi(iday).ent_ratio_qt = ent_ratio_qt;
        clouds.oraclespdi(iday).AF = ql./ql_adb_prof;
    end
    %% pcasp validity testing
    a_t = clouds.oraclespdi(iday).a_t;
    a_ntot = clouds.oraclespdi(iday).a_ntot;
    ai_ntot = clouds.oraclespdi(iday).a_ntot_ex;

    [cmt, cmt_ipdi, cmt_ipcasp] = intersect(t, a_t);
%         figure
%         plot(s_lwc(cmt_ipdi),ai_ntot(cmt_ipcasp),'.')
%         xlabel('LWC (g/m^3)')
%         ylabel('interstitial (PCASP) aerosol # conc')
%         set(gca,'fontsize',16)
end
%%
% save('oraclespdi_flight_basics.mat', 'oraclespdi_flight_basics')
save('clouds.mat','clouds', '-v7.3')

% %%
% function txt = customDataCursorUpdateFcn(~, event)
% pos = event.Position;
% txt = {sprintf('X: %.5f', pos(1)), sprintf('Y: %.5f', pos(2))};
% end % <- may not be needed