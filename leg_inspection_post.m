clear
close all
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
load postpdi_flight_basics.mat
load clouds.mat

%%

thet_BL=[285;285;285;286;286;286;287;286;286;287;287;286;287;286;286;286];
thet_FT=[306;305;312;305;312;310;304;307;306;306;310;306;311;313;308;312];

qt_BL = [9.8;8.3;8.7;8.5;9.0;8.8;9.8;10.2;8.8;8.8;8.8;9.4;10.1;9.5;9.5;9.7];
qt_FT = [3.0;3.0;1.7;1.5;1.8;1.5;6.5;2.5 ;1.6;1.6;2.5;2.0;2.0 ;2.8;2.0;2.0];

%%
close all
iday=16;
figure('Position',[1732 209 794 630])
yyaxis right
plot(clouds.postpdi(iday).s_ap,'.')
yyaxis left
plot(clouds.postpdi(iday).s_thet)
grid('on')

%%

ti = {[65654 69695 73729]; [63827 67257 71462 75374]; [7194 11371 15285];...
    [63938 66426 70931]; [7536 11911 16168]; [6806 10713 14882];...
    [66075 69652 73252 76821]; [65185 68929 72239 76062];...
    [62258 65476 69202 73128 76991]; [5940 9020 11186 17866];...
    [4425 7770 11145 15282]; [4386 8027 11876 15603]; [4541 6971 11279 14936];...
    [4773 7400 11367 14984]; [66356 69475 73502]; [62798 66081 70506]};

tf = {[69043 73016 76760]; [65880 70666 74660 78408]; [10697 14700 16855]; ...
    [65810 69955 74560]; [10500 14865 17988]; [9993 14195 19740]; ...
    [68986 72691 76197 79565]; [67726 71575 75387 79208];...
    [64910 68519 72487 76241 79242]; [8483 10588 13736 19670];...
    [6628 10433 14325 18000]; [6500 11339 15156 18987]; [6748 10798 14659 19577];...
    [7020 10806 14448 18224]; [68924 72890 76685]; [65499 69307 72721]};

z_CB = {[250 210 150]; [180 250 290 290]; [235 210 140]; [360 710 550]; [200 160 280];...
    [310 260 180]; [40 85 70 100]; [30 80 100 40]; [30 340 370 330 280];...
    [360 400 370 200]; [290 520 380 430]; [320 240 240 330]; [60 95 100 75];...
    [50 100 180 220]; [130 100 140]; [150 130 150]};
z_CT = {[560 555 540]; [500 540 530 500]; [480 480 460]; [910 1050 960]; [490 440 550];...
    [670 630 550]; [330 350 350 370]; [400 400 400 270]; [640 690 670 630 530];...
    [540 570 585 480]; [660 760 790 740]; [650 500 580 610]; [460 500 520 520];...
    [375 410 480 480]; [440 460 450]; [450 425 450]};
%%
close all
for iday = 1:length(clouds.postpdi)
    
    t = floor(clouds.postpdi(iday).s_t);
    lwc = clouds.postpdi(iday).s_lwc_pdi;
    z = clouds.postpdi(iday).s_ap;
    s_ntot = clouds.postpdi(iday).s_ntot_pdi;
    ta = clouds.postpdi(iday).s_ta;
    thet = clouds.postpdi(iday).s_thet;
    thete = clouds.postpdi(iday).s_thete;
    wz = clouds.postpdi(iday).s_wz;
    s_qt = clouds.postpdi(iday).s_qt;
    
%     thet_BL(iday,1) = nanmean(thet(z<nanmean(z_CB{iday})));
%     thet_FT(iday,1) = nanmean(thet(z<nanmean(z_CT{iday})*1.3 & z>nanmean(z_CT{iday})*1.2));
    
    T_BL(iday,1) = nanmean(ta(z<nanmean(z_CB{iday})));
    T_FB(iday,1) = nanmean(ta(z<nanmean(z_CT{iday})*1.3 & z>nanmean(z_CT{iday})*1.2));
    
%     qt_BL(iday,1) = nanmean(s_qt(z<nanmean(z_CB{iday})));
%     qt_FT(iday,1) = nanmean(s_qt(z<nanmean(z_CT{iday})*1.3 & z>nanmean(z_CT{iday})*1.2));
    
    ent_ratio_T = (thet-thet_BL(iday))./(thet_FT(iday)-thet_BL(iday));
    ent_ratio_qt = (s_qt-qt_BL(iday))./(qt_FT(iday)-qt_BL(iday));
    
    
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
    %     figure
    %     plot(thet,z)
    %     xline(thet_BL(iday));
    %     xline(thet_FB(iday));
    %     plot(t,lwc,'.')
    %     yyaxis right
    %     plot(t,z)
    
    
    ti_idx{iday} = arrayfun(@(x) findInSorted(t, ti{iday}(x)), 1:length(ti{iday}));
    tf_idx{iday} = arrayfun(@(x) findInSorted(t, tf{iday}(x)), 1:length(tf{iday}));
    %     hold on
    %     scatter(ti{iday}, z(ti_idx{iday}), 'g')
    %     scatter(tf{iday}, z(tf_idx{iday}), 'b')
    %     for ileg = 1:length(ti{iday})
    %         z_leg = z(ti_idx{iday}(ileg):tf_idx{iday}(ileg));
    %         lwc_leg = lwc(ti_idx{iday}(ileg):tf_idx{iday}(ileg));
    %         t_leg = t(ti_idx{iday}(ileg):tf_idx{iday}(ileg));
    %         thet_leg = thet(ti_idx{iday}(ileg):tf_idx{iday}(ileg));
    %         wz_leg = wz(ti_idx{iday}(ileg):tf_idx{iday}(ileg));
    %         s_ntot_leg = s_ntot(ti_idx{iday}(ileg):tf_idx{iday}(ileg));
    %
    %         figure
    %         scatter(thet_leg,z_leg,[],t_leg, '.')
    %         title(['day ' num2str(iday) ', leg' num2str(ileg)])
    %
    %         colorbar
    %         CBline = refline(0, z_CB{iday}(ileg));
    %         CTline = refline(0, z_CT{iday}(ileg));
    %
    %         CBline.Color = 'r'; CBline.LineWidth = 1;
    %         CTline.Color = 'g'; CTline.LineWidth = 1;
    %
    %     end
    
    postpdi_flight_basics(iday).ti = ti{iday};
    postpdi_flight_basics(iday).tf = tf{iday};
    postpdi_flight_basics(iday).ti_idx = ti_idx{iday};
    postpdi_flight_basics(iday).tf_idx = tf_idx{iday};
    postpdi_flight_basics(iday).z_CB = z_CB{iday};
    postpdi_flight_basics(iday).z_CT = z_CT{iday};
    postpdi_flight_basics(iday).thet_BL = thet_BL(iday);
    postpdi_flight_basics(iday).thet_FB = thet_FT(iday);
    postpdi_flight_basics(iday).T_BL = T_BL(iday);
    postpdi_flight_basics(iday).T_FB = T_FB(iday);
    postpdi_flight_basics(iday).qt_BL = qt_BL(iday);
    postpdi_flight_basics(iday).qt_FB = qt_FT(iday);
    
    clouds.postpdi(iday).ent_ratio_T = ent_ratio_T;
    clouds.postpdi(iday).ent_ratio_qt = ent_ratio_qt;
    
    %% pcasp validity testing
%     a_t = clouds.postpdi(iday).a_t;
%     a_ntot = clouds.postpdi(iday).a_ntot;
%     ai_ntot = clouds.postpdi(iday).a_ntot_ex;
%     
%     [cmt, cmt_ipdi, cmt_ipcasp] = intersect(t, a_t);
%     figure
%     plot(lwc(cmt_ipdi),ai_ntot(cmt_ipcasp),'.')
end

save('postpdi_flight_basics.mat', 'postpdi_flight_basics')
save('clouds.mat','clouds', '-v7.3')