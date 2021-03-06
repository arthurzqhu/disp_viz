clearvars -except clouds
cd '~/MEGAsync/grad/research/aerosol_reldisp/datasets/'
if ~exist('clouds','var') load clouds.mat, end
load oraclespdi_flight_basics.mat
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

%%
for iday = 1:length(clouds.oraclespdi)
   %%
   %     ti{iday} = oraclespdi_flight_basics(iday).ti;
   %     tf{iday} = oraclespdi_flight_basics(iday).tf;
   %     z_CB{iday} = oraclespdi_flight_basics(iday).z_CB;
   %     z_CT{iday} = oraclespdi_flight_basics(iday).z_CT;
   
   nleg = length(ti{iday});
   
   t = clouds.oraclespdi(iday).s_t;
   lwc = clouds.oraclespdi(iday).s_lwc_pdi;
   z = clouds.oraclespdi(iday).s_ap;
   s_lwc = clouds.oraclespdi(iday).s_lwc_pdi;
   s_ap = clouds.oraclespdi(iday).s_ap;
   s_qt = clouds.oraclespdi(iday).s_qt;
   s_thet = clouds.oraclespdi(iday).s_thet;
   Tk = clouds.oraclespdi(iday).s_ta+273.15;
   mr = clouds.oraclespdi(iday).s_mr/1e3;
   p = clouds.oraclespdi(iday).s_p;
   ql = clouds.oraclespdi(iday).s_ql;
   
   %     thet_BL(iday,1) = mean(thet(z<mean(z_CB{iday})));
   %     thet_FT(iday,1) = mean(thet(z<mean(z_CT{iday})*1.3 & z>mean(z_CT{iday})*1.2));
   
   
   qt_BL(iday,1) = mean(s_qt(z<mean(z_CB{iday})));
   qt_FB(iday,1) = mean(s_qt(z<mean(z_CT{iday})*1.3 & z>mean(z_CT{iday})*1.2));
   
   ent_ratio_qt = (s_qt-qt_BL(iday))./(qt_FB(iday)-qt_BL(iday));
   
   ti_idx{iday} = arrayfun(@(x) findInSorted(t, ti{iday}(x)), 1:nleg);
   tf_idx{iday} = arrayfun(@(x) findInSorted(t, tf{iday}(x)), 1:nleg);
   %
%    for ileg=1:length(ti{iday})
%       figure('Position',[1773 440 560 420])
%       plot(lwc(ti_idx{iday}(ileg):tf_idx{iday}(ileg)),...
%          z(ti_idx{iday}(ileg):tf_idx{iday}(ileg)),'.')
%       
%       try
%          yline(z_CB{iday}(ileg),'--','linewidth',3)
%          yline(z_CT{iday}(ileg),'--','linewidth',3)
%       catch
%       end
%       
%       figure
%       plot(t,lwc,'.')
%       yyaxis right
%       plot(t,z)
%    end
   
   
   %     dcm = datacursormode(gcf);
   %
   %     set(dcm, 'UpdateFcn', @customDataCursorUpdateFcn, 'Enable', 'On');
   % Here's the function that specifies 5 decimal places
   %
   
   
   %
   low_lwc_idx = lwc<0.01;
   
   normAC = nan(size(s_ap));
   %
   for ileg = 1:nleg
      leg_ti = ti{iday}(ileg);
      leg_tf = tf{iday}(ileg);
      leg_CB = z_CB{iday}(ileg);
      leg_CT = z_CT{iday}(ileg);
      
      if ileg == 1
%          a_be4_1st_cld = t<leg_ti;
%          normAC(a_be4_1st_cld) = (s_ap(a_be4_1st_cld)-leg_CB)/(leg_CT-leg_CB);
      end
      
      a_incloud_idx = t>=leg_ti & t<=leg_tf;
      normAC(a_incloud_idx) = (s_ap(a_incloud_idx)-leg_CB)/(leg_CT-leg_CB);
      
      if ileg > 1
%          prev_int_idx = t>tf{iday}(ileg-1) & t<leg_ti;
%          normAC(prev_int_idx) = (s_ap(prev_int_idx)-leg_CB)/(leg_CT-leg_CB);
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
   
   thet_BL(iday) = nanmean(s_thet(normAC<0.05));
   thet_FT(iday) = nanmean(s_thet(normAC>1.01));
   
   ent_ratio_T = (s_thet-thet_BL(iday))./(thet_FT(iday)-thet_BL(iday));
   
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
   clouds.oraclespdi(iday).normAC = normAC;
   clouds.oraclespdi(iday).ent_ratio_T = ent_ratio_T;
   clouds.oraclespdi(iday).ent_ratio_qt = ent_ratio_qt;
   
   if nleg>0
      AF = ql./ql_adb_prof;
      AF(AF>1)=nan;
      clouds.oraclespdi(iday).AF = AF;
      clouds.oraclespdi(iday).ql_adb_prof = ql_adb_prof;
   else
      clouds.oraclespdi(iday).AF = ql*nan;
      clouds.oraclespdi(iday).ql_adb_prof = ql*nan;
   end
   
%    figure
%    plot(t,s_lwc)
%    yyaxis right
%    plot(t,normAC)
   
end
%%
% save('oraclespdi_flight_basics.mat', 'oraclespdi_flight_basics')
% save('clouds.mat','clouds', '-v7.3')

%%