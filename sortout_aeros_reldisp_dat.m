clear
cd '~/Box/grad/research/datasets/PDI Data/'
load clouds.mat

%%
close all
campaigns = {'vocalspdi','postpdi','oraclespdi','gomaccspdi', 'masecas'};

for c = 1:length(campaigns)
    camp=campaigns{c};
    
    % load the time the aircraft entering the cloud and exiting the cloud
    fb = load([camp,'_flight_basics.mat']);
    fbvar = [camp,'_flight_basics'];
    
    ndays = length(clouds.(camp));
    days_analyzed = 1:ndays;
    
    % remove the days that have incomplete flights in vocals
    if c==1
        days_analyzed(days_analyzed==10 | days_analyzed==13 | days_analyzed==11)=[];
    end
    
    for iday = days_analyzed
        cloudlegs_i = fb.(fbvar)(iday).ti;
        cloudlegs_f = fb.(fbvar)(iday).tf;
        
        s_t = clouds.(camp)(iday).s_t;
        ccn_a = clouds.(camp)(iday).ccn_a;
        ccn_b = clouds.(camp)(iday).ccn_b;
        a_ntot = clouds.(camp)(iday).a_ntot;
        
        if c == 5
            a_t = s_t; 
            ccn_t = s_t;
            s_disp = clouds.(camp)(iday).s_disp_cas;
            s_ntot = clouds.(camp)(iday).s_ntot_cas;
            s_lwc = clouds.(camp)(iday).s_lwc_cas;
        else
            a_t = clouds.(camp)(iday).a_t;
            ccn_t = clouds.(camp)(iday).ccn_t;
            s_disp = clouds.(camp)(iday).s_disp_pdi;
            s_ntot = clouds.(camp)(iday).s_ntot_pdi;
            s_lwc = clouds.(camp)(iday).s_lwc_pdi;
        end
        
        allinst_commontime = mintersect(floor(s_t), floor(a_t), floor(ccn_t));
            
        cmt_ipdi = find(ismember(floor(s_t), allinst_commontime));
        cmt_ipcasp = find(ismember(floor(a_t), allinst_commontime));
        cmt_iccn = find(ismember(floor(ccn_t), allinst_commontime));
        
        cm_t = s_t(cmt_ipdi);
        cm_lwc = s_lwc(cmt_ipdi);
        cm_disp = s_disp(cmt_ipdi);
        cm_s_ntot = s_ntot(cmt_ipdi);
        cm_ccn_a = ccn_a(cmt_iccn);
        cm_ccn_b = ccn_b(cmt_iccn);
        cm_a_ntot = a_ntot(cmt_ipcasp);
        
        % each campaign can be optimized to include as many datapoints as
        % possible. also deal with other campaign-specific variables.
        switch camp
            case {'vocalspdi','postpdi', 'masecas'}
                cm_z = clouds.(camp)(iday).s_ap(cmt_ipdi);
                s_ap = clouds.(camp)(iday).s_ap;
                cm_s_ntot = s_ntot;
                cm_disp = s_disp;
            case 'oraclespdi'
                cm_z = clouds.(camp)(iday).hskp_z(cmt_ipcasp);
            case 'gomaccspdi'
                cm_z = clouds.(camp)(iday).a_z(cmt_ipcasp);
                cm_ccn_a = ccn_a;
                cm_ccn_b = ccn_b;
                cm_a_ntot = a_ntot;
        end
        
        if ~isempty(cloudlegs_i)
            for ileg = 1:length(cloudlegs_i)
                ti = cloudlegs_i(ileg);
                tf = cloudlegs_f(ileg);
                ti_idx = findInSorted(cm_t, ti);
                tf_idx = findInSorted(cm_t, tf);
                
                if ti_idx<=0 || isnan(ti_idx) || isnan(tf_idx)
                    continue
                end

                % to sample some extra distance below the cloud top, in case
                % the flight doesnt only go from low to high. does not apply to
                % mase due to its G-1 flight path
                if c~=5
                    if cm_z(ti_idx) < cm_z(tf_idx)
                        ti = ti - 600;
                    else
                        tf = tf + 600;
                    end
                end
                
                z_min = min(cm_z(cm_t < tf & cm_t > ti));
                z_max = max(cm_z(cm_t < tf & cm_t > ti));
                z_max_sampled = (z_min+z_max)/2;
                
                try
                    aerCMS = @(x) calcMeanSampsize(x, cm_t < tf & cm_t > ti & cm_z < z_max_sampled);
                    drpCMS = @(x) calcMeanSampsize(x, cm_t < cloudlegs_f(ileg) & cm_t > cloudlegs_i(ileg) & cm_s_ntot > 25 & cm_z < z_max);
                    
                    % include as many datapoints as possible for vocals, post, 
                    % and mase
                    if ismember(c,[1,2,5])
                        drpCMS = @(x) calcMeanSampsize(x, s_t < tf & s_t > ti & s_ntot > 25 & s_ap < z_max_sampled);
                    end
                    
                    [clouds.(camp)(iday).a_ntot_CB(ileg), clouds.(camp)(iday).a_ntot_CB_sampsize(ileg)] = ...
                        aerCMS(cm_a_ntot);
                    [clouds.(camp)(iday).s_ntot_CB(ileg), clouds.(camp)(iday).s_ntot_CB_sampsize(ileg)] = ...
                        drpCMS(cm_s_ntot);
                    [clouds.(camp)(iday).ccn_a_CB(ileg), clouds.(camp)(iday).ccn_a_CB_sampsize(ileg)] = ...
                        aerCMS(cm_ccn_a);
                    [clouds.(camp)(iday).ccn_b_CB(ileg), clouds.(camp)(iday).ccn_b_CB_sampsize(ileg)] = ...
                        aerCMS(cm_ccn_b);
                    [clouds.(camp)(iday).reldisp_CB(ileg), clouds.(camp)(iday).reldisp_CB_sampsize(ileg)] = ...
                        drpCMS(cm_disp);
                catch
                    [clouds.(camp)(iday).a_ntot_CB(ileg), clouds.(camp)(iday).a_ntot_CB_sampsize(ileg), ...
                        clouds.(camp)(iday).s_ntot_CB(ileg), clouds.(camp)(iday).s_ntot_CB_sampsize(ileg), ...
                        clouds.(camp)(iday).ccn_a_CB(ileg), clouds.(camp)(iday).ccn_a_CB_sampsize(ileg), ...
                        clouds.(camp)(iday).ccn_b_CB(ileg), clouds.(camp)(iday).ccn_b_CB_sampsize(ileg), ...
                        clouds.(camp)(iday).reldisp_CB(ileg), clouds.(camp)(iday).reldisp_CB_sampsize(ileg)] = ...
                        deal(nan);
                end
            end
        else
            [clouds.(camp)(iday).a_ntot_CB, clouds.(camp)(iday).a_ntot_CB_sampsize, ...
                clouds.(camp)(iday).s_ntot_CB, clouds.(camp)(iday).s_ntot_CB_sampsize, ...
                clouds.(camp)(iday).ccn_a_CB, clouds.(camp)(iday).ccn_a_CB_sampsize, ...
                clouds.(camp)(iday).ccn_b_CB, clouds.(camp)(iday).ccn_b_CB_sampsize, ...
                clouds.(camp)(iday).reldisp_CB, clouds.(camp)(iday).reldisp_CB_sampsize] = ...
                deal(nan);
        end
    end
    
end

%%
save('clouds.mat','clouds')

finishingTaskSound