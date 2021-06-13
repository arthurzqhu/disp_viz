clear
cd '~/Box/grad/research/datasets/PDI Data/'
load clouds.mat

%%
close all

campaigns={'vocals','post','oracles','gomaccs','mase'};
ctitle={'VOCALS','POST','ORACLES','GoMACCS','MASE'};
drp_instr={'pdi','pdi','pdi','pdi','cas'};

c = 1;
camp = [campaigns{c} drp_instr{c}];

fb = load([camp,'_flight_basics.mat']);
fbvar = [camp,'_flight_basics'];
    
epsvar = ['s_disp_' drp_instr{c}];
Nvar = ['s_ntot_' drp_instr{c}];
    
nday = length(clouds.(camp));

days_analyzed = 1:nday;

if c==1
    days_analyzed(days_analyzed==10 | days_analyzed==13 | days_analyzed==11)=[];
elseif c==3
    days_analyzed(ismember(days_analyzed, [2,3,5,14,15,16,17,18]))=[];
end
%%
close all
layers = 4;
layer_seg = linspace(0, 1.2, layers+1);
% layer_seg = [0 0.3 0.8 1.2];

max_pay_gap = 3;
samp_frac = 0.5;
% rand_sampled = zeros(length(days_analyzed), 1);
allb = [];
allpred = [];
allrsq = [];
all_mean_theta = [];
all_mean_ap = [];
% all_mean_

for iday = days_analyzed
    %%
    s_t = floor(clouds.(camp)(iday).s_t);
%     a_t = clouds.(camp)(iday).a_t;
    
%     [~, cmt_ipdi{c,iday}, cmt_ipcasp{c,iday}] = intersect(s_t,a_t);

    reldisp_raw = clouds.(camp)(iday).(epsvar);
    s_ntot_raw = clouds.(camp)(iday).(Nvar);
    s_ap_raw = clouds.(camp)(iday).s_ap;
    normAC_raw = clouds.(camp)(iday).normAC;
    s_actfrac_raw = clouds.(camp)(iday).s_actfrac;
    s_rh_raw = clouds.(camp)(iday).s_rh;
    s_thet_raw = clouds.(camp)(iday).s_thet;
    s_wz_raw = clouds.(camp)(iday).s_wz;
    
    cloudlegs_i = fb.(fbvar)(iday).ti;
    cloudlegs_f = fb.(fbvar)(iday).tf;
    %%
    if ~isempty(cloudlegs_i)
        %%
        for ileg = 1:length(cloudlegs_i)
            ti = cloudlegs_i(ileg);
            tf = cloudlegs_f(ileg);
            for ilayer = 1:length(layer_seg)-1
                %%
                filt_crit = s_t > ti & s_t < tf & normAC_raw > layer_seg(ilayer) & ...
                    normAC_raw < layer_seg(ilayer+1) & s_ntot_raw > 25;
                num_vdp{iday,ileg}(ilayer) = sum(filt_crit);
                idx_vdp{iday,ileg,ilayer} = find(filt_crit);
            end
            
            % subsampling
            vdp_idx = s_t > ti & s_t < tf & s_ntot_raw > 25 & normAC_raw > layer_seg(1) & ...
                normAC_raw < layer_seg(end);
            total_vdp = sum(vdp_idx);
            max_per_bin = round(total_vdp/layers);
            sampsize_per_bin = round(num_vdp{iday,ileg}*samp_frac);
            % set to nan before downsampling
            sampsize_per_bin(sampsize_per_bin<0.05*max_per_bin) = nan; 
            sampsize_per_bin(sampsize_per_bin > max_pay_gap*min(sampsize_per_bin)) = ...
                round(max_pay_gap*min(sampsize_per_bin));
            sampsize_per_bin(isnan(sampsize_per_bin)) = 0; % to avoid problems afterwards

            rand_sampled{iday,ileg} = [];
            % random selection from the subsampled dataset
            for ilayer = 1:length(layer_seg)-1
                % random sampling index within the layer

                rand_samp_idx_layer = randsample(idx_vdp{iday,ileg,ilayer}, sampsize_per_bin(ilayer));

                rand_sampled{iday,ileg} = [rand_sampled{iday,ileg};rand_samp_idx_layer];
            end
            reldisp = reldisp_raw(rand_sampled{iday,ileg});
%             a_ntot = a_ntot_raw(rand_sampled{iday,ileg});
%             a_ntot_int = a_ntot_ex_raw(rand_sampled{iday,ileg});
            s_ntot = s_ntot_raw(rand_sampled{iday,ileg}); 
            s_actfrac = s_actfrac_raw(rand_sampled{iday,ileg});
            s_thet = s_thet_raw(rand_sampled{iday,ileg});
            s_ap = s_ap_raw(rand_sampled{iday,ileg});
            s_ap(s_ap>4000)=nan;
            normAC = normAC_raw(rand_sampled{iday,ileg});
            s_rh = s_rh_raw(rand_sampled{iday,ileg});
            s_ss = (s_rh-100)/100;
            s_ss(s_ss<0) = nan;
            
            %     nbins = 30;
            color = s_thet;    
            x1 = s_actfrac-nanmean(s_actfrac);
            x2 = normAC-nanmean(normAC);
            x3 = s_thet-nanmean(s_thet);
            y = reldisp;
            
            all_mean_theta = [all_mean_theta nanmean(s_thet)];
            all_mean_ap = [all_mean_ap nanmean(s_ap)];
    % 
            X = [ones(size(x1)) x1 x2 x3];
            b = regress(y,X);

            y_hat = X*b;
            figure
            scatter(y_hat, y, [], color, '.');
            refline(1,0)
            axis equal
            
            pred = 1-sum((y_hat-y).^2)/sum((y-mean(y)).^2);
            
            stats = regstats(y, X, 'linear');
            rsq_raw = stats.rsquare;
            allb = [allb;b'];
            allpred = [allpred;pred];
            allrsq = [allrsq;rsq_raw];
        end
        
        
        
    end

end

meanb{c} = nanmean(allb)';
reldispb{c} = std(allb)'./abs(mean(allb))';
meanpred{c} = nanmean(allpred)';
meanrsq{c} = nanmean(allrsq)';
disp(meanpred{c})
disp(meanrsq{c})
disp(reldispb{c})
% meanrsq = nanmean(rsq_raw)
% bar(num_vdp(iday,:)); hold on
% plot(sampsize_per_bin,'LineWidth',2)