close all
clearvars -except clouds
if ~exist('clouds','var') load clouds.mat, end
campaigns={'vocalspdi','masepdi','postpdi','oraclespdi','gomaccspdi'};
camp_proper = {'VOCALS','MASE','POST','ORACLES','GoMACCS'};
nc=length(campaigns);

%% PCR
demean = @(x) (x - nanmean(x));
zscore = @(x) (x - nanmean(x))./nanstd(x);
[reldisp,lwc,actfrac,ntot_aer,ntot_pdi,ntot_pcasp,AF,NH,ETqt,ETT,w]=deal([]);

for c=4%1:nc
   camp=campaigns{c};
   
   for iday=1:length(clouds.(camp))

      reldisp_raw=clouds.(camp)(iday).s_disp_pdi;
      vidx=~isnan(reldisp_raw);
      reldisp=[reldisp;reldisp_raw(vidx)];

      lwc=[lwc;clouds.(camp)(iday).s_lwc_pdi(vidx)];
      actfrac=[actfrac;clouds.(camp)(iday).s_actfrac(vidx)];
      ntot_pdi=[ntot_pdi;clouds.(camp)(iday).s_ntot_pdi(vidx)];
      ntot_aer=[ntot_aer;clouds.(camp)(iday).s_ntot_aer(vidx)];
      ntot_pcasp=[ntot_pcasp;clouds.(camp)(iday).s_ntot_pcasp(vidx)];
      AF=[AF;clouds.(camp)(iday).AF(vidx)];
      NH=[NH;clouds.(camp)(iday).normAC(vidx)];
      ETqt=[ETqt;clouds.(camp)(iday).ent_ratio_qt(vidx)];
      ETT=[ETT;clouds.(camp)(iday).ent_ratio_T(vidx)];
%       w=[w;clouds.(camp)(iday).s_wz(vidx)];

   end
end

YY=reldisp;
XX=[lwc actfrac AF NH];

% XX(XX<0)=0;
% XXlog=log10(XX);

XX_demean=zscore(XX);

[PCAcoeff,~,PCAlatt]=pca(XX_demean);
PCAexp=PCAlatt/sum(PCAlatt);
% nfeat=3;
nfeat=size(XX,2);

XX_new=[ones(size(XX,1),1) XX*PCAcoeff(:,1:nfeat)];
[b,bint]=regress(reldisp,XX_new);
Yhat=XX_new*b;
plot(Yhat,reldisp,'.')
xlim([0 1])
ylim([0 1])
refline(1,0)