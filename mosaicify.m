function avgd=mosaicify(raw, rat)
% convert raw data into segment-averaged data

nraw=length(raw);
navg=floor(nraw/rat);

avgd=zeros(navg,1);

for iavg=1:navg
   idxs=(iavg-1)*rat+1;
   idxe=iavg*rat;
   
   avgd(iavg)=nanmean(raw(idxs:idxe));
end

end