function avgd=mosaicify(raw, rat)
% convert raw data into segment-averaged data

if rat==1
   avgd=raw;
   return
end

nraw=max(size(raw)); % the longer dimension is the length by default
navg=floor(nraw/rat);

avgd=zeros(navg,min(size(raw)));

for iavg=1:navg
   idxs=(iavg-1)*rat+1;
   idxe=iavg*rat;
   
   avgd(iavg,:)=nanmean(raw(idxs:idxe,:));
end

end