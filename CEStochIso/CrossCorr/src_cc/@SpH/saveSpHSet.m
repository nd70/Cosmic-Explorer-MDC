function vSpH=saveSpHSet(vSpH,terminate)

try
  terminate;
catch
  terminate = true;
end

numOfNewSegs=vSpH.out.counter-vSpH.setOffset;
if numOfNewSegs>0
  Sky=vSpH.out.data;
  save(vSpH.currentFilename,'Sky');
  clear Sky;
  if ~terminate
    vSpH.out.data={};
    vSpH.setOffset=vSpH.setOffset+numOfNewSegs;
    vSpH.setNumber=vSpH.setNumber+1;
    vSpH.currentFilename=[vSpH.setPrefix,num2str(vSpH.setNumber),vSpH.suffix];
  end
end

if and(terminate,vSpH.out.counter>0)
  Sky=vSpH.out.meta;
  save([vSpH.prefix,vSpH.suffix],'Sky');
  clear Sky;
end
