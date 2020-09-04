function [scale_factor] = determine_scale_factor(ifos,scale_factors,ifopair)
for mm=1:length(ifos)
  for kk=1:length(ifos)
    if strcmp(ifopair,[ifos{mm},ifos{kk}])
      scale_factor = scale_factors(mm)*scale_factors(kk);
    end
  end
end

try 
   scale_factor;
catch 
   scale_factor = 1;
   fprintf('Scale factor was not set. Setting it to 1...\n');
end
