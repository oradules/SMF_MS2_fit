function str=strrephi(str,idel,npars,nparinv)
for i=nparinv+1:npars
   if ~ismember(i,idel)   
        str=strrep(str,['k',num2str(i)],['k',num2str(i),'hi']);
   end 
end
str=strrep(str,'k1hi0','k10hi');
str=strrep(str,'L1','L1hi');
str=strrep(str,'L2','L2hi');
str=strrep(str,'L3','L3hi');
str=strrep(str,'L4','L4hi');
str=strrep(str,'L5','L5hi');