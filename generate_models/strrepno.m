function str=strrepno(str,idel,npars,nparinv)
for i=nparinv+1:npars
   if ~ismember(i,idel)   
        str=strrep(str,['k',num2str(i)],['k',num2str(i),'no']);
   end 
end
str=strrep(str,'k1no0','k10no');
str=strrep(str,'L1','L1no');
str=strrep(str,'L2','L2no');
str=strrep(str,'L3','L3no');
str=strrep(str,'L4','L4no');
str=strrep(str,'L5','L5no');