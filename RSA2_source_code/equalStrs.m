function n = equalStrs(str1,str2)

while ~strcmp(str1,str2)
    str1=str1(1:end-1); str2=str2(1:end-1);
end

n=length(str1);
    
