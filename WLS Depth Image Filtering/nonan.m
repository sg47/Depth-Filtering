function [ out_XYZ ] = nonan( XYZ )
[ N_Scan, N_Point]=size(XYZ);
for i=1:N_Point
    for j=1:N_Scan
     if (isnan(XYZ(j,i))) XYZ(j,i)=0;
     end
    end
end
out_XYZ=XYZ;
end

