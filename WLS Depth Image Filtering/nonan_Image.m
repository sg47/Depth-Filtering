function [ out ] = nonan_Image( img )
out = zeros(size(img));
for k=1:size(img,3)
    for i=1:size(img,1)
        for j=1:size(img,2)
            if ( isnan(img(i,j,k))) out(i,j,k)=0; 
            else out(i,j,k) = img(i,j,k); end
        end
    end
end

end

