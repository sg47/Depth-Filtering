	

    % Bilateral Filter for Matlab/Octave
     
    % Author: Barak Itkin
    % Author Website: http://barak-itkin.blogspot.com/
    % Created: 2011-09-08
     
    % This work is licensed under a Creative Commons Attribution-ShareAlike
    % 3.0 Unported License.
    % See http://creativecommons.org/licenses/by-sa/3.0/
     
    % The bilateral filter is a 'Smart Blur' filter that avoids smoothing
    % edges. It does so by considering the intensity differences between the
    % pixels - if the difference between the 'center pixel' and another
    % pixel is too high, then that pixel get's less weight in the weight
    % matrix of the blur. In addition to considering the intensity
    % difference, we also consider the actual distance (as in a regular
    % Gaussian blur).
    % The bilateral filter is most commonly used as a noise reduction filter
    % due to it's ability to smooth large areas without destroying edges.
    % For more information, see http://en.wikipedia.org/wiki/Bilateral_filter
    %
    % The equation of the bilateral filter is
    %
    %         (       dx ^ 2       )       (         dI ^2        )
    % F = exp (- ----------------- ) * exp (- ------------------- )
    %         (  sigma_spatial ^ 2 )       (  sigma_Intensity ^ 2 )
    %     ~~~~~~~~~~~~~~~~~~~~~~~~~~
    %     This is a guassian filter!
    %
    % dx - The 'geometric' distance between the 'center pixel' and the pixel
    %      to sample
    % dI - The difference between the intensity of the 'center pixel' and
    %      the pixel to sample
    % sigma_spatial and sigma_Intesity are constants. Higher values mean
    % that we 'tolerate more' higher value of the distances dx and dI.
     
    % Usage:
    %   im         - A GRAY image
    %   filterSize - The size of the filter (the lookup area would be
    %                filterSize X filterSize). This is a real positive
    %                integer. This MUST BE AN ODD (NON EVEN) REAL NUMBER!
    %   sSpatial   - The sigma_spatial of the bilateral filter equation.
    %                Typically, this is half of the filterSize.
    %   SIntensity - The sigma_Intensity of the bilateral filter equation.
    %
    % Return:
    %   ret        - The resulting image.
     
    function [ ret ] = bilateral (im, filterSize, sSpatial, sIntensity)
     
    if (not (nargin == 4))
            error ('Expected exactly 4 arguments!')
    end
     
    if (not (ismatrix(im)) || not (ndims (im) == 2))
            error ('Argument 1 (im) must be a GRAY image')
    end
     
    if (not (isscalar (filterSize)) || filterSize <= 0 || not (isreal (filterSize)) || not (rem (filterSize, 2) == 1))
            error ('Argument 2 (filterSize) must be a positive real odd integer!')
    end
     
    if (not (isscalar (sSpatial)) || not (isreal (sSpatial)))
            error ('Argument 3 (sSpatial) must be a real number!')
    end
     
    if (not (isscalar (sIntensity)) || not (isreal (sIntensity)))
            error ('Argument 4 (sIntensity must be a real number!')
    end
     
    % Compute the Gaussian filter part of the Bilateral filter
    gauss = fspecial ('gaussian', filterSize, sSpatial);
     
    % The radius of the filter
    radius = floor (filterSize / 2);
     
    % The original image dimensions
    [height,width] = size (im);
     
    % The resulting image
    ret = zeros (size (im));
     
    % We will now compute the result by iterating over the result image,
    % computing one pixel at a time. For the pixel at (i,j) we want to look
    % at:
    %
    %      im((i - radius):(i + radius), (j - radius):(j + radius))
    %
    % But, we we must also remember to restrict ourselves to remain inside
    % the image boundries
    for i=1:height
     
            % The top part of the lookup area, clamped to the image
            ymin  = max (i - radius, 1);
            % How many rows were outside the image, on the top?
            dymin = ymin - (i - radius);
           
            % The bottom part of the lookup area, clamped to the image
            ymax  = min (i + radius, height);
            % How many rows were outside the image, on the bottom?
            dymax = (i + radius) - ymax;
           
            for j=1:width
           
                    % The left part of the lookup area, clamped to the image
                    xmin = max (j - radius, 1);
                    % How many columns were outside the image, on the left?
                    dxmin = xmin - (j - radius);
                   
                    % The right part of the lookup area, clamped to the image
                    xmax = min (j + radius, width);
                    % How many columns were outside the image, on the right?
                    dxmax = (j + radius) - xmax;
                   
                    % The actual area of the image we will look at
                    area = im (ymin:ymax, xmin:xmax);
                   
                    % The center pixel
                    center = im (i,j);
     
                    % The left expression in the bilateral filter equation
                    % We take only the relevant parts of the matrix of the
                    % Gaussian weights - we use dxmin, dxmax, dymin, dymax to
                    % ignore the parts that are outside the image
                    expS = gauss((1+dymin):(filterSize-dymax),(1+dxmin):(filterSize-dxmax));
     
                    % The intensity difference from the center pixel
                    dI = area - center;
                    % The right expression in the bilateral filter equation
                    expI = exp((-dI .* dI) / (sIntensity * sIntensity));
                   
                    % The bilater filter (weights matrix)
                    F = expI .* expS;
                    % Normalized bilateral filter
                    Fnormalized = F / sum(F(:));
                   
                    % Multiply the area by the filter
                    temp = area .* Fnormalized;
                   
                    % The resulting pixel is the sum of all the pixels in
                    % the area, according to the weights of the filter
                    ret(i,j) = sum (temp(:));
            end
    end

