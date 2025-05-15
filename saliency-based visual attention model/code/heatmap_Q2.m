
A=ETdata.pos(:,4:5) + 0.5;
A=round(A);

backimag = img; 
X = A(:,1); 
Y = A(:,2); 
W = ones(size(X)); 
heatmap = mat2cell(hot,256,[1 1 1]);
sigma = 16;

[dimY, dimX, ~] = size(backimag);

% Create "mask"
mask = zeros(dimY, dimX);
for k = 1:1:length(X)  
    if(Y(k)<=768)&&(Y(k)>0)&&(X(k)<=1024)&&(X(k)>0)
        mask(Y(k), X(k)) = mask(Y(k), X(k)) + W(k);
    end
end

% Filter using a gaussian kernel
mask = imgaussfilt(mask, sigma);

% Normalise total mass of heatmap
% Here we equal the total mass with the one of one rectangle out of a 8x8 grid
% Decrease (increase) the constant to make the heatmap hotter (colder)
normConstant = 8;
normMask = dimX * dimY / normConstant^2 / sum(mask , [1 2]);
mask = mask * normMask;

% Colour the background image with the heatmap
newImage = backimag;
for rgbInd = 1:3
    a = floor(mask*255) + 1;
    temp = find(a>256);
    a(temp) = 256;
    thisHeat = heatmap{rgbInd}( a );
    newImage(:,:,rgbInd) = (newImage(:,:,rgbInd) + uint8(thisHeat*255));
end
figure; imshow(newImage)
