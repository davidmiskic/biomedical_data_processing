function MHdetector( record )
% prepare filter
threshold = 100;
im = imread(record);
sigma=5;
filterdim = ceil(sigma*4*sqrt(2));
kernel = fspecial('log', filterdim, sigma);
kernel = kernel - mean(kernel(:));
kernel = kernel/sum(kernel(:));
%mesh(kernel);

% filter image
dslice = im2double(im);
fslice = conv2(dslice, kernel, 'same');

% zero crossing
% 3x3 window, check for sign change in each diagonal out of 4 - if signs at each end differ then it is an edge
[h, w] = size(im);
zero_detect = zeros(h,w);
for i=1:h
    for j=1:w
        if i > 1 && j > 1 && i < h && j < w
            if (~isequal(sign(fslice(i-1,j-1)), sign(fslice(i+1,j+1))) && abs(fslice(i-1,j-1)-fslice(i+1,j+1)) > threshold) || ...
               (~isequal(sign(fslice(i-1,j+1)), sign(fslice(i+1,j-1))) && abs(fslice(i-1,j+1)-fslice(i+1,j-1)) > threshold) ||...
               (~isequal(sign(fslice(i-1,j)), sign(fslice(i+1,j))) && abs(fslice(i-1,j)-fslice(i+1,j)) > threshold) || ...
               (~isequal(sign(fslice(i,j+1)), sign(fslice(i,j-1))) && abs(fslice(i,j+1)-fslice(i,j-1)) > threshold)
                zero_detect(i,j) = 1;
            end
        end
    end
end

% binimage = imbinarize(fslice, 'adaptive');
binimage = imbinarize(zero_detect, 'global');

% matlabim = edge(im, 'log');
binimagebefore = binimage(:,:);
% edge linking 
max_line = 50;
for i=1:h
    for j=1:w
       if binimage(i,j)==1 && binimage(i, j+1) == 0
           gap_end = j;
           while binimage(i, gap_end)==0 && gap_end < w
               gap_end = gap_end + 1;
           end
           if gap_end < max_line && binimage(i, gap_end) == 1
               disp("GAP FILL");
               while gap_end ~= j
                    binimage(i, gap_end) = 1;
                    gap_end = gap_end - 1;
               end
           end   
       end
    end
end

for j=1:w
    for i=1:h
       if binimage(i,j)==1 && binimage(i+1, j) == 0
           gap_end = i;
           while binimage(gap_end, j)==0  && gap_end < h
               gap_end = gap_end + 1;
           end
           if gap_end < max_line 
               while gap_end ~= i && binimage(gap_end, j) == 1
                   disp("GAP FILL");
                   binimage(gap_end, j) = 1;
                   gap_end = gap_end - 1;
               end
           end   
       end
    end
end
name = split(record, ".");
imname = name(1)+"out_image.png";
imwrite(binimage, imname);
figure;
  subplot(2,2,1);imshow(im);title('Original image');
  subplot(2,2,2);imshow(binimagebefore);title('Binarized image');
  subplot(2,2,3);imshow(binimage);title('Final image');
  subplot(2,2,4);imshow(edge(im, 'log'));title('MATLAB auto image');

saveas(gcf, name(1) + "out.png");
