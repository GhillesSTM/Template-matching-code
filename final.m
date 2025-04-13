function mask_defects_in_stm_image()
    % Load the STM image and template from .txt files
       image = dlmread('AFC-ZrTe3- Hf Doped-91.8k_0023.txt');
    template = dlmread('EDD-ZrTe3- Hf Doped-91.8k_0023.txt');
    
% Get the size of the template
[template_height, template_width] = size(template);

% Display the dimensions
disp(['Template size: ', num2str(template_height), 'x', num2str(template_width)]);
    % Perform template matching
    correlation_output = normxcorr2(template, image);
    threshold = 0.55;
    [yPeak, xPeak] = find(correlation_output >= threshold);

    % Create a mask for the detected defects
    mask = zeros(size(image));
    template_size = size(template);
    used_region = zeros(size(image)); % To track used regions
    
    for i = 1:length(yPeak)
        y = yPeak(i) - template_size(1) + 1;
        x = xPeak(i) - template_size(2) + 1;
        if y > 0 && x > 0 && y + template_size(1) - 1 <= size(image, 1) && x + template_size(2) - 1 <= size(image, 2)
            % Check if the region is already used
            if sum(used_region(y:y + template_size(1) - 1, x:x + template_size(2) - 1), 'all') == 0
                mask(y:y + template_size(1) - 1, x:x + template_size(2) - 1) = 1;
                used_region(y:y + template_size(1) - 1, x:x + template_size(2) - 1) = 1;
 
        end
        end
    end
  raw_adjusted = imadjust(mat2gray(image), [0.5 0.95], []); % Enhance contrast
    % Remove small regions or noise
    min_size = 1;  % Minimum size of a defect (area)
    mask = bwareaopen(mask, min_size);
    masked_image = raw_adjusted;
    masked_image(mask == 1) = NaN;  % Masked defects set to NaN or another appropriate value

    % Find the connected components in the mask
    CC = bwconncomp(mask);
     
    % Find the centroids of the defects
    properties = regionprops(CC, 'Centroid');
    centroids = cat(1, properties.Centroid);

    % Display results
   
    % imshow(raw_adjusted);
    
    figure;
    imshow(image, []);
    colormap('gray');
    title('Original Image');
    
    figure;
   % imshow(template, []);
    title('Template');
    
    figure;
    imshow(raw_adjusted, []);
    title('Masked Image');
%}
    hold on;
      % Overlay the transparent mask
    red = cat(3, zeros(size(mask)), ones(size(mask)), zeros(size(mask))); % Create a red color mask
    h = imshow(red);
    set(h, 'AlphaData', mask * 0.2); % Adjust the transparency of the mask
    % Plot the centroids on the masked image
    if ~isempty(centroids)
        plot(centroids(:, 1), centroids(:, 2), 'r.', 'MarkerSize', 15);
        title('Masked Image with Defect Centers');
    end
    
    hold off;
%}
    % Display the centroids coordinates
    disp('Centroids of the detected defects:');
    disp(centroids);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the number of defects
    num_defects = size(centroids, 1);
    disp('Number of defects:');
    disp(num_defects);

    % Calculate the density of defects
    density = (num_defects / 17920) * 100;
    disp('Density of defects:');
    disp(density);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create a matrix to represent the center of defects in matrix shape
    % Create a matrix to represent the clusters
    PixNum1 = size(mask, 1);
    PixNum2 = size(mask, 2);
   % Initialize cluster matrix
cluster_matrix = zeros(PixNum1, PixNum2);
%%%%%%%%% CIRCLE AROUND CENTER OF DEFECTS %%%%%%%%
%
% Define the radius around each centroid
radius =1; % You can adjust this value as needed

% Set elements corresponding to cluster centroids to 1 within the radius
for i = 1:size(centroids, 1)
    x = round(centroids(i, 1));
    y = round(centroids(i, 2));
    
    % Check if the coordinates are within the matrix bounds
    for dx = -radius:radius
        for dy = -radius:radius
            if sqrt(dx^2 + dy^2) <= radius
                nx = x + dx;
                ny = y + dy;
                if nx >= 1 && nx <= PixNum2 && ny >= 1 && ny <= PixNum1
                    cluster_matrix(ny, nx) = 1;
                end
            end
        end
    end
end

%{
%%%%%%%%%% RECTANGLE TO CAPTURE ALL THE DEFECTS%%%%%%
% Define the half-size of the region around each centroid for X and Y
x_half_size = 5; % Half-width for the X direction
y_half_size = 39;  % Half-height for the Y direction

% Set elements corresponding to defect centroids within the specified X and Y ranges
for i = 1:size(centroids, 1)
    x = round(centroids(i, 1));
    y = round(centroids(i, 2));
    
    % Check if the coordinates are within the matrix bounds
    for dx = -x_half_size:x_half_size
        for dy = -y_half_size:y_half_size
            nx = x + dx;
            ny = y + dy;
            if nx >= 1 && nx <= PixNum2 && ny >= 1 && ny <= PixNum1
                cluster_matrix(ny, nx) = 1;
            end
        end
    end
end
%}
% Plotting
hold on;
axis equal;

% Assuming 'find_maxima' is a function to find maxima in a matrix
def = double(find_maxima(cluster_matrix, 100000, 21, 1));
plot(def(:, 1), def(:, 2), 'r.', 'MarkerSize', 5);
%plot(CDW(:, 1), CDW(:, 2), 'g.', 'MarkerSize', 10);

title('center-defects+CDWmax');

axis([0 PixNum2 0 PixNum1]);
hold off;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Cross correlation with the CDW %%%%%%
    CDW = dlmread('CDW-image-ZrTe3- Hf Doped-91.8k_0023.txt');
    figure;
    imshow(CDW);
    C = normxcorr2(CDW,cluster_matrix);

    M = max(C, [], "all"); % Max cross correlation
    disp('Max cross-correlation defects with CDW:');
    disp(M);
    disp('Number of defects:');
    disp(num_defects);
    disp('Density of defects:');
    disp(density);
    display('threshold%:');
    display(threshold * 100);
    % Display the cropped cross-correlation map
    % figure;
    % imagesc(C);
    colorbar;
    axis image;
    caxis([0 max(C(:))]); % Adjust the color scale if needed
    title('Cropped Cross-Correlation Map');
end
