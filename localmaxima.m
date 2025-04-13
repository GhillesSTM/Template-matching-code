I1 = dlmread('RIEBD-AFC-ZrTe3 - Hf Dopped-9k_0892 (1).txt');

I1 = flipud(I1);
PixNum1 = size(I1, 1);
PixNum2 = size(I1, 2);


% Define the size of the square to search for maxima
r = 11;

% Find the high-intensity pixels and remove edge points

EBD = double(find_maxima(I1, 1000000, r, 1));


% Define the size of the rectangle to search for center within
dx = 80; % Set dx size
dy = 80; % Set dy size

% Initialize variables to store centroid points
centroid_EBD = [];

% Iterate over the image with steps of dx and dy
for x = 1:dx:PixNum2
    for y = 1:dy:PixNum1
        % Calculate the boundaries of the current segment
        x1 = max(1, x);
        x2 = min(PixNum2, x + dx - 1);
        y1 = max(1, y);
        y2 = min(PixNum1, y + dy - 1);

        % Extract the region within the rectangle
        rect_region = I1(y1:y2, x1:x2);

        % Find indices (ii, jj) of non-zero elements in rect_region
        [ii, jj, values] = find(rect_region);

        % Calculate the centroid of non-zero values
        if ~isempty(values)
            sum_v = sum(values);
            sum_iv = sum(ii .* values);
            sum_jv = sum(jj .* values);

            % Compute the centroid for non-zero values
            centroid_x = sum_jv / sum_v + x1 - 1;
            centroid_y = sum_iv / sum_v + y1 - 1;
            centroid = [centroid_x, centroid_y];

            % Store the centroid
            centroid_EBD = [centroid_EBD; centroid];
        end
    end
end


%I2 = dlmread('point dark defects filtered-ZrTe3-fullimage- Hf Doped-91.8k_0023.txt');
I2 = flipud(I2);
PixNum12 = size(I2, 1);
PixNum22 = size(I2, 2);



%Plot the  points
figure;
hold on;
axis equal;
plot(EBD(:, 1), EBD(:, 2), 'k.', 'MarkerSize', 5);
plot(centroid_EBD(:, 1), centroid_EBD(:, 2), 'g.', 'MarkerSize', 10);
%plot(EDD(:, 1), EDD(:, 2), 'k.', 'MarkerSize', 30);
%plot(CENTER_EDD(:, 1), CENTER_EDD(:, 2), 'r.', 'MarkerSize', 10);
axis([0 PixNum2 - 1 0 PixNum1 - 1]);
hold off;

%{
Rotate the plot by 20 degrees
angle_degrees = 15;
angle_radians = deg2rad(angle_degrees);
rotation_matrix = [cos(angle_radians), -sin(angle_radians); sin(angle_radians), cos(angle_radians)];
rotated_bright_p = bright_p * rotation_matrix;
rotated_center_points = center_points * rotation_matrix;
% Plot the rotated points
%plot(rotated_bright_p(:, 1), rotated_bright_p(:, 2), 'k.', 'MarkerSize', 5);
plot(rotated_center_points(:, 1), rotated_center_points(:, 2), 'r.', 'MarkerSize', 10);
plot(rotated_center_points1(:, 1), rotated_center_points1(:, 2), 'g.', 'MarkerSize', 10);
%}