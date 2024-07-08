function [] = Plot_density(x, y, Nx, Ny, varargin)

%% Default values
xmin = min(x, [], 'all');
xmax = max(x, [], 'all');
ymin = min(y, [], 'all');
ymax = max(y, [], 'all');
color = [0, 0, 0];
linecolor = [1, 1, 1];
linewdith = 2;

%% Check that x and y are the same size
if ~all(size(x) == size(y))
    error('ERROR: Plot_density requires that ''x'' and ''y'' be the same size!')
end

%% Check ''varargin''
n_varargin = length(varargin);
n_varargin_name_and_value = floor(n_varargin/2);
if 2*n_varargin_name_and_value ~= n_varargin
    error('ERROR: varargin inputs should have an even total number of inputs: each name must has a value!')
end

for i = 1:n_varargin_name_and_value
   
    switch varargin{1+2*(i-1)}
        case 'xlim'
            if length(varargin{2*i}) == 2
                xmin = varargin{2*i}(1);
                xmax = varargin{2*i}(2);
            else
                error('ERROR: ''xlim'' must be a 2-element array!')
            end
        case 'ylim'
            if length(varargin{2*i}) == 2
                ymin = varargin{2*i}(1);
                ymax = varargin{2*i}(2);
            else
                error('ERROR: ''ylim'' must be a 2-element array!')
            end
        case 'color'
            color = varargin{2*i};
        case 'linecolor'
            linecolor = varargin{2*i};
        case 'linewidth'
            linewdith = varargin{2*i};
        otherwise
            if isstring(varargin{1+2*(i-1)}) || ischar(varargin{1+2*(i-1)})
                error(['ERROR: unsupported varargin name: ''', varargin{1+2*(i-1)}, '''', newline, 'Acceptable varagin names are ''xlim'', ''ylim'', ''color'', ''linecolor'', and ''linewidth'''])
            else
                error('ERROR: odd-numbered varargin inputs must be strings or characters!')
            end
    end
    
end

%% Dimensions and coordinates of pixels
dx = (xmax - xmin)/Nx;
dy = (ymax - ymin)/Ny; 
xcenter = linspace(xmin+dx/2, xmax-dx/2, Nx);
xcenter = repmat(xcenter, Ny, 1); %x-coordinate of center of square
ycenter = linspace(ymin+dy/2, ymax-dy/2, Ny)';
ycenter = repmat(ycenter, 1, Nx); %y-coordinate of center of square
dx_square = dx * [-1, 1, 1, -1, -1] / 2;
dy_square = dy * [1, 1, -1, -1, 1] / 2;

%% Count the number of x-y pairs within each square
count = nan(Ny, Nx); %number of x-y pairs within each square
for i = 1:Ny
    for j = 1:Nx
        count(i,j) = sum(inpolygon(x, y, xcenter(i,j) + dx_square, ycenter(i,j) + dy_square), 'all');
    end
end

%% Normalize count
count = count/max(count, [], 'all');

%% Plot density
set(gca, 'Layer', 'top')
hold on
for i = 1:Ny
    for j = 1:Nx
        if count(i,j) > 0
            fill(xcenter(i,j) + dx_square, ycenter(i,j) + dy_square, color, 'FaceAlpha', count(i,j), 'LineStyle', 'none')
            plot(xcenter(i,j) + dx_square, ycenter(i,j) + dy_square, '-', 'LineWidth', linewdith, 'color', linecolor)
        end
    end
end
hold off

end

