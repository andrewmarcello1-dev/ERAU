clear; clc; close all

%% Coefficients
% Axis 1
Beta  = 2*sqrt(3)/3 - 1;
Kappa = -(2*sqrt(3)/3 + 1);
R     = -4/(sqrt(2)*((3/2 - sqrt(3))*Beta + (sqrt(3) - 3/2)*Kappa - 2*sqrt(3)));

% Axis 2
Beta2 = -(2*sqrt(3)/3+1);
Kappa2 = 2*sqrt(3)/3-1;
R2 = -4/(sqrt(2)*(Beta2*(-3/2-sqrt(3))+Kappa2*(-sqrt(3)+3/2)));

fprintf('==========Coefficients (Axis 1)==========\n')
fprintf('Beta (a): %f\n',Beta)
fprintf('Kappa (b): %f\n',Kappa)
fprintf('R (c): %f\n\n',R)

fprintf('==========Coefficients (Axis 2)==========\n')
fprintf('Beta (a): %f\n',Beta)
fprintf('Kappa (b): %f\n',Kappa)
fprintf('R (c): %f\n\n',R)

fprintf('Displaying Data...\n\n')

%% Shared angle grid
N     = 41;
alpha = linspace(-20, 20, N);
i0    = round((N+1)/2);   % index where alpha ~ 0

% Sequences:
idx1 = [1:N, N-1:-1:i0];  % Phase 1: -20 -> 20 -> 0
idx2 = [i0:N, N-1:-1:1];  % Phase 2:  0  -> 20 -> -20

%% Geometry for both axes
% Axis 1
c1 = R * tand(alpha);  a1 = Beta  * c1;  b1 = Kappa  * c1;
P1a = zeros(N,3); P2a = P1a; P3a = P1a;
P1a(:,1) = -1; P1a(:,3) = a1;
P2a(:,1) =  1; P2a(:,3) = b1;
P3a(:,2) =  sqrt(3); P3a(:,3) = c1;

% Axis 2
c2 = R2 * tand(alpha); a2 = Beta2 * c2; b2 = Kappa2 * c2;
P1b = zeros(N,3); P2b = P1b; P3b = P1b;
P1b(:,1) = -1; P1b(:,3) = a2;
P2b(:,1) =  1; P2b(:,3) = b2;
P3b(:,2) =  sqrt(3); P3b(:,3) = c2;

% Combined axis limits (stable view)
allXYZ = [P1a; P2a; P3a; P1b; P2b; P3b];
xl = [min(allXYZ(:,1)) max(allXYZ(:,1))];
yl = [min(allXYZ(:,2)) max(allXYZ(:,2))];
zl = [min(allXYZ(:,3)) max(allXYZ(:,3))];

%% Figure
figure('Color','w','Position',[100 100 800 600]);
ax = gca; view(3); grid on; axis(ax,'equal');
xlabel('x'); ylabel('y'); zlabel('z'); hold on;

% Both dashed reference lines (persist throughout)
xline = linspace(xl(1)-1, xl(2)+1, 500);
y_line1 =  xline + sqrt(3)/2;   % Axis 1 line
y_line2 = -xline + sqrt(3)/2;   % Axis 2 line
zline   = zeros(size(xline));
plot3(xline, y_line1, zline, 'k--', 'LineWidth', 1.2);
plot3(xline, y_line2, zline, 'k--', 'LineWidth', 1.2);

% Fixed reference point
plot3(0, sqrt(3)/2, 0, 'o', 'MarkerSize', 8, ...
      'MarkerEdgeColor','k', 'MarkerFaceColor',[1 0 0]);

% Final axis limits (ensure lines are visible)
yl_final = [min([yl(1), min(y_line1), min(y_line2)]) , max([yl(2), max(y_line1), max(y_line2)])];
xlim(xl); ylim(yl_final); zlim(zl);

% Pre-build line point arrays for distance calc
Qline1 = [xline(:),  y_line1(:), zline(:)];
Qline2 = [xline(:),  y_line2(:), zline(:)];

% Triangle patch (reused)
V0 = [P1a(1,:); P2a(1,:); P3a(1,:)];
htri = patch('Faces',[1 2 3],'Vertices',V0, ...
             'FaceColor',[0.30 0.65 1.00],'FaceAlpha',0.6, ...
             'EdgeColor','k','LineWidth',1.5);

%% GIF setup
gifFile   = 'Tilt_Both_Axes_inPlane.gif';
delayTime = 0.05;
isFirst   = true;  % first frame flag

%% Phase 1: use Axis 1 triangle, sweep -20->20->0
for k = 1:numel(idx1)
    r = idx1(k);

    V = [P1a(r,:); P2a(r,:); P3a(r,:)];
    set(htri, 'Vertices', V);

    % Distances to BOTH lines (report both)
    n  = cross(V(2,:)-V(1,:), V(3,:)-V(1,:));
    nn = norm(n);
    if nn > 0
        d1 = max( abs((Qline1 - V(1,:)) * n.') / nn );  % to y = x + sqrt(3)/2
        d2 = max( abs((Qline2 - V(1,:)) * n.') / nn );  % to y = -x + sqrt(3)/2
    else
        d1 = NaN; d2 = NaN;
    end

    title(sprintf('Phase 1 (Axis 1) | \\alpha = %.2f^\\circ  |  max dist: line1=%.4f, line2=%.4f', ...
          alpha(r), d1, d2));

    drawnow;
    writeGifFrame(gifFile, delayTime, isFirst);  % first call initializes palette
    isFirst = false;
end

%% Phase 2: use Axis 2 triangle, sweep 0->20->-20
for k = 1:numel(idx2)
    r = idx2(k);

    V = [P1b(r,:); P2b(r,:); P3b(r,:)];
    set(htri, 'Vertices', V);

    % Distances to BOTH lines
    n  = cross(V(2,:)-V(1,:), V(3,:)-V(1,:));
    nn = norm(n);
    if nn > 0
        d1 = max( abs((Qline1 - V(1,:)) * n.') / nn );
        d2 = max( abs((Qline2 - V(1,:)) * n.') / nn );
    else
        d1 = NaN; d2 = NaN;
    end

    title(sprintf('Phase 2 (Axis 2) | \\alpha = %.2f^\\circ  |  max dist: line1=%.4f, line2=%.4f', ...
          alpha(r), d1, d2));

    drawnow;
    writeGifFrame(gifFile, delayTime, false);   % append frames
end

disp(['Saved GIF: ', gifFile]);

%% ---- local function (place after the main script code) ----
function writeGifFrame(gifFile, delayTime, firstFrame)
%WRITEGIFFRAME Capture current figure and append to a GIF with fixed colormap.
%   writeGifFrame(gifFile, delayTime, firstFrame)
%   - gifFile:   output GIF filename
%   - delayTime: seconds per frame
%   - firstFrame (logical): true for the very first frame; false afterwards
persistent map
frame = getframe(gcf);
im    = frame2im(frame);
if firstFrame || isempty(map)
    [A,map] = rgb2ind(im, 256);  % initialize fixed palette
    imwrite(A, map, gifFile, 'gif', 'LoopCount', Inf, 'DelayTime', delayTime);
else
    A = rgb2ind(im, map);        % reuse palette to avoid color shift
    imwrite(A, map, gifFile, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
end
end
