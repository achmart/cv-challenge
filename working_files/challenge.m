%% Computer Vision Challenge

% Groupnumber:
group_number = 40;

% Groupmembers:
% members = {'Max Mustermann', 'Johannes Daten'};
members = {'Martin Achtner', 'Ridon Arifi', 'Daniel Ehrhardt', 'Lukas Fichtner', 'Thomas Holzner'};

% Email-Adress (from Moodle!):
% mail = {'ga99abc@tum.de', 'daten.hannes@tum.de'};
mail = {'martin.achtner@tum.de','ridon.arifi@tum.de','ga53coc@tum.de','lukas.fichtner@tum.de','thomas.holzner@tum.de'};

addpath('RectifKitE')
USEGUI = 0; % 1: use a gui, 0 : do not use gui
DEBUG = 1; % Intermediate figures will be shown

%% Load images
if USEGUI == 1
    % call GUI and save handles
    gui_handle = challengeGUI;
    gui_figure_handle = gui_handle.FreeViewpointGUIUIFigure;
    % wait for GUI to be closed or button press
    uiwait(gui_figure_handle);
    % terminate programm if GUI was closed
    if ~isgraphics(gui_figure_handle)
        return
    end
    img1 = gui_handle.Image1;
    img2 = gui_handle.Image2;
    p = gui_handle.p;
    img_src = gui_handle.camera;
else
    img1 = imread('./img/L2.JPG');
    img2 = imread('./img/R2.JPG');
    p = 0.70;
    img_src = 2;
end

%% Free Viewpoint Rendering
% start execution timer -> tic;
tic
%img_src = 1;
output_image = free_viewpoint(img1, img2, p, img_src, DEBUG);
toc

% stop execution timer -> toc;
elapsed_time = toc-tic;

%% Display Output
% Display Virtual View
figure
imshow(output_image);
title('Output image');