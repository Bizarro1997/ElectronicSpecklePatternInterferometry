function   max_tip = plot3D(num)
% Demo to create a surface from a gray scale image and have the coloration of the surface taken from a different RGB image.
%clc;    % Clear the command window.
%clear all;
%close all;
workspace;  % Make sure the workspace panel is showing.
format short g;
format compact;
fontSize = 15;
fprintf('Beginning to run %s.m ...\n', mfilename);
%===============================================================================
% Get the name of the demo image the user wants to use.
% Let's let the user select from a list of all the demo images that ship with the Image Processing Toolbox.
folder = './image2/plot3D'; % Determine where demo folder is (works with all versions).
% Demo images have extensions of TIF, PNG, and JPG.  Get a list of all of them.

baseFileName = ['ipSystematicErrorCorrected_124-', num2str(num),'.png'];

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Get the full filename, with path prepended.
fullFileName = fullfile(folder, baseFileName);
img = imread(fullFileName);
rgbImage = img;

textureImageBaseName = ['newNormfinal_124-', num2str(num),'.png'];
textureImageName = fullfile(folder, textureImageBaseName);
peaksImage = imread(textureImageName);
peaksImage = 21*0.543/4*double(peaksImage)/196;
max_tip = max(max(peaksImage));

surf(peaksImage, ...
'FaceColor', 'texturemap',...
'EdgeColor', 'none',...
'Cdata', rgbImage);
view(3);
axis ij;
xlabel('X', 'FontSize', fontSize);
ylabel('Y', 'FontSize', fontSize);
zlabel('Z', 'FontSize', fontSize);
%if numberOfColorChannels == 1
%colorbar;
%end
%
zlim([0,22*0.543/4]);

% Create zlabel
zlabel('{\mu}m ','FontSize',15);

% Create ylabel
ylabel('[mm]','FontSize',15);

% Create xlabel
xlabel('[mm]','FontSize',15);

xlim(axes1,[1 1280]);
ylim(axes1,[1 960]);
zlim(axes1,[0 2.9865]);
view(axes1,[-37.5 30]);
grid(axes1,'on');
axis(axes1,'ij');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'Colormap',...
[0 0 0;0.00392156862745098 0.00392156862745098 0.00392156862745098;0.00784313725490196 0.00784313725490196 0.00784313725490196;0.0117647058823529 0.0117647058823529 0.0117647058823529;0.0156862745098039 0.0156862745098039 0.0156862745098039;0.0196078431372549 0.0196078431372549 0.0196078431372549;0.0235294117647059 0.0235294117647059 0.0235294117647059;0.0274509803921569 0.0274509803921569 0.0274509803921569;0.0313725490196078 0.0313725490196078 0.0313725490196078;0.0352941176470588 0.0352941176470588 0.0352941176470588;0.0392156862745098 0.0392156862745098 0.0392156862745098;0.0431372549019608 0.0431372549019608 0.0431372549019608;0.0470588235294118 0.0470588235294118 0.0470588235294118;0.0509803921568627 0.0509803921568627 0.0509803921568627;0.0549019607843137 0.0549019607843137 0.0549019607843137;0.0588235294117647 0.0588235294117647 0.0588235294117647;0.0627450980392157 0.0627450980392157 0.0627450980392157;0.0666666666666667 0.0666666666666667 0.0666666666666667;0.0705882352941176 0.0705882352941176 0.0705882352941176;0.0745098039215686 0.0745098039215686 0.0745098039215686;0.0784313725490196 0.0784313725490196 0.0784313725490196;0.0823529411764706 0.0823529411764706 0.0823529411764706;0.0862745098039216 0.0862745098039216 0.0862745098039216;0.0901960784313725 0.0901960784313725 0.0901960784313725;0.0941176470588235 0.0941176470588235 0.0941176470588235;0.0980392156862745 0.0980392156862745 0.0980392156862745;0.101960784313725 0.101960784313725 0.101960784313725;0.105882352941176 0.105882352941176 0.105882352941176;0.109803921568627 0.109803921568627 0.109803921568627;0.113725490196078 0.113725490196078 0.113725490196078;0.117647058823529 0.117647058823529 0.117647058823529;0.12156862745098 0.12156862745098 0.12156862745098;0.125490196078431 0.125490196078431 0.125490196078431;0.129411764705882 0.129411764705882 0.129411764705882;0.133333333333333 0.133333333333333 0.133333333333333;0.137254901960784 0.137254901960784 0.137254901960784;0.141176470588235 0.141176470588235 0.141176470588235;0.145098039215686 0.145098039215686 0.145098039215686;0.149019607843137 0.149019607843137 0.149019607843137;0.152941176470588 0.152941176470588 0.152941176470588;0.156862745098039 0.156862745098039 0.156862745098039;0.16078431372549 0.16078431372549 0.16078431372549;0.164705882352941 0.164705882352941 0.164705882352941;0.168627450980392 0.168627450980392 0.168627450980392;0.172549019607843 0.172549019607843 0.172549019607843;0.176470588235294 0.176470588235294 0.176470588235294;0.180392156862745 0.180392156862745 0.180392156862745;0.184313725490196 0.184313725490196 0.184313725490196;0.188235294117647 0.188235294117647 0.188235294117647;0.192156862745098 0.192156862745098 0.192156862745098;0.196078431372549 0.196078431372549 0.196078431372549;0.2 0.2 0.2;0.203921568627451 0.203921568627451 0.203921568627451;0.207843137254902 0.207843137254902 0.207843137254902;0.211764705882353 0.211764705882353 0.211764705882353;0.215686274509804 0.215686274509804 0.215686274509804;0.219607843137255 0.219607843137255 0.219607843137255;0.223529411764706 0.223529411764706 0.223529411764706;0.227450980392157 0.227450980392157 0.227450980392157;0.231372549019608 0.231372549019608 0.231372549019608;0.235294117647059 0.235294117647059 0.235294117647059;0.23921568627451 0.23921568627451 0.23921568627451;0.243137254901961 0.243137254901961 0.243137254901961;0.247058823529412 0.247058823529412 0.247058823529412;0.250980392156863 0.250980392156863 0.250980392156863;0.254901960784314 0.254901960784314 0.254901960784314;0.258823529411765 0.258823529411765 0.258823529411765;0.262745098039216 0.262745098039216 0.262745098039216;0.266666666666667 0.266666666666667 0.266666666666667;0.270588235294118 0.270588235294118 0.270588235294118;0.274509803921569 0.274509803921569 0.274509803921569;0.27843137254902 0.27843137254902 0.27843137254902;0.282352941176471 0.282352941176471 0.282352941176471;0.286274509803922 0.286274509803922 0.286274509803922;0.290196078431373 0.290196078431373 0.290196078431373;0.294117647058824 0.294117647058824 0.294117647058824;0.298039215686275 0.298039215686275 0.298039215686275;0.301960784313725 0.301960784313725 0.301960784313725;0.305882352941176 0.305882352941176 0.305882352941176;0.309803921568627 0.309803921568627 0.309803921568627;0.313725490196078 0.313725490196078 0.313725490196078;0.317647058823529 0.317647058823529 0.317647058823529;0.32156862745098 0.32156862745098 0.32156862745098;0.325490196078431 0.325490196078431 0.325490196078431;0.329411764705882 0.329411764705882 0.329411764705882;0.333333333333333 0.333333333333333 0.333333333333333;0.337254901960784 0.337254901960784 0.337254901960784;0.341176470588235 0.341176470588235 0.341176470588235;0.345098039215686 0.345098039215686 0.345098039215686;0.349019607843137 0.349019607843137 0.349019607843137;0.352941176470588 0.352941176470588 0.352941176470588;0.356862745098039 0.356862745098039 0.356862745098039;0.36078431372549 0.36078431372549 0.36078431372549;0.364705882352941 0.364705882352941 0.364705882352941;0.368627450980392 0.368627450980392 0.368627450980392;0.372549019607843 0.372549019607843 0.372549019607843;0.376470588235294 0.376470588235294 0.376470588235294;0.380392156862745 0.380392156862745 0.380392156862745;0.384313725490196 0.384313725490196 0.384313725490196;0.388235294117647 0.388235294117647 0.388235294117647;0.392156862745098 0.392156862745098 0.392156862745098;0.396078431372549 0.396078431372549 0.396078431372549;0.4 0.4 0.4;0.403921568627451 0.403921568627451 0.403921568627451;0.407843137254902 0.407843137254902 0.407843137254902;0.411764705882353 0.411764705882353 0.411764705882353;0.415686274509804 0.415686274509804 0.415686274509804;0.419607843137255 0.419607843137255 0.419607843137255;0.423529411764706 0.423529411764706 0.423529411764706;0.427450980392157 0.427450980392157 0.427450980392157;0.431372549019608 0.431372549019608 0.431372549019608;0.435294117647059 0.435294117647059 0.435294117647059;0.43921568627451 0.43921568627451 0.43921568627451;0.443137254901961 0.443137254901961 0.443137254901961;0.447058823529412 0.447058823529412 0.447058823529412;0.450980392156863 0.450980392156863 0.450980392156863;0.454901960784314 0.454901960784314 0.454901960784314;0.458823529411765 0.458823529411765 0.458823529411765;0.462745098039216 0.462745098039216 0.462745098039216;0.466666666666667 0.466666666666667 0.466666666666667;0.470588235294118 0.470588235294118 0.470588235294118;0.474509803921569 0.474509803921569 0.474509803921569;0.47843137254902 0.47843137254902 0.47843137254902;0.482352941176471 0.482352941176471 0.482352941176471;0.486274509803922 0.486274509803922 0.486274509803922;0.490196078431373 0.490196078431373 0.490196078431373;0.494117647058824 0.494117647058824 0.494117647058824;0.498039215686275 0.498039215686275 0.498039215686275;0.501960784313725 0.501960784313725 0.501960784313725;0.505882352941176 0.505882352941176 0.505882352941176;0.509803921568627 0.509803921568627 0.509803921568627;0.513725490196078 0.513725490196078 0.513725490196078;0.517647058823529 0.517647058823529 0.517647058823529;0.52156862745098 0.52156862745098 0.52156862745098;0.525490196078431 0.525490196078431 0.525490196078431;0.529411764705882 0.529411764705882 0.529411764705882;0.533333333333333 0.533333333333333 0.533333333333333;0.537254901960784 0.537254901960784 0.537254901960784;0.541176470588235 0.541176470588235 0.541176470588235;0.545098039215686 0.545098039215686 0.545098039215686;0.549019607843137 0.549019607843137 0.549019607843137;0.552941176470588 0.552941176470588 0.552941176470588;0.556862745098039 0.556862745098039 0.556862745098039;0.56078431372549 0.56078431372549 0.56078431372549;0.564705882352941 0.564705882352941 0.564705882352941;0.568627450980392 0.568627450980392 0.568627450980392;0.572549019607843 0.572549019607843 0.572549019607843;0.576470588235294 0.576470588235294 0.576470588235294;0.580392156862745 0.580392156862745 0.580392156862745;0.584313725490196 0.584313725490196 0.584313725490196;0.588235294117647 0.588235294117647 0.588235294117647;0.592156862745098 0.592156862745098 0.592156862745098;0.596078431372549 0.596078431372549 0.596078431372549;0.6 0.6 0.6;0.603921568627451 0.603921568627451 0.603921568627451;0.607843137254902 0.607843137254902 0.607843137254902;0.611764705882353 0.611764705882353 0.611764705882353;0.615686274509804 0.615686274509804 0.615686274509804;0.619607843137255 0.619607843137255 0.619607843137255;0.623529411764706 0.623529411764706 0.623529411764706;0.627450980392157 0.627450980392157 0.627450980392157;0.631372549019608 0.631372549019608 0.631372549019608;0.635294117647059 0.635294117647059 0.635294117647059;0.63921568627451 0.63921568627451 0.63921568627451;0.643137254901961 0.643137254901961 0.643137254901961;0.647058823529412 0.647058823529412 0.647058823529412;0.650980392156863 0.650980392156863 0.650980392156863;0.654901960784314 0.654901960784314 0.654901960784314;0.658823529411765 0.658823529411765 0.658823529411765;0.662745098039216 0.662745098039216 0.662745098039216;0.666666666666667 0.666666666666667 0.666666666666667;0.670588235294118 0.670588235294118 0.670588235294118;0.674509803921569 0.674509803921569 0.674509803921569;0.67843137254902 0.67843137254902 0.67843137254902;0.682352941176471 0.682352941176471 0.682352941176471;0.686274509803922 0.686274509803922 0.686274509803922;0.690196078431373 0.690196078431373 0.690196078431373;0.694117647058824 0.694117647058824 0.694117647058824;0.698039215686274 0.698039215686274 0.698039215686274;0.701960784313725 0.701960784313725 0.701960784313725;0.705882352941177 0.705882352941177 0.705882352941177;0.709803921568627 0.709803921568627 0.709803921568627;0.713725490196078 0.713725490196078 0.713725490196078;0.717647058823529 0.717647058823529 0.717647058823529;0.72156862745098 0.72156862745098 0.72156862745098;0.725490196078431 0.725490196078431 0.725490196078431;0.729411764705882 0.729411764705882 0.729411764705882;0.733333333333333 0.733333333333333 0.733333333333333;0.737254901960784 0.737254901960784 0.737254901960784;0.741176470588235 0.741176470588235 0.741176470588235;0.745098039215686 0.745098039215686 0.745098039215686;0.749019607843137 0.749019607843137 0.749019607843137;0.752941176470588 0.752941176470588 0.752941176470588;0.756862745098039 0.756862745098039 0.756862745098039;0.76078431372549 0.76078431372549 0.76078431372549;0.764705882352941 0.764705882352941 0.764705882352941;0.768627450980392 0.768627450980392 0.768627450980392;0.772549019607843 0.772549019607843 0.772549019607843;0.776470588235294 0.776470588235294 0.776470588235294;0.780392156862745 0.780392156862745 0.780392156862745;0.784313725490196 0.784313725490196 0.784313725490196;0.788235294117647 0.788235294117647 0.788235294117647;0.792156862745098 0.792156862745098 0.792156862745098;0.796078431372549 0.796078431372549 0.796078431372549;0.8 0.8 0.8;0.803921568627451 0.803921568627451 0.803921568627451;0.807843137254902 0.807843137254902 0.807843137254902;0.811764705882353 0.811764705882353 0.811764705882353;0.815686274509804 0.815686274509804 0.815686274509804;0.819607843137255 0.819607843137255 0.819607843137255;0.823529411764706 0.823529411764706 0.823529411764706;0.827450980392157 0.827450980392157 0.827450980392157;0.831372549019608 0.831372549019608 0.831372549019608;0.835294117647059 0.835294117647059 0.835294117647059;0.83921568627451 0.83921568627451 0.83921568627451;0.843137254901961 0.843137254901961 0.843137254901961;0.847058823529412 0.847058823529412 0.847058823529412;0.850980392156863 0.850980392156863 0.850980392156863;0.854901960784314 0.854901960784314 0.854901960784314;0.858823529411765 0.858823529411765 0.858823529411765;0.862745098039216 0.862745098039216 0.862745098039216;0.866666666666667 0.866666666666667 0.866666666666667;0.870588235294118 0.870588235294118 0.870588235294118;0.874509803921569 0.874509803921569 0.874509803921569;0.87843137254902 0.87843137254902 0.87843137254902;0.882352941176471 0.882352941176471 0.882352941176471;0.886274509803922 0.886274509803922 0.886274509803922;0.890196078431372 0.890196078431372 0.890196078431372;0.894117647058824 0.894117647058824 0.894117647058824;0.898039215686275 0.898039215686275 0.898039215686275;0.901960784313726 0.901960784313726 0.901960784313726;0.905882352941176 0.905882352941176 0.905882352941176;0.909803921568627 0.909803921568627 0.909803921568627;0.913725490196078 0.913725490196078 0.913725490196078;0.917647058823529 0.917647058823529 0.917647058823529;0.92156862745098 0.92156862745098 0.92156862745098;0.925490196078431 0.925490196078431 0.925490196078431;0.929411764705882 0.929411764705882 0.929411764705882;0.933333333333333 0.933333333333333 0.933333333333333;0.937254901960784 0.937254901960784 0.937254901960784;0.941176470588235 0.941176470588235 0.941176470588235;0.945098039215686 0.945098039215686 0.945098039215686;0.949019607843137 0.949019607843137 0.949019607843137;0.952941176470588 0.952941176470588 0.952941176470588;0.956862745098039 0.956862745098039 0.956862745098039;0.96078431372549 0.96078431372549 0.96078431372549;0.964705882352941 0.964705882352941 0.964705882352941;0.968627450980392 0.968627450980392 0.968627450980392;0.972549019607843 0.972549019607843 0.972549019607843;0.976470588235294 0.976470588235294 0.976470588235294;0.980392156862745 0.980392156862745 0.980392156862745;0.984313725490196 0.984313725490196 0.984313725490196;0.988235294117647 0.988235294117647 0.988235294117647;0.992156862745098 0.992156862745098 0.992156862745098;0.996078431372549 0.996078431372549 0.996078431372549;1 1 1],...
'XTickLabel',{'1.4','2.8','4.2','5.6','7.0','8.4'},'YTickLabel',...
{'1.4','2.8','4.2','5.6'});
grid on

end