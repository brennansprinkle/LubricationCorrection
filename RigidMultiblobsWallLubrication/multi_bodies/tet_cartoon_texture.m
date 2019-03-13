clc
clf
%close all

tet = [0.0	0.0	4.449489742783178
1.732050807568877	0.0	2.0
-0.866025403784439	-1.500000000000000	2.0
-0.866025403784439	1.500000000000000	2.0];

Ry = @(theta) [cos(theta) 0 sin(theta);0 1 0; -sin(theta) 0 cos(theta)];



[x,y,z] = sphere(400);

% [m,n] = size(x);
% theta = -pi/7;
% for i = 1:m
%     for j = 1:n
%         X = [x(i,j); y(i,j); z(i,j)];
%         RX = Ry(theta)*X;
%         x(i,j) = RX(1);
%         y(i,j) = RX(2);
%         z(i,j) = RX(3);
%     end
% end

% th = 2*pi*rand(400,1);
% u = 2*rand(400,1) - 1;
% 
% [U,T] = meshgrid(,sort(th));
% x = sqrt(1 - U.^2).*cos(T);
% y = sqrt(1 - U.^2).*sin(T);
% z = U;


% derp = dlmread('Structures/shell_N_2562_Rg_1_Rh_1_0113.vertex');
% % derp = dlmread('Structures/shell_N_642_Rg_1_Rh_1_0239.vertex');
% derp(1,:) = [];
% x = derp(:,1);
% y = derp(:,2);
% z = derp(:,3);
% DT_A= DelaunayTri(x,y,z); 
% [tri, XA]= freeBoundary(DT_A); 
% h = trisurf(tri,XA(:,1),XA(:,2),XA(:,3));
% xdata = get(h,'xdata');
% ydata = get(h,'ydata');
% zdata = get(h,'zdata');

for k = 1:4
   [I,map] = imread(['./number_pictures/' num2str(k) '.png']);
    I(I == 0) = 255;
    I(I < 250) = 0;
    I = [0*I+255 I 0*I+255];
    I = [0*I+255 I 0*I+255];
    I = [0*I+255; I; 0*I+255];
    I = I - uint8(randi(30,size(I)));
%     I = imresize(I,5*size(x));
   %((abs(z)<0.99) | ((abs(z)<0.9995) & (mod(fix(atan2(y,x)*1000),5) == 0)))
   r_texture = (1 + 0.005*randn(size(x)).*((abs(z)<(1-1*abs(randn(size(z)))))));
   
   
   h = warp(r_texture.*x+tet(k,1),r_texture.*y+tet(k,2),r_texture.*z+tet(k,3),fliplr(flipud(I)));
   set(h,'FaceLighting','gouraud', ...
    'AmbientStrength',0.1, ...
    'DiffuseStrength',0.9, ... 
    'Clipping','off',...
    'BackFaceLighting','lit', ...
    'SpecularStrength',0.3, ...
    'SpecularColorReflectance',0.7, ...
    'SpecularExponent',1)
   hold all
end
[X,Y] = meshgrid(-6:0.1:4,-5:0.1:4);
surf(X,Y,0*X + 0.05*rand(size(X)),...
    'EdgeColor','none', ...
    'AlphaData',(X+2*min(X(:))).^3+(Y+2*min(Y(:))).^3,...
    'FaceAlpha','flat',...
    'FaceColor',[33 182 57]/255, ... %[33 182 57]/255
    'FaceLighting','gouraud', ...
    'AmbientStrength',0.8, ...
    'DiffuseStrength',0.9, ... 
    'Clipping','off',...
    'BackFaceLighting','lit', ...
    'SpecularStrength',0.9, ...
    'SpecularColorReflectance',0.3, ...
    'SpecularExponent',20)

hold all
l1 = light('Position',10*[4 -2 8], ...
    'Style','local', ...
    'Color',0.8*[0.60000 0.10000 0.590000]); %


l4 = light('Position',10*[4 4 1], ...
    'Style','local', ...
    'Color',0.9*[1.000000 0.900000 0.100000]); %

hold all
h = arrow3D([1.732050807568877	0.0	2.0], [-2.5 0 0], 'r', .75);
set(h,'EdgeColor','none', ...
    'FaceColor',0.99*[1 0.2 .4], ... %0.5*[0 0.7 1]
    'FaceLighting','gouraud', ...
    'AmbientStrength',0.2, ...
    'DiffuseStrength',0.2, ... 
    'Clipping','off',...
    'BackFaceLighting','lit', ...
    'SpecularStrength',0.9, ...
    'SpecularColorReflectance',0.5, ...
    'SpecularExponent',.2)

for i = 1:6
if(i<4)
    [X, Y, Z] = cylinder2P(0.05, 50,tet(1,:),tet(i+1,:));
elseif(i == 4)
    [X, Y, Z] = cylinder2P(0.05, 50,tet(2,:),tet(4,:));
elseif(i == 5)
    [X, Y, Z] = cylinder2P(0.05, 50,tet(2,:),tet(3,:));
elseif(i == 6)
    [X, Y, Z] = cylinder2P(0.05, 50,tet(3,:),tet(4,:));
end
h = surf(X,Y,Z);
set(h,'EdgeColor','none', ...
    'FaceColor',0.2*[0.5 0.5 0.5], ... 0.4*[0.8 0.3 0.9]
    'facealpha',0.4,...
    'FaceLighting','gouraud', ...
    'AmbientStrength',0.9, ...
    'DiffuseStrength',0.9, ... 
    'Clipping','off',...
    'BackFaceLighting','lit', ...
    'SpecularStrength',0.9, ...
    'SpecularColorReflectance',0.1, ...
    'SpecularExponent',5)
hold all
end

for i = 2:4
    r1 = tet(i,:);
    r2 = r1;
    r2(3) = 0;
    [X, Y, Z] = cylinder2P(0.05, 50,r1,r2);
h = surf(X,Y,Z);
set(h,'EdgeColor','none', ...
    'FaceColor',0.2*[0.5 0.5 0.5], ... 0.4*[0.8 0.3 0.9]
    'facealpha',0.4,...
    'FaceLighting','gouraud', ...
    'AmbientStrength',0.9, ...
    'DiffuseStrength',0.9, ... 
    'Clipping','off',...
    'BackFaceLighting','lit', ...
    'SpecularStrength',0.9, ...
    'SpecularColorReflectance',0.1, ...
    'SpecularExponent',5)
hold all
end

annotation(gcf,'textarrow',[0.458854166666667 0.475],...
    [0.290172595520422 0.412384716732543],'String',{'$$\epsilon$$'},...
    'Interpreter','latex',...
    'FontSize',40);
annotation(gcf,'arrow',[0.4453125 0.388541666666667],...
    [0.276688453159042 0.316859989494261]);


daspect([1 1 1])
xlim(2*[-1 1])
ylim(2*[-1 1])
zlim([-0.471304456776432 4.528695543223568])
view([57 37])
set(gca, 'Position', [0 0 1 1]);
zoom(1.1)
set(gca,'Visible','off')
box off
grid off

print('-dpng','-r300',['/home/bs162/LubricationCode/Tetra_Texture_fuzz_Cartoon.png'])