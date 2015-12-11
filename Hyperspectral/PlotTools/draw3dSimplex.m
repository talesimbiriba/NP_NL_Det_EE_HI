function draw3dSimplex(M, figureHandle, M2, simplexColor )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


d = size(M,2);

figure(figureHandle);
if nargin<4
    simplexColor = 'k';
end

if (d==3)
    if nargin>=3
        mm = M'*M2;
    else
        mm = M'*M;
    end
    scatter3(mm(:,1),mm(:,2),mm(:,3),['o',simplexColor],'fill')

    hold on
    plot3([mm(1,1) mm(2,1)],[mm(1,2) mm(2,2)],[mm(1,3) mm(2,3)],['-',simplexColor],'LineWidth',1);
    plot3([mm(2,1) mm(3,1)],[mm(2,2) mm(3,2)],[mm(2,3) mm(3,3)],['-',simplexColor],'LineWidth',1);
    plot3([mm(1,1) mm(3,1)],[mm(1,2) mm(3,2)],[mm(1,3) mm(3,3)],['-',simplexColor],'LineWidth',1);
end

if d==2
    disp('Not implemented yet!')
end
if d==4
    disp('Not implemented yet!')
end
end

