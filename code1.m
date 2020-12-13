% Clear workspace
clear
close(findall(0,'Type','figure'));
clc;
%% Input 
nC = 2; % number of curves
nS = 30; % number of subdivisions for each curve
ng = 30;  % number of subdivisions for grid
nB = nC*nS - 2*(nC-1); % Number of vectors on boundary
%% Curves
t = linspace(0,pi,nS);
x = linspace(-1,1,nS);
c1 = [cos(t)',sin(t)'];
c2 = [x',zeros(nS,1)];
Curves = {c1,c2};
pol = [c1;c2(2:end-1,:)]; % convex hull
%% Create mesh
bd = getBd(Curves,nS);   % bounding box 
line1 = linspace(bd(1),bd(2),ng);
line2 = linspace(bd(3),bd(4),ng);
[X,Y] = meshgrid(line1,line2);
inside = getGridInside(X,Y,pol);
total = [pol;inside];  % total points 
meshTri = delaunay(total(:,1),total(:,2));
mesh = createMesh(meshTri,total);
%% get representative VF and Solve
[vfrBd,tan] = createVfrBd(Curves,nB);
mef = FEM(meshTri,total,nB);
ux = mef.Solver(vfrBd(:,1));
uy = mef.Solver(vfrBd(:,2));
u = [ux,uy];
[ptsSing,triSing,idxPre,ptSing] = getSingularPts(mesh,total,u);
figure
set(gca,'XColor', 'none','YColor','none')
hold on
%% main loop
for i = 1:length(ptSing)
    ptS = ptSing(i);
    ptS.color = i*[0.25,0.4,0.2];
    ptS.triSing = triSing;
    ptS.initialize();
    ptS.integrate();
end
keyboard;
% %% Divide to transfinite
% patch = cell(1,4);
% bigCurves = cell.empty;
% sideCurves = cell.empty;
% %Get little parts
% for i = 1:length(ptSing)
%     littleCurves = cell.empty;
%     ptS = ptSing(i);
%     ptBound = ptS.ptBound;
%     ptBoundEdge = ptS.ptBoundEdge;
%     for j = 1: size(ptBound,1)
%         pt = ptBound(j,:);
%         inter = ptBoundEdge(j,1);
%         sideCurve = {whichCurve(pt,ptS.streams)};
%         littleCurves(end+1) = sideCurve;
%         sideCurves(end+1) = sideCurve; 
%         idC = floor(inter(1)/nS) + 1;
%         curve = Curves{idC};
%         id = inter - (idC - 1)*nS - 2*(idC - 1);
%         newC1 = [curve(1:id,:);pt];
%         newC2 = [pt;curve(id+1:end,:)];
%         if size(newC1,1) < size(newC2,1)
%             littleCurves(end+1) = {newC1};
%             bigCurves(end+1) = {newC2};
%         else
%             littleCurves(end+1) = {newC2};
%             bigCurves(end+1) = {newC1};            
%         end
%     end
%     patch(i,1) = littleCurves;
% end
% %Internal patchs
% for i = length(sideCurves)
%     c = sideCurves(i,:);
%     ptEnd = c(end,:);
%     intercept = {whichCurve(pt,bigCurves)};
%     
%     
% end

%% plot
%1
plot(c1(:,1),c1(:,2),'Color','b');
plot(c2(:,1),c2(:,2),'Color','r');
%2
quiver(pol(:,1),pol(:,2),0.1*tan(:,1),0.1*tan(:,2),'AutoScale','off','Color','g')
%3
crossBd = createCrossVF(vfrBd);
%quiver(pol(:,1),pol(:,2),0.025*crossBd(:,1),0.025*crossBd(:,2),'AutoScale','off','Color','g')
%quiver(pol(:,1),pol(:,2),0.025*crossBd(:,3),0.025*crossBd(:,4),'AutoScale','off','Color','g')
%quiver(pol(:,1),pol(:,2),0.025*crossBd(:,5),0.025*crossBd(:,6),'AutoScale','off','Color','g')
%quiver(pol(:,1),pol(:,2),0.025*crossBd(:,7),0.025*crossBd(:,8),'AutoScale','off','Color','g')
%4
%quiver(pol(:,1),pol(:,2),0.1*vfrBd(:,1),0.1*vfrBd(:,2),'AutoScale','off','Color','g')
%5
%plot(inside(:,1),inside(:,2),'x','Color','y');
%6
%triplot(meshTri,total(:,1),total(:,2),'Color','r');
%7
%quiver(total(:,1),total(:,2),0.1*u(:,1),0.1*u(:,2),'AutoScale','off','Color','g')
%8
%plot(ptsSing(:,1),ptsSing(:,2),'.','Color','cyan','MarkerSize',20)
%9
% crossBd = createCrossVF(u);
% quiver(total(:,1),total(:,2),0.025*crossBd(:,1),0.025*crossBd(:,2),'AutoScale','off','Color','g')
% quiver(total(:,1),total(:,2),0.025*crossBd(:,3),0.025*crossBd(:,4),'AutoScale','off','Color','g')
% quiver(total(:,1),total(:,2),0.025*crossBd(:,5),0.025*crossBd(:,6),'AutoScale','off','Color','g')
% quiver(total(:,1),total(:,2),0.025*crossBd(:,7),0.025*crossBd(:,8),'AutoScale','off','Color','g')
axis equal
ax = gca;
% 

%exportgraphics(ax,'myplot11.png','Resolution',1000) 



%% Functions
function out = getBd(Curves,n)
nt = length(Curves)*n;
x = zeros(nt,1);
j =1;
for i = 1: length(Curves)
    c = Curves{i};
    x(j:j+n-1,1) = c(:,1)';
    y(j:j+n-1,1) = c(:,2)';
    j = j + n;
end
xmax = max(x);
xmin = min(x);
ymax = max(y);
ymin = min(y);
out = [xmin,xmax,ymin,ymax];
end
function out = getGridInside(X,Y,pol)
out = [];
n = size(X,2);
for j = 1:n
    for i = 1:n
        pt = [X(i,j),Y(i,j)];
        if isInside(pol,pt)
            out(end+1,1:2) = pt;
        end
    end
end
end
function out = isInside(pol,pt)
n = size(pol,1);
ni = 0;
x = pt(1);
y = pt(2);
for i = 1:n
    p1 = pol(i,:);
    p2 = pol(mod(i,n)+1,:);
    sign12p = det([p2-p1;pt-p1]);
    if (p1(2) == p2(2))
        continue;
    elseif (p1(2) >= y && p2(2) >= y)
        continue;
        % teste se o ponto está acima da aresta
    elseif (p1(2) <= y && p2(2) <= y)
        continue;
        % teste se o ponto está a esquerda da aresta
    elseif ((y >= p1(2) && sign12p < 0) || (y <= p1(2) && sign12p > 0))
        continue;
        % testa interseção no ponto superior
    elseif (p2(2) == y)
        continue;
    end
    ni = ni + 1;
end
out = false;
if (mod(ni,2) > 0 )
    out = true;
end
end
function [VfrBd,tan] = createVfrBd(curves,nB)
VfrBd = zeros(nB,2);
tan = zeros(nB,2);
nc = length(curves);
count = 2;
count2 = 1;
bd = zeros(2*nc,2);
nv = zeros(nc,1);
for i = 1:nc
    c = curves{i};
    nv(i) = size(c,1);
    for j = 2: (nv(i)-1)
        pt1 = c(j,:);
        pt2 = c(j+1,:);
        VfrBd(count,:) = VFrpr(pt1,pt2);
        tan(count,:) = (pt2 - pt1)/norm((pt2 - pt1));
        count = count + 1;
    end
    count = count + 1;
    bd(count2,:) = VFrpr(c(1,:),c(2,:));
    bd(count2+1,:) = VFrpr(c(end-1,:),c(end,:));
    count2 = count2 + 2;
end
count1 = 1;
count2 = 1;
  for i = 1:nc
      VfrBd(count1,:) = (bd(count2,:) + bd(count2+2,:))/2;
      count1 = count1 + nv(i)-1;
      count2 = count2 + 1;
  end
end
function out = VFrpr(pt1,pt2)
   v = (pt2 - pt1);
   teta = mod(mod(atan2(v(2),v(1)),2*pi),2*pi);
   tetaDr = mod(4*teta,2*pi); 
   out = [cos(tetaDr),sin(tetaDr)];
end
function [ptsSing,triSing,idxPre,ptSing] = getSingularPts(mesh,pts,u)
n = length(mesh);
ptsSing = [];
triSing = [];
idxPre = [];
ptSing = Singular.empty;
for i = 1:n
    nodes = mesh(i).nodes;
    v = u(nodes,:);
    idx = idxPoincare(v);
    if (idx == 0)
        continue
    end
    ptSingular = Singular();
    Pe = [pts(nodes,:)';ones(1,3)];
    M = [v(1,:);v(2,:);v(3,:)];
    A = 2*M*M';
    A(end+1,1:3) = ones(1,3);
    A(1:3,end+1) = ones(3,1);
    b = [zeros(3,1);1];
    lam = A\b;
    pt = Pe*lam(1:3);
    ptsSing(end+1,1:2) = (pt(1:2))';
    triSing(end+1,1:3) = nodes;
    idxPre(end+1) = idx;
    ptSingular.pt = (pt(1:2))';
    ptSingular.tri = mesh(i);
    ptSingular.idxPoincare = idx; 
    ptSingular.u = u;
    ptSing(end+1) = ptSingular; 
end
end
function int = idxPoincare(v)
v1 = v(1,:);
v2 = v(2,:);
v3 = v(3,:);
teta1 = mod(mod(atan2(v1(2),v1(1)),2*pi),2*pi);
teta2 = mod(mod(atan2(v2(2),v2(1)),2*pi),2*pi);
teta3 = mod(mod(atan2(v3(2),v3(1)),2*pi),2*pi);
int = (deltapi(teta1,teta2) + deltapi(teta2,teta3) + deltapi(teta3,teta1))/(2*pi);
int = round(int);
end
function out = deltapi(tetai,tetaj)
out = mod((tetaj - tetai + pi),2*pi) - pi;
end
function mesh = createMesh(meshTri,pts)
   n = size(meshTri,1);
   mesh(n) = triangle;
   for i = 1:n
       nodes = meshTri(i,:);
       mesh(i).nodes = nodes; 
       mesh(i).pts = pts(nodes,:);
       mesh(i).triadj(3) = triangle;
       edge1 = [nodes(2),nodes(1)];
       edge2 = [nodes(2),nodes(3)];
       edge3 = [nodes(3),nodes(1)];
       idx1 = find(sum(ismember(meshTri,edge1),2) == 2); 
       idx2 = find(sum(ismember(meshTri,edge2),2) == 2); 
       idx3 = find(sum(ismember(meshTri,edge3),2) == 2); 
       if (length(idx1) == 2)
           t1 = triangle;
           t1.nodes = meshTri(idx1(find(idx1 ~= i)),:);
           t1.pts = pts(t1.nodes,:);
           mesh(i).triadj(1) = t1;
       end
       if (length(idx2) == 2)
           t2 = triangle;
           t2.nodes = meshTri(idx2(find(idx2 ~= i)),:);
           t2.pts = pts(t2.nodes,:);
           mesh(i).triadj(2) = t2;
       end
       if (length(idx3) == 2)
           t3 = triangle;
           t3.nodes = meshTri(idx3(find(idx3 ~= i)),:);
           t3.pts = pts(t3.nodes,:);
           mesh(i).triadj(3) = t3;
       end
   end
   for i = 1:n
       t = mesh(i);
       t.mesh = mesh;
   end
end
function crossVF = createCrossVF(u)
 teta = mod(mod(atan2(u(:,2),u(:,1)),2*pi),2*pi);
 teta = teta/4;
 v0 = [cos(teta),sin(teta)];
 v1 = [cos(teta + pi/2),sin(teta + pi/2)];
 v2 = [cos(teta + 2*pi/2),sin(teta + 2*pi/2)];
 v3 = [cos(teta + 3*pi/2),sin(teta + 3*pi/2)];
 crossVF = [v0,v1,v2,v3];
end
function out = whichCurve(pt,curves)
out = [];
bool = false;
for i = 1:length(curves)
    c = curves{i};
    if pt == c(1,:)
        out = c; 
        bool = true;
    elseif pt == c(end,:)
        out = c; 
        bool = false; 
    end
end
end