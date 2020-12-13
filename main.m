% Clear workspace
% clear
% close(findall(0,'Type','figure'));
% clc;
%% Input
nC = 2; % number of curves
nS = 200; % number of subdivisions for each curve
ng = 60;  % number of subdivisions for grid
nB = nC*nS - 2*(nC-1); % Number of vectors on boundary
%% Curves
t = linspace(0,2.3561944901923448,nS);
% x = linspace(0,1,nS);
radius = 6;
c1 = radius * [cos(t)',sin(t)'] + [3,4];
c2 = bezier([-1.2426406871192857, 8.2426406871192857],[9, 4],[-2.7207293666026877, 5.1017274472168905],[4.6055662188099804, 8.0153550863723613],nS);
Curves = {c1,c2};
% edgesX = 1:nB - 1;
% edgesY = 2:nB;
% edges = [edgesX', edgesY'];
pol = [c1;c2(2:end-1,:)]; % convex hull
% bool = [1;zeros(nS-2,1);1;zeros(nS-2,1)];
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

%% main loop
for i = 1:length(ptSing)
    ptS = ptSing(i);
    ptS.color = i*[0.25,0.4,0.2];
    ptS.triSing = triSing;
    ptS.initialize();
    ptS.integrate();
end
% [p1,p2,p3,p4] = transfinito(ptSing,Curves,nS,nC,pol);
%% plot
%1
figure
set(gca,'XColor', 'none','YColor','none')
hold on

axis equal
ax = gca;
plot(c1(:,1),c1(:,2),'Color','b');
plot(c2(:,1),c2(:,2),'Color','r');
%exportgraphics(ax,'myplot17.png','Resolution',1000)
%2
%quiver(pol(:,1),pol(:,2),0.1*tan(:,1),0.1*tan(:,2),'AutoScale','off','Color','g')
%3
crossBd = createCrossVF(vfrBd);
%quiver(pol(:,1),pol(:,2),0.025*crossBd(:,1),0.025*crossBd(:,2),'AutoScale','off','Color','g')
%quiver(pol(:,1),pol(:,2),0.025*crossBd(:,3),0.025*crossBd(:,4),'AutoScale','off','Color','g')
%quiver(pol(:,1),pol(:,2),0.025*crossBd(:,5),0.025*crossBd(:,6),'AutoScale','off','Color','g')
%quiver(pol(:,1),pol(:,2),0.025*crossBd(:,7),0.025*crossBd(:,8),'AutoScale','off','Color','g')
%4
quiver(pol(:,1),pol(:,2),0.1*vfrBd(:,1),0.1*vfrBd(:,2),'AutoScale','off','Color','g')
%5
%plot(inside(:,1),inside(:,2),'x','Color','y');
%6
%triplot(meshTri,total(:,1),total(:,2),'Color','r');
%7
quiver(total(:,1),total(:,2),0.1*u(:,1),0.1*u(:,2),'AutoScale','off','Color','g')
%8
plot(ptsSing(:,1),ptsSing(:,2),'.','Color','cyan','MarkerSize',20)
%9
% crossBd = createCrossVF(u);
% quiver(total(:,1),total(:,2),0.025*crossBd(:,1),0.025*crossBd(:,2),'AutoScale','off','Color','g')
% quiver(total(:,1),total(:,2),0.025*crossBd(:,3),0.025*crossBd(:,4),'AutoScale','off','Color','g')
% quiver(total(:,1),total(:,2),0.025*crossBd(:,5),0.025*crossBd(:,6),'AutoScale','off','Color','g')
% quiver(total(:,1),total(:,2),0.025*crossBd(:,7),0.025*crossBd(:,8),'AutoScale','off','Color','g')

%plot esqueleto

% p1.plot();
% %exportgraphics(ax,'myplot12.png','Resolution',1000)
% p2.plot();
% %exportgraphics(ax,'myplot13.png','Resolution',1000)
% p3.plot();
% %exportgraphics(ax,'myplot14.png','Resolution',1000)
% p4.plot();
%exportgraphics(ax,'myplot15.png','Resolution',1000)
%
% for i = 1:length(internalCurves)
%     ci = internalCurves{i};
%     plot(ci(:,1),ci(:,2),'color','red');
% end

%exportgraphics(ax,'myplot16.png','Resolution',1000)

%exportgraphics(ax,'myplot11.png','Resolution',1000)



%% Functions
function out = getBd(Curves,n)
nt = length(Curves)*n;
x = zeros(nt,1);
y = zeros(nt,1);
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
    if i == 13
        a=1;
    end
    idx = idxPoincare(v);
    if (idx == 0 )
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
function out = whichSing(id,ptSing)
for i=1:length(ptSing)
    ptS = ptSing(i);
    if ismember(id,ptS.ptBoundEdge(:,1))
        [~,loc] = ismember(id,ptS.ptBoundEdge(:,1));
        out = ptS.ptBound(loc,:);
        return
    elseif ismember(id,ptS.ptBoundEdge(:,2))
        [~,loc] = ismember(id,ptS.ptBoundEdge(:,2));
        out = ptS.ptBound(loc,:);
        return
    end
end
end
function [p1,p2,p3,p4] = transfinito(ptSing,Curves,nS,nC,pol)
%% Divide to transfinite
%Get little parts
internalCurves = cell.empty;
boundaryCurves = cell.empty;
mainCurves = cell.empty;
%Internal curves
for i = 1:length(ptSing)
    ptS = ptSing(i);
    ptBound = ptS.ptBound;
    ptBoundEdge = ptS.ptBoundEdge;
    for j = 1: size(ptBound,1)
        pt = ptBound(j,:);
        inter = ptBoundEdge(j,1);
        sideCurve = {whichCurve(pt,ptS.streams)};
        boundaryCurves(end+1) = sideCurve;
        internalCurves(end+1) = sideCurve;
    end
    for j = 1:length(ptS.streams)
        c = ptS.streams{j};
        cbd2 = boundaryCurves{1};
        cbd3 = boundaryCurves{2};
        if c(end,:) ~= cbd2(end,:) & c(end,:) ~= cbd3(end,:)
            mainCurves(end+1) = {c};
            break;
        end
    end
    boundaryCurves = cell.empty;
end
cmain1 = mainCurves{1};
cmain2 = mainCurves{2};
N = [size(cmain1,1),size(cmain2,1)];
N1 = size(cmain1,1);
N2 = size(cmain2,1);
Nmax = max(N1,N2);
def = Nmax - N2;
pt2 = cmain2(end,:);
teste = cmain1 - pt2;
norm = vecnorm(teste,2,2);
[~,id] = min(norm);
mergedCurve = (cmain1(id:end,:) + flip(cmain2((id-def):end,:) ) )/2;
mainCurve = [cmain1(1:id,:);mergedCurve;flip(cmain2(1:(id-def),:))];
internalCurves(end+1) = {mainCurve};
% External curves
externalCurves = cell.empty;
n = 1;
ptBoundEdgeParent = [ptSing(:).ptBoundEdge];
ptBoundEdge = ptBoundEdgeParent(:);
indexBef = 1;
for i = 1:length(Curves)
    newN = ( n + nS -1 - nC*(i-1)) ;
    ptBoundEdge = [ptBoundEdge ; n ; newN];
    ptBoundEdge = sort(ptBoundEdge);
    [~,index] = ismember(newN,ptBoundEdge);
    indexCurves = ptBoundEdge(indexBef:index); %será sempre um número par
    for j = 1:2:size(indexCurves)
        newCurve = pol(indexCurves(j):indexCurves(j+1),:);
        if ismember(indexCurves(j),ptBoundEdgeParent)
            pt = whichSing(indexCurves(j),ptSing);
            newCurve = [pt ; newCurve];
            if i > 1
                if (indexCurves(j) == n )
                    c = Curves{i};
                    pt = c(1,:);
                    newCurve = [pt ; newCurve];
                end
            end
        elseif ismember(indexCurves(j+1),ptBoundEdgeParent)
            pt = whichSing(indexCurves(j+1),ptSing);
            newCurve = [newCurve ; pt];
            if i > 1
                if (indexCurves(j+1) == newN  )
                    c = Curves{i};
                    pt = c(end,:);
                    newCurve = [newCurve;pt];
                end
            end
        end
        externalCurves(end+1) = {newCurve};
    end
    n = ( n + nS -1 - nC*(i-1))  + 1;
    indexBef = index+1;
end
patch1 = {internalCurves{3},flip(externalCurves{1}),internalCurves{4},externalCurves{6}};
p1 = Projetor(patch1);
patch2 = {internalCurves{5},flip(externalCurves{2}),internalCurves{2},internalCurves{4}};
p2 = Projetor(patch2);
patch3 = {externalCurves{3},internalCurves{1},flip(internalCurves{2}),externalCurves{4}};
p3 = Projetor(patch3);
patch4 = {externalCurves{5},flip(internalCurves{5}),flip(internalCurves{1}),flip(internalCurves{3})};
p4 = Projetor(patch4);
end