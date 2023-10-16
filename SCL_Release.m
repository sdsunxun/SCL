%SCL3 利用系统函数进行极小值过滤,动态参数

function [cout,curve,curveindex,cd] = SCL(varargin)



% Parse the input parameters
[I,High,Low,Gap_size,EP,ChordLength,Th] = parse_inputs(varargin{:});
if size(I,3)==3
    I=rgb2gray(I); % Transform RGB image to a Gray one. 
end

% I         : the input grayscale image
% High      : the high threshold for the Canny edge detector in [0,1] (default = 0.7)
% Low       : the low threshold for the Canny edge detector in [0,1] (default = 0.2)
% Gap_size  : The maximum gap between the successive edge-points. If the
%           gap between the sucessive points is more than Gap_size, then
%           they would be in two different edges (default = 1 pixel)
% EP        : A flag to indicate whether the end-points of a curve be
%           detected as corners (default = 0)

% Detect edges
BW = edge(I,'canny',[Low,High]);  % [Low,High] low and high sensitivity thresholds
% BW = edge(I,'canny'); % If sensitivity thresholds are not specified Canny edge detector selects itself
% output:  binary edge-image BW of the same size as I, with 1's where the function finds edges in I and 0's elsewhere

% Extract curves from the edge-image
[curve,curve_start,curve_end,curve_mode,curve_num,TJ,img1] = extract_curve(BW,Gap_size);  
% curve         : MATLAB cell data structure where each cell is a 2D array containing
%               pixel locations (x and y values)
% curve_start   : starting points of extracted curves
% curve_end     : ending points of extracted curves
% curve_mode    : two curve modes - 'line' and 'loop'. If the both ends of
%               a curve are at maximum 25 square pixels (default) away, then the
%               curve is a loop curve, otherwise a line curve
% curve_num     : number of extracted curves
% TJ            : the T-junction found in the edge-extraction process
% img1          : output image containing the extracted edges

[sizex sizey] = size(I);
if size(curve{1})>0
    % Detect corners on the extracted edges
%     tic
    [corner_out curveindex cd2] = getcorner(curve,curve_mode,curve_start,curve_num,sizex,sizey,I,ChordLength,Th); 
    
    % corner_out	: n by 2 matrix containing the positions of the
    %               detected corners, where n is the number of detected
    %               corners
    % index         : MATLAB cell data structure where each cell is an 1D
    %               column matrix contaning the edge pixel numbers (in curve) where
    %               the corners are detected
    % Sig           : the sigma values used to smooth the curves
    % cd2           : cpda curvature values of the detected corners
    
    % Update the T-junctions
    [corner_final cd3] = Refine_TJunctions(corner_out,TJ,cd2,curve, curve_num, curve_start, curve_end, curve_mode,EP);
    % corner_final  : n by 2 matrix containing the positions of the
    %               detected corners, where n is the number of detected
    %               corners
    % cd3           : cpda curvature values of the detected corners
%     time = toc
    
%     img=img1; % to show corners on the edge image
%     %img=I; % to show corners on the input image
%     for i=1:size(corner_final,1)
%         img=mark(img,corner_final(i,1),corner_final(i,2),7);
%     end
% 
%     marked_img=img;
%     figure(); imshow(marked_img);
%     figure(); imshow(BW); hold on; plot(corner_final(:,2),corner_final(:,1),'r*');
     %imwrite(marked_img,'acclaim002_cpda_corner.bmp','BMP');
    cout = corner_final;
    cd = cd3;
else
    cout = [];
    marked_img = [];
    cd = [];
end

here = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [corners curveindex cd] = getcorner(curve,curve_mode,curve_start,curve_num,sizex,sizey,img,ChordLength,Th);
    corners = [];
    cor = []; % candidate corners
    cd = [];
    curveindex={};
    T = Th; % define the curvature threshold   %阈值修改
    sig = 3.0;
    [gau W] = makeGFilter(sig);
    for i=1:curve_num;
        C = []; C3 = [];
        x = curve{i}(:,2);
        y = curve{i}(:,1);
        curveLen = size(x,1);    
        [xs ys W] = smoothing(x,y,curveLen,curve_mode(i,:),gau,W); % smooth the curve with Gaussian kernel
        if size(xs,1)>1   
            if curve_mode(i,:)=='loop'
                xs1=[xs(curveLen-W+1:curveLen);xs;xs(1:W)];
                ys1=[ys(curveLen-W+1:curveLen);ys;ys(1:W)];
            else %expand the ends to gaussian window
                xs1=[ones(W,1)*2*xs(1)-xs(W+1:-1:2);xs;ones(W,1)*2*xs(curveLen)-xs(curveLen-1:-1:curveLen-W)];
                ys1=[ones(W,1)*2*ys(1)-ys(W+1:-1:2);ys;ones(W,1)*2*ys(curveLen)-ys(curveLen-1:-1:curveLen-W)];
            end   
            xs = xs1;
            ys = ys1;
            L = curveLen+2*W;   
            tic
            %曲率计算开始
            C = DirectChordLength(xs,ys,ChordLength,L,curve_mode(i,:),W);    %计算曲率修改
             %曲率计算结束
            time=toc;
            
            
%             figure; imshow(img); hold on; plot(xs,ys,'r-');   
%             hold on;
%             plot(1:size(C),T,'r-');
            L = curveLen;
    
    
%     figure()
%     plot(C,'K') 
%     axis([0,L,0,2]);
    
    % Find curvature local maxima as corner candidates
            [C_extremum, Pos_extremum]=findpeaks(-C);
            C_extremum=-C_extremum;
            Pos_extremum2=[];
            for j=1:1:size(Pos_extremum,1) % if the maxima is less than local minima, remove it as flase corner
                if (C_extremum(j) < T)
                   Pos_extremum2=[Pos_extremum2,Pos_extremum(j)];
                end
            end
            
            extremum=Pos_extremum2;
            extremum=extremum(find(extremum>0 & extremum<=curveLen)); % find corners which are not endpoints of the curve             
            index{i} = extremum';
            Sig(i,1) = sig;
            n = size(extremum,2); 
            curveindex{i}=extremum;
            for j = 1:n
                corners = [corners; curve{i}(extremum(j),:)];
                cd = [cd; C(extremum(j))];
            end    
            
            if curve_mode(i,:)=='loop'& size(corners,1)>0   
                compare_corner=corners-ones(size(corners,1),1)*curve{i}(1,:);
                compare_corner=compare_corner.^2;
                compare_corner=compare_corner(:,1)+compare_corner(:,2);
                if min(compare_corner)>100       % Add end points far from detected corners, i.e. outside of 5 by 5 neighbor                                                
                    if(C(1)<T)
                        corners = [corners; curve{i}(1,:)];
                        cd = [cd; C(1)];
                    end
                end
            end
            
%             hold on; plot(corners(:,2),corners(:,1),'b*');
%             figure;plot(C,'-');   

        end 
    here = 1;
    end 

here = 1;



function Cd = DirectChordLength(xs,ys,chordLen,curveLen,curve_mode,W)

Cd=[];

pixel_distance=[];

for i=1:curveLen-1
    pixel_distance=[pixel_distance;pdist([xs(i),ys(i);xs(i+1),ys(i+1)],'euclidean')];
end

for i=W+1:curveLen-W
    
    
    left_index=2;
    right_index=2;
    
    RootPoint=[xs(i) ys(i)];
    

    try
    while(pdist([xs(i),ys(i);xs(i-left_index),ys(i-left_index)],'euclidean')<chordLen)
        left_index=left_index+1;
    end
    catch exec
        left_index=left_index-1;
    end

    PointA=[xs(i-left_index+1) ys(i-left_index+1)];
    distA=pdist([RootPoint;PointA],'euclidean');
    

    PointB=[xs(i-left_index) ys(i-left_index)]; 
    distB=pdist([RootPoint;PointB],'euclidean');
    LeftPoint=getPointWhereDistAndTwoPoint(RootPoint,PointA,PointB,chordLen);
    
    

    try
    while(pdist([xs(i),ys(i);xs(i+right_index),ys(i+right_index)],'euclidean')<chordLen)
        right_index=right_index+1;
    end
    catch exec
        right_index=right_index-1;
    end
    PointA=[xs(i+right_index-1) ys(i+right_index-1)];
    distA=pdist([RootPoint;PointA],'euclidean');
    PointB=[xs(i+right_index) ys(i+right_index)];
    distB=pdist([RootPoint;PointB],'euclidean');
    RightPoint=getPointWhereDistAndTwoPoint(RootPoint,PointA,PointB,chordLen);
 
   MiddlePoint=[(LeftPoint(1,1)+RightPoint(1,1))/2 (LeftPoint(1,2)+RightPoint(1,2))/2]; 
   
%    figure; plot(xs,ys,'b-.'); hold on; plot(LeftPoint(1,1),LeftPoint(1,2),'r*');hold on; plot(RightPoint(1,1),RightPoint(1,2),'r*');hold on; plot(MiddlePoint(1,1),MiddlePoint(1,2),'r*');
   

   Cd=[Cd;pdist([LeftPoint;RightPoint],'euclidean')];
   
end


function [point]=getPointWhereDistAndTwoPoint(root_point,a_point,b_point,chord_len)

x0=root_point(1,1);
y0=root_point(1,2);
x1=a_point(1,1);
y1=a_point(1,2);
x2=b_point(1,1);
y2=b_point(1,2);
d=chord_len;
x = -((x2-x1)*sqrt(((-x1^2)+2*x0*x1-x0^2+d^2)*y2^2+(((2*x1-2*x0)*x2-2*x0*x1+2*x0^2-2*d^2)*y1+((2*x0-2*x1)*x2+2*x1^2-2*x0*x1)*y0)*y2+((-x2^2)+2*x0*x2-x0^2+d^2)*y1^2+(2*x2^2+((-2*x1)-2*x0)*x2+2*x0*x1)*y0*y1+((-x2^2)+2*x1*x2-x1^2)*y0^2+d^2*x2^2-2*d^2*x1*x2+d^2*x1^2)-x1*y2^2+((x2+x1)*y1+(x1-x2)*y0)*y2-x2*y1^2+(x2-x1)*y0*y1-x0*x2^2+2*x0*x1*x2-x0*x1^2)/(y2^2-2*y1*y2+y1^2+x2^2-2*x1*x2+x1^2);
y = -((y2-y1)*sqrt(((-x1^2)+2*x0*x1-x0^2+d^2)*y2^2+(((2*x1-2*x0)*x2-2*x0*x1+2*x0^2-2*d^2)*y1+((2*x0-2*x1)*x2+2*x1^2-2*x0*x1)*y0)*y2+((-x2^2)+2*x0*x2-x0^2+d^2)*y1^2+(2*x2^2+((-2*x1)-2*x0)*x2+2*x0*x1)*y0*y1+((-x2^2)+2*x1*x2-x1^2)*y0^2+d^2*x2^2-2*d^2*x1*x2+d^2*x1^2)-y0*y2^2+(2*y0*y1+(x1-x0)*x2-x1^2+x0*x1)*y2-y0*y1^2+((-x2^2)+(x1+x0)*x2-x0*x1)*y1)/(y2^2-2*y1*y2+y1^2+x2^2-2*x1*x2+x1^2);

if(pdist([[x y];[x1 y1]],'euclidean')+pdist([[x y];[x2 y2]],'euclidean')==pdist([[x1 y1];[x2 y2]],'euclidean'))
    point=[x y]; 
else
    x = ((x2-x1)*sqrt(((-x1^2)+2*x0*x1-x0^2+d^2)*y2^2+(((2*x1-2*x0)*x2-2*x0*x1+2*x0^2-2*d^2)*y1+((2*x0-2*x1)*x2+2*x1^2-2*x0*x1)*y0)*y2+((-x2^2)+2*x0*x2-x0^2+d^2)*y1^2+(2*x2^2+((-2*x1)-2*x0)*x2+2*x0*x1)*y0*y1+((-x2^2)+2*x1*x2-x1^2)*y0^2+d^2*x2^2-2*d^2*x1*x2+d^2*x1^2)+x1*y2^2+(((-x2)-x1)*y1+(x2-x1)*y0)*y2+x2*y1^2+(x1-x2)*y0*y1+x0*x2^2-2*x0*x1*x2+x0*x1^2)/(y2^2-2*y1*y2+y1^2+x2^2-2*x1*x2+x1^2);
    y = ((y2-y1)*sqrt(((-x1^2)+2*x0*x1-x0^2+d^2)*y2^2+(((2*x1-2*x0)*x2-2*x0*x1+2*x0^2-2*d^2)*y1+((2*x0-2*x1)*x2+2*x1^2-2*x0*x1)*y0)*y2+((-x2^2)+2*x0*x2-x0^2+d^2)*y1^2+(2*x2^2+((-2*x1)-2*x0)*x2+2*x0*x1)*y0*y1+((-x2^2)+2*x1*x2-x1^2)*y0^2+d^2*x2^2-2*d^2*x1*x2+d^2*x1^2)+y0*y2^2+((-2*y0*y1)+(x0-x1)*x2+x1^2-x0*x1)*y2+y0*y1^2+(x2^2+((-x1)-x0)*x2+x0*x1)*y1)/(y2^2-2*y1*y2+y1^2+x2^2-2*x1*x2+x1^2);
    point=[x y]; 
end

    
    
function [xse yse] = enlarge(xs,ys,CL,curve_mode);
%CL = chord length
L = size(xs,1);
if curve_mode=='loop' % wrap around the curve by CL pixles at both ends
    xse = [xs(L-CL+1:L);xs;xs(1:CL)];
    yse = [ys(L-CL+1:L);ys;ys(1:CL)];
else % extend each line curve by CL pixels at both ends
    xse = [ones(CL,1)*2*xs(1)-xs(CL+1:-1:2);xs;ones(CL,1)*2*xs(L)-xs(L-1:-1:L-CL)];
    yse = [ones(CL,1)*2*ys(1)-ys(CL+1:-1:2);ys;ones(CL,1)*2*ys(L)-ys(L-1:-1:L-CL)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%
function [xs ys W] = smoothing(x,y,L,curve_mode,gau,W);

if L>W
    if curve_mode=='loop' % wrap around the curve by W pixles at both ends
        x1 = [x(L-W+1:L);x;x(1:W)];
        y1 = [y(L-W+1:L);y;y(1:W)];
    else % extend each line curve by W pixels at both ends
        x1 = [ones(W,1)*2*x(1)-x(W+1:-1:2);x;ones(W,1)*2*x(L)-x(L-1:-1:L-W)];
        y1 = [ones(W,1)*2*y(1)-y(W+1:-1:2);y;ones(W,1)*2*y(L)-y(L-1:-1:L-W)];
    end
    
    xx=conv(x1,gau);
    xs=xx(2*W+1:L+2*W);
    yy=conv(y1,gau);
    ys=yy(2*W+1:L+2*W);    
else
    xs = [];
    ys = [];    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% extract curves from input edge-image
function [curve,curve_start,curve_end,curve_mode,cur_num,TJ,img]=extract_curve(BW,Gap_size)
%   Function to extract curves from binary edge map, if the endpoint of a
%   contour is nearly connected to another endpoint, fill the gap and continue
%   the extraction. The default gap size is 1 pixel
[L,W]=size(BW);
BW1=zeros(L+2*Gap_size,W+2*Gap_size);
BW_edge=zeros(L,W);
BW1(Gap_size+1:Gap_size+L,Gap_size+1:Gap_size+W)=BW;
[r,c]=find(BW1==1); %returns indices of non-zero elements
cur_num=0;

while size(r,1)>0 %when number of rows > 0
    point=[r(1),c(1)];
    cur=point;
    BW1(point(1),point(2))=0; %make the pixel black
    [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1); 
                               %find if any pixel around the current point is an edge pixel
    while size(I,1)>0 %if number of row > 0
        dist=(I-Gap_size-1).^2+(J-Gap_size-1).^2;
        [min_dist,index]=min(dist);
        p=point+[I(index),J(index)];
        point = p-Gap_size-1; % next is the current point
        cur=[cur;point]; %add point to curve 
        BW1(point(1),point(2))=0;%make the pixel black
        [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
                                %find if any pixel around the current point 
                                %is an edge pixel
    end
    
    % Extract edge towards another direction
    point=[r(1),c(1)];
    BW1(point(1),point(2))=0;
    [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    
    while size(I,1)>0
        dist=(I-Gap_size-1).^2+(J-Gap_size-1).^2;
        [min_dist,index]=min(dist);
        point=point+[I(index),J(index)]-Gap_size-1;
        cur=[point;cur];
        BW1(point(1),point(2))=0;
        [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    end
        
    if size(cur,1)>(size(BW,1)+size(BW,2))/25 % for 512 by 512 image, choose curve if its length > 40
        cur_num=cur_num+1;                    % One can change this value to control the length of the extracted edges
        curve{cur_num}=cur-Gap_size;
    end
    [r,c]=find(BW1==1);
    
end

for i=1:cur_num
    curve_start(i,:)=curve{i}(1,:);
    curve_end(i,:)=curve{i}(size(curve{i},1),:);
    if (curve_start(i,1)-curve_end(i,1))^2+...
        (curve_start(i,2)-curve_end(i,2))^2<=25  %if curve's ends are within sqrt(32) pixels
        curve_mode(i,:)='loop';
    else
        curve_mode(i,:)='line';
    end
    BW_edge(curve{i}(:,1)+(curve{i}(:,2)-1)*L)=1;
end
%%%
if cur_num>0
    TJ = find_TJunctions(curve, cur_num, Gap_size+1); % if a contour goes just outsize of ends, i.e., outside of gapsize we note a T-junction there
else
    curve{1} = [];
    curve_start = [];
    curve_end = [];
    curve_mode = [];
    cur_num = [];
    TJ = [];    
end
%%%
img=~BW_edge;
%%%%%%%%%%%%%%%%%%%%%%%%%%%


% find T-junctions within (gap by gap) neighborhood, where gap = Gap_size +
% 1; edges were continued (see edge_extraction procedure) when ends are within (Gap_size by Gap_size)
function TJ = find_TJunctions(curve, cur_num, gap); % finds T-Junctions in planar curves
TJ = [];
for i = 1:cur_num
    cur = curve{i};
    szi = size(cur,1);
    for j = 1:cur_num
        if i ~= j
            temp_cur = curve{j};
            compare_send = temp_cur - ones(size(temp_cur, 1),1)* cur(1,:);
            compare_send = compare_send.^2;
            compare_send = compare_send(:,1)+compare_send(:,2);
            if min(compare_send)<=gap*gap       % Add curve strat-points as T-junctions using a (gap by gap) neighborhood
                TJ = [TJ; cur(1,:)];
            end
            
            compare_eend = temp_cur - ones(size(temp_cur, 1),1)* cur(szi,:);
            compare_eend = compare_eend.^2;
            compare_eend = compare_eend(:,1)+compare_eend(:,2);
            if min(compare_eend)<=gap*gap       % Add end-points T-junctions using a (gap by gap) neighborhood
                TJ = [TJ; cur(szi,:)];
            end
        end
    end
end
%%%

% Compare T-junctions with obtained corners and add T-junctions to corners
% which are far away (outside a 5 by 5 neighborhood) from detected corners
function [corner_final c3] = Refine_TJunctions(corner_out,TJ,c2,curve, curve_num, curve_start, curve_end, curve_mode,EP);
%corner_final = corner_out;
c3=c2;

%%%%% add end points
if EP
    corner_num = size(corner_out,1);
    for i=1:curve_num
            if size(curve{i},1)>0 & curve_mode(i,:)=='line'

                % Start point compare with detected corners
                compare_corner=corner_out-ones(size(corner_out,1),1)*curve_start(i,:);
                compare_corner=compare_corner.^2;
                compare_corner=compare_corner(:,1)+compare_corner(:,2);
                if min(compare_corner)>100       % Add end points far from detected corners 
                    corner_num=corner_num+1;
                    corner_out(corner_num,:)=curve_start(i,:);
                    c3 = [c3;8];
                end

                % End point compare with detected corners
                compare_corner=corner_out-ones(size(corner_out,1),1)*curve_end(i,:);
                compare_corner=compare_corner.^2;
                compare_corner=compare_corner(:,1)+compare_corner(:,2);
                if min(compare_corner)>100
                    corner_num=corner_num+1;
                    corner_out(corner_num,:)=curve_end(i,:);
                    c3 = [c3;9];
                end
            end
    end
end
%%%%%%%%%%%%%%%5

%%%%%Add T-junctions
corner_final = corner_out;
for i=1:size(TJ,1)
    % T-junctions compared with detected corners
    if size(corner_final)>0
        compare_corner=corner_final-ones(size(corner_final,1),1)*TJ(i,:);
        compare_corner=compare_corner.^2;
        compare_corner=compare_corner(:,1)+compare_corner(:,2);
        if min(compare_corner)>100       % Add end points far from detected corners, i.e. outside of 5 by 5 neighbor
            corner_final = [corner_final; TJ(i,:)];
            c3 = [c3;10];
        end
    else
        corner_final = [corner_final; TJ(i,:)];
        c3 = [c3;10];
    end
end

corner_last=[];
corner_final2=corner_final;
flag(1:size(corner_final2,1),1)=1;
for k=1:size(corner_final2,1)
    if(flag(k)==0)
        continue;
    end
    comp_corner=corner_final2(k,:);
    compare_corner=corner_final2-ones(size(corner_final2,1),1)*comp_corner;
    compare_corner=compare_corner.^2;
    compare_corner=compare_corner(:,1)+compare_corner(:,2);
    mini_index=find(compare_corner<50&flag==1&compare_corner>0);
    if(size(mini_index,1)>0)
        for m=1:size(mini_index)
            corner2=corner_final2(mini_index(m),:);
            valid_corner=[floor((comp_corner(1,1)+corner2(1,1))/2),floor((comp_corner(1,2)+corner2(1,2))/2),];
            flag(k)=0;
            flag(mini_index(m))=0;
            corner_last=[corner_last;valid_corner];
        end
    else
        corner_last=[corner_last;comp_corner];
    end
end
corner_final=corner_last;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% show corners into the output images or into the edge-image
function img1=mark(img,x,y,w)
[M,N,C]=size(img);
img1=img;
if isa(img,'logical')
    img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)=...
        (img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)<1);
    img1(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:)=...
        img(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:);
else
    img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)=...
        (img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)<128)*255;
    img1(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:)=...
        img(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
function ang=curve_tangent(cur,center) % center is always the position of the corresponding extrema in cur

for i=1:2
    if i==1
        curve=cur(center:-1:1,:);
    else
        curve=cur(center:size(cur,1),:);
    end
    L=size(curve,1);
    
    if L>3
        if sum(curve(1,:)~=curve(L,:))~=0 % if not collinear
            M=ceil(L/2);
            x1=curve(1,1);
            y1=curve(1,2);
            x2=curve(M,1);
            y2=curve(M,2);
            x3=curve(L,1);
            y3=curve(L,2);
        else
            M1=ceil(L/3);
            M2=ceil(2*L/3);
            x1=curve(1,1);
            y1=curve(1,2);
            x2=curve(M1,1);
            y2=curve(M1,2);
            x3=curve(M2,1);
            y3=curve(M2,2);
        end
        
        if abs((x1-x2)*(y1-y3)-(x1-x3)*(y1-y2))<1e-8  % straight line
            tangent_direction=angle(complex(curve(L,1)-curve(1,1),curve(L,2)-curve(1,2)));
        else
            % Fit a circle 
            x0 = 1/2*(-y1*x2^2+y3*x2^2-y3*y1^2-y3*x1^2-y2*y3^2+x3^2*y1+y2*y1^2-y2*x3^2-y2^2*y1+y2*x1^2+y3^2*y1+y2^2*y3)/(-y1*x2+y1*x3+y3*x2+x1*y2-x1*y3-x3*y2);
            y0 = -1/2*(x1^2*x2-x1^2*x3+y1^2*x2-y1^2*x3+x1*x3^2-x1*x2^2-x3^2*x2-y3^2*x2+x3*y2^2+x1*y3^2-x1*y2^2+x3*x2^2)/(-y1*x2+y1*x3+y3*x2+x1*y2-x1*y3-x3*y2);
            % R = (x0-x1)^2+(y0-y1)^2;

            radius_direction=angle(complex(x0-x1,y0-y1));
            if radius_direction<0
                radius_direction = 2*pi-abs(radius_direction);
            end
            
            adjacent_direction=angle(complex(x2-x1,y2-y1));
            
            if adjacent_direction<0
                adjacent_direction = 2*pi-abs(adjacent_direction);
            end
            
            tangent_direction=sign(sin(adjacent_direction-radius_direction))*pi/2+radius_direction;
            if tangent_direction<0
                tangent_direction = 2*pi-abs(tangent_direction);
            elseif tangent_direction>2*pi
                tangent_direction = tangent_direction-2*pi;
            end
        end
    
    else % very short line
        tangent_direction=angle(complex(curve(L,1)-curve(1,1),curve(L,2)-curve(1,2)));
    end
    direction(i)=tangent_direction*180/pi;
end
ang=abs(direction(1)-direction(2));
%%%%%%%%%%%%%%%%%%5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parses the inputs into input_parameters
function [I,Hg,Lo,Gap_size,EP,ChordLength,Th] = parse_inputs(varargin);

error(nargchk(0,5,nargin));
Para=[1,0,0.35,0,5,9.86]; %Default experience value;
if nargin>2
    I=varargin{1};
    for i=2:nargin
        if size(varargin{i},1)>0
            Para(i-1)=varargin{i};
        end
    end
end

if nargin==1
    I=varargin{1};
end

if nargin==2
    I=varargin{1};
    args=varargin{2};
    Para=args';
end
    
if nargin==0 | size(I,1)==0
    [fname,dire]=uigetfile('*.bmp;*.jpg;*.gif','Open the image to be detected');
    I=imread([dire,fname]);
end


Gap_size = Para(1);
EP = Para(2);
Hg = Para(3); % high edge detection threshold
Lo = Para(4); % low edge detection threshold
ChordLength= Para(5);
Th=Para(6);
%%%%%%%%%%%%%%%%%%%%%%%%

function [G W] = makeGFilter(sig);

GaussianDieOff = .0001; 
pw = 1:100;

ssq = sig*sig;
W = max(find(exp(-(pw.*pw)/(2*ssq))>GaussianDieOff));
if isempty(W)
    W = 1;  
end
t = (-W:W);
gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq); 
G=gau/sum(gau);
