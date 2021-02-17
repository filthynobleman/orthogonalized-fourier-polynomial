function drawaw(S,mtl)
%% drawawobj
% Draw S with material M

if isfield(S,'f3');p3=patch('Vertices',S.v','Faces',S.f3');end
if isfield(S,'f4');p4=patch('Vertices',S.v','Faces',S.f4');end
if isfield(S,'f5');p5=patch('Vertices',S.v','Faces',S.f5');end
if isfield(S,'f6');p6=patch('Vertices',S.v','Faces',S.f6');end

if nargin==1
    if isfield(S,'f3') 
        p3.FaceColor=[0 1 1];
        p3.EdgeColor=[0 0 0];
    end
    if isfield(S,'f4') 
        p4.FaceColor=[0 1 0];
        p4.EdgeColor=[0 0 0];
    end
    if isfield(S,'f5') 
        p5.FaceColor=[1 0 0];
        p5.EdgeColor=[0 0 0];
    end
    if isfield(S,'f6') 
        p6.FaceColor=[1 1 0];
        p6.EdgeColor=[0 0 0];
    end
else
    if isfield(S,'f3') 
        for ii=1:length(S.umat3)
            mtlnum=S.umat3(ii);
            for jj=1:length(mtl)
                if strcmp(mtl(jj).name,S.usemtl(mtlnum-1));
                    break;
                end
            end
            fvcd3(ii,:)=mtl(jj).Kd';
            %fvcd(ii,:)=rand(1,3);
        end
        p3.FaceVertexCData=fvcd3;
        p3.FaceColor='flat';
        p3.EdgeColor='none'; 
    end
    if isfield(S,'f4') 
        for ii=1:length(S.umat4)
            mtlnum=S.umat4(ii);
            for jj=1:length(mtl)
                if strcmp(mtl(jj).name,S.usemtl(mtlnum-1));
                    break;
                end
            end
            fvcd4(ii,:)=mtl(jj).Kd';
            %fvcd(ii,:)=rand(1,3);
        end
        p4.FaceVertexCData=fvcd4;
        p4.FaceColor='flat';
        p4.EdgeColor='none'; 
    end
    if isfield(S,'f5') 
        for ii=1:length(S.umat5)
            mtlnum=S.umat5(ii);
            for jj=1:length(mtl)
                if strcmp(mtl(jj).name,S.usemtl(mtlnum-1));
                    break;
                end
            end
            fvcd5(ii,:)=mtl(jj).Kd';
            %fvcd(ii,:)=rand(1,3);
        end
        p5.FaceVertexCData=fvcd5;
        p5.FaceColor='flat';
        p5.EdgeColor='none'; 
    end
    if isfield(S,'f6') 
        for ii=1:length(S.umat6)
            mtlnum=S.umat6(ii);
            for jj=1:length(mtl)
                if strcmp(mtl(jj).name,S.usemtl(mtlnum-1));
                    break;
                end
            end
            fvcd6(ii,:)=mtl(jj).Kd';
            %fvcd(ii,:)=rand(1,3);
        end
        p6.FaceVertexCData=fvcd6;
        p6.FaceColor='flat';
        p6.EdgeColor='none'; 
    end
end
    

axis('equal')