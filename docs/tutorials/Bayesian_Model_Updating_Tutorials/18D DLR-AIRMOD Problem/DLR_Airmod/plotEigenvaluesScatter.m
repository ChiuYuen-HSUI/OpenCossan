function plotEigenvaluesScatter(varargin)

MMEigs=varargin{1};
ExpEigs=varargin{2};


assert(size(MMEigs,2)==size(ExpEigs,2),'Matrices must have the same number of columns.')

Neigs = size(MMEigs,2);


if nargin==3
    FEEigs=varargin{3};
    assert(size(FEEigs,2)==size(ExpEigs,2),'Matrices must have the same number of columns.')
    
    [lh1,ah1]=gplotmatrix([MMEigs(:,1:end-1);ExpEigs(:,1:end-1);FEEigs(:,1:end-1)],...
    [MMEigs(:,2:end);ExpEigs(:,2:end);FEEigs(:,2:end)],...
    [zeros(length(MMEigs),1);ones(length(ExpEigs),1);2*ones(length(FEEigs),1)],...
    'brg','...',[3,3,3]);

else
   
% for i= 1:Neigs
%     for j= i+1:Neigs
%         subplot(Neigs-1,Neigs-1,(j-2)*(Neigs-1) + i)
%         hold on
%         scatter(FEEigs(:,i),FEEigs(:,j),'.b')
%         scatter(ExpEigs(:,i),ExpEigs(:,j),'.r')
%         hold off
%     end
% end

[lh1,ah1]=gplotmatrix([MMEigs(:,1:end-1);ExpEigs(:,1:end-1)],...
    [MMEigs(:,2:end);ExpEigs(:,2:end)],...
    [zeros(length(MMEigs),1);ones(length(ExpEigs),1)],...
    'br','..',[8,8]);

end

for i= 1:Neigs-1
    for j= i+1:Neigs-1
        delete(ah1(i,j))
    end
end

end