

function [radgrad, axialgrad] = gradCentralDifference_adb(varb, delx1, delx2)

nx=size(varb,1); ny=size(varb,2);
radgrad = zeros(nx,ny);
axialgrad = zeros(nx,ny);

    % radial gradient
    radgrad(1,:) = (1/delx2) * (varb(2,:)-varb(1,:));
    for i=2:nx-1
        radgrad(i,:) = (1/(2*delx2)) * (varb(i+1,:)-varb(i-1,:));
%         % take care of solid boundaries
%         if isnan(varb(i+1)) == 1
%             radgrad(i,:) = (1/(delx2)) * (varb(i,:)-varb(i-1,:));
%         end
%         if isnan(varb(i-1)) == 1
%             radgrad(i,:) = (1/(delx2)) * (varb(i+1,:)-varb(i,:));
%         end
        
    end
    radgrad(nx,:) = (1/delx2) * (varb(nx,:)-varb(nx-1,:));
    
    % axial gradient
    axialgrad(:,1) = (1/delx1) * (varb(:,2)-varb(:,1));
    for j=2:ny-1
        axialgrad(:,j) = (1/(2*delx1)) * (varb(:,j+1)-varb(:,j-1));
    end
    axialgrad(:,ny) = (1/delx1) * (varb(:,ny)-varb(:,ny-1));
end