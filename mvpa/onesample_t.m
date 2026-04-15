function onesample_t = onesample_t(X)
    % Compute t-map on X
    % X: [nSubjects x nVoxels]
    % tvec: [1 x nVoxels]
    n  = size(X,1);
    mu = mean(X,1);
    sd = std(X,0,1);           % sample std (dof = n-1)
    se = sd ./ sqrt(n);
    onesample_t = zeros(1,size(X,2));
    nz = se>0 & isfinite(se);
    onesample_t(nz) = mu(nz) ./ se(nz);
    % if se==0 leave t=0
end
