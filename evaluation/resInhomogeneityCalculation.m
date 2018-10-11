function [residuals, originals] = resInhomogeneityCalculation(estField, targetField)
if(~isequal(size(estField),size(targetField)))
    error('Input field sizes must be equal.')
end

nROI    = length(estField);
const   = 10^(floor(log10(std(targetField))) - 3);

for j=1:100
    k = const;
    while(100*sum(abs((targetField-estField)-mean(targetField-estField))< k)/nROI < j)
        k= k + const;
    end
    ppmCorr(j) = k;
    
    k = const;
    while(100*sum(abs((targetField)-mean(targetField))< k)/nROI < j)
        k= k + const;
    end
    ppmBB(j) = k;
end

residuals = [ppmCorr(95) ppmCorr(90) ppmCorr(80)];
originals = [ppmBB(95) ppmBB(90) ppmBB(80)];
end