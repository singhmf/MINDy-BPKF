function[X]=NanToZero(X)
%% Replaces nans with zeros
X(isnan(X))=0;
end