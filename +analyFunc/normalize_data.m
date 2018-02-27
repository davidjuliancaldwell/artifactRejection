function [ Z ] = normalize_data( data,ref )
%NORMALIZE_DATA normalize the data relative to a reference period, based
%off of normalize PLV from felix 

%   need reference data and 

m_ref1=mean(ref,2);
s_ref1=std(ref')';
Z=(data-m_ref1*ones(1,size(data,2)))./(s_ref1*ones(1,size(data,2)));


end

