function s = vectorize(S,mask)
Ssz = size(S);
Sdims = length(Ssz);
masksz = size(mask);
if Sdims > 2 && all(size(mask) == Ssz(1:end-1))
    s = permute(S,[Sdims,1:Sdims-1]);
    s = s(:,mask)';
elseif ismatrix(S) && Ssz(1) == sum(mask(:))
    s = zeros([Ssz(2),masksz], "like",S);
    s(:,mask) = S';
    s = permute(s,[2:length(masksz) + 1, 1]);
end
end