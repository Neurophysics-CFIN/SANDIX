function varargout = Get_uniformly_distributed_SMEX_prior(Mtrain,lb,ub)
if nargout ~= length(lb)
    error('Syntax error: bounds and number of output arguments inconsistent')
end
varargout = cell(length(lb));
X = rand(Mtrain,length(lb)) .* (ub(:)' - lb(:)') + lb(:)';
for j = 1:length(lb)
    varargout{j} = X(:,j);
end