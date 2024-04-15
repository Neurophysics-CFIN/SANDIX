function S = powderaverage_finitePulses(pars,info,use_ODF)
fs = pars(1);
Ds = pars(2);
Db = pars(3);
rs = pars(4);

if nargin<3
    if isfield(info,"use_ODF")
        use_ODF = info.use_ODF;
    else
        use_ODF = true;
    end
end

b = info.b;
Delta = info.Delta;
delta = info.delta;
if numel(delta)==1
    delta = delta*ones(length(b),1);
    Delta = Delta*ones(length(b),1);
end

persistent W bias

if use_ODF
    fb = 1-fs;
    rb = rs*fs/(fb+1e-4);
    
    % can solve differential equation for all eps simultanously
    [eps,w] = legpts(20,[0 1]);
    D = zeros(2*length(w),1);
    D(1:2:end) = Ds*eps.^2;
    D(2:2:end) = Db;
    D = diag(D);
    R = zeros(2*length(w));
    for j = 1:2:2*length(w)
        R(j:j+1,j:j+1) = [-rs  rb
                           rs -rb];
    end
    f = zeros(2*length(w),1);
    f(1:2:end) = fs;
    f(2:2:end) = fb;
    
    q = sqrt(b./(Delta-delta/3));
    
    S = zeros(length(b),length(w));
    for i = 1:length(b)
        St = fun.karger.finitePulses(f,D,R,delta(i),Delta(i),q(i));
        S(i,:) = St(1:2:end)+St(2:2:end);
    end
    S = S*w';
else
    if isempty(W)
        loaded = load("combined_net_5_b100.mat");
        net = loaded.net;
        for n = 1:5
            W{n}{1} = net.IW{1+(n-1)*4};
            bias{n}{1} = net.b{1+(n-1)*4};
            W{n}{2} = net.LW{2+(n-1)*4,1+(n-1)*4};
            bias{n}{2} = net.b{2+(n-1)*4};
            W{n}{3} = net.LW{3+(n-1)*4,2+(n-1)*4};
            W{n}{4} = net.IW{3+(n-1)*4};
            bias{n}{3} = net.b{3+(n-1)*4};
            W{n}{5} = net.LW{4+(n-1)*4,3+(n-1)*4};
            bias{n}{4} = net.b{4+(n-1)*4};
        end
        W{6}{1} = net.IW{22};
        bias{6}{1} = net.b{22};
        W{6}{2} = net.LW{23,22};
        bias{6}{2} = net.b{23};
    end
    X = [fs*ones(length(b),1) Ds*b Db*b rs*(Delta-delta/3) rs*delta]';
    X = (X-[0.5 15 15 25 5]')./[0.5 15 15 25 5]';
    Sp = zeros(length(b),numel(W));
    weights = softmax(W{end}{2}*tansig(W{end}{1}*X+bias{end}{1})+bias{end}{2});
    for n = 1:numel(W)-1
        Sp(:,n) = weights(n,:).*(W{n}{5}*tansig(W{n}{3}*tansig(W{n}{2}*tansig(W{n}{1}*X+bias{n}{1})+bias{n}{2})+bias{n}{3}+W{n}{4}*X)+bias{n}{4});
    end
    S = 0.5+0.5*sum(Sp,2);
end