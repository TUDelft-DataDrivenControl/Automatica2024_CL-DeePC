function hankel_matrix = data2hankel(s,data,varargin)
% creates hankel matrix with s x N blocks
%   s    = number of block rows
%   data = data for which to construct block-hankel matrix
%   varargin{1} = number of columns (can truncate data set, or returns
%                 error if there is insufficient data)
%   varargin{2} = use last available data? Bool, default=True

% finding matrix shape
data_size = size(data);
if data_size(2) < data_size(1)
    data = data.';
    data_size = size(data);
end

% find sizes of data
n_dim = data_size(1); % number of inputs/outputs/states/etc
N_bar = data_size(2); % number of data points in data

% define number of columns (N) and select requested part of data
if isempty(varargin) % use all data
    N = N_bar-s+1;      % N_bar = s+N-1

else % select data to use
    N = varargin{1}; % number of columns -> deterimines amount of data

    if N_bar < N+s-1 % check for sufficient data, otherwise raise error
        error(['Not enough data to construct block-Hankel matrix with ', ...
               num2str(s),' block rows and ',num2str(N) ' columns.\n%s'], ...
               [num2str(N+s-1-N_bar),' additional samples needed.'])
    
    else % sufficient data available
        N_bar = N+s-1; % adjust amount of data to use

        % determine whether to use most recent (UseLast=true)
        %  or oldest data (UseLast=false)
        if length(varargin) > 1
            UseLast = varargin{2};
        else
            UseLast = true; % Default
        end
        
        % select data to use:
        if UseLast % use most recent data
            data = data(:,end-(N_bar-1):end);
        else       % use oldest data
            data = data(:,1:N_bar);
        end
    end
end

% create hankel matrix
hankel_matrix = zeros(n_dim*s,N);

% fill hankel matrix
rows = 1:n_dim;
for i = 1:s
    hankel_matrix(rows,:) = data(:,i:i+N-1);
    rows = rows+n_dim;
end

end

