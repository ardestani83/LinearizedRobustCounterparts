% ROME_VAR Describes a variable in ROME
%
%
% Modification History: 
% 1. Joel 
% 

classdef rome_var
%     % CONSTANTS
%     % ----------
%     properties (Constant = true, Hidden = true)
%         % Constants representing Conic Constraints
%         NO_CONE  = -1;  % No Cone Restrictions
%         ZERO     =  0;  % Restricted to 0
%         NNOC     =  1;  % Nonnegative Cone
%         SOC      =  2;  % Second Order Cone
%         
%         % Constant representing Continuity Constraints
%         CONTINUOUS  = 0; % Continuous Variable
%         INTEGER     = 1; % Integer Variable
%         BINARY      = 2; % Binary Variable
%     end
    % PRIVATE PROPERTIES
    % ------------------
    % Note to Joel: All properties are currently public. To be made
    % protected on beta release    
    properties % (SetAccess = protected, GetAccess = protected) % TO BE REMOVED LATER
        % Affine Mapping (should usually be sparse) from the global model
        % variable to this user variable. 
        % AffineMap has the form [b | A];
        BiAffineMap;             % Stores the Affine Mapping from global table into current variable
    end
    properties(Dependent = true)
        NumMappedVars;  % Number of mapped variables (excluding constant)
    end
    properties
        NumUnmappedVars = 0;     % Number of variables (before first one) that are unmapped        
        NumMappedRandVars   = 0; % Excludes constant
        NumUnmappedRandVars = 0; % Offset into uncertainty table
        
        Continuity = rome_constants.CONTINUOUS;% Refer to Constants
        Cone = rome_constants.NO_CONE;         % Refer to Constants
        DiagMult = NaN;                        % Diagonal Multiplier (Used to denote diagonal matrices)
        
        % For SQ optimization
        NonLinearity = '';               % Type of NonLinearity        
        ExtraData = [];                  
    end
    
    properties (Hidden = true, SetAccess = public)
        % Mean = 0;       % Mean of variable (deprecated as at v1.0.8)
        Covar= [];      % Covariance Matrix
        FDev = [];      % Forward Deviation
        BDev = [];      % Backward Deviation
    end
    
    properties (Dependent = true)
        TotalSize;      % Number of output variables (scalar integer)
        IsDiag;         % Flag to check if AffineMap is diagonal (for optimization)
        IsRand;         % Flag to check if the variable is pure uncertain
        IsCertain;      % Flag to check if the variable is pure certain
        IsLDR;          % Flag to check if the variable is an LDR
        IsConst;        % Flag to check if the variable is a constant
    end
    properties %(SetAccess = protected)
        Size = [];  % Vector representing size
    end
    
    % METHODS
    % --------
    methods
        % Constructor
        % ------------
        function obj = rome_var(varargin)
            % implement a copy-constructor
%             if(nargin == 1 && isa(varargin{1}, 'rome_var'))
%                 old_obj = varargin{1};
%                 obj.BiAffineMap = old_obj.BiAffineMap;
%                 obj.NumUnmappedVars = old_obj.NumUnmappedVars;
%                 obj.NumMappedRandVars  = old_obj.NumMappedRandVars;
%                 obj.NumUnmappedRandVars = old_obj.NumUnmappedRandVars;
%                 obj.Continuity = old_obj.Continuity;
%                 obj.Cone       = old_obj.Cone;
%                 obj.DiagMult   = old_obj.DiagMult;
%                 obj.Size       = old_obj.Size;
%                 return;
%             end
            
            % initialize number of size arguments
            num_size_args = 0;
            
            % check if there were optional arguments
            for ii = 1:nargin
                if(~ischar(varargin{ii}))
                    num_size_args = ii;
                else break;
                end
            end
            
            % parse optional arguments
            if(num_size_args ~= nargin)
                for ii = (num_size_args + 1):2:nargin
                    if(~ischar(varargin{ii}) || ~isnumeric(varargin{ii+1}))
                        error('rome_var:ctor:InvalidOps', 'Optional Arguments must be entered in pairs');
                    else
                        eval(sprintf('obj.%s = %g;', varargin{ii}, varargin{ii+1}));
                    end
                end
            end
            
            % if number of size args = 0, we create a single (scalar)
            if(num_size_args == 0)
                obj.Size = [1, 1];
           
            % if number of size args == 1, special case, will make dim = 2
            % instead
            elseif(num_size_args == 1)
                dim_size = varargin{1};
                if(isscalar(dim_size))
                    obj.Size = [floor(dim_size), 1];
                else
                    obj.Size = floor(dim_size);
                end
                
            % if number of size args > 1, each argument represents the
            % size of the variable in each dimension respectively.
            elseif(num_size_args > 1)
                obj.Size = zeros(1, num_size_args);
                
                % iterate over each argument and allocate memory
                for ii = 1:num_size_args
                    dim_size = varargin{ii};
                    if(isscalar(dim_size))
                        obj.Size(ii) = floor(dim_size);
                    else
                        warning('rome_var:ctor:InvalidArg', ...
                            'Constructor requires scalar integer arguments');
                    end 
                end % end for each argument
            end
        end
        
        % PROPERTY ACCESS
        function val = get.TotalSize(obj)
            val = prod(obj.Size);   % total size formed from product of individual dimensions
        end        
        
        function val = get.NumMappedVars(obj)
            val = (size(obj.BiAffineMap, 2) ./ (obj.NumMappedRandVars + 1)) - 1;   % Note: change to not include affine term
        end 
        
        function val = get.IsDiag(obj)
            val = all(~isnan(obj.DiagMult));
        end
        
        function val = get.IsCertain(obj)
            val = (obj.NumMappedRandVars == 0);
        end
        
        function val = get.IsRand(obj)
            val = (obj.NumMappedVars == 0);
        end
        
        function val = get.IsLDR(obj)
            val = ~(obj.IsRand || obj.IsCertain);
        end
        function val = get.IsConst(obj)
            val = obj.IsRand && obj.IsCertain;
        end
        
        % might want to take this out at some point
        function obj = set.BiAffineMap(obj, A)
            obj.BiAffineMap = A;
            obj.DiagMult = NaN;
        end
        
        % B) COVARIANCE
        % --------------
        % Setting
        function obj = set.Covar(obj, val)
            % Debugging
%             disp('Setting Covar');
            
            % Error check uncertainness
            if(~obj.IsRand)
                error('rome_var:Covar:ReqRandDiag', 'Object must be uncertain to set covariance');
            end
            
            % expand scalar
            if(isscalar(val))
                val = repmat(val, obj.TotalSize, 1);
            end
            
            % expand vector
            if(isvector(val))
                if(obj.TotalSize ~= numel(val))
                    error('rome_var:Covar:WrongInputSize', 'For vector input, supplied covariance must be the same size as the object');
                end
                val = spdiags(val(:), 0, obj.TotalSize, obj.TotalSize); % notice vectorization to make into column vector
            end
            
            % error check size
            if(obj.TotalSize ~= sqrt(numel(val)))
                error('rome_var:Covar:WrongInputSize', 'Supplied covariance must a square matrix with each side having the same size as the object');
            end
                        
            % add to global model (vectorize for convenience)
%             rome_get_current_model().set_cov(obj.NumUnmappedRandVars, val, obj.BiAffineMap(:, 2:end));
            rome_get_current_model().set_cov(obj(:), val);
        end

        % Getting
        function [cov_mat, cov_mix] = get.Covar(obj)
            % Debugging
%             disp('Getting Covar');

            % temp:
            if(~obj.IsDiag)
                error('rome_var:Covar:DirectGet', 'Cannot get covariance directly');
            end

            % Error check uncertainness
            if(~obj.IsRand)
                error('rome_var:Covar:ReqRand', 'Object must be uncertain to get covariance');
            end

            % get from global variable
            h = rome_get_current_model();
            [cov_mat, cov_mix] = h.get_cov(obj);

            % insert value and return
%             val = obj.BiAffineMap(:, 2:end) * cov_mat * (obj.BiAffineMap(:, 2:end))';
        end
        
        % C) Forward Deviation
        % --------------
        % Setting
        function obj = set.FDev(obj, val)
%              % Error check diagonality and uncertainness
%             if(~obj.IsDiag || ~obj.IsRand)
%                 error('rome_var:FDev:ReqRandDiag', 'Object must be uncertain and diagonal to set Forward Deviation');
%             end
            
             % Error check uncertainness
            if(~obj.IsRand)
                error('rome_var:FDev:ReqRand', 'Object must be uncertain to set Forward Deviation');
            end
            
            % expand scalars
            if(isscalar(val))
                val = repmat(val, obj.TotalSize, 1);
            end
            
            % add to global model (vectorize for convenience)
            % rome_get_current_model().set_dirdev(obj.NumUnmappedRandVars, val(:), 1);
            rome_get_current_model().set_dirdev(obj(:), val(:), 1);            
        end
        % Getting
        function val = get.FDev(obj)
            % Debugging
%             disp('Getting FDev');

            % Error check uncertainness
            if(~obj.IsRand)
                error('rome_var:FDev:ReqRand', 'Object must be uncertain to get Forward Deviation');
            end

            % get from global variable
            h = rome_get_current_model();
            primitive_pdev = h.get_dirdev(obj, 1);

            % potential bug: how does directional deviation change with
            % affine composion?            
            
            % Just substitute for now
            val = primitive_pdev;
        end
        
        % D) Backward Deviation
        % ----------------------
        % Setting
        function obj = set.BDev(obj, val)
            % Debugging
%             disp('Setting BDev');
            
%              % Error check diagonality and uncertainness
%             if(~obj.IsDiag || ~obj.IsRand)
%                 error('rome_var:BDev:ReqRandDiag', 'Object must be uncertain and diagonal to set Backward Deviation');
%             end
            
             % Error check uncertainness
            if(~obj.IsRand)
                error('rome_var:BDev:ReqRand', 'Object must be uncertain to set Backward Deviation');
            end
            
            % expand scalars
            if(isscalar(val))
                val = repmat(val, obj.TotalSize, 1);
            end
            
            % add to global model (vectorize for convenience)
            % rome_get_current_model().set_dirdev(obj.NumUnmappedRandVars, val(:), -1);
            rome_get_current_model().set_dirdev(obj(:), val(:), -1);
        end
        % Getting
        function val = get.BDev(obj)
            % Debugging
%             disp('Getting BDev');
            
            % Error check uncertainness
            if(~obj.IsRand)
                error('rome_var:BDev:ReqRand', 'Object must be uncertain to get Backward Deviation');
            end

            % get from global variable
            h = rome_get_current_model();
            primitive_qdev = h.get_dirdev(obj, -1);

            % potential bug: how does directional deviation change with
            % affine composion?            
            
            % Just substitute for now
            val = primitive_qdev;
        end
    end
end


% ROME: Copyright (C) 2009 by Joel Goh and Melvyn Sim
% See the file COPYING.txt for full copyright information.
