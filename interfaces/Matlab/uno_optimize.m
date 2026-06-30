% UNO_OPTIMIZE - UNO solver
%   Call the UNO solver and solve the specified optimization model.
%
%   Syntax
%     result = UNO_OPTIMIZE(model)
%     result = UNO_OPTIMIZE(model,options)
%     result = UNO_OPTIMIZE(model,options,callbacks)
%
%   Input Arguments
%     model - Optimization model
%       structure
%     options - UNO options, as returned by <a href="matlab:helpPopup('uno_options')">uno_options</a>
%       structure
%     callbacks - UNO solver callbacks
%       structure
%
%   Output Arguments
%     result - UNO result
%       structure
%
%   Examples
%     <a href="matlab:open([fileparts(which('uno')) filesep 'example/example_options.m'])">example_options</a>
%     <a href="matlab:open([fileparts(which('uno')) filesep 'example/example_hs015.m'])">example_hs015</a>
%     <a href="matlab:open([fileparts(which('uno')) filesep 'example/example_polak5.m'])">example_polak5</a>
%
%   See also <a href="matlab:helpPopup('uno_options')">uno_options</a> <a href="matlab:helpPopup('uno')">uno</a>
%
% Copyright (c) 2025 Stefano Lovato and Charlie Vanaret
% Licensed under the MIT license. See LICENSE file in the project directory for details.