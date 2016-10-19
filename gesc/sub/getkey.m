function [ch, tim] = getkey(N,nonascii)

% GETKEY - get a keypress
%   CH = GETKEY waits for a single keypress and returns the ASCII code. It
%   accepts all ascii characters, including backspace (8), space (32),
%   enter (13), etc, that can be typed on the keyboard.
%   Non-ascii keys (ctrl, alt, ..) return a NaN. CH is a double.
%
%   CH = GETKEY(N) waits for N keypresses and returns their ASCII codes.
%   GETKEY(1) is the same as GETKEY without arguments.
%
%   GETKEY('non-ascii') or GETKEY(N,'non-ascii') uses non-documented
%   matlab features to return a string describing the key pressed.
%   In this way keys like ctrl, alt, tab etc. can also distinguished.
%   The return is a string (when N = 1) or a cell array of strings.
%
%   [CH,T] = GETKEY(...) also returns the time between the start of the
%   function and each keypress. This is, however, not that accurate.
%
%   This function is kind of a workaround for getch in C. It uses a modal,
%   but non-visible window, which does show up in the taskbar.
%   C-language keywords: KBHIT, KEYPRESS, GETKEY, GETCH
%
%   Examples:
%
%    fprintf('\nPress any key: ') ;
%    ch = getkey ;
%    fprintf('%c\n',ch) ;
%
%    fprintf('\nPress the Ctrl-key within 3 presses: ') ;
%    ch = getkey(3,'non-ascii')
%    if ismemmber('control', ch),
%      fprintf('OK\n') ;
%    else
%      fprintf(' ... wrong keys ...\n') ;
%    end
%
%  See also INPUT, UIWAIT
%           GETKEYWAIT (File Exchange)

% for Matlab 6.5 and higher
% version 2.0 (jun 2012)
% author : Jos van der Geest
% email  : jos@jasen.nl
%
% History
% 1.0 2005 - creation
% 1.1 dec 2006 - modified lay-out and help
% 1.2 apr 2009 - tested for more recent MatLab releases
% 1.3 jan 2012 - modified a few properties, included check is figure still
%            exists (after comment on FEX by Andrew).
% 2.0 jun 2012 - added functionality to accept multiple key presses

t00 = tic ; % start time of this function

% check the input arguments
error(nargchk(0,2,nargin))
switch nargin
    case 0
        nonascii = '' ;
        N = 1 ;
    case 1
        if ischar(N),
            nonascii = N ;
            N = 1 ;
        else
            nonascii = '' ;
        end
end

if numel(N) ~= 1 || ~isnumeric(N) || N < 1 || fix(N) ~= N
    error('N should be a positive integer scalar.') ;
end

% Determine the callback string to use
if strcmpi(nonascii,'non-ascii'),
    % non-ascii characters are accepted
    nonascii = true ;
    callstr = 'set(gcbf,''Userdata'',get(gcbf,''Currentkey'')) ; uiresume ' ;
elseif isempty(nonascii)
    nonascii = false ;
    % only standard ascii characters are accepted
    callstr = 'set(gcbf,''Userdata'',double(get(gcbf,''Currentcharacter''))) ; uiresume ' ;
else
    error('String argument should be the string ''non-ascii''') ;
end

% Set up the figure
% May be the position property  should be individually tweaked to avoid visibility
fh = figure(...
    'name','Press a key', ...
    'keypressfcn',callstr, ...
    'windowstyle','modal',...
    'numbertitle','off', ...
    'position',[0 0  1 1],...
    'userdata','timeout') ;
try
    ch = cell(1,N) ;
    tim = zeros(1,N) ;
    
    % loop to get N keypresses
    for k=1:N
        % Wait for something to happen, usually a key press so uiresume is
        % executed
        uiwait ;
        tim(k) = toc(t00) ; % get the time of the key press
        ch{k} = get(fh,'Userdata') ;  % and the key itself
        if isempty(ch{k}),
            if nonascii
                ch{k} = NaN ;
            else
                ch{k} = '' ;
            end
        end
    end
    if ~nonascii
        ch = [ch{:}] ;
    else
        if N==1
            ch = ch{1} ; % return as a string
        end
        % return as a cell array of strings
    end
 catch
    % Something went wrong, return empty matrices.
    ch = [] ;
    tim = [] ;
end

% clean up the figure, if it still exists
if ishandle(fh)
    delete(fh) ;
end

% 
% 
% Copyright (c) 2005-prese, Jos van der Geest
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
