function [sources,ialp,gamma,pow,wts] = champagne_wrapper(sess,t1,t2,t3,t4,vcs,init,ax)
%CHAMPAGNE_WRAPPER  Wrapper for champagne analysis.
%   CHAMPAGNE_WRAPPER(SESS,T1,T2,T3,T4,VCS) is a wrapper for run_champ_code
%   which automatically creates an output directory.  SESS is the full path
%   of the session file (i.e. '/data/path/to/sess_listen_av1.mat').  T1 and
%   T2 are the start and end timepoints for the baseline; T3 and T4 are the
%   start and end timepoints for the signal to be analyzed.
%
%CN 3/2012

if nargin < 1, error('First argument must be path to session file.'); end
if nargin < 2, t1 = -100; end
if nargin < 3, t2 = -5; end
if nargin < 4, t3 = 25; end
if nargin < 5, t4 = 250; end
if nargin < 6, vcs = 2; end
if nargin < 7, init = 1; end
if nargin < 8, ax = [0 0]; end

nuts = load(sess);
[pathstr,sessname] = fileparts(sess);

trials=size(nuts.meg.data,3);
curr_dir = pwd;
output_dir = fullfile(pathstr,sprintf('%s_pre%dto%dms_post%dto%dms_vcs%d',sessname,t1,t2,t3,t4,vcs));

bSave = savecheck(output_dir,'dir');

if ~bSave
    [sources,ialp,gamma,pow,wts] = deal([]);
    disp('Directory exists: Champagne not run.')
    return
end

mkdir(output_dir)
cd(output_dir)
if length(trials)==1, trials=1:trials; end    % use this on local server

[sources, ialp, gamma, pow, wts] = run_champ_code(nuts,trials,t1,t2,t3,t4,init,vcs,ax(1),ax(2));

cd(curr_dir)