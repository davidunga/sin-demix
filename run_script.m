function ret = run_script(mfile, opts)

arguments
    mfile
    opts.allow_figs = 0
end

assert(exist(mfile,"file"));

pre_existing = struct();
pre_existing.figs = allchild(0);
pre_existing.opts = opts;
pre_existing.mfile = mfile;
clear mfile opts;
assert(length(who())==1);

run(pre_existing.mfile);

if ~pre_existing.opts.allow_figs
    close(setdiff(allchild(0),pre_existing.figs));
end

clear pre_existing;
ret = vars2struct(who());