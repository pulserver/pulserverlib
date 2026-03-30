function run_generators()
this_dir = fileparts(mfilename('fullpath'));
scripts_dir = fullfile(this_dir, 'scripts');
pulseq_matlab_dir = fullfile(this_dir, 'pulseq', 'matlab');

addpath(this_dir);
addpath(scripts_dir);

if isfolder(pulseq_matlab_dir) || isdir(pulseq_matlab_dir)
    addpath(pulseq_matlab_dir);
else
    error('Pulseq MATLAB library not found at: %s', pulseq_matlab_dir);
end

scripts = dir(fullfile(scripts_dir, 'generate_*.m'));
names = {scripts.name};
[~, script_order] = sort(names);

for idx = 1:numel(script_order)
    script_name = scripts(script_order(idx)).name;
    script_path = fullfile(scripts_dir, script_name);
    fprintf('Running %s\n', script_name);
    run_generator_script(script_path);
end
end

function run_generator_script(script_path)
run(script_path);
end

run_generators();