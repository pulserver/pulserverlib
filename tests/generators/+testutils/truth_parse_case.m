function truth = truth_parse_case(base)
%TRUTH_PARSE_CASE Parse all TruthBuilder outputs for one test case.
%
% base can be:
%   - basename (e.g. 'gre_2d_1sl_1avg')
%   - full path without suffix
%   - .seq path
%   - one of the exported artifact paths

    [data_dir, base_name] = normalize_case_path(base);

    truth = struct();
    truth.base_name = base_name;
    truth.data_dir = data_dir;

    truth.paths = struct( ...
        'meta', fullfile(data_dir, [base_name '_meta.txt']), ...
        'tr_waveform', fullfile(data_dir, [base_name '_tr_waveform.bin']), ...
        'segment_def', fullfile(data_dir, [base_name '_segment_def.bin']), ...
        'freqmod_def', fullfile(data_dir, [base_name '_freqmod_def.bin']), ...
        'scan_table', fullfile(data_dir, [base_name '_scan_table.bin']), ...
        'seq', fullfile(data_dir, [base_name '.seq']));

    required = {'meta', 'tr_waveform', 'segment_def', 'freqmod_def', 'scan_table'};
    for i = 1:numel(required)
        key = required{i};
        p = truth.paths.(key);
        if ~isfile(p)
            error('truth:io', 'Missing required truth artifact: %s', p);
        end
    end

    truth.meta = testutils.truth_parse_meta(truth.paths.meta);
    truth.tr_waveforms = testutils.truth_parse_tr_waveform_bin(truth.paths.tr_waveform);
    truth.segment_def = testutils.truth_parse_segment_def_bin(truth.paths.segment_def);
    truth.freqmod_def = testutils.truth_parse_freqmod_def_bin(truth.paths.freqmod_def);
    truth.scan_table = testutils.truth_parse_scan_table_bin(truth.paths.scan_table);

    truth.validation = local_validate(truth);
end

function [data_dir, base_name] = normalize_case_path(base)
    if nargin < 1 || isempty(base)
        error('truth:input', 'A base case name or path is required.');
    end

    if isstring(base)
        base = char(base);
    end

    [in_dir, in_name, in_ext] = fileparts(base);
    if isempty(in_dir)
        data_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data');
    else
        data_dir = in_dir;
    end

    if strcmp(in_ext, '.seq')
        base_name = in_name;
        return;
    end

    suffixes = {'_meta', '_tr_waveform', '_segment_def', '_freqmod_def', '_scan_table'};
    base_name = in_name;
    for i = 1:numel(suffixes)
        sfx = suffixes{i};
        if endsWith(base_name, sfx)
            base_name = base_name(1:(length(base_name) - length(sfx)));
            return;
        end
    end

    if isempty(in_ext)
        base_name = in_name;
    else
        base_name = [in_name in_ext];
    end
end

function validation = local_validate(truth)
    warnings = {};

    if truth.meta.num_canonical_trs ~= truth.tr_waveforms.num_trs
        warnings{end + 1} = sprintf( ...
            'num_canonical_trs mismatch: meta=%d, tr_waveforms=%d', ...
            truth.meta.num_canonical_trs, truth.tr_waveforms.num_trs); %#ok<AGROW>
    end

    if truth.meta.num_segments ~= truth.segment_def.num_segments
        warnings{end + 1} = sprintf( ...
            'num_segments mismatch: meta=%d, segment_def=%d', ...
            truth.meta.num_segments, truth.segment_def.num_segments); %#ok<AGROW>
    end

    nseg = min(numel(truth.meta.segment_num_blocks), truth.segment_def.num_segments);
    for s = 1:nseg
        expected = truth.meta.segment_num_blocks(s);
        got = truth.segment_def.segments(s).num_blocks;
        if expected ~= got
            warnings{end + 1} = sprintf( ...
                'segment %d block-count mismatch: meta=%d, segment_def=%d', ...
                s - 1, expected, got); %#ok<AGROW>
        end
    end

    bad_fmod = [];
    for i = 1:truth.scan_table.num_entries
        fid = truth.scan_table.entries(i).freq_mod_id;
        if fid < 0 || fid > truth.freqmod_def.num_defs
            bad_fmod(end + 1) = i; %#ok<AGROW>
        end
    end
    if ~isempty(bad_fmod)
        warnings{end + 1} = sprintf( ...
            'scan table has invalid freq_mod_id at rows: %s', ...
            mat2str(bad_fmod)); %#ok<AGROW>
    end

    validation = struct();
    validation.ok = isempty(warnings);
    validation.warnings = warnings;
end
