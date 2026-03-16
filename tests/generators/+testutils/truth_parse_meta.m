function meta = truth_parse_meta(path)
%TRUTH_PARSE_META Parse TruthBuilder _meta.txt file.

    fid = fopen(path, 'r');
    if fid < 0
        error('truth:io', 'Failed to open meta file: %s', path);
    end
    cleanup = onCleanup(@() fclose(fid));

    lines = textscan(fid, '%s %f', 'Delimiter', {' ', '\t'}, 'MultipleDelimsAsOne', true);
    if isempty(lines) || numel(lines) < 2
        error('truth:format', 'Malformed meta file: %s', path);
    end

    keys = lines{1};
    vals = lines{2};

    meta = struct();
    meta.num_unique_adcs = 0;
    meta.adc_samples = [];
    meta.adc_dwell_ns = [];
    meta.max_b1_subseq = 0;
    meta.tr_duration_us = 0;
    meta.num_segments = 0;
    meta.segment_num_blocks = [];
    meta.num_canonical_trs = 0;

    for i = 1:numel(keys)
        key = keys{i};
        val = vals(i);

        if strcmp(key, 'num_unique_adcs')
            meta.num_unique_adcs = round(val);
        elseif strcmp(key, 'max_b1_subseq')
            meta.max_b1_subseq = round(val);
        elseif strcmp(key, 'tr_duration_us')
            meta.tr_duration_us = round(val);
        elseif strcmp(key, 'num_segments')
            meta.num_segments = round(val);
        elseif strcmp(key, 'num_canonical_trs')
            meta.num_canonical_trs = round(val);
        else
            tok = regexp(key, '^adc_(\d+)_(samples|dwell_ns)$', 'tokens', 'once');
            if ~isempty(tok)
                idx = str2double(tok{1}) + 1;
                suffix = tok{2};
                if strcmp(suffix, 'samples')
                    meta.adc_samples(idx) = round(val);
                else
                    meta.adc_dwell_ns(idx) = round(val);
                end
                continue;
            end

            tok = regexp(key, '^segment_(\d+)_num_blocks$', 'tokens', 'once');
            if ~isempty(tok)
                idx = str2double(tok{1}) + 1;
                meta.segment_num_blocks(idx) = round(val);
            end
        end
    end

    if meta.num_unique_adcs ~= numel(meta.adc_samples)
        error('truth:consistency', ...
            'Meta mismatch in %s: num_unique_adcs=%d but adc_samples has %d entries.', ...
            path, meta.num_unique_adcs, numel(meta.adc_samples));
    end

    if meta.num_unique_adcs ~= numel(meta.adc_dwell_ns)
        error('truth:consistency', ...
            'Meta mismatch in %s: num_unique_adcs=%d but adc_dwell_ns has %d entries.', ...
            path, meta.num_unique_adcs, numel(meta.adc_dwell_ns));
    end

    if meta.num_segments ~= numel(meta.segment_num_blocks)
        error('truth:consistency', ...
            'Meta mismatch in %s: num_segments=%d but segment_num_blocks has %d entries.', ...
            path, meta.num_segments, numel(meta.segment_num_blocks));
    end
end
