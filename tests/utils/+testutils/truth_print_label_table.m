function truth_print_label_table(truth)
%TRUTH_PRINT_LABEL_TABLE Print ordered ADC label sequence to the terminal.
%
%   truth_print_label_table(truth)
%
%   For every ADC readout row the following is printed:
%     ADC# | ave=AVE | rot=rot_id | traj=traj_val | n_samp=NS | center=CI
%
%   ave      — value of the "AVE" or "REP" label at this ADC row.
%   rot_id   — 1-based index of the rotation matrix at the corresponding
%              scan-table row (first-appearance ordering).
%   traj_val — value of the "SEG" or "SHT" label (= shot / trajectory
%              index); 0 if the label is absent.
%   n_samp   — number of ADC samples for this readout (from meta).
%   center   — 1-based sample index corresponding to k-space centre
%              (adc_kzero_us / dwell, rounded to nearest integer + 1).

    lbs  = truth.label_state;
    meta = truth.meta;
    st   = truth.scan_table;
    sdef = truth.segment_def;

    if lbs.adc_rows == 0
        fprintf('\n=== Label table: %s — no ADC rows ===\n', truth.base_name);
        return;
    end

    % ── Find label columns ───────────────────────────────────────────
    col_ave  = find_label_col(lbs, 'AVE', 'REP');
    col_traj = find_label_col(lbs, 'SEG', 'SHT');

    % ── Build unique rotation-matrix index (first-appearance) ────────
    n_scan = st.num_entries;
    rot_ids = zeros(n_scan, 1);
    rot_list = zeros(0, 9);
    for i = 1:n_scan
        rm = double(st.entries(i).rotmat);
        found = false;
        for r = 1:size(rot_list, 1)
            if max(abs(rot_list(r,:) - rm)) < 1e-5
                rot_ids(i) = r;
                found = true;
                break;
            end
        end
        if ~found
            rot_list(end+1, :) = rm; %#ok<AGROW>
            rot_ids(i) = size(rot_list, 1);
        end
    end

    % ── Build map: adc_def_id (0-based) → (n_samples, center_us) ────
    % Scan all segment blocks that have an ADC event.
    adc_n_samples  = zeros(meta.num_unique_adcs, 1);
    adc_center_us  = zeros(meta.num_unique_adcs, 1);
    adc_seen       = false(meta.num_unique_adcs, 1);
    for s = 1:sdef.num_segments
        seg = sdef.segments(s);
        for b = 1:seg.num_blocks
            blk = seg.blocks(b);
            if ~blk.has_adc
                continue;
            end
            id1 = blk.adc_def_id + 1;   % 0-based → 1-based index
            if id1 < 1 || id1 > meta.num_unique_adcs
                continue;
            end
            if ~adc_seen(id1)
                adc_n_samples(id1) = meta.adc_samples(id1);
                adc_center_us(id1) = blk.adc_kzero_us;
                adc_seen(id1) = true;
            end
        end
    end

    % ── Header ──────────────────────────────────────────────────────
    fprintf('\n=== Label table: %s (%d ADC rows) ===\n', ...
        truth.base_name, lbs.adc_rows);
    fprintf('  %6s  %5s  %5s  %6s  %7s  %7s\n', ...
        'ADC#', 'ave', 'rot', 'traj', 'n_samp', 'center');
    sep = repmat('-', 1, 52);
    fprintf('%s\n', sep);

    % ── Row loop ────────────────────────────────────────────────────
    for a = 1:lbs.adc_rows
        scan_row = lbs.adc_scan_row_idx(a) + 1;   % 0-based → 1-based

        % Rotation id from scan table
        if scan_row >= 1 && scan_row <= n_scan
            rot_id = rot_ids(scan_row);
        else
            rot_id = 0;
        end

        % Label values
        if col_ave > 0
            ave_val = lbs.adc_states(a, col_ave);
        else
            ave_val = 0;
        end
        if col_traj > 0
            traj_val = lbs.adc_states(a, col_traj);
        else
            traj_val = 0;
        end

        % ADC def id from scan table entry (0-based encoded as adc_flag > 0)
        % adc_flag stores the 1-based ADC definition ID in the scan table; 0 = no ADC.
        adc_def_id1 = 1;   % default to first ADC def
        if scan_row >= 1 && scan_row <= n_scan
            flag = double(st.entries(scan_row).adc_flag);
            if flag > 0 && flag <= meta.num_unique_adcs
                adc_def_id1 = flag;
            end
        end
        n_samp = adc_n_samples(adc_def_id1);
        dwell_us = meta.adc_dwell_ns(adc_def_id1) * 1e-3;
        if dwell_us > 0
            center_idx = round(adc_center_us(adc_def_id1) / dwell_us) + 1;
        else
            center_idx = 0;
        end

        fprintf('  %6d  %5d  %5d  %6d  %7d  %7d\n', ...
            a, ave_val, rot_id, traj_val, n_samp, center_idx);
    end
    fprintf('%s\n', sep);
end

% ── Helpers ────────────────────────────────────────────────────────────
function col = find_label_col(lbs, varargin)
    % Find first matching label column (case-insensitive) from candidate names.
    col = -1;
    for k = 1:lbs.num_labels
        for j = 1:numel(varargin)
            if strcmpi(lbs.label_names{k}, varargin{j})
                col = k;
                return;
            end
        end
    end
end
