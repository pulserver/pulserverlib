classdef TruthBuilder < handle
% TRUTHBUILDER  Automate ground-truth export for pulseqlib C-test sequences.
%
%   The user constructs an mr.Sequence object, then hands it to TruthBuilder.
%   The builder derives all bookkeeping (peak RF, canonical TR scale, segment
%   energy, freq-mod definitions, scan table) from the sequence blocks and a
%   small set of user-supplied hints.
%
%   Usage:
%       sys = mr.opts(...);
%       seq = mr.Sequence(sys);
%       % ... build sequence ...
%
%       tb = TruthBuilder(seq, sys);
%       tb.setBlocksPerTR(4);
%       tb.setSegments([4]);               % unique segment sizes
%       tb.setSegmentOrder([1]);           % segment instances per TR
%       tb.setNumAverages(3);
%       tb.anchorPoints.adc = [0.5 0.0];   % optional, per ADC definition
%       tb.setBaseRotation(eye(3));        % optional, default eye(3)
%       tb.export(out_dir, 'gre_2d_1sl_1avg');

    properties
        anchorPoints = struct('adc', [])
    end

    % ----- public properties (read-only after export) -----
    properties (SetAccess = private)
        seq         % mr.Sequence object
        sys         % mr.opts system struct

        % User-supplied hints
        num_blocks_in_tr  = 0
        segment_sizes     = []   % [S1, S2, ...] blocks in each unique segment definition
        segment_reps      = []   % derived from segment_order
        segment_order     = []   % expanded segment instances per TR (1-based segment IDs)
        num_averages      = 1
        base_rot          = eye(3)
        canonical_mode    = 'tr'  % 'tr' or 'pass_expanded'
        multipass_info    = struct('enabled', false, 'pass_starts', [], 'pass_len', 0, 'once_flags', [])

        % Derived quantities (computed by prepare())
        peakRF              = 0
        num_canonical_trs   = 0    % number of unique canonical TR groups
        tr_group_labels     = []   % (num_imaging_trs, 1) 0-based group IDs
        canonical_scales    = {}   % cell array of (nbt, 3) per group
        canonical_seqs      = {}   % cell array of mr.Sequence per group
        TR                  = 0
        tr_times_list       = {}   % cell array of time vectors per group
        tr_waveforms        = {}   % cell array of (N_i, 3) per group
        segment_data        = []   % struct for binary export
        unique_adcs         = []   % struct array: (numSamples, dwell)
        fmod_defs           = {}
        fmod_types          = []
        fmod_durations      = []   % blockDuration of each def (for scan-table matching)
        fmod_kinds          = []   % 0=RF, 1=ADC (same as fmod_types)
        fmod_def_tr_group_ids = [] % [n] TR-group ID associated with each def
        fmod_waveforms_cell = {}   % cell array of waveform cell arrays for dedup
        fmod_build_mode     = 'full_collection'
        scan_table          = []
        rotmat_table        = []
        freq_mod_table      = []
        scan_src_block_idx  = []
        scan_src_tr_start_idx = []
        scan_src_tr_group_id = []
        fmod_plan_truth     = struct()
        label_names         = {}
        label_states_per_block = int32([])
        label_states_per_scan = int32([])
        label_states_per_adc = int32([])
        label_adc_scan_rows = int32([])
        label_adc_value_min = int32([])
        label_adc_value_max = int32([])

        prepared          = false
    end

    methods
        % ---- constructor ----
        function obj = TruthBuilder(seq, sys)
            obj.seq = seq;
            obj.sys = sys;
        end

        % ---- setters ----
        function setBlocksPerTR(obj, n)
            obj.num_blocks_in_tr = n;
            obj.prepared = false;
        end

        function setSegments(obj, sizes)
        % SETSEGMENTS  Define unique segment definitions.
        %   sizes: row vector of block counts per unique segment definition
            assert(all(sizes > 0), ...
                'segment sizes must be positive');
            obj.segment_sizes = sizes(:)';
            obj.segment_reps  = ones(size(obj.segment_sizes));
            obj.segment_order = 1:length(obj.segment_sizes);
            obj.prepared = false;
        end

        function setSegmentOrder(obj, order)
        % SETSEGMENTORDER  Override segment instance order within a TR.
        %   order: row vector of 1-based segment IDs, e.g. [1 2 3 3 3 2]
            assert(~isempty(obj.segment_sizes), ...
                'Must call setSegments before setSegmentOrder');

            ord = order(:)';
            num_seg = length(obj.segment_sizes);
            assert(all(ord >= 1 & ord <= num_seg & floor(ord) == ord), ...
                'segment order must contain valid 1-based segment IDs');

            obj.segment_order = ord;

            % Keep reps consistent with effective order.
            reps = zeros(1, num_seg);
            for s = 1:num_seg
                reps(s) = sum(ord == s);
            end
            obj.segment_reps = reps;
            obj.prepared = false;
        end

        function setNumAverages(obj, n)
            obj.num_averages = n;
            obj.prepared = false;
        end

        function setBaseRotation(obj, R)
            obj.base_rot = R;
            obj.prepared = false;
        end

        function setFreqModBuildMode(obj, mode)
        % SETFREQMODBUILDMODE  Select freq-mod deduplication scope.
        %   mode: 'full_collection' (default) or 'tr_scoped'
            is_string_scalar = false;
            if exist('isstring', 'builtin') || exist('isstring', 'file')
                is_string_scalar = isstring(mode) && isscalar(mode);
            end
            if ~ischar(mode) && ~is_string_scalar
                error('freq-mod build mode must be a string');
            end
            mode = strtrim(char(mode));
            if ~strcmp(mode, 'full_collection') && ~strcmp(mode, 'tr_scoped')
                error('Unsupported freq-mod build mode: %s', mode);
            end
            obj.fmod_build_mode = mode;
            obj.prepared = false;
        end

        % ---- main entry point ----
        function export(obj, out_dir, base_name)
        % EXPORT  Compute all derived quantities and write binary files.
            obj.validate();
            obj.prepared = false;
            obj.prepare();

            if ~exist(out_dir, 'dir')
                mkdir(out_dir);
            end

            % Write .seq file
            % Ensure RequiredExtensions is a char so Pulseq write() strfind
            % does not fail on empty [] in Octave.
            if isempty(obj.seq.getDefinition('RequiredExtensions'))
                obj.seq.setDefinition('RequiredExtensions', '');
            end
            obj.seq.write(fullfile(out_dir, [base_name '.seq']));

            % Meta text
            obj.exportMeta(fullfile(out_dir, [base_name '_meta.txt']));

            % TR waveform binary
            obj.exportTrWaveform(fullfile(out_dir, [base_name '_tr_waveform.bin']));

            % Segment definition binary
            obj.exportSegmentDef(fullfile(out_dir, [base_name '_segment_def.bin']));

            % Freq-mod definition binary
            obj.exportFreqModDefs(fullfile(out_dir, [base_name '_freqmod_def.bin']));

            % Scan table binary
            obj.exportScanTable(fullfile(out_dir, [base_name '_scan_table.bin']));

            % Supplemental label-state truth for metadata validation
            obj.exportLabelState(fullfile(out_dir, [base_name '_label_state.bin']));

            % Supplemental plan-level freq-mod truth (rotation-aware projected library)
            obj.exportFreqModPlan(fullfile(out_dir, [base_name '_freqmod_plan.bin']));
        end
    end

    methods (Access = private)
        % ---- validation ----
        function validate(obj)
            assert(obj.num_blocks_in_tr > 0, ...
                'Must call setBlocksPerTR before export');
            assert(~isempty(obj.segment_sizes), ...
                'Must call setSegments before export');
            assert(~isempty(obj.segment_order), ...
                'Segment order is empty');

            blocks_in_pattern = sum(obj.segment_sizes(obj.segment_order));
            assert(blocks_in_pattern == obj.num_blocks_in_tr, ...
                'Segment topology does not match blocks-per-TR');
        end

        % ---- prepare all derived data ----
        function prepare(obj)
            if obj.prepared, return; end

            obj.discoverUniqueAdcs();
            obj.validateAnchorPoints();
            obj.determineCanonicalMode();
            obj.computePeakRFAndCanonicalScales();
            obj.buildCanonicalTRs();
            obj.buildFreqModDefs();
            obj.buildSegmentData();
            obj.buildScanTableData();
            obj.buildFreqModPlanTruth();

            % Store definitions on the sequence object.
            obj.seq.setDefinition('PeakRF', obj.peakRF);
            obj.seq.setDefinition('TR', obj.TR);
            obj.seq.setDefinition('TotalDuration', sum(obj.seq.blockDurations));

            obj.prepared = true;
        end

        % ---- ADC discovery ----
        function discoverUniqueAdcs(obj)
        % DISCOVERUNIQUEADCS  Scan all blocks, collect unique (numSamples, dwell) pairs.
            adcs = zeros(0, 2);  % each row: [numSamples, dwell]
            for b = 1:length(obj.seq.blockEvents)
                block = obj.seq.getBlock(b);
                if isfield(block, 'adc') && ~isempty(block.adc)
                    key = [block.adc.numSamples, block.adc.dwell];
                    if isempty(adcs) || ~any(adcs(:,1) == key(1) & abs(adcs(:,2) - key(2)) < 1e-15)
                        adcs(end+1, :) = key; %#ok<AGROW>
                    end
                end
            end
            obj.unique_adcs = adcs;
        end

        function validateAnchorPoints(obj)
            if ~isstruct(obj.anchorPoints)
                error('anchorPoints must be a struct');
            end
            if ~isfield(obj.anchorPoints, 'adc') || isempty(obj.anchorPoints.adc)
                return;
            end

            adc = obj.anchorPoints.adc(:)';
            if any(~isfinite(adc)) || any(adc < 0) || any(adc > 1)
                error('anchorPoints.adc values must be in [0, 1]');
            end

            num_adc_defs = size(obj.unique_adcs, 1);
            if length(adc) ~= num_adc_defs
                error('anchorPoints.adc must have one value per unique ADC definition');
            end

            obj.anchorPoints.adc = adc;
        end

        function determineCanonicalMode(obj)
        % DETERMINECANONICALMODE  Choose TR-based or pass-expanded canonical mode.
        %
        %   Two modes:
        %   1) Degenerate / standard TR: canonical_mode = 'tr', multipass_info disabled.
        %      Used when: no once flags, pass_len != nbt, or every pass is purely
        %      imaging (no once!=0 blocks).
        %   2) Non-degenerate (pass_expanded): canonical_mode = 'pass_expanded',
        %      multipass_info enabled. Used when: each pass is exactly nbt blocks AND
        %      at least one block has once==1 OR once==2 (structurally different from
        %      the main imaging TR). Applies to both single-pass (e.g. bSSFP 1sl) and
        %      multi-pass (e.g. bSSFP 3sl) sequences.
        %      Canonical TR duration = seg_cseq.duration (full expanded pass including
        %      all averages).  Canonical TR waveform = full single pass (one per-slice
        %      pass from the block table, without average expansion).
            nbt    = obj.num_blocks_in_tr;
            nblocks = length(obj.seq.blockEvents);
            once_flags = obj.getOnceFlags();

            obj.canonical_mode = 'tr';
            obj.multipass_info = struct('enabled', false, ...
                'pass_starts', [], 'pass_len', 0, 'once_flags', once_flags);

            if isempty(once_flags) || nblocks == 0
                return;
            end

            starts = find(once_flags == 1 & [true, once_flags(1:end-1) ~= 1]);
            if isempty(starts)
                return;
            end

            pass_lens = diff([starts, nblocks + 1]);
            if any(pass_lens ~= pass_lens(1))
                return;
            end

            pass_len = pass_lens(1);

            % Pass must be exactly one TR unit.  Sequences whose "pass" spans many
            % TR-unit widths (e.g. GRE with dummy TRs) are left in standard 'tr' mode.
            if pass_len ~= nbt
                return;
            end

            has_main    = true;
            has_nonmain = false;
            for p = 1:length(starts)
                i0 = starts(p);
                i1 = i0 + pass_len - 1;
                st = once_flags(i0:i1);
                if ~any(st == 0)
                    has_main = false;
                    break;
                end
                if any(st == 1) || any(st == 2)
                    has_nonmain = true;
                end
            end

            if ~has_main || ~has_nonmain
                return;
            end

            % Non-degenerate detected: use pass_expanded for both single and
            % multi-pass sequences.
            obj.canonical_mode = 'pass_expanded';
            obj.multipass_info = struct('enabled', true, ...
                'pass_starts', starts(:)', ...
                'pass_len', pass_len, ...
                'once_flags', once_flags);
        end

        % ---- TR group discovery ----
        function discoverTRGroups(obj)
        % DISCOVERTGROUPS  Fingerprint imaging TRs by gradient shape;
        %   group TRs with identical shot-index patterns.
        %   Mirrors C library's pulseqlib__find_unique_shot_trs().
            if obj.multipass_info.enabled
                np = length(obj.multipass_info.pass_starts);
                obj.num_canonical_trs = 1;
                obj.tr_group_labels = zeros(np, 1);
                return;
            end

            nbt = obj.num_blocks_in_tr;
            num_dummy = obj.findNumDummyBlocks();
            num_imaging_blocks = length(obj.seq.blockEvents) - num_dummy;
            num_trs = num_imaging_blocks / nbt;
            assert(num_trs == floor(num_trs) && num_trs > 0, ...
                'Imaging block count is not a multiple of blocks-per-TR');

            % Shape cache: cell array of fingerprint vectors (one per unique shape).
            shapes = {};

            % Build per-TR fingerprint: one shape ID per (position, axis).
            fp_matrix = zeros(num_trs, nbt * 3);

            for tr_i = 1:num_trs
                base_blk = num_dummy + (tr_i - 1) * nbt;
                for pos = 1:nbt
                    block = obj.seq.getBlock(base_blk + pos);
                    ax_names = {'gx', 'gy', 'gz'};
                    for a = 1:3
                        axn = ax_names{a};
                        if isfield(block, axn) && ~isempty(block.(axn))
                            [sid, shapes] = testutils.TruthBuilder.matchOrAddShape(block.(axn), shapes);
                        else
                            sid = 0;
                        end
                        fp_matrix(tr_i, (pos - 1) * 3 + a) = sid;
                    end
                end
            end

            % Deduplicate fingerprints.
            % unique(...,'rows','stable') with 3 outputs is not available in
            % Octave < 9, so we replicate it manually.
            [unique_rows, ia] = unique(fp_matrix, 'rows', 'stable');
            ic = zeros(num_trs, 1);
            for k = 1:num_trs
                ic(k) = find(all(bsxfun(@eq, unique_rows, fp_matrix(k,:)), 2), 1);
            end

            obj.num_canonical_trs = length(ia);
            obj.tr_group_labels = ic(:) - 1;  % 0-based to match C convention
        end

        % ---- Phase 1: peak RF + per-group canonical scale + segment energy ----
        function computePeakRFAndCanonicalScales(obj)
            obj.discoverTRGroups();

            if obj.multipass_info.enabled
                pass_starts = obj.multipass_info.pass_starts;
                pass_len = obj.multipass_info.pass_len;
                num_total_blocks = length(obj.seq.blockEvents);

                % Global peak RF across all blocks.
                peak_rf = 0;
                for blk = 1:num_total_blocks
                    block = obj.seq.getBlock(blk);
                    if isfield(block, 'rf') && ~isempty(block.rf)
                        pk = max(abs(block.rf.signal));
                        if pk > peak_rf, peak_rf = pk; end
                    end
                end
                obj.peakRF = peak_rf;

                % Choose worst-energy pass as canonical representative.
                best_energy = -inf;
                best_idx = pass_starts(1);
                for p = 1:length(pass_starts)
                    e = 0;
                    base_blk = pass_starts(p);
                    for pos = 1:pass_len
                        block = obj.seq.getBlock(base_blk + pos - 1);
                        ax_names = {'gx', 'gy', 'gz'};
                        for a = 1:3
                            axn = ax_names{a};
                            if isfield(block, axn) && ~isempty(block.(axn))
                                e = e + testutils.TruthBuilder.gradEnergy(block.(axn));
                            end
                        end
                    end
                    if e > best_energy
                        best_energy = e;
                        best_idx = base_blk;
                    end
                end

                obj.canonical_scales = {[]};
                obj.segment_data = struct();
                obj.segment_data.max_seg_energy_idx = best_idx;
                obj.segment_data.repr_energy_indices = best_idx;
                obj.segment_data.first_group_indices = best_idx;
                return;
            end

            nbt = obj.num_blocks_in_tr;
            num_total_blocks = length(obj.seq.blockEvents);
            num_dummy = obj.findNumDummyBlocks();
            num_imaging_blocks = num_total_blocks - num_dummy;
            num_trs = num_imaging_blocks / nbt;
            M = obj.num_canonical_trs;

            % Global peak RF (across all blocks including dummies).
            peak_rf = 0;
            for blk = 1:num_total_blocks
                block = obj.seq.getBlock(blk);
                if isfield(block, 'rf') && ~isempty(block.rf)
                    pk = max(abs(block.rf.signal));
                    if pk > peak_rf, peak_rf = pk; end
                end
            end
            obj.peakRF = peak_rf;

            % Per-group: canonical scale + max-energy representative.
            obj.canonical_scales = cell(M, 1);
            repr_energy_indices = zeros(M, 1);
            first_group_indices = zeros(M, 1);
            overall_best_energy = 0;
            overall_best_idx = num_dummy + 1;
            shapes = {};

            ax_names = {'gx', 'gy', 'gz'};

            for g = 1:M
                ref_shape_ids = zeros(nbt, 3);
                ref_grads = cell(nbt, 3);
                can_scale = zeros(nbt, 3);
                best_energy = 0;
                best_idx = 0;

                for tr_i = 1:num_trs
                    if obj.tr_group_labels(tr_i) ~= (g - 1)
                        continue;
                    end

                    base_blk = num_dummy + (tr_i - 1) * nbt;
                    tr_energy = 0;

                    if first_group_indices(g) == 0
                        first_group_indices(g) = base_blk + 1;
                    end

                    for pos = 1:nbt
                        block = obj.seq.getBlock(base_blk + pos);

                        for a = 1:3
                            axn = ax_names{a};
                            if isfield(block, axn) && ~isempty(block.(axn))
                                grad = block.(axn);
                                amp = testutils.TruthBuilder.gradPeakAmpSigned(grad);
                                [sid, shapes] = testutils.TruthBuilder.matchOrAddShape(grad, shapes);
                                tr_energy = tr_energy + testutils.TruthBuilder.gradEnergy(grad);

                                if ref_shape_ids(pos, a) == 0
                                    ref_shape_ids(pos, a) = sid;
                                    ref_grads{pos, a} = grad;
                                    can_scale(pos, a) = amp;
                                else
                                    same_state = false;
                                    if strcmp(grad.type, 'trap') || strcmp(grad.type, 'trapezoid')
                                        same_state = (sid == ref_shape_ids(pos, a));
                                    else
                                        same_state = testutils.TruthBuilder.gradArbStateMatches(ref_grads{pos, a}, grad);
                                    end
                                    if same_state && abs(amp) > abs(can_scale(pos, a))
                                        can_scale(pos, a) = sign(can_scale(pos, a)) * abs(amp);
                                    end
                                end
                            elseif ref_shape_ids(pos, a) == 0
                                can_scale(pos, a) = 0;
                            end
                        end
                    end

                    if tr_energy > best_energy
                        best_energy = tr_energy;
                        best_idx = base_blk + 1;
                    end
                end

                obj.canonical_scales{g} = can_scale;
                repr_energy_indices(g) = best_idx;

                if best_energy > overall_best_energy
                    overall_best_energy = best_energy;
                    overall_best_idx = best_idx;
                end
            end

            % Store the overall best-energy index for segment def export.
            obj.segment_data = struct();
            obj.segment_data.max_seg_energy_idx = overall_best_idx;
            obj.segment_data.repr_energy_indices = repr_energy_indices;
            obj.segment_data.first_group_indices = first_group_indices;
        end

        % ---- Phase 2: canonical TR waveforms (one per group) ----
        function buildCanonicalTRs(obj)
            if obj.multipass_info.enabled
                obj.canonical_seqs = cell(1, 1);
                obj.tr_waveforms   = cell(1, 1);
                obj.tr_times_list  = cell(1, 1);

                pass_len   = obj.multipass_info.pass_len;
                once_flags = obj.multipass_info.once_flags;
                pass_start = obj.segment_data.first_group_indices(1);

                % ---- Build segment canonical sequence (full pass, all averages) ----
                % Used by buildSegmentData (segment block count, timing, etc.).
                % Includes prep on first avg, imaging on all avgs, cooldown on last avg.
                seg_cseq = mr.Sequence(obj.sys);
                for avg = 1:obj.num_averages
                    for pos = 1:pass_len
                        blk_idx = pass_start + pos - 1;
                        once = once_flags(blk_idx);
                        keep = (once == 0) || ...
                               (once == 1 && avg == 1) || ...
                               (once == 2 && avg == obj.num_averages);
                        if ~keep
                            continue;
                        end

                        block = obj.seq.getBlock(blk_idx);
                        args = {};

                        if isfield(block, 'rf') && ~isempty(block.rf)
                            args{end+1} = block.rf; %#ok<AGROW>
                        end

                        ax_names = {'gx', 'gy', 'gz'};
                        for a = 1:3
                            axn = ax_names{a};
                            if isfield(block, axn) && ~isempty(block.(axn))
                                args{end+1} = block.(axn); %#ok<AGROW>
                            end
                        end

                        if isfield(block, 'adc') && ~isempty(block.adc)
                            args{end+1} = block.adc; %#ok<AGROW>
                        end

                        if isempty(args)
                            seg_cseq.addBlock(mr.makeDelay(block.blockDuration));
                        else
                            seg_cseq.addBlock(args{:});
                        end
                    end
                end
                obj.canonical_seqs{1} = seg_cseq;

                % ---- Build canonical waveform: full SINGLE pass (not avg-expanded) ----
                % Mirrors pulseqlib_get_tr_gradient_waveforms which renders block_table
                % [0 .. pass_len-1] for non-degenerate sequences — i.e. the prep,
                % one imaging loop, and cooldown, all in one slice-execution.
                % This is the same for all avg counts since gradients are invariant
                % to RF phase (which is the only thing that changes between averages).
                wf_cseq = mr.Sequence(obj.sys);
                for pos = 1:pass_len
                    block = obj.seq.getBlock(pass_start + pos - 1);
                    args = {};
                    if isfield(block, 'rf') && ~isempty(block.rf)
                        args{end+1} = block.rf; %#ok<AGROW>
                    end
                    ax_names = {'gx', 'gy', 'gz'};
                    for a = 1:3
                        axn = ax_names{a};
                        if isfield(block, axn) && ~isempty(block.(axn))
                            args{end+1} = block.(axn); %#ok<AGROW>
                        end
                    end
                    if isfield(block, 'adc') && ~isempty(block.adc)
                        args{end+1} = block.adc; %#ok<AGROW>
                    end
                    if isempty(args)
                        wf_cseq.addBlock(mr.makeDelay(block.blockDuration));
                    else
                        wf_cseq.addBlock(args{:});
                    end
                end

                wave_data = wf_cseq.waveforms_and_times(false);
                raster = 0.5 * obj.sys.gradRasterTime;
                times = 0.0 : raster : wf_cseq.duration;
                samples = zeros(length(times), 3);
                for c = 1:3
                    if c <= length(wave_data) && ~isempty(wave_data{c})
                        t_raw = wave_data{c}(1, :);
                        w_raw = wave_data{c}(2, :);
                        [t_u, iu] = unique(t_raw, 'last');
                        w_u = w_raw(iu);
                        samples(:, c) = interp1(t_u, w_u, times, 'linear', 0);
                    end
                end

                obj.tr_times_list{1} = times(:) * 1e6;
                obj.tr_waveforms{1}  = samples;

                % ---- obj.TR = canonical TR duration (written to tr_duration_us) ----
                % Full expanded pass (prep + N_avg*imaging + cooldown) regardless of
                % num_passes. Matches C library: total_scan_dur / num_passes.
                obj.TR = seg_cseq.duration;
                return;
            end

            nbt = obj.num_blocks_in_tr;
            M = obj.num_canonical_trs;
            repr = obj.segment_data.first_group_indices;

            obj.canonical_seqs = cell(M, 1);
            obj.tr_waveforms   = cell(M, 1);
            obj.tr_times_list  = cell(M, 1);

            for g = 1:M
                max_idx = repr(g);
                can_scale = obj.canonical_scales{g};

                cseq = mr.Sequence(obj.sys);
                tr_dur = 0;

                for pos = 1:nbt
                    block = obj.seq.getBlock(max_idx + pos - 1);
                    args = {};

                    % RF: same shape in every TR, use as-is
                    if isfield(block, 'rf') && ~isempty(block.rf)
                        args{end+1} = block.rf; %#ok<AGROW>
                    end

                    % Gradients: scale to canonical amplitude
                    ax_names = {'gx', 'gy', 'gz'};
                    for a = 1:3
                        axn = ax_names{a};
                        if isfield(block, axn) && ~isempty(block.(axn))
                            ref_amp = testutils.TruthBuilder.gradPeakAmpSigned(block.(axn));
                            can_amp = can_scale(pos, a);
                            if ref_amp ~= 0
                                scale = can_amp / ref_amp;
                            else
                                scale = 0;
                            end
                            args{end+1} = mr.scaleGrad(block.(axn), scale); %#ok<AGROW>
                        end
                    end

                    % ADC: include if present
                    if isfield(block, 'adc') && ~isempty(block.adc)
                        args{end+1} = block.adc; %#ok<AGROW>
                    end

                    if isempty(args)
                        % Delay-only blocks have no RF/grad/ADC events; emit explicit delay.
                        cseq.addBlock(mr.makeDelay(block.blockDuration));
                    else
                        cseq.addBlock(args{:});
                    end
                    tr_dur = tr_dur + block.blockDuration;
                end

                obj.canonical_seqs{g} = cseq;

                % Resample to half-gradient raster (matches C library).
                wave_data = cseq.waveforms_and_times(false);
                raster = 0.5 * obj.sys.gradRasterTime;
                times = 0.0 : raster : cseq.duration;
                samples = zeros(length(times), 3);
                for c = 1:3
                    if c <= length(wave_data) && ~isempty(wave_data{c})
                        t_raw = wave_data{c}(1, :);
                        w_raw = wave_data{c}(2, :);
                        % interp1 requires monotonic unique sample times.
                        % Keep the last value on duplicated timestamps.
                        [t_u, iu] = unique(t_raw, 'last');
                        w_u = w_raw(iu);
                        samples(:, c) = interp1(t_u, w_u, times, 'linear', 0);
                    end
                end
                obj.tr_times_list{g} = times(:) * 1e6;
                obj.tr_waveforms{g}  = samples;
            end

            % TR duration (same for all groups — same block structure).
            first_idx = obj.segment_data.repr_energy_indices(1);
            obj.TR = 0;
            for pos = 1:nbt
                block = obj.seq.getBlock(first_idx + pos - 1);
                obj.TR = obj.TR + block.blockDuration;
            end
        end

        % ---- Phase 3: segment definition ----
        function buildSegmentData(obj)
            nbt = obj.num_blocks_in_tr;
            if obj.multipass_info.enabled
                tr_starts = obj.multipass_info.pass_starts;
                num_trs = length(tr_starts);
            else
                num_total_blocks = length(obj.seq.blockEvents);
                num_dummy = obj.findNumDummyBlocks();
                num_imaging_blocks = num_total_blocks - num_dummy;
                num_trs = num_imaging_blocks / nbt;
                tr_starts = zeros(1, num_trs);
                for tr_i = 1:num_trs
                    tr_starts(tr_i) = num_dummy + (tr_i - 1) * nbt + 1;
                end
            end
            num_seg = length(obj.segment_sizes);
            seg_order = obj.segment_order;

            obj.segment_data.num_segments = num_seg;
            obj.segment_data.segments = cell(1, num_seg);

            % In pass-expanded mode, canonical TR already includes ONCE filtering
            % over all averages; reuse it so segment blocks match exported TR.
            if obj.multipass_info.enabled && obj.num_averages > 1 && ...
               num_seg == 1 && all(seg_order == 1)
                cseq = obj.canonical_seqs{1};
                n_canon_blocks = length(cseq.blockEvents);

                obj.segment_sizes(1) = n_canon_blocks;

                seg_blocks = cell(1, n_canon_blocks);
                block_start = 0.0;
                for b = 1:n_canon_blocks
                    block = cseq.getBlock(b);
                    seg_blocks{b} = obj.extractBlockData(block, block_start);
                    block_start = block_start + block.blockDuration;
                end

                [rf_adc_gap, adc_adc_gap] = obj.computeSegmentGaps(seg_blocks);

                seg = struct();
                seg.blocks = seg_blocks;
                seg.rf_adc_gap_us = rf_adc_gap;
                seg.adc_adc_gap_us = adc_adc_gap;
                obj.segment_data.segments{1} = seg;
                return;
            end

            % Compute start offsets of each segment instance in one TR.
            num_inst = length(seg_order);
            inst_start = zeros(1, num_inst);
            cum_offset = 0;
            for k = 1:num_inst
                inst_start(k) = cum_offset;
                cum_offset = cum_offset + obj.segment_sizes(seg_order(k));
            end
            assert(cum_offset == nbt, ...
                'Expanded segment order does not span full TR block count');

            for s = 1:num_seg
                inst_ids = find(seg_order == s);
                assert(~isempty(inst_ids), ...
                    'Every segment definition must appear at least once in segment order');

                best_energy = -inf;
                best_start = -1;

                for tr_i = 1:num_trs
                    tr_base = tr_starts(tr_i);

                    for ii = 1:length(inst_ids)
                        inst = inst_ids(ii);
                        seg_offset = inst_start(inst);
                        e = 0;

                        for b = 1:obj.segment_sizes(s)
                            block = obj.seq.getBlock(tr_base + seg_offset + b - 1);
                            ax_names = {'gx', 'gy', 'gz'};
                            for a = 1:3
                                axn = ax_names{a};
                                if isfield(block, axn) && ~isempty(block.(axn))
                                    e = e + testutils.TruthBuilder.gradEnergy(block.(axn));
                                end
                            end
                        end

                        % Use a relative tolerance to avoid picking a later
                        % instance over the first when energies are equal up
                        % to floating-point noise (e.g. conjugate spiral
                        % interleaves where sum(w^2) differs by ~1 ULP).
                        fprintf('[DEBUG] seg%d TR%d inst%d blk_base=%d energy=%.15g\n', s, tr_i, ii, tr_base + seg_offset, e);
                        if e > best_energy + 1e-9 * abs(e)
                            best_energy = e;
                            best_start = tr_base + seg_offset;
                        end
                    end
                end

                assert(best_start > 0, 'Could not find representative segment instance');

                seg_blocks = cell(1, obj.segment_sizes(s));
                block_start = 0.0;

                for b = 1:obj.segment_sizes(s)
                    block = obj.seq.getBlock(best_start + b - 1);
                    seg_blocks{b} = obj.extractBlockData(block, block_start);
                    block_start = block_start + block.blockDuration;
                end

                % Compute segment-level gaps.
                [rf_adc_gap, adc_adc_gap] = obj.computeSegmentGaps(seg_blocks);

                seg = struct();
                seg.blocks = seg_blocks;
                seg.rf_adc_gap_us = rf_adc_gap;
                seg.adc_adc_gap_us = adc_adc_gap;
                obj.segment_data.segments{s} = seg;
            end
        end

        % ---- Phase 4: frequency modulation definitions ----
        function buildFreqModDefs(obj)
            num_total_blocks = length(obj.seq.blockEvents);
            num_dummy_blocks = obj.findNumDummyBlocks();

            % Collect ALL unique freq-mod block types.
            % Deduplication: compare actual waveforms (not just peak amplitudes).
            % This correctly distinguishes rotated spirals from identical ones.
            defs = {};
            types = [];
            def_durations = [];   % blockDuration of each def (for matching)
            def_kinds     = [];   % 0=RF, 1=ADC
            def_waveforms_cell = {};  % per-def gradient waveform samples for dedup

            for blk = 1:num_total_blocks
                block = obj.seq.getBlock(blk);

                if isfield(block, 'rf') && ~isempty(block.rf)
                    rf_start = block.rf.delay;
                    rf_end   = block.rf.delay + block.rf.t(end);
                    if obj.anyGradNonzeroInWindow(block, rf_start, rf_end)
                        dur = block.blockDuration;
                        gwaves = testutils.TruthBuilder.blockGradWaveformsInWindow(block, rf_start, rf_end);
                        already = false;
                        for k = 1:length(def_durations)
                            if def_kinds(k) == 0 && abs(def_durations(k) - dur) < 1e-9
                                % Compare actual waveforms with tolerance
                                if testutils.TruthBuilder.waveformsAlmostEqual(def_waveforms_cell{k}, gwaves)
                                    already = true;
                                    break;
                                end
                            end
                        end
                        if ~already
                            rf_active_start = block.rf.delay;
                            rf_active_end   = block.rf.delay + block.rf.t(end);
                            rf_isodelay     = block.rf.t(end) - mr.calcRfCenter(block.rf);
                            defs{end+1} = testutils.TruthBuilder.buildFreqModDefinition( ...
                                block, rf_active_start, rf_active_end, rf_isodelay, ...
                                obj.sys.gradRasterTime, obj.sys.rfRasterTime); %#ok<AGROW>
                            types(end+1) = 0;         %#ok<AGROW>
                            def_durations(end+1) = dur; %#ok<AGROW>
                            def_kinds(end+1) = 0;      %#ok<AGROW>
                            def_waveforms_cell{end+1} = gwaves; %#ok<AGROW>
                        end
                    end
                end

                if isfield(block, 'adc') && ~isempty(block.adc) && blk > num_dummy_blocks
                    adc_start = block.adc.delay;
                    adc_end   = block.adc.delay + block.adc.numSamples * block.adc.dwell;
                    if obj.anyGradNonzeroInWindow(block, adc_start, adc_end)
                        dur = block.blockDuration;
                        gwaves = testutils.TruthBuilder.blockGradWaveformsInWindow(block, adc_start, adc_end);
                        already = false;
                        for k = 1:length(def_durations)
                            if def_kinds(k) == 1 && abs(def_durations(k) - dur) < 1e-9
                                % Compare actual waveforms with tolerance
                                if testutils.TruthBuilder.waveformsAlmostEqual(def_waveforms_cell{k}, gwaves)
                                    already = true;
                                    break;
                                end
                            end
                        end
                        if ~already
                            adc_dur = block.adc.numSamples * block.adc.dwell;
                            adc_active_start = block.adc.delay;
                            adc_active_end   = adc_active_start + adc_dur;
                            key = [block.adc.numSamples, block.adc.dwell];
                            match = find(obj.unique_adcs(:,1) == key(1) & ...
                                         abs(obj.unique_adcs(:,2) - key(2)) < 1e-15, 1);
                            adc_ref_time     = obj.getADCAnchorFraction(match - 1) * adc_dur;
                            defs{end+1} = testutils.TruthBuilder.buildFreqModDefinition( ...
                                block, adc_active_start, adc_active_end, adc_ref_time, ...
                                obj.sys.gradRasterTime, obj.sys.adcRasterTime); %#ok<AGROW>
                            types(end+1) = 1;           %#ok<AGROW>
                            def_durations(end+1) = dur;  %#ok<AGROW>
                            def_kinds(end+1) = 1;        %#ok<AGROW>
                            def_waveforms_cell{end+1} = gwaves; %#ok<AGROW>
                        end
                    end
                end
            end

            obj.fmod_defs              = defs;
            obj.fmod_types             = types;
            obj.fmod_durations         = def_durations;
            obj.fmod_kinds             = def_kinds;
            obj.fmod_waveforms_cell    = def_waveforms_cell;
        end

        % ---- Phase 5: scan table ----
        function buildScanTableData(obj)
            num_blocks_per_pass = length(obj.seq.blockEvents);
            num_cols = 15;
            max_entries = obj.num_averages * num_blocks_per_pass;

            % Per-block mapping within one TR/pass:
            %   segment id (0-based) and block index within that segment (0-based).
            % For pass-expanded cases, one TR can contain repeated segment-order patterns.
            pattern_len = 0;
            for si = 1:length(obj.segment_order)
                pattern_len = pattern_len + obj.segment_sizes(obj.segment_order(si));
            end
            assert(pattern_len > 0, 'invalid segment_order/segment_sizes pattern');

            seg_id_pattern = zeros(1, pattern_len);
            seg_blk_pattern = zeros(1, pattern_len);
            pat_cursor = 1;
            for si = 1:length(obj.segment_order)
                seg_id_1b = obj.segment_order(si);
                seg_sz = obj.segment_sizes(seg_id_1b);
                for bi = 1:seg_sz
                    seg_id_pattern(pat_cursor) = seg_id_1b - 1;
                    seg_blk_pattern(pat_cursor) = bi - 1;
                    pat_cursor = pat_cursor + 1;
                end
            end

            seg_id_map = zeros(1, obj.num_blocks_in_tr);
            seg_blk_map = zeros(1, obj.num_blocks_in_tr);
            for k = 1:obj.num_blocks_in_tr
                p = mod(k - 1, pattern_len) + 1;
                seg_id_map(k) = seg_id_pattern(p);
                seg_blk_map(k) = seg_blk_pattern(p);
            end

            st  = zeros(max_entries, num_cols);
            rot = zeros(max_entries, 9);
            fmt = zeros(max_entries, 1);
            src_blk = zeros(max_entries, 1);
            src_tr_start = zeros(max_entries, 1);
            src_tr_group = zeros(max_entries, 1);
            act = 1;

            ppm_to_hz = 1e-6 * obj.sys.gamma * obj.sys.B0;
            num_dummy_blocks = obj.findNumDummyBlocks();

            obj.ensureLabelTimeline();
            block_label_states = obj.label_states_per_block;
            nlabels = size(block_label_states, 2);
            once_col = obj.findLabelColumn('ONCE');
            norot_col = obj.findLabelColumn('NOROT');
            label_scan = zeros(max_entries, nlabels, 'int32');

            once  = 0;
            norot_active = 0;

            % Build ordered traversal list: (block_index, play_prep, play_cool).
            % For multi-pass sequences (e.g. bSSFP where each slice is a pass),
            % the C library expands averages WITHIN each pass (outer=pass,
            % inner=avg).  For standard single-pass sequences the order is
            % outer=avg, inner=block (identical to before).
            trav_b    = zeros(1, max_entries, 'int32');
            trav_prep = false(1, max_entries);
            trav_cool = false(1, max_entries);
            trav_len  = 0;
            if obj.multipass_info.enabled && obj.num_averages > 1
                mp_starts  = obj.multipass_info.pass_starts;
                mp_plen    = obj.multipass_info.pass_len;
                for p = 1:length(mp_starts)
                    ps = mp_starts(p);
                    for avg = 1:obj.num_averages
                        for bi = 0:mp_plen-1
                            trav_len = trav_len + 1;
                            trav_b(trav_len)    = ps + bi;
                            trav_prep(trav_len) = (avg == 1);
                            trav_cool(trav_len) = (avg == obj.num_averages);
                        end
                    end
                end
            else
                for avg = 1:obj.num_averages
                    for b = 1:num_blocks_per_pass
                        trav_len = trav_len + 1;
                        trav_b(trav_len)    = b;
                        trav_prep(trav_len) = (avg == 1);
                        trav_cool(trav_len) = (avg == obj.num_averages);
                    end
                end
            end

            for li = 1:trav_len
                b         = trav_b(li);
                play_prep = trav_prep(li);
                play_cool = trav_cool(li);
                block = obj.seq.getBlock(b);

                    if once_col > 0
                        once = double(block_label_states(b, once_col));
                    else
                        once = 0;
                    end
                    if norot_col > 0
                        norot_active = double(block_label_states(b, norot_col));
                    else
                        norot_active = 0;
                    end

                    % ONCE filter: once==1 → first avg of this pass only;
                    %              once==2 → last avg of this pass only.
                    if once == 0 || (once == 1 && play_prep) || (once == 2 && play_cool)
                        tr_local_idx = mod(b - 1, obj.num_blocks_in_tr) + 1;

                        % Prepended debug columns:
                        %   1) TR start flag (1 at start of each TR/pass)
                        %   2) segment id (0-based)
                        %   3) block index within segment (0-based)
                        st(act, 1) = double(tr_local_idx == 1);
                        st(act, 2) = seg_id_map(tr_local_idx);
                        st(act, 3) = seg_blk_map(tr_local_idx);

                        % RF
                        if isfield(block, 'rf') && ~isempty(block.rf)
                            st(act, 4) = max(abs(block.rf.signal));
                            st(act, 5) = block.rf.phaseOffset + ppm_to_hz * block.rf.phasePPM;
                            st(act, 6) = block.rf.freqOffset  + ppm_to_hz * block.rf.freqPPM;
                            rf_start = block.rf.delay;
                            rf_end   = block.rf.delay + block.rf.t(end);
                            if obj.anyGradNonzeroInWindow(block, rf_start, rf_end)
                                dur = block.blockDuration;
                                gwaves = testutils.TruthBuilder.blockGradWaveformsInWindow(block, rf_start, rf_end);
                                for ki = 1:length(obj.fmod_durations)
                                    if obj.fmod_kinds(ki) == 0 && abs(obj.fmod_durations(ki) - dur) < 1e-9
                                        if testutils.TruthBuilder.waveformsAlmostEqual(obj.fmod_waveforms_cell{ki}, gwaves)
                                            fmt(act) = ki;
                                            break;
                                        end
                                    end
                                end
                            end
                        end

                        % Gradients
                        if isfield(block, 'gx') && ~isempty(block.gx)
                            st(act, 7) = testutils.TruthBuilder.gradPeakAmpSigned(block.gx);
                        end
                        if isfield(block, 'gy') && ~isempty(block.gy)
                            st(act, 8) = testutils.TruthBuilder.gradPeakAmpSigned(block.gy);
                        end
                        if isfield(block, 'gz') && ~isempty(block.gz)
                            st(act, 9) = testutils.TruthBuilder.gradPeakAmpSigned(block.gz);
                        end

                        % ADC
                        if isfield(block, 'adc') && ~isempty(block.adc)
                            st(act, 10) = 1;
                            st(act, 11) = block.adc.phaseOffset + ppm_to_hz * block.adc.phasePPM;
                            st(act, 12) = block.adc.freqOffset  + ppm_to_hz * block.adc.freqPPM;
                            adc_start = block.adc.delay;
                            adc_end   = block.adc.delay + block.adc.numSamples * block.adc.dwell;
                            if obj.anyGradNonzeroInWindow(block, adc_start, adc_end)
                                dur = block.blockDuration;
                                gwaves = testutils.TruthBuilder.blockGradWaveformsInWindow(block, adc_start, adc_end);
                                for ki = 1:length(obj.fmod_durations)
                                    if obj.fmod_kinds(ki) == 1 && abs(obj.fmod_durations(ki) - dur) < 1e-9
                                        if testutils.TruthBuilder.waveformsAlmostEqual(obj.fmod_waveforms_cell{ki}, gwaves)
                                            fmt(act) = ki;
                                            break;
                                        end
                                    end
                                end
                            end
                        end

                        % Triggers
                        if isfield(block, 'trig') && ~isempty(block.trig)
                            for t = 1:length(block.trig)
                                if strcmp(block.trig(t).type, 'output')
                                    st(act, 13) = 1;
                                end
                                if strcmp(block.trig(t).type, 'trigger')
                                    st(act, 14) = 1;
                                end
                            end
                        end

                        % Block duration in µs (column 15)
                        st(act, 15) = round(block.blockDuration * 1e6);

                        % Rotation
                        if isfield(block, 'rotation')
                            rotmat = mr.aux.quat.toRotMat(block.rotation.rotQuaternion);
                        else
                            rotmat = eye(3);
                        end
                        if norot_active == 1
                            act_rotmat = rotmat;
                        else
                            act_rotmat = obj.base_rot * rotmat;
                        end
                        rot(act, :) = reshape(act_rotmat', 1, 9);
                        src_blk(act) = b;
                        src_tr_start(act) = floor((b - 1) / obj.num_blocks_in_tr) * obj.num_blocks_in_tr + 1;
                        src_tr_group(act) = floor((b - 1) / obj.num_blocks_in_tr);
                        if nlabels > 0
                            label_scan(act, :) = block_label_states(b, :);
                        end

                        act = act + 1;
                    end
            end

            n = act - 1;
            obj.scan_table    = st(1:n, :);
            obj.rotmat_table  = rot(1:n, :);
            obj.freq_mod_table = fmt(1:n);
            obj.scan_src_block_idx = src_blk(1:n);
            obj.scan_src_tr_start_idx = src_tr_start(1:n);
            obj.scan_src_tr_group_id = src_tr_group(1:n);

            if nlabels > 0
                obj.label_states_per_scan = label_scan(1:n, :);
                adc_rows = find(obj.scan_table(:, 10) ~= 0);
                obj.label_adc_scan_rows = int32(adc_rows(:) - 1);  % 0-based rows
                if isempty(adc_rows)
                    obj.label_states_per_adc = zeros(0, nlabels, 'int32');
                    obj.label_adc_value_min = zeros(1, nlabels, 'int32');
                    obj.label_adc_value_max = zeros(1, nlabels, 'int32');
                else
                    obj.label_states_per_adc = obj.label_states_per_scan(adc_rows, :);
                    obj.label_adc_value_min = min(obj.label_states_per_adc, [], 1);
                    obj.label_adc_value_max = max(obj.label_states_per_adc, [], 1);
                end
            else
                obj.label_states_per_scan = zeros(n, 0, 'int32');
                adc_rows = find(obj.scan_table(:, 10) ~= 0);
                obj.label_states_per_adc = zeros(length(adc_rows), 0, 'int32');
                obj.label_adc_scan_rows = int32(adc_rows(:) - 1);  % 0-based rows
                obj.label_adc_value_min = zeros(1, 0, 'int32');
                obj.label_adc_value_max = zeros(1, 0, 'int32');
            end
        end

        function buildFreqModPlanTruth(obj)
        % BUILDFREQMODPLANTRUTH  Build supplemental plan-level projected freq-mod truth.
            nscan = length(obj.freq_mod_table);
            oblique = [1, 2, 3];
            oblique = oblique / norm(oblique);
            probes = [1, 0, 0; 0, 1, 0; 0, 0, 1; oblique];

            plan_map = -1 * ones(nscan, 1);
            plans = struct('def_id', {}, 'rot', {}, 'scan_row', {});

            for i = 1:nscan
                def_id = obj.freq_mod_table(i);
                if def_id <= 0
                    continue;
                end
                tr_scope_id = 0;
                rot = reshape(obj.rotmat_table(i, :), [3, 3])';
                if strcmp(obj.fmod_build_mode, 'tr_scoped')
                    tr_scope_id = obj.scan_src_tr_start_idx(i);
                end

                found = 0;
                for p = 1:length(plans)
                    if plans(p).def_id == def_id && plans(p).tr_scope_id == tr_scope_id ...
                            && max(abs(plans(p).rot(:) - rot(:))) < 1e-7
                        plan_map(i) = p - 1;
                        found = 1;
                        break;
                    end
                end

                if ~found
                    plans(end+1).def_id = def_id; %#ok<AGROW>
                    plans(end).rot = rot;
                    plans(end).scan_row = i;
                    plans(end).tr_scope_id = tr_scope_id;
                    plan_map(i) = length(plans) - 1;
                end
            end

            num_plans = length(plans);
            max_samples = 0;
            for p = 1:num_plans
                fm = obj.fmod_defs{plans(p).def_id};
                max_samples = max(max_samples, fm.num_samples);
            end

            plan_defs = zeros(num_plans, 1);
            plan_num_samples = zeros(num_plans, 1);
            plan_rot = zeros(num_plans, 9);
            plan_phase_active = zeros(num_plans, size(probes, 1));
            plan_phase_extra = zeros(num_plans, size(probes, 1));
            plan_phase_total = zeros(num_plans, size(probes, 1));
            plan_waveforms = zeros(num_plans, size(probes, 1), max_samples);

            for p = 1:num_plans
                def_id = plans(p).def_id;
                fm = obj.fmod_defs{def_id};
                plan_defs(p) = def_id;
                plan_num_samples(p) = fm.num_samples;
                plan_rot(p, :) = reshape(plans(p).rot', [1, 9]);

                for q = 1:size(probes, 1)
                    d = probes(q, :).';
                    u = plans(p).rot' * d;
                    w = fm.waveform * u;

                    phase_active = dot(fm.ref_integral, u.');
                    phase_extra = 0.0;

                    plan_phase_active(p, q) = phase_active;
                    plan_phase_extra(p, q) = phase_extra;
                    plan_phase_total(p, q) = phase_active + phase_extra;
                    plan_waveforms(p, q, 1:fm.num_samples) = w;
                end
            end

            obj.fmod_plan_truth = struct();
            obj.fmod_plan_truth.probes = probes;
            obj.fmod_plan_truth.num_plans = num_plans;
            obj.fmod_plan_truth.max_samples = max_samples;
            obj.fmod_plan_truth.plan_defs = plan_defs;
            obj.fmod_plan_truth.plan_num_samples = plan_num_samples;
            obj.fmod_plan_truth.plan_rot = plan_rot;
            obj.fmod_plan_truth.plan_phase_active = plan_phase_active;
            obj.fmod_plan_truth.plan_phase_extra = plan_phase_extra;
            obj.fmod_plan_truth.plan_phase_total = plan_phase_total;
            obj.fmod_plan_truth.plan_waveforms = plan_waveforms;
            obj.fmod_plan_truth.scan_to_plan = plan_map;
            obj.fmod_plan_truth.mode = obj.fmod_build_mode;
        end

        function [ref_abs_s, ok] = getFreqModReferenceTime(obj, block, def_id)
            ok = true;
            ref_abs_s = 0;

            if obj.fmod_types(def_id) == 0
                if isfield(block, 'rf') && ~isempty(block.rf)
                    rf_iso_delay = block.rf.t(end) - mr.calcRfCenter(block.rf);
                    ref_abs_s = block.rf.delay + rf_iso_delay;
                else
                    ok = false;
                end
            else
                if isfield(block, 'adc') && ~isempty(block.adc)
                    key = [block.adc.numSamples, block.adc.dwell];
                    match = find(obj.unique_adcs(:,1) == key(1) & ...
                                 abs(obj.unique_adcs(:,2) - key(2)) < 1e-15, 1);
                    if isempty(match)
                        ok = false;
                        return;
                    end
                    anchor = obj.getADCAnchorFraction(match - 1);
                    adc_dur = block.adc.numSamples * block.adc.dwell;
                    ref_abs_s = block.adc.delay + anchor * adc_dur;
                else
                    ok = false;
                end
            end
        end

        function phase = integrateInactivePhaseFromTrStart(obj, tr_start_idx, block_idx, ref_abs_s, active_axes, u)
            phase = 0;
            if tr_start_idx <= 0 || block_idx <= 0
                return;
            end

            area = zeros(1, 3);
            for bb = tr_start_idx:block_idx
                blk = obj.seq.getBlock(bb);
                if bb < block_idx
                    t0 = 0;
                    t1 = blk.blockDuration;
                else
                    t0 = 0;
                    t1 = max(0, min(ref_abs_s, blk.blockDuration));
                end

                if t1 <= t0
                    continue;
                end

                ax_names = {'gx', 'gy', 'gz'};
                for ch = 1:3
                    if active_axes(ch)
                        continue;
                    end
                    axn = ax_names{ch};
                    if isfield(blk, axn) && ~isempty(blk.(axn))
                        area(ch) = area(ch) + testutils.TruthBuilder.integrateGradientInterval(blk.(axn), t0, t1);
                    end
                end
            end

            phase = 2 * pi * dot(area, u.');
        end

        function exportFreqModPlan(obj, path)
            fid = fopen(path, 'wb');
            if fid < 0, error('Failed to open %s', path); end

            p = obj.fmod_plan_truth;
            if isempty(fieldnames(p))
                fwrite(fid, int32(0), 'int32');
                fwrite(fid, int32(0), 'int32');
                fwrite(fid, int32(0), 'int32');
                fwrite(fid, int32(0), 'int32');
                fclose(fid);
                return;
            end

            nprobe = size(p.probes, 1);
            fwrite(fid, int32(nprobe), 'int32');
            fwrite(fid, single(p.probes(:)), 'float32');

            fwrite(fid, int32(p.num_plans), 'int32');
            fwrite(fid, int32(p.max_samples), 'int32');

            for i = 1:p.num_plans
                fwrite(fid, int32(p.plan_defs(i)), 'int32');
                fwrite(fid, int32(p.plan_num_samples(i)), 'int32');
                fwrite(fid, single(p.plan_rot(i, :)), 'float32');
                fwrite(fid, single(p.plan_phase_active(i, :)), 'float32');
                fwrite(fid, single(p.plan_phase_extra(i, :)), 'float32');
                fwrite(fid, single(p.plan_phase_total(i, :)), 'float32');
                for q = 1:nprobe
                    wf = squeeze(p.plan_waveforms(i, q, :));
                    fwrite(fid, single(wf), 'float32');
                end
            end

            fwrite(fid, int32(length(p.scan_to_plan)), 'int32');
            fwrite(fid, int32(p.scan_to_plan), 'int32');

            fclose(fid);
        end

        function once_flags = getOnceFlags(obj)
        % GETONCEFLAGS  ONCE state per block after applying label updates.
            obj.ensureLabelTimeline();
            once_col = obj.findLabelColumn('ONCE');
            if once_col <= 0
                once_flags = zeros(1, length(obj.seq.blockEvents));
            else
                once_flags = double(obj.label_states_per_block(:, once_col)).';
            end
        end

        function ensureLabelTimeline(obj)
            nblocks = length(obj.seq.blockEvents);
            if size(obj.label_states_per_block, 1) == nblocks
                return;
            end

            names = {};
            states = zeros(nblocks, 0, 'int32');
            curr = zeros(1, 0, 'int32');

            for b = 1:nblocks
                block = obj.seq.getBlock(b);
                if isfield(block, 'label') && ~isempty(block.label)
                    labels_in_block = {block.label.label};
                    [ordered_labels, ~] = unique(labels_in_block, 'stable');

                    % Pulseq label semantics are sticky; apply all SETs then INCs per label.
                    for li = 1:length(ordered_labels)
                        name = ordered_labels{li};
                        idx = find(strcmp(names, name), 1);
                        if isempty(idx)
                            names{end+1} = name; %#ok<AGROW>
                            states(:, end+1) = int32(0); %#ok<AGROW>
                            curr(1, end+1) = int32(0); %#ok<AGROW>
                            idx = numel(names);
                        end

                        set_mask = strcmp({block.label.label}, name) & strcmp({block.label.type}, 'labelset');
                        if any(set_mask)
                            set_vals = [block.label(set_mask).value];
                            curr(idx) = int32(round(set_vals(end)));
                        end

                        inc_mask = strcmp({block.label.label}, name) & strcmp({block.label.type}, 'labelinc');
                        if any(inc_mask)
                            inc_vals = [block.label(inc_mask).value];
                            curr(idx) = curr(idx) + int32(round(sum(inc_vals)));
                        end
                    end
                end

                if ~isempty(curr)
                    states(b, :) = curr;
                end
            end

            obj.label_names = names;
            obj.label_states_per_block = states;
        end

        function idx = findLabelColumn(obj, label_name)
            idx = find(strcmp(obj.label_names, label_name), 1);
            if isempty(idx)
                idx = 0;
            end
        end

        % ---- helper: find number of dummy blocks via ONCE label ----
        function n = findNumDummyBlocks(obj)
            % Walk states until first transition from ONCE=1 region to ONCE=0.
            n = 0;
            once_flags = obj.getOnceFlags();
            if isempty(once_flags)
                return;
            end

            in_prep = false;
            for b = 1:length(once_flags)
                if once_flags(b) == 1
                    in_prep = true;
                elseif once_flags(b) == 0 && in_prep
                    n = b - 1;
                    return;
                end
            end
            n = 0;
        end

        % ---- helper: find fundamental imaging TR period within a pass ----
        function period = findImagingTRPeriod(obj, pass_start, pass_len)
        % FINDIMAGINGTRPERIOD  Shortest repeating period among once==0 blocks.
        %
        %   Mirrors the C library's first_repeating_segment() logic applied to
        %   the imaging region of one pass.  Fingerprint per block:
        %   (blockDuration_us, hasRF, hasADC) — sufficient for typical MR
        %   sequences.  Returns the full imaging block count if no shorter
        %   period is found.
            once_flags = obj.getOnceFlags();

            % Collect once==0 block indices in this pass.
            img_blks = zeros(1, pass_len);
            n = 0;
            for pos = 1:pass_len
                idx = pass_start + pos - 1;
                if once_flags(idx) == 0
                    n = n + 1;
                    img_blks(n) = idx;
                end
            end
            img_blks = img_blks(1:n);

            if n == 0
                period = 1;
                return;
            end

            % Fingerprint: (blockDuration_us, hasRF, hasADC).
            fp = zeros(n, 3);
            for i = 1:n
                blk = obj.seq.getBlock(img_blks(i));
                fp(i, 1) = round(blk.blockDuration * 1e6);
                fp(i, 2) = double(isfield(blk, 'rf')  && ~isempty(blk.rf));
                fp(i, 3) = double(isfield(blk, 'adc') && ~isempty(blk.adc));
            end

            % Find shortest repeating period (compatible with Octave).
            for l = 1:floor(n / 2)
                match = true;
                for i = 1:l
                    if any(fp(i, :) ~= fp(i + l, :))
                        match = false;
                        break;
                    end
                end
                if match
                    period = l;
                    return;
                end
            end
            period = n;
        end

        % ---- helper: extract block data for segment def ----
        function bd = extractBlockData(obj, block, block_start)
            bd = struct();

            % RF
            if isfield(block, 'rf') && ~isempty(block.rf)
                bd.has_rf = true;
                bd.rf_delay = block.rf.delay;
                rf_samples = block.rf.signal;
                bd.rf_amp = max(abs(rf_samples));
                if bd.rf_amp > 0
                    rf_samples = rf_samples / bd.rf_amp;
                end
                if ~any(imag(rf_samples))
                    bd.rf_rho = real(rf_samples);
                    bd.rf_theta = [];
                else
                    bd.rf_rho = abs(rf_samples);
                    bd.rf_theta = angle(rf_samples);
                end
                rf_time = block.rf.t;
                dt = rf_time(2) - rf_time(1);
                if length(unique(diff(rf_time))) == 1 && dt == obj.sys.rfRasterTime
                    bd.rf_time = [];
                else
                    bd.rf_time = rf_time;
                end
            else
                bd.has_rf = false;
                bd.rf_delay = 0;
                bd.rf_rho = [];
                bd.rf_theta = [];
                bd.rf_amp = 0;
                bd.rf_time = [];
            end

            % Gradients
            ax_names = {'gx', 'gy', 'gz'};
            for a = 1:3
                axn = ax_names{a};
                if isfield(block, axn) && ~isempty(block.(axn))
                    grad = block.(axn);
                    bd.([axn '_delay']) = grad.delay;
                    if strcmp(grad.type, 'trap')
                        bd.([axn '_amp']) = grad.amplitude;
                        if grad.flatTime == 0
                            bd.([axn '_wave']) = [0, 1, 0];
                            bd.([axn '_time']) = cumsum([0, grad.riseTime, grad.fallTime]);
                        else
                            bd.([axn '_wave']) = [0, 1, 1, 0];
                            bd.([axn '_time']) = cumsum([0, grad.riseTime, grad.flatTime, grad.fallTime]);
                        end
                    else
                        w = grad.waveform;
                        bd.([axn '_amp']) = testutils.TruthBuilder.gradPeakAmpSigned(grad);
                        if bd.([axn '_amp']) ~= 0
                            bd.([axn '_wave']) = w / bd.([axn '_amp']);
                        else
                            bd.([axn '_wave']) = w;
                        end
                        t = grad.tt;
                        dt = t(2) - t(1);
                        if length(unique(diff(t))) == 1 && dt == obj.sys.gradRasterTime
                            bd.([axn '_time']) = [];
                        else
                            bd.([axn '_time']) = t;
                        end
                    end
                else
                    bd.([axn '_delay']) = 0;
                    bd.([axn '_wave']) = [];
                    bd.([axn '_amp']) = 0;
                    bd.([axn '_time']) = [];
                end
            end

            % ADC
            has_grad = ~isempty(bd.gx_wave) || ~isempty(bd.gy_wave) || ~isempty(bd.gz_wave);
            if isfield(block, 'adc') && ~isempty(block.adc)
                bd.has_adc = 1;
                bd.adc_delay = block.adc.delay;
                % Look up 0-based unique ADC index.
                key = [block.adc.numSamples, block.adc.dwell];
                match = find(obj.unique_adcs(:,1) == key(1) & ...
                             abs(obj.unique_adcs(:,2) - key(2)) < 1e-15, 1);
                bd.adc_id = match - 1;  % 0-based
            else
                bd.has_adc = 0;
                bd.adc_delay = 0;
                bd.adc_id = -1;
            end

            % Rotation: derive from labels or block rotation presence.
            bd.rotate = isfield(block, 'rotation') && ~isempty(block.rotation);

            % Digital output
            bd.has_digital_out = 0;
            bd.digital_out_delay = 0;
            bd.digital_out_duration = 0;
            if isfield(block, 'trig') && ~isempty(block.trig)
                for t = 1:length(block.trig)
                    if strcmp(block.trig(t).type, 'output')
                        bd.has_digital_out = 1;
                        bd.digital_out_delay = block.trig(t).delay;
                        bd.digital_out_duration = block.trig(t).duration;
                    end
                end
            end

            % Frequency modulation
            bd.has_freq_mod = false;
            bd.num_freq_mod_samples = 0;
            bd.freq_mod_def_id = -1;  % 0-based; -1 = none
            if bd.has_rf && has_grad
                rf_ws = block.rf.delay;
                rf_we = block.rf.delay + block.rf.t(end);
                if obj.anyGradNonzeroInWindow(block, rf_ws, rf_we)
                    bd.has_freq_mod = true;
                    bd.num_freq_mod_samples = round(block.blockDuration / obj.sys.rfRasterTime);
                    dur = block.blockDuration;
                    gwaves = testutils.TruthBuilder.blockGradWaveformsInWindow(block, rf_ws, rf_we);
                    for ki = 1:length(obj.fmod_durations)
                        if obj.fmod_kinds(ki) == 0 && abs(obj.fmod_durations(ki) - dur) < 1e-9
                            if testutils.TruthBuilder.waveformsAlmostEqual(obj.fmod_waveforms_cell{ki}, gwaves)
                                bd.freq_mod_def_id = ki - 1;  % 0-based
                                break;
                            end
                        end
                    end
                end
            end
            if bd.has_adc && has_grad
                adc_ws = block.adc.delay;
                adc_we = block.adc.delay + block.adc.numSamples * block.adc.dwell;
                if obj.anyGradNonzeroInWindow(block, adc_ws, adc_we)
                    bd.has_freq_mod = true;
                    bd.num_freq_mod_samples = round(block.blockDuration / obj.sys.adcRasterTime);
                    dur = block.blockDuration;
                    gwaves = testutils.TruthBuilder.blockGradWaveformsInWindow(block, adc_ws, adc_we);
                    for ki = 1:length(obj.fmod_durations)
                        if obj.fmod_kinds(ki) == 1 && abs(obj.fmod_durations(ki) - dur) < 1e-9
                            if testutils.TruthBuilder.waveformsAlmostEqual(obj.fmod_waveforms_cell{ki}, gwaves)
                                bd.freq_mod_def_id = ki - 1;  % 0-based
                                break;
                            end
                        end
                    end
                end
            end

            % Anchor times (relative to segment start, in us)
            if bd.has_rf
                rf_iso = mr.calcRfCenter(block.rf);
                bd.rf_isocenter_us = (block_start + block.rf.delay + rf_iso) * 1e6;
                bd.rf_start_us = (block_start + block.rf.delay) * 1e6;
                bd.rf_end_us   = (block_start + block.rf.delay + block.rf.t(end)) * 1e6;
            else
                bd.rf_isocenter_us = -1;
                bd.rf_start_us = -1;
                bd.rf_end_us   = -1;
            end
            if bd.has_adc
                adc_dur_s = block.adc.numSamples * block.adc.dwell;
                adc_anchor = obj.getADCAnchorFraction(bd.adc_id);
                bd.adc_kzero_us = (block_start + block.adc.delay + adc_anchor * adc_dur_s) * 1e6;
                bd.adc_start_us = (block_start + block.adc.delay) * 1e6;
                bd.adc_end_us   = (block_start + block.adc.delay + adc_dur_s) * 1e6;
            else
                bd.adc_kzero_us = -1;
                bd.adc_start_us = -1;
                bd.adc_end_us   = -1;
            end
        end

        % ---- helper: compute segment-level RF->ADC and ADC->ADC gaps ----
        function [rf_adc_gap, adc_adc_gap] = computeSegmentGaps(~, seg_blocks)
            rf_starts = [];
            rf_ends = [];
            adc_starts = [];
            adc_ends = [];
            for b = 1:length(seg_blocks)
                bd = seg_blocks{b};
                if bd.rf_end_us >= 0
                    rf_starts = [rf_starts, bd.rf_start_us]; %#ok<AGROW>
                    rf_ends = [rf_ends, bd.rf_end_us]; %#ok<AGROW>
                end
                if bd.adc_start_us >= 0
                    adc_starts = [adc_starts, bd.adc_start_us]; %#ok<AGROW>
                    adc_ends   = [adc_ends,   bd.adc_end_us];   %#ok<AGROW>
                end
            end

            rf_adc_gap = -1;
            for r = 1:length(rf_ends)
                candidates = adc_starts(adc_starts >= rf_ends(r));
                if ~isempty(candidates)
                    gap = min(candidates) - rf_ends(r);
                    if rf_adc_gap < 0 || gap < rf_adc_gap
                        rf_adc_gap = gap;
                    end
                end
            end

            adc_adc_gap = -1;
            if length(adc_starts) >= 2
                sorted_starts = sort(adc_starts);
                sorted_ends   = sort(adc_ends);
                for a = 2:length(sorted_starts)
                    % skip this ADC pair if any RF event starts between them
                    interleaved = any(rf_starts > sorted_ends(a-1) & ...
                                      rf_starts < sorted_starts(a));
                    if interleaved
                        continue;
                    end
                    gap = sorted_starts(a) - sorted_ends(a-1);
                    if adc_adc_gap < 0 || gap < adc_adc_gap
                        adc_adc_gap = gap;
                    end
                end
            end
        end

        function anchor = getADCAnchorFraction(obj, adc_id)
            if adc_id < 0
                anchor = 0.5;
                return;
            end

            anchor = 0.5;
            if ~isempty(obj.anchorPoints.adc)
                idx = adc_id + 1;
                if idx > length(obj.anchorPoints.adc)
                    error('Missing anchorPoints.adc entry for ADC definition %d', adc_id);
                end
                anchor = obj.anchorPoints.adc(idx);
            end
        end

        % ---- helper: check if any gradient in block is nonzero in window ----
        function result = anyGradNonzeroInWindow(~, block, wstart, wend)
            result = false;
            ax_names = {'gx', 'gy', 'gz'};
            for a = 1:3
                axn = ax_names{a};
                if isfield(block, axn) && ~isempty(block.(axn))
                    if testutils.TruthBuilder.gradNonzeroInWindow(block.(axn), wstart, wend)
                        result = true;
                        return;
                    end
                end
            end
        end

        % ---- export: meta text ----
        function exportMeta(obj, path)
            fid = fopen(path, 'w');
            if fid < 0, error('Failed to open %s', path); end

            num_adcs = size(obj.unique_adcs, 1);
            fprintf(fid, 'num_unique_adcs %d\n', num_adcs);
            for a = 1:num_adcs
                fprintf(fid, 'adc_%d_samples %d\n', a - 1, obj.unique_adcs(a, 1));
                fprintf(fid, 'adc_%d_dwell_ns %d\n', a - 1, round(obj.unique_adcs(a, 2) * 1e9));
                fprintf(fid, 'adc_%d_anchor %.9g\n', a - 1, obj.getADCAnchorFraction(a - 1));
            end
            fprintf(fid, 'max_b1_subseq %d\n', 0);
            fprintf(fid, 'tr_duration_us %d\n', round(obj.TR * 1e6));
            fprintf(fid, 'num_segments %d\n', length(obj.segment_sizes));
            for s = 1:length(obj.segment_sizes)
                fprintf(fid, 'segment_%d_num_blocks %d\n', s - 1, obj.segment_sizes(s));
            end
            fprintf(fid, 'num_canonical_trs %d\n', obj.num_canonical_trs);
            fprintf(fid, 'canonical_mode %s\n', obj.canonical_mode);
            fprintf(fid, 'fmod_build_mode %s\n', obj.fmod_build_mode);
            if ~isempty(obj.tr_times_list) && ~isempty(obj.tr_times_list{1})
                fprintf(fid, 'canonical_duration_us %d\n', round(obj.tr_times_list{1}(end)));
            end
            if ~isempty(fieldnames(obj.fmod_plan_truth))
                fprintf(fid, 'num_freqmod_plan_probes %d\n', size(obj.fmod_plan_truth.probes, 1));
                fprintf(fid, 'num_freqmod_plan_entries %d\n', obj.fmod_plan_truth.num_plans);
            end
            fprintf(fid, 'num_labels %d\n', numel(obj.label_names));
            fprintf(fid, 'num_label_scan_rows %d\n', size(obj.label_states_per_scan, 1));
            fprintf(fid, 'num_label_adc_rows %d\n', size(obj.label_states_per_adc, 1));
            fprintf(fid, 'segment_order');
            for k = 1:length(obj.segment_order)
                fprintf(fid, ' %d', obj.segment_order(k) - 1);  % 0-based
            end
            fprintf(fid, '\n');

            fclose(fid);
        end

        % ---- export: TR waveform binary ----
        function exportTrWaveform(obj, path)
            fid = fopen(path, 'w');
            if fid < 0, error('Failed to open %s', path); end

            M = obj.num_canonical_trs;
            fwrite(fid, int32(M), 'int32');

            for g = 1:M
                times = obj.tr_times_list{g};
                wf    = obj.tr_waveforms{g};
                N = length(times);
                fwrite(fid, int32(N), 'int32');
                fwrite(fid, single(times), 'float32');
                fwrite(fid, single(wf(:, 1)), 'float32');
                fwrite(fid, single(wf(:, 2)), 'float32');
                fwrite(fid, single(wf(:, 3)), 'float32');
            end

            fclose(fid);
        end

        % ---- export: segment definition binary ----
        function exportSegmentDef(obj, path)
            fid = fopen(path, 'wb');
            if fid < 0, error('Failed to open %s', path); end

            fwrite(fid, obj.segment_data.num_segments, 'int32');

            for s = 1:obj.segment_data.num_segments
                seg = obj.segment_data.segments{s};
                blocks = seg.blocks;
                fwrite(fid, length(blocks), 'int32');

                for b = 1:length(blocks)
                    bd = blocks{b};

                    flags = 0;
                    flags = flags + (bd.has_rf * (2^0));
                    flags = flags + (~isempty(bd.gx_wave) * (2^1));
                    flags = flags + (~isempty(bd.gy_wave) * (2^2));
                    flags = flags + (~isempty(bd.gz_wave) * (2^3));
                    flags = flags + (bd.has_adc * (2^4));
                    flags = flags + (bd.rotate * (2^5));
                    flags = flags + (bd.has_digital_out * (2^6));
                    flags = flags + (bd.has_freq_mod * (2^7));
                    fwrite(fid, uint8(flags), 'uint8');

                    % RF
                    fwrite(fid, single(bd.rf_delay), 'float32');
                    fwrite(fid, single(bd.rf_amp), 'float32');
                    fwrite(fid, single(obj.sys.rfRasterTime * 1e6), 'float32');  % rf_raster_us
                    if bd.has_rf && ~isempty(bd.rf_rho)
                        fwrite(fid, int32(length(bd.rf_rho)), 'int32');
                        fwrite(fid, single(bd.rf_rho), 'float32');
                        % rf_time_s: sample times in seconds (0-based).
                        % Empty bd.rf_time → uniform at rfRasterTime;
                        % non-empty → non-uniform / extended (e.g. hard pulse).
                        if isempty(bd.rf_time)
                            t_s = (0:length(bd.rf_rho)-1).' * obj.sys.rfRasterTime;
                        else
                            t_s = bd.rf_time(:);
                        end
                        fwrite(fid, single(t_s), 'float32');
                    else
                        fwrite(fid, int32(0), 'int32');
                    end

                    % Gradients
                    ax_names = {'gx', 'gy', 'gz'};
                    for a = 1:3
                        axn = ax_names{a};
                        wave = bd.([axn '_wave']);
                        delay = bd.([axn '_delay']);
                        amp = bd.([axn '_amp']);
                        fwrite(fid, single(delay), 'float32');
                        fwrite(fid, single(amp), 'float32');
                        if ~isempty(wave)
                            fwrite(fid, int32(length(wave)), 'int32');
                            fwrite(fid, single(wave), 'float32');
                            % grad_time_s: local knot/sample times in s (0-based).
                            % For traps: cumsum([0, riseTime, flatTime, fallTime]).
                            % For arb:   grad.tt (uniform at gradRasterTime if _time==[]).
                            t_s = bd.([axn '_time']);
                            if isempty(t_s)
                                t_s = (0:length(wave)-1).' * obj.sys.gradRasterTime;
                            end
                            fwrite(fid, single(t_s(:)), 'float32');
                        else
                            fwrite(fid, int32(0), 'int32');
                        end
                    end

                    % ADC
                    fwrite(fid, single(bd.adc_delay), 'float32');
                    fwrite(fid, int32(bd.adc_id), 'int32');  % 0-based ADC def index

                    % Digital output
                    fwrite(fid, single(bd.digital_out_delay), 'float32');
                    fwrite(fid, single(bd.digital_out_duration), 'float32');

                    % Freq-mod
                    fwrite(fid, int32(bd.num_freq_mod_samples), 'int32');
                    fwrite(fid, int32(bd.freq_mod_def_id), 'int32');  % 0-based def index

                    % Anchors
                    fwrite(fid, single(bd.rf_isocenter_us), 'float32');
                    fwrite(fid, single(bd.adc_kzero_us), 'float32');
                end

                % Segment-level gaps
                fwrite(fid, single(seg.rf_adc_gap_us), 'float32');
                fwrite(fid, single(seg.adc_adc_gap_us), 'float32');
            end

            fclose(fid);
        end

        % ---- export: freq-mod definitions binary ----
        function exportFreqModDefs(obj, path)
            fid = fopen(path, 'wb');
            if fid < 0, error('Failed to open %s', path); end

            nd = length(obj.fmod_defs);
            fwrite(fid, int32(nd), 'int32');

            for d = 1:nd
                fm = obj.fmod_defs{d};
                fwrite(fid, int32(obj.fmod_types(d)), 'int32');
                fwrite(fid, int32(fm.num_samples), 'int32');
                fwrite(fid, single(fm.raster_us), 'float32');
                fwrite(fid, single(fm.duration_us), 'float32');
                fwrite(fid, single(fm.ref_time_us), 'float32');
                fwrite(fid, single(fm.ref_integral), 'float32');
                fwrite(fid, single(fm.waveform(:, 1)), 'float32');
                fwrite(fid, single(fm.waveform(:, 2)), 'float32');
                fwrite(fid, single(fm.waveform(:, 3)), 'float32');
            end

            fclose(fid);
        end

        % ---- export: scan table binary ----
        function exportScanTable(obj, path)
            fid = fopen(path, 'wb');
            if fid < 0, error('Failed to open %s', path); end

            n = size(obj.scan_table, 1);
            fwrite(fid, int32(n), 'int32');

            for i = 1:n
                fwrite(fid, int32(obj.scan_table(i, 1)),  'int32');
                fwrite(fid, int32(obj.scan_table(i, 2)),  'int32');
                fwrite(fid, int32(obj.scan_table(i, 3)),  'int32');
                fwrite(fid, single(obj.scan_table(i, 4)), 'float32');
                fwrite(fid, single(obj.scan_table(i, 5)), 'float32');
                fwrite(fid, single(obj.scan_table(i, 6)), 'float32');
                fwrite(fid, single(obj.scan_table(i, 7)), 'float32');
                fwrite(fid, single(obj.scan_table(i, 8)), 'float32');
                fwrite(fid, single(obj.scan_table(i, 9)), 'float32');
                fwrite(fid, int32(obj.scan_table(i, 10)), 'int32');
                fwrite(fid, single(obj.scan_table(i, 11)), 'float32');
                fwrite(fid, single(obj.scan_table(i, 12)), 'float32');
                fwrite(fid, int32(obj.scan_table(i, 13)), 'int32');
                fwrite(fid, int32(obj.scan_table(i, 14)), 'int32');
                fwrite(fid, single(obj.rotmat_table(i, :)), 'float32');
                fwrite(fid, int32(obj.freq_mod_table(i)), 'int32');
                fwrite(fid, int32(obj.scan_table(i, 15)), 'int32');
            end

            fclose(fid);
        end

        function exportLabelState(obj, path)
            fid = fopen(path, 'wb');
            if fid < 0, error('Failed to open %s', path); end

            nlabels = numel(obj.label_names);
            fwrite(fid, int32(nlabels), 'int32');
            for i = 1:nlabels
                name = obj.label_names{i};
                bytes = uint8(name);
                fwrite(fid, int32(numel(bytes)), 'int32');
                fwrite(fid, bytes, 'uint8');
            end

            nscan = size(obj.label_states_per_scan, 1);
            nadc = size(obj.label_states_per_adc, 1);
            fwrite(fid, int32(nscan), 'int32');
            fwrite(fid, int32(nadc), 'int32');

            if nlabels > 0 && nscan > 0
                fwrite(fid, int32(obj.label_states_per_scan.'), 'int32');
            end

            if nadc > 0
                fwrite(fid, int32(obj.label_adc_scan_rows(:)), 'int32');
            end

            if nlabels > 0 && nadc > 0
                fwrite(fid, int32(obj.label_states_per_adc.'), 'int32');
                fwrite(fid, int32(obj.label_adc_value_min(:)), 'int32');
                fwrite(fid, int32(obj.label_adc_value_max(:)), 'int32');
            elseif nlabels > 0
                fwrite(fid, int32(zeros(nlabels, 1)), 'int32');
                fwrite(fid, int32(zeros(nlabels, 1)), 'int32');
            end

            fclose(fid);
        end
    end

    methods (Static)
        function e = gradEnergy(g)
        % GRADENERGY  Gradient energy: integral of amplitude^2 over time.
        %   Uses the trapezoidal rule, matching pulseqlib__trapz_real_* in C.
            if isfield(g, 'amplitude')
                e = (g.amplitude)^2 / 3 * g.riseTime ...
                  + (g.amplitude)^2 * g.flatTime ...
                  + (g.amplitude)^2 / 3 * g.fallTime;
            elseif isfield(g, 'waveform') && ~isempty(g.waveform)
                w2 = g.waveform .^ 2;
                e = sum(0.5 * (w2(1:end-1) + w2(2:end)) .* diff(g.tt));
            else
                e = 0;
            end
        end

        function s = gradPeakAmpSigned(grad)
        % GRADPEAKAMPSIGNED  Return the amplitude stored in the sequence file.
        %   Matches the value that pulseq writes to the [GRADIENTS]/[TRAP]
        %   sections of the .seq file, and that the C library reads back:
        %
        %   Trapezoid: grad.amplitude (already signed; can be negative).
        %
        %   Arbitrary: max(abs(waveform)) * sign(first_non_zero_sample).
        %     This reproduces the formula in pulseq Sequence.m addBlock() /
        %     registerGradEvent() (see lines ~643-646 of Sequence.m):
        %       amplitude = max(abs(event.waveform));
        %       amplitude = amplitude * sign(fnz);   % fnz = first non-zero
        %     The waveform returned by getBlock() is amplitude*normalized_shape,
        %     so dividing back gives the normalized shape; the seq file stores
        %     this signed amplitude together with the normalized shape.
            if strcmp(grad.type, 'trap') || strcmp(grad.type, 'trapezoid')
                s = grad.amplitude;
            else
                w = grad.waveform;
                peak = max(abs(w));
                if peak == 0
                    s = 0;
                else
                    fnz_idx = find(w ~= 0, 1, 'first');
                    s = peak * sign(w(fnz_idx));
                end
            end
        end

        function amp = gradPeakAmp(grad)
        % GRADPEAKAMP  Return the peak magnitude of a gradient.
            if strcmp(grad.type, 'trap') || strcmp(grad.type, 'trapezoid')
                amp = grad.amplitude;
            else
                amp = max(abs(grad.waveform));
            end
        end

        function sig = blockGradSig(block)
        % BLOCKGRADSIG  Return [|peak_gx|, |peak_gy|, |peak_gz|] for a block.
            sig = [0, 0, 0];
            axes = {'gx', 'gy', 'gz'};
            for ch = 1:3
                ax = axes{ch};
                if isfield(block, ax) && ~isempty(block.(ax))
                    sig(ch) = abs(testutils.TruthBuilder.gradPeakAmp(block.(ax)));
                end
            end
        end

        function sig = blockGradSigInWindow(block, wstart, wend)
        % BLOCKGRADSIGWINDOW  Peak |gradient| per axis within [wstart, wend].
            sig = [0, 0, 0];
            axes = {'gx', 'gy', 'gz'};
            for ch = 1:3
                ax = axes{ch};
                if isfield(block, ax) && ~isempty(block.(ax))
                    [t, w] = testutils.TruthBuilder.gradToKnots(block.(ax));
                    pts = [wstart; t(t > wstart & t < wend); wend];
                    vals = interp1(t, w, pts, 'linear', 0);
                    sig(ch) = max(abs(vals));
                end
            end
        end

        function waveforms = blockGradWaveformsInWindow(block, wstart, wend)
        % BLOCKGRADWAVEFORMSINWINDOW  Extract gradient waveforms in active window.
        %   Uses knot interpolation so windows fully inside flat segments are preserved.
        %   Returns cell array {gx_samples, gy_samples, gz_samples}.
            waveforms = {};
            axes = {'gx', 'gy', 'gz'};
            for ch = 1:3
                ax = axes{ch};
                if isfield(block, ax) && ~isempty(block.(ax))
                    grad = block.(ax);
                    [t, w] = testutils.TruthBuilder.gradToKnots(grad);
                    [t, iu] = unique(t(:), 'last');
                    w = w(iu);

                    if wend <= wstart || numel(t) < 2
                        waveforms{ch} = [];
                        continue;
                    end

                    tk = t(t > wstart & t < wend);
                    pts = [wstart; tk; wend];
                    vals = interp1(t, w, pts, 'linear', 0);

                    if any(abs(vals) > 1.0)
                        waveforms{ch} = vals;
                    else
                        waveforms{ch} = [];
                    end
                else
                    waveforms{ch} = [];
                end
            end
        end

        function match = waveformsAlmostEqual(w1, w2, tol)
        % WAVEFORMSALMOSTEQUAL  Compare two cell-array gradient waveforms.
        %   w1, w2: cell arrays {gx, gy, gz} with variable-length samples
        %   tol: relative tolerance (default 1e-6)
            if nargin < 3
                tol = 1e-6;
            end
            match = true;
            for ch = 1:3
                if isempty(w1{ch}) && isempty(w2{ch})
                    continue;
                end
                if isempty(w1{ch}) || isempty(w2{ch})
                    match = false;
                    return;
                end
                % Quick checks: length and peak
                if length(w1{ch}) ~= length(w2{ch})
                    match = false;
                    return;
                end
                if abs(max(abs(w1{ch})) - max(abs(w2{ch}))) > tol * max(abs(w1{ch}))
                    match = false;
                    return;
                end
                % Full allclose comparison
                max_abs = max(abs(w1{ch}), abs(w2{ch}));
                if any(abs(w1{ch} - w2{ch}) > tol * max(max_abs, 1))
                    match = false;
                    return;
                end
            end
        end

        function result = gradNonzeroInWindow(grad, wstart, wend)
        % GRADNONZEROINWINDOW  Check if gradient has nonzero samples in [wstart, wend].
            result = false;
            if isempty(grad), return; end

            if wend <= wstart
                return;
            end

            [t, w] = testutils.TruthBuilder.gradToKnots(grad);
            [t, iu] = unique(t(:), 'last');
            w = w(iu);
            if numel(t) < 2
                return;
            end

            tk = t(t > wstart & t < wend);
            pts = [wstart; tk; wend];
            vals = interp1(t, w, pts, 'linear', 0);
            result = any(abs(vals) > 1.0);
        end

        function [t, w] = gradToKnots(grad)
        % GRADTOKNOTS  Convert Pulseq gradient to time/amplitude knots.
            if strcmp(grad.type, 'trap') || strcmp(grad.type, 'trapezoid')
                d = grad.delay;
                r = grad.riseTime;
                f = grad.flatTime;
                l = grad.fallTime;
                a = grad.amplitude;
                if f > 0
                    t = [d; d+r; d+r+f; d+r+f+l];
                    w = [0; a; a; 0];
                else
                    t = [d; d+r; d+r+l];
                    w = [0; a; 0];
                end
            else
                t = grad.delay + grad.tt(:);
                w = grad.waveform(:);
            end
        end

        function area = integrateGradientInterval(grad, t0, t1)
        % INTEGRATEGRADIENTINTERVAL  Integrate gradient waveform over [t0, t1] in seconds.
            area = 0;
            if t1 <= t0
                return;
            end

            [t, w] = testutils.TruthBuilder.gradToKnots(grad);
            [t, iu] = unique(t(:), 'last');
            w = w(iu);

            if isempty(t) || numel(t) < 2
                return;
            end

            tk = t(t > t0 & t < t1);
            tq = [t0; tk; t1];
            wq = interp1(t, w, tq, 'linear', 0);
            area = trapz(tq, wq);
        end

        function fmod = buildFreqModDefinition(block, active_start_s, active_end_s, ref_time_s, grad_raster_s, target_raster_s)
        % BUILDFREQMODDEFINITION  Build freq-mod base definition (static).
            active_dur_s   = active_end_s - active_start_s;
            grad_raster_us = grad_raster_s * 1e6;
            active_dur_us  = active_dur_s * 1e6;
            ref_time_us    = ref_time_s * 1e6;

            num_samples = floor(active_dur_us / grad_raster_us) + 1;
            if num_samples < 2, num_samples = 2; end

            uniform_t = (0:num_samples-1)' * grad_raster_s;
            axes = {'gx', 'gy', 'gz'};
            waveform = zeros(num_samples, 3);
            ref_integral = zeros(1, 3);

            for ch = 1:3
                ax = axes{ch};
                if isfield(block, ax) && ~isempty(block.(ax))
                    [raw_t, raw_w] = testutils.TruthBuilder.gradToKnots(block.(ax));
                    raw_t = raw_t - active_start_s;
                    if raw_t(1) > 0
                        raw_t = [0; raw_t(:)]; %#ok<AGROW>
                        raw_w = [0; raw_w(:)]; %#ok<AGROW>
                    end
                    if raw_t(end) < active_dur_s
                        raw_t = [raw_t(:); active_dur_s]; %#ok<AGROW>
                        raw_w = [raw_w(:); 0];            %#ok<AGROW>
                    end
                    waveform(:, ch) = interp1(raw_t(:), raw_w(:), uniform_t, 'linear', 0);
                end

                ref_sample = floor(ref_time_us / grad_raster_us);
                ref_sample = max(0, min(ref_sample, num_samples - 1));
                if ref_sample > 0
                    ref_integral(ch) = 2 * pi * 1e-6 * ...
                        trapz(waveform(1:ref_sample+1, ch)) * grad_raster_us;
                else
                    ref_integral(ch) = 0;
                end
            end

            fmod.num_samples   = num_samples;
            fmod.raster_us     = grad_raster_us;
            fmod.duration_us   = active_dur_us;
            fmod.ref_time_us   = ref_time_us;
            fmod.ref_integral  = ref_integral;
            fmod.waveform      = waveform;

            if target_raster_s > 0 && target_raster_s < grad_raster_s - 1e-9
                target_raster_us = target_raster_s * 1e6;
                fine_num = floor(active_dur_us / target_raster_us) + 1;
                if fine_num < 2, fine_num = 2; end

                fine_waveform = zeros(fine_num, 3);
                for j = 1:fine_num
                    orig_idx = min(floor((j-1) * target_raster_us / grad_raster_us) + 1, num_samples);
                    fine_waveform(j, :) = waveform(orig_idx, :);
                end

                for ch = 1:3
                    ref_sample = floor(ref_time_us / target_raster_us);
                    ref_sample = max(0, min(ref_sample, fine_num - 1));
                    if ref_sample > 0
                        ref_integral(ch) = 2 * pi * 1e-6 * ...
                            trapz(fine_waveform(1:ref_sample+1, ch)) * target_raster_us;
                    else
                        ref_integral(ch) = 0;
                    end
                end

                fmod.num_samples  = fine_num;
                fmod.raster_us    = target_raster_us;
                fmod.ref_integral = ref_integral;
                fmod.waveform     = fine_waveform;
            end
        end

        function [sid, shapes] = matchOrAddShape(grad, shapes)
        % MATCHORADDSHAPE  Return shape ID for a gradient, adding to cache if new.
        %   For trapezoids: fingerprint = (riseTime, flatTime, fallTime)
        %   For arbitrary:  fingerprint = normalized waveform (waveform / peak)
            if strcmp(grad.type, 'trap') || strcmp(grad.type, 'trapezoid')
                fp = [grad.riseTime, grad.flatTime, grad.fallTime];
            else
                w = grad.waveform(:)';
                pk = max(abs(w));
                if pk > 0
                    fp = w / pk;
                else
                    fp = w;
                end
            end

            for i = 1:length(shapes)
                if length(fp) == length(shapes{i}) && ...
                        max(abs(fp - shapes{i})) < 1e-8
                    sid = i;
                    return;
                end
            end

            sid = length(shapes) + 1;
            shapes{sid} = fp;
        end

        function tf = gradArbStateMatches(ref_grad, cand_grad)
        % GRADARBSTATEMATCHES  True when arbitrary gradients share the same raw state.
            tol = 1e-9;
            if isempty(ref_grad) || isempty(cand_grad)
                tf = false;
                return;
            end
            if ~(strcmp(ref_grad.type, cand_grad.type) && ...
                    ~(strcmp(ref_grad.type, 'trap') || strcmp(ref_grad.type, 'trapezoid')))
                tf = false;
                return;
            end
            if length(ref_grad.waveform) ~= length(cand_grad.waveform) || ...
                    length(ref_grad.tt) ~= length(cand_grad.tt)
                tf = false;
                return;
            end
            tf = max(abs(ref_grad.waveform(:) - cand_grad.waveform(:))) < 1e-8 && ...
                 max(abs(ref_grad.tt(:) - cand_grad.tt(:))) < tol && ...
                 abs(ref_grad.delay - cand_grad.delay) < tol;
        end
    end
end
