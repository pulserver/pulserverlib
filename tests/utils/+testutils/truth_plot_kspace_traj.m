function truth_plot_kspace_traj(truth)
%TRUTH_PLOT_KSPACE_TRAJ Plot non-Cartesian k-space trajectories.
%
%   truth_plot_kspace_traj(truth)
%
%   Computes the k-space trajectory for each unique (ADC-def, gradient-
%   amplitude) combination found in the segment definitions.  Straight-line
%   trajectories (Cartesian readout lines, radial spokes) are detected
%   automatically and excluded.  If every trajectory is a straight line, no
%   figure is created and an informational message is printed.
%
%   K-space is expressed in cycles/m (Hzs/m, standard Pulseq convention).
%   Sample axis is centred at 0 (= k-space centre, i.e. adc_kzero_us).
%
%   Requires MATLAB R2017b or later (vecnorm) and R2018b or later (sgtitle).

    meta = truth.meta;
    sdef = truth.segment_def;

    if meta.num_unique_adcs == 0
        return;
    end

    % ---- Collect unique (adc_def_id 0-based, grad_amp[3]) blocks ----------
    % Two blocks are considered identical when they have the same ADC
    % definition ID and identical (bitwise, after double conversion)
    % gradient-amplitude triple.  This matches the library's keying.
    key_list = zeros(0, 4);   % rows: [adc_def_id0, amp_x, amp_y, amp_z]
    blk_list = {};

    for s = 1:sdef.num_segments
        seg = sdef.segments(s);
        for b = 1:seg.num_blocks
            blk = seg.blocks(b);
            if ~blk.has_adc
                continue;
            end
            key = [double(blk.adc_def_id), double(blk.grad_amp(:)')];
            % Exact-match deduplication (same binary representation)
            is_dup = false;
            for ki = 1:size(key_list, 1)
                if key_list(ki,1) == key(1) && all(key_list(ki,2:4) == key(2:4))
                    is_dup = true;
                    break;
                end
            end
            if ~is_dup
                key_list(end+1, :) = key; %#ok<AGROW>
                blk_list{end+1}    = blk; %#ok<AGROW>
            end
        end
    end

    if isempty(key_list)
        return;
    end

    % ---- Compute k-space trajectory per unique block ------------------------------------
    n_keys      = size(key_list, 1);
    kspace      = cell(n_keys, 3);   % kspace{i, ax}: n_samp x 1 (cycles/m)
    n_samp_arr  = zeros(n_keys, 1);  % sample count per key
    center_samp_arr = zeros(n_keys, 1); % 0-based centre sample per key
    is_trivial  = true(n_keys, 1);   % straight-line flag

    K_TOL_REL = 1e-3;   % relative deviation threshold for straight-line test

    for i = 1:n_keys
        blk     = blk_list{i};
        adc_id0 = double(blk.adc_def_id);
        adc_id1 = adc_id0 + 1;
        if adc_id1 < 1 || adc_id1 > meta.num_unique_adcs
            continue;
        end
        dwell_s  = double(meta.adc_dwell_ns(adc_id1)) * 1e-9;
        n_samp   = double(meta.adc_samples(adc_id1));
        if n_samp < 2 || dwell_s <= 0
            continue;
        end
        dwell_us     = dwell_s * 1e6;
        center_samp0 = max(0, min(n_samp-1, ...
            round(double(blk.adc_kzero_us) / dwell_us)));   % 0-based

        n_samp_arr(i)      = n_samp;
        center_samp_arr(i) = center_samp0;

        % Absolute time (within the block, seconds) of each ADC sample.
        % Sample 0-based index j is at:
        %   t_j = adc_kzero_us*1e-6 + (j - center_samp0) * dwell_s
        j_vec = (0:n_samp-1)';
        t_s   = double(blk.adc_kzero_us)*1e-6 + (j_vec - center_samp0) * dwell_s;

        k_traj = zeros(n_samp, 3);
        for ax = 1:3
            gt = double(blk.grad_time_s{ax}(:));
            if numel(gt) < 2
                kspace{i, ax} = zeros(n_samp, 1);
                continue;
            end
            gw     = double(blk.grad_wave{ax}(:)) .* double(blk.grad_amp(ax));
            g_at_t = interp1(gt, gw, t_s, 'linear', 0);
            k_ax   = cumsum(g_at_t) * dwell_s;           % cycles/m
            k_ax   = k_ax - k_ax(center_samp0 + 1);      % zero at k-centre
            kspace{i, ax}  = k_ax;
            k_traj(:, ax)  = k_ax;
        end

        % Straight-line detection in 3-D k-space.
        % A trajectory is a straight line when every sample lies on the
        % chord from the first to the last k-space point within a relative
        % tolerance K_TOL_REL.
        k_first = k_traj(1, :);
        k_last  = k_traj(end, :);
        t_lin   = j_vec / max(n_samp-1, 1);
        k_line  = k_first + t_lin .* (k_last - k_first);
        deviation = max(vecnorm(k_traj - k_line, 2, 2));
        k_scale   = max(vecnorm(k_traj, 2, 2)) + eps;
        if deviation > K_TOL_REL * k_scale
            is_trivial(i) = false;
        end
    end

    % ---- If everything is straight-line, skip --------------------------------------------------
    if all(is_trivial)
        fprintf('truth_plot_kspace_traj: all %d trajectory(ies) are straight lines (Cartesian/radial) - no figure created.\n', n_keys);
        return;
    end

    % ---- Plot non-trivial trajectories ----------------------------------------------------------------
    non_trivial = find(~is_trivial);
    n_plot      = numel(non_trivial);

    figure('Name', sprintf('K-space trajectories: %s', truth.base_name), ...
           'NumberTitle', 'off');

    ax_names = {'kx', 'ky', 'kz'};

    for row = 1:3
        subplot(3, 1, row);
        hold on;
        for jj = 1:n_plot
            i        = non_trivial(jj);
            blk      = blk_list{i};
            adc_id1  = double(blk.adc_def_id) + 1;
            n_samp   = n_samp_arr(i);
            center_s = center_samp_arr(i);
            samp_ax  = (0:n_samp-1) - center_s;   % centre at 0

            k_ax = kspace{i, row};
            if isempty(k_ax)
                k_ax = zeros(1, n_samp);
            end
            plot(samp_ax, k_ax(:)', ...
                'DisplayName', sprintf('traj %d  (adc=%d)', jj, adc_id1-1));
        end
        hold off;
        xlabel('Sample (0 = k-centre)');
        ylabel('cycles/m');
        title(ax_names{row});
        if n_plot > 1
            legend('Location', 'best');
        end
        grid on;
    end

    sgtitle(sprintf('K-space trajectories - %s', truth.base_name), ...
            'Interpreter', 'none');
end
