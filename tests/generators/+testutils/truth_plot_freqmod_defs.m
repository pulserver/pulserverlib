function fig = truth_plot_freqmod_defs(base_or_truth, varargin)
%TRUTH_PLOT_FREQMOD_DEFS Plot frequency modulation definitions.

    p = inputParser;
    addParameter(p, 'def_idx', [], @(x) isempty(x) || isnumeric(x));
    parse(p, varargin{:});

    truth = coerce_truth(base_or_truth);

    if truth.freqmod_def.num_defs == 0
        error('truth:empty', 'No frequency modulation definitions found for %s', truth.base_name);
    end

    if isempty(p.Results.def_idx)
        ids = 1:truth.freqmod_def.num_defs;
    else
        ids = p.Results.def_idx(:).';
    end

    fig = figure('Name', sprintf('Freqmod defs: %s', truth.base_name), 'Color', 'w');
    tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    labels = {'Gx', 'Gy', 'Gz'};
    hold_state = gobjects(1, 3);
    for ax = 1:3
        hold_state(ax) = nexttile;
        hold(hold_state(ax), 'on');
        ylabel(hold_state(ax), labels{ax});
        grid(hold_state(ax), 'on');
    end
    xlabel(hold_state(3), 'Time (us)');

    lgd = {};
    for i = 1:numel(ids)
        d = ids(i);
        if d < 1 || d > truth.freqmod_def.num_defs
            error('truth:input', 'def_idx %d out of range', d);
        end
        def = truth.freqmod_def.defs(d);
        t = (0:double(def.num_samples)-1).' * double(def.raster_us);
        kind = 'RF';
        if def.type == 1
            kind = 'ADC';
        end
        lgd{end + 1} = sprintf('Def %d (%s)', d, kind); %#ok<AGROW>

        for ax = 1:3
            plot(hold_state(ax), t, double(def.waveform(:, ax)), 'LineWidth', 1.2);
        end
    end

    legend(hold_state(1), lgd, 'Location', 'eastoutside');
    title(tl, sprintf('Frequency Modulation Definitions (%s)', truth.base_name), 'Interpreter', 'none');
end

function truth = coerce_truth(base_or_truth)
    if isstruct(base_or_truth)
        truth = base_or_truth;
    else
        truth = truth_parse_case(base_or_truth);
    end
end
