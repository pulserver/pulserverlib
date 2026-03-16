function truth = truth_plot_case(base, varargin)
%TRUTH_PLOT_CASE Convenience entrypoint for report and plots.
%
% truth_plot_case('gre_2d_1sl_1avg')
% truth_plot_case('gre_2d_1sl_1avg', 'show_report', true)

    p = inputParser;
    addParameter(p, 'show_report', true, @islogical);
    addParameter(p, 'plot_tr', true, @islogical);
    addParameter(p, 'plot_segments', true, @islogical);
    addParameter(p, 'plot_freqmod', true, @islogical);
    parse(p, varargin{:});

    truth = testutils.truth_parse_case(base);

    if p.Results.show_report
        testutils.truth_report_case(base);
    end
    if p.Results.plot_tr
        testutils.truth_plot_tr_waveforms(truth);
    end
    if p.Results.plot_segments
        testutils.truth_plot_segments(truth);
    end
    if p.Results.plot_freqmod && truth.freqmod_def.num_defs > 0
        testutils.truth_plot_freqmod_defs(truth);
    end
end
