function [gxCells, gyCells, adc, k] = makeTestSpiral(sys, interleaves, nSamples, fov)
%MAKETESTSPIRAL  Generate rotated spiral readout gradients and ADC event.
%
%   [gxCells, gyCells, adc, k] = makeTestSpiral(sys, interleaves, nSamples, fov)
%
%   Returns cell arrays of per-interleave gradient events (gx, gy), a
%   single shared ADC event, and the k-space trajectory (2 x nSamples).

    deltak = 1 / fov;
    nTurns = nSamples / (2 * pi);

    t = linspace(0, 1, nSamples);
    r   = deltak * nTurns * t;
    phi = 2 * pi * nTurns * t;

    k = [r .* cos(phi); r .* sin(phi)];

    dt = sys.gradRasterTime;
    [g, ~] = mr.traj2grad(k, 'RasterTime', dt);

    adc = mr.makeAdc(nSamples, 'Dwell', dt);

    gxCells = cell(interleaves, 1);
    gyCells = cell(interleaves, 1);

    for i = 1:interleaves
        angle = 2 * pi * (i - 1) / interleaves;
        R = [cos(angle), -sin(angle); sin(angle), cos(angle)];
        grot = R * g;
        gx = grot(1, :); gx = [gx(1), gx, gx(end)];
        gy = grot(2, :); gy = [gy(1), gy, gy(end)];
        gxCells{i} = mr.makeArbitraryGrad('x', gx, 'system', sys, 'first', gx(1), 'last', gx(end));
        gyCells{i} = mr.makeArbitraryGrad('y', gy, 'system', sys, 'first', gy(1), 'last', gy(end));
    end
end
