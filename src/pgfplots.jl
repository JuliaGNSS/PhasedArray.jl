
function plot(pattern::Pattern)
    @pgf PolarAxis({ymin = 0, ymax = pattern.max_el * 180 / π, "colormap/viridis", "axis on top", "colorbar style={title=Amplification (dB)}", xticklabel = "{\$\\pgfmathparse{Mod(90-\\tick,360)}\\pgfmathprintnumber{\\pgfmathresult}\$}", yticklabel = "{\$\\pgfmathparse{90-\\tick}\\pgfmathprintnumber{\\pgfmathresult}\$}",}, #colorbar
        Plot3({surf, shader="interp"},
            Coordinates(pattern.azs .* 180 ./ π, pattern.els .* 180 ./ π, pattern.values')
        )
    )
end

function plot(pattern::Pattern3D)
    @pgf Axis({"colormap/viridis"}, #colorbar
        Plot3({surf, "z buffer=sort", "point meta=sqrt(x^2+y^2+z^2)"},
            Coordinates(vec(pattern.X), vec(pattern.Y), vec(pattern.Z))
        )
    )
end
