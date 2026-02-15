function volt_amper_char(
    I::Float64,
    G::Float64,
    d::Float64
)
    return 49022 * I ^ (-0.5) * G ^ (0.75) * d ^ (-0.5)
end

ch_diameter = 16e-3
println(ch_diameter)