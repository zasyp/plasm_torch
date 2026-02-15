function volt_amper_char(
    I::Float64,
    G::Float64,
    d::Float64
)
    return 49022 * I ^ (-0.5) * G ^ (0.75) * d ^ (-0.5)
end
P = 3e5
h = 6.877e7
ρ0 = 1.25
η1 = 0.75
P_plasm = P * η1
G = P_plasm/h
ch_diameter = 16e-3

# Подбор тока для достижения заданной мощности
I_range = range(0, 7000, 1000)
P_target = 3e5
best_I   = 0.0
best_U   = 0.0
min_diff = 1e5

for I in I_range
    global  min_diff, best_I, best_U
    if I == 0 continue end
    local U = volt_amper_char(I, G, ch_diameter)
    power   = U * I
    diff    = abs(power - P_target)
    
    if diff < min_diff
        min_diff = diff
        best_I   = I
        best_U   = U
    end
end

println("I = ", round(best_I, digits=1), " A")
println("U = ", round(best_U, digits=1), " V")
println("P = ", round(best_U * best_I), " W")


d = diameter_fvach(best_U, best_I, G)
w0 = 4*G/(ρ0 * π * d ^ 2)
E_avg = 281.8 * best_I * (best_I^2 / w0)^(-0.59) * (G / d)^(-0.53)
ΔU_A = 17
ΔU_K = 25

# Расчеты с L_СУД
L_SUD = (best_U - ΔU_A - ΔU_K) / E_avg
L_zsh = 2.5d
L_chan = L_zsh - 0.5L_zsh
ΔL = d
L_anode = L_zsh + ΔL
j_anode = best_I / (π * d * L_anode)
j_cathode = 1e8
d_cathode = sqrt((4best_I) / (π * j_cathode))
println("d (расчетный) = ", round(d, digits=4), " m")
println("E_avg = ", round(E_avg, digits=2), " V/m")
println("L_SUD = ", round(L_SUD, digits=3), " m")
println("L_anode = ", round(L_anode, digits=3), " m")
println("j_anode = ", round(j_anode, digits=2), " A/m^2")
println("d_cathode = ", round(d_cathode, digits=4), " m")

# Расчеты с L_ФИК
L_FIK = (best_U - ΔU_A - ΔU_K) / E_avg
L_UST = L_FIK
d_anode = 2d
L_0 = (d_anode - d) / (2 * tan(0.226892803))
L_anode = 2L_0 + ΔL
L_chan = L_UST - L_0
L_overall = L_chan + L_anode

# Газодинамические расчеты
a = sqrt(1.4 * 8.31 * 5000/0.028)
w = 4*G/(0.07*π*d^2)
d_critical = sqrt(4G/(0.07*π*a))
a_0 = sqrt(1.4 * 8.31 * 293/0.028)
F_hole = G/(ρ0*a_0)
d_hole = sqrt(4*F_hole/(3*π))
d_hole = 0.002

# Ресурс электродов
g_anode  = 1e-10
ρ_anode = 8960
δ_anode = 0.005
τ_anode = (ρ_anode * δ_anode) / (3600 * g_anode * 1.5 * j_anode)

g_cathode  = 1e-12
ρ_cathode = 19300
δ_cathode = 0.005
τ_cathode = (ρ_cathode * δ_cathode) / (3600 * g_cathode * 1.5 * j_cathode)
println("Длина анода, м:", L_anode)
println("Длина канала, м:", L_chan)
println("Длина всей системы, м:", L_overall)
println("Скорость истечения, м/с:", w)
println("Критический диаметр, м:", d_critical)
println("Диаметр отверстия, м:", d_hole)

println("Ресурс анода = ", round(τ_anode, digits=2), " часов")
println("Ресурс катода = ", round(τ_cathode, digits=2), " часов")
