function volt_amper_char(
    I::Float64,
    G::Float64,
    d::Float64
)
    return 49022 * I ^ (-0.5) * G ^ (0.75) * d ^ (-0.5)
end
p = 1e5
P = 3e5
h = 6.877e7
ρ0 = 1.25
η1 = 0.8943
P_plasm = P * η1
G = P_plasm/h
ch_diameter = 16e-3
R = 8.31/0.028
T = 5000
ρ = p / (R * T)

# Подбор тока для достижения заданной мощности
I_range = range(0, 7000, 1000)
P_target = 3e5 * η1
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

function diameter_fvach(
    U::Float64,
    I::Float64,
    G::Float64
)
    return ((49022 * G ^ (0.75))/(U * I ^ (0.5)))^2
end

d = diameter_fvach(best_U, best_I, G)
println("Диаметр из ВАХ, мм: ", d)
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
a = sqrt(1.4 * R * T)
w = 4*G/(0.07*π*d^2)
d_critical = sqrt(4G/(0.07*π*a))
a_0 = sqrt(1.4 * R * 293)
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


# Тепловые потоки в элементы конструкции
U_ae = 6
U_ak = 4
Q_at = U_ae * best_I
Q_kt = U_ak * best_I

function Nusselt(Re, Pr, L_SUD, d)
    Nu = 0.28 * (Re ^ (0.5)) * (Pr ^ (0.33)) * ((L_SUD / d) ^ (-0.5))
    return Nu
end

function Stanton(h, h0, Re, L_anode, d)
    St = exp(0.0156 * h / h₀ - 5.3) - 
        real(exp(0.0402 * h / h₀ - 16)) * 
        (0.341 * (d / L_A - 0.1)^0.241 + 0.63)
    return St
end

function convective_heat_Stanton(St, h, hw, ρ, w)
    return St * (h - hw) * ρ * w
end

function convective_heat_Nusselt(α, Tsm, Tw, Nu, λ, d)
    α = Nu * λ / d
    return α * (Tsm - Tw)
end
μ = 1.3e-4
k = 1.0
cp = 2500
Pr = cp * μ / k
Re = ρ * w * d / μ
Nu = Nusselt(Re, Pr, L_SUD, ch_diameter)
println("Pr = ", Pr)
println("Re = ", Re)
println("Nu = ", Nu)
convective_heat_Nusselt_calc = convective_heat_Nusselt(Nu, T, 293, Nu, k, d)

Q_conv_anode = convective_heat_Nusselt_calc * π * d_anode * L_anode
Q_conv_channel = convective_heat_Nusselt_calc * π * ch_diameter * L_chan
println("Тепловой поток на аноде, Вт:", Q_conv_anode)
println("Тепловой поток в канале, Вт:", Q_conv_channel)

Q_total_anode = Q_at + Q_conv_anode
Q_total_channel = Q_kt + Q_conv_channel
println("Общий тепловой поток на аноде, Вт:", Q_total_anode)
println("Общий тепловой поток в канале/катоде, Вт:", Q_total_channel)

Q_total = Q_total_anode + Q_total_channel
println("Общий тепловой поток, Вт:", Q_total)
η_thermal = (P - Q_total) / P
println("Тепловой КПД, %:", round(η_thermal * 100, digits=2))


# Теплонапряженность элементов конструкции
q_at = best_I * U_ae / (π * d_anode * L_anode)
q_at_max = 1.75q_at
q_a_max = q_at_max + convective_heat_Nusselt_calc

println("Максимальное значение суммарного удельного теплового потока в анод, Вт/м^2: ", q_a_max)
f = (d_anode + 2δ_anode) / d_anode 
q_a_cooling = q_a_max / f

cp_water = 4200
# Расчет охлаждения
Δ = 0.001
ΔT_cooling = 30
G_a_cooling = Q_total_anode / (cp_water * ΔT_cooling)
w_water = G_a_cooling / (1000 * π * (d_anode + 2δ_anode + Δ) * Δ)
println("Расход воды на охлаждение анода, кг/с: ", G_a_cooling)

Re_water = 1000 * w_water * (2Δ) / μ
