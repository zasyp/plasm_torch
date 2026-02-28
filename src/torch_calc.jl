using Plots

function volt_amper_char(I::Float64, G::Float64, d::Float64)
    return 49022 * I ^ (-0.5) * G ^ (0.75) * d ^ (-0.5)
end

function volt_amper_char_alt1(I::Float64, j_k::Float64)
    r_k = sqrt(I / j_k) * 1000  # в мм (La в мм)
    La = 130.0
    ρ = 0.233
    
    a = 3.2 * r_k
    b = 2.2 * r_k
    m = 1.0 / (5 * r_k)

    term1 = m * (1 - a^2) / (a*b) * (1 - a^2 * b * exp(m*La)) * log(a)
    term2 = (a^2 * m * La) / a^2 * (log(a) * exp(m*La) + a^2 / a^2)
    V = I * ρ * (term1 + term2)
    return V
end
function diameter_fvach(U::Float64, I::Float64, G::Float64)
    return ((49022 * G ^ (0.75)) / (U * I ^ (0.5)))^2
end

function Nusselt(Re, Pr, L_SUD, d)
    Nu = 0.28 * (Re ^ (0.5)) * (Pr ^ (0.33)) * ((L_SUD / d) ^ (-0.5))
    return Nu
end

function convective_heat_Nusselt(Nu, Tsm, Tw, λ, d)
    α = Nu * λ / d
    return α * (Tsm - Tw)
end

j_cathode = 1e8

p = 1e5
P = 1e5
h = 6.877e7
ρ0 = 1.25
η1 = 0.94
P_plasm = P * η1
G = P_plasm / h
G = 0.0035
println("Расход газа, кг/с: ", G)
ch_diameter = 10e-3
R = 8.31 / 0.028
T = 5000
ρ = p / (R * T)
P_target = P_plasm

I_range = range(50, 180, 1000)
best_I   = 0.0
best_U   = 0.0
min_diff = 1e5

for I in I_range
    if I == 0 continue end
    local U = volt_amper_char(I, G, ch_diameter)
    local power = U * I
    local diff  = abs(power - P_target)
    global min_diff, best_I, best_U
    if diff < min_diff
        min_diff = diff
        best_I   = I
        best_U   = U
    end
end

println("I = ", round(best_I, digits=1), " A")
println("U = ", round(best_U, digits=1), " V")
println("P = ", round(best_U * best_I), " W")

d = ch_diameter
println("Диаметр из ВАХ, мм: ", round(d * 1000, digits=2))

w0 = 4 * G / (ρ0 * π * d ^ 2)
println("Скорость истечения, м/с: ", round(w0, digits=2))
E_avg = 281.8 * best_I * (best_I^2 / w0)^(-0.59) * (G / d)^(-0.53)
ΔU_A = 17
ΔU_K = 25
L_SUD = (best_U - ΔU_A - ΔU_K) / E_avg
L_zsh = 2.5 * d
L_chan = L_zsh - 0.5 * L_zsh
ΔL = d
L_anode = L_zsh + ΔL
j_anode = best_I / (π * d * L_anode)
d_cathode = sqrt((4 * best_I) / (π * j_cathode))

println("d (расчетный) = ", round(d, digits=4), " m")
println("E_avg = ", round(E_avg, digits=2), " V/m")
println("L_SUD = ", round(L_SUD, digits=3), " m")
println("L_anode = ", round(L_anode, digits=3), " m")
println("j_anode = ", round(j_anode, digits=2), " A/m^2")
println("d_cathode = ", round(d_cathode, digits=4), " m")

L_FIK = (best_U - ΔU_A - ΔU_K) / E_avg
L_UST = L_FIK
d_anode = 2 * d
L_0 = (d_anode - d) / (2 * tan(0.226892803))
println("L0, мм: ", L_0)
L_anode = 2 * L_0 + ΔL
L_chan = L_UST - L_0
L_overall = L_chan + L_anode
j_anode = best_I / (π * d_anode * L_anode)

μ = 1.3e-4
a = sqrt(1.4 * R * T / μ)
println("Скорость звука, м/с: ", a)
w = 4 * G / (0.07 * π * d^2)
println("Скорость истечения для критического диаметра, м/с: ", round(w, digits=2))
d_critical = sqrt(4 * G / (0.07 * π * a))
a_0 = sqrt(1.4 * R * 293)
println("Скорость звука при 293K, м/с: ", round(a_0, digits=2))
F_hole = G / (ρ0 * a_0)
println("Площадь отверстий, м2: ", F_hole)
d_hole = sqrt(4 * F_hole / (3 * π))
d_hole = 0.002

g_anode = 1e-9
ρ_anode = 8960
δ_anode = 0.007
τ_anode = (ρ_anode * δ_anode) / (3600 * g_anode * 1.5 * j_anode)
g_cathode = 1e-12
ρ_cathode = 19300
δ_cathode = 0.005
τ_cathode = (ρ_cathode * δ_cathode) / (3600 * g_cathode * 1.5 * j_cathode)

println("Длина анода, м: ", round(L_anode, digits=3))
println("Длина канала, м: ", round(L_chan, digits=3))
println("Длина всей системы, м: ", round(L_overall, digits=3))
println("Скорость истечения, м/с: ", round(w, digits=2))
println("Критический диаметр, м: ", round(d_critical, digits=4))
println("Диаметр отверстия, м: ", d_hole)
println("Ресурс анода = ", round(τ_anode, digits=2), " часов")
println("Ресурс катода = ", round(τ_cathode, digits=2), " часов")

U_ae = 6
U_ak = 6
Q_at = U_ae * best_I
Q_kt = U_ak * best_I

k = 1.0
cp = 2500
Pr = cp * μ / k
Re = ρ * w * d / μ
Nu = Nusselt(Re, Pr, L_SUD, ch_diameter)
println("Pr = ", round(Pr, digits=2))
println("Re = ", round(Re, digits=0))
println("Nu = ", round(Nu, digits=2))

convective_heat_Nusselt_calc = convective_heat_Nusselt(Nu, T, 293, 0.3, d)
println("Конвективный тепловой поток по Нуссельту, Вт/м^2: ", round(convective_heat_Nusselt_calc, digits=0))
Q_conv_anode = convective_heat_Nusselt_calc * π * d_anode * L_anode
Q_conv_channel = convective_heat_Nusselt_calc * π * ch_diameter * L_chan
println("Тепловой поток на аноде, Вт: ", round(Q_conv_anode, digits=0))
println("Тепловой поток в канале, Вт: ", round(Q_conv_channel, digits=0))
Q_total_anode = Q_at + Q_conv_anode
Q_total_channel = Q_kt + Q_conv_channel
println("Общий тепловой поток на аноде, Вт: ", round(Q_total_anode, digits=0))
println("Общий тепловой поток в канале/катоде, Вт: ", round(Q_total_channel, digits=0))
Q_total = Q_total_anode + Q_total_channel
println("Общий тепловой поток, Вт: ", round(Q_total, digits=0))
η_thermal = (P - Q_total) / P
println("Тепловой КПД, %: ", round(η_thermal * 100, digits=2))

q_at = best_I * U_ae / (π * d_anode * L_anode)
println("Удельный тепловой поток в анод, Вт/м^2: ", round(q_at, digits=0))
q_at_max = 1.75 * q_at
q_a_max = q_at_max + convective_heat_Nusselt_calc
println("Максимальное значение суммарного удельного теплового потока в анод, Вт/м^2: ", round(q_a_max, digits=0))
f = (d_anode + 2 * δ_anode) / d_anode
println("f:", f)
F_anode = π * d_anode * L_anode
F_cooling_anode = f * F_anode
println("Охлаждаемая площадь анода: ", round(F_cooling_anode, digits=4))
q_a_cooling = q_a_max / f
println("Тепловой поток в систему охлаждения: ", round(q_a_cooling, digits=0))

# Охлаждение
cp_water = 4180.0
λ_water = 0.62
μ_water = 8.9e-4
ρ_water = 990.0
λ_cu = 380.0
T_water_in = 293.0
ΔT_cooling = 30.0
T_water_out = T_water_in + ΔT_cooling
p_water = 6e5


p_MPa = p_water / 1e6  
# T_boil_C = 179.47 + 84.45 * log10(p_MPa)  
# T_boil = T_boil_C + 273.15  
# T_crit = T_boil + 20

Δ = 0.003
d_g = 2Δ
Q_remove = q_a_cooling * F_cooling_anode
G_water = (Q_remove / (cp_water * ΔT_cooling))
println("Расход воды на анод, кг/с: ", round(G_water, digits=4))
w_water = G_water / (1000 * π * (d_anode + 2δ_anode + Δ) * Δ)
println("Скорость воды, м/с: ",w_water)
Re_water = 1000 * w_water * d_g / μ_water
println("Re воды = ", round(Re_water))
Pr_water = μ_water * cp_water / λ_water
Nu_water = 0.023 * Re_water^0.8 * Pr_water^0.4
α_water = Nu_water * λ_water / d_g
println("α воды = ", round(α_water), " Вт/м²·К")
T_wv = 323 - (q_a_cooling * d_anode / 2λ_water) * log((d_anode + 2δ_anode) / (d_anode))
println("T_wv, К: ", round(T_wv, digits=1))
T_sr = (T_water_in + T_water_out) / 2
println("T_sr, К: ", round(T_sr, digits=1))
T_boil = 6.5997 * (p_water)^(0.2391) +273.15
ΔT_wall_water = q_a_cooling / α_water
println("ΔT_wall_water для q_a_cooling, К: ", round(ΔT_wall_water, digits=1))
T_crit = T_boil + 20
println("Температура кипения, К: ", T_boil)
T_ned = T_boil - T_sr
println("Недогрев воды, К: ", T_ned)
q_w = q_a_cooling * (d_anode/ (d+2δ_anode))
println("Тепловой поток, отводимый водой, Вт: ", round(q_w, digits=0))
q_cr = 4q_w
w_water = w_water * (q_cr / (q_a_cooling*(1+0.0078*T_ned))) ^ 2
println("Скорость воды для q_cr, м/с: ", round(w_water, digits=2))
δ_water = G_water / (π * (d_anode + 2δ_anode) * w_water * 1000)
println("Толщина слоя воды для q_cr, м: ", round(δ_water, digits=8))
x = Δ / δ_water
println("Отношение Δ/δ_water: ", round(x, digits=7))
G_new = G_water * x
println("Требуемый расход воды для q_cr, кг/с: ", round(G_new, digits=4))
w_new_water = G_new / (1000 * π * (d_anode + 2δ_anode + Δ) * Δ)
println("Скорость воды для q_cr, м/с: ", round(w_new_water, digits=2))
Re_new_water = 1000 * w_new_water * d_g / μ_water
println("Re для q_cr = ", round(Re_new_water))
Pr_new_water = μ_water * cp_water / λ_water
println("Pr для q_cr = ", round(Pr_new_water, digits=2))
Nu_new_water = 0.023 * Re_new_water^0.8 * Pr_new_water^0.4
println("Nu для q_cr = ", round(Nu_new_water, digits=2))
α_new_water = Nu_new_water * λ_water / d_g
println("α для q_cr = ", round(α_new_water), " Вт/м²·К")
ΔT_wall_water = q_a_cooling / α_new_water
println("ΔT_wall_water для q_cr, К: ", round(ΔT_wall_water, digits=1))
r1 = d_anode / 2
r2 = r1 + δ_anode
ΔT_cu = (q_cr*d_anode / 2λ_cu) * log((d_anode + 0.5δ_anode) / 0.5d_anode)
println("ΔT в стенке для q_cr, К: ", round(ΔT_cu, digits=1))
T_wall = ΔT_wall_water + ΔT_cu
println("T стенки для q_cr, К: ", round(T_wall, digits=1))
println("T допустимая, К: ", round(T_crit, digits=1))
println(335.6/1085)


# # Переменные для результатов
# Δ_solution = 0.0
# w_water_solution = 0.0
# Re_water_solution = 0.0
# α_water_solution = 0.0
# ΔT_wall_water_solution = 0.0
# ΔT_cu_solution = 0.0
# T_wall_solution = 0.0
# solution_found = false
# r1 = D_inner / 2
# r2 = r1 + δ
# ΔT_cu = (q_a_cooling * r1 / λ_cu) * log(r2 / r1)
# T_wall = T_water_out + ΔT_wall_water + ΔT_cu

# println(T_wall)
# for Δ in 0.001:0.0001:0.03          # ← изменён диапазон
#     D_outer = D_inner + 2 * δ + 2 * Δ
#     F_flow = π / 4 * (D_outer^2 - (D_inner + 2 * δ)^2)
#     w_water = G_water / (ρ_water * F_flow)
#     d_h = 2 * Δ
#     Re_water = ρ_water * w_water * d_h / μ_water
#     if Re_water < Re_min || w_water > w_max
#         continue
#     end
#     Pr_water = μ_water * cp_water / λ_water
#     Nu_water = 0.023 * Re_water^0.8 * Pr_water^0.4
#     α_water = Nu_water * λ_water / d_h
#     ΔT_wall_water = q_a_cooling / α_water
#     r1 = D_inner / 2
#     r2 = r1 + δ
#     ΔT_cu = (q_a_max * r1 / λ_cu) * log(r2 / r1)
#     T_wall = T_water_out + ΔT_wall_water + ΔT_cu
#     if T_wall < T_crit
#         println("----- РЕШЕНИЕ НАЙДЕНО -----")
#         println("Зазор Δ = ", round(Δ * 1000, digits=2), " мм")
#         println("Скорость воды = ", round(w_water, digits=2), " м/с")
#         println("Re воды = ", round(Re_water))
#         println("α воды = ", round(α_water), " Вт/м²·К")
#         println("T стенки = ", round(T_wall, digits=1), " К")
#         println("T допустимая = ", round(T_crit, digits=1), " К")
#         println("ΔT_wall_water", ΔT_wall_water)
        
#         global Δ_solution = Δ
#         global w_water_solution = w_water
#         global Re_water_solution = Re_water
#         global α_water_solution = α_water
#         global ΔT_wall_water_solution = ΔT_wall_water
#         global ΔT_cu_solution = ΔT_cu
#         global T_wall_solution = T_wall
#         global solution_found = true
#         break
#     end
# end

# if !solution_found
#     println("Решение не найдено.")
# end

# Построение ВАХ + отметка рабочей точки
U_values = Float64[]
I_values = collect(range(100, 7000, length=200))
for I in I_values
    push!(U_values, volt_amper_char(I, G, ch_diameter))
end

p = plot(I_values, U_values, label="ВАХ", xlabel="I, A", ylabel="U, V", linewidth=2)
scatter!(p, [best_I], [best_U], label="Рабочая точка", markersize=6, color=:red)
savefig(p, "vah.png")
println("ВАХ с рабочей точкой сохранена в vah.png")

# Сбор всех ключевых параметров в словарь
params = Dict(
    "I, A"                      => round(best_I, digits=1),
    "U, V"                      => round(best_U, digits=1),
    "P, W"                      => round(best_U * best_I),
    "d из ВАХ, мм"              => round(d * 1000, digits=2),
    "d расчетный, м"            => round(d, digits=4),
    "E_avg, В/м"                => round(E_avg, digits=2),
    "L_SUD, м"                  => round(L_SUD, digits=3),
    "L_zsh, м"                  => round(L_zsh, digits=3),
    "L_chan (СУД), м"           => round(L_chan, digits=3),  # для СУД
    "L_anode (уступ), м"        => round(L_anode, digits=3),
    "L_chan (уступ), м"         => round(L_chan, digits=3),
    "L_overall, м"              => round(L_overall, digits=3),
    "d_anode, м"                => round(d_anode, digits=4),
    "j_anode, А/м²"             => round(j_anode, digits=0),
    "d_cathode, м"              => round(d_cathode, digits=4),
    "w, м/с"                    => round(w, digits=2),
    "d_critical, м"             => round(d_critical, digits=4),
    "d_hole, м"                 => d_hole,
    "Ресурс анода, ч"           => round(τ_anode, digits=2),
    "Ресурс катода, ч"          => round(τ_cathode, digits=2),
    "Pr"                        => round(Pr, digits=2),
    "Re"                        => round(Re, digits=0),
    "Nu"                        => round(Nu, digits=2),
    "Q_conv_anode, Вт"          => round(Q_conv_anode, digits=0),
    "Q_conv_channel, Вт"        => round(Q_conv_channel, digits=0),
    "Q_at, Вт"                  => round(Q_at, digits=0),
    "Q_kt, Вт"                  => round(Q_kt, digits=0),
    "Q_total_anode, Вт"         => round(Q_total_anode, digits=0),
    "Q_total_channel, Вт"       => round(Q_total_channel, digits=0),
    "Q_total, Вт"               => round(Q_total, digits=0),
    "Тепловой КПД, %"           => round(η_thermal * 100, digits=2),
    "q_at, Вт/м²"               => round(q_at, digits=0),
    "q_at_max, Вт/м²"           => round(q_at_max, digits=0),
    "q_a_max, Вт/м²"            => round(q_a_max, digits=0),
    "F_anode, м²"               => round(F_anode, digits=6),
    "F_cooling_anode, м²"       => round(F_cooling_anode, digits=6),
    "q_a_cooling, Вт/м²"        => round(q_a_cooling, digits=0),
    "G_water, кг/с"             => round(G_water, digits=4)
)

# Добавляем параметры охлаждения, только если решение найдено
if solution_found
    params["Δ охлаждения, мм"]  = round(Δ_solution * 1000, digits=2)
    params["w_water, м/с"]      = round(w_water_solution, digits=2)
    params["Re_water"]          = round(Re_water_solution, digits=0)
    params["α_water, Вт/м²·К"]  = round(α_water_solution, digits=0)
    params["ΔT стенка-вода, К"] = round(ΔT_wall_water_solution, digits=1)
    params["ΔT в стенке, К"]    = round(ΔT_cu_solution, digits=1)
    params["T_wall, К"]         = round(T_wall_solution, digits=1)
    params["T_crit, К"]         = round(T_crit, digits=1)
end

# Запись в файл
open("расчет_плазмотрона.txt", "w") do io
    println(io, "Расчёт плазмотрона")
    println(io, "")
    for (key, val) in sort(collect(params))
        println(io, rpad(key, 28), ": ", val)
    end
end

println("Все параметры записаны в расчет_плазмотрона.txt")