#=
A control to generate a sinthetic case 
=#
linear_interpol = (t, t₁, t₂, v₁, v₂) -> (v₂ - v₁)*(t - t₁)/(t₂ - t₁) + v₁

# Definimos un control lineal por pedazos 
function control_pieces(t)
    # control grande γₑ = 1/7 (con 1/14 funciona feo)
    t₀ = 0.; t₁ = 100.; t₂ = 130.; t₃ = 200.; t₄ = 280.; t₅ = 300.; t₆ = 330.; t₇ = 400.
    v₀ = 1.6e-2; v₁ = 7e-3; v₂ = 1.4e-2; v₃ = 8e-3; v₄ = 1.2e-2; v₅ = 1e-2; 
    
    if t >= t₀ && t < t₁
        linear_interpol(t, t₀, t₁, v₀, v₁)
    elseif t >= t₁ && t < t₂
        linear_interpol(t, t₁, t₂, v₁, v₂)
    elseif t >= t₂ && t < t₃
        linear_interpol(t, t₂, t₃, v₂, v₃)
    elseif t >= t₃ && t < t₄
        v₃
    elseif t >= t₄ && t < t₅
        linear_interpol(t, t₄, t₅, v₃, v₄)
    elseif t >= t₅ && t < t₆
        linear_interpol(t, t₅, t₆, v₄, v₅)
    elseif t >= t₆ 
        v₅
    end
end
